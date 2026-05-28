//! Bayesian consensus SNP/indel caller from mpileup genotype likelihoods.
//!
//! Implements `bcftools call -c` (consensus Bayesian model) — the original
//! Li 2011 diploid caller that reads Phred-scaled PL values from a VCF
//! produced by `bcftools mpileup` or `rsomics-vcf-mpileup` and outputs only
//! variant sites (by default) or all sites (`--all-sites`).
//!
//! Algorithm (consensus model, `-c`):
//!   For a biallelic site with REF/ALT and PL vector [RR, RA, AA]:
//!     1. Convert PL to likelihoods: L[GT] = 10^(-PL[GT]/10)
//!     2. Apply diploid Hardy-Weinberg prior using theta:
//!        prior[RR] = (1 - AF)^2
//!        prior[RA] = 2 * AF * (1 - AF)
//!        prior[AA] = AF^2
//!        where AF = theta / (1 + theta) with default theta = 1e-3
//!     3. Best GT = argmax(L[GT] * prior[GT])
//!     4. Emit the site if best GT is RA or AA
//!
//! Multi-allelic sites use the same logic over all PL entries.
//!
//! Source reference (MIT-licensed, clean-room algorithm read):
//!   bcftools 1.21 call.c, vcfcall.c (consensus caller `-c`)
//!   Li, H. (2011). A statistical framework for SNP calling, mutation
//!   discovery, association mapping and population genetical parameter
//!   estimation from sequencing data. Bioinformatics 27(21):2987-2993.
//!   DOI: 10.1093/bioinformatics/btr509

#![allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss
)]

use std::io::{BufRead, Write};

use rsomics_common::{Result, RsomicsError};

/// Default site-frequency prior (bcftools `-P theta`, default 1e-3).
pub const DEFAULT_THETA: f64 = 1e-3;
/// Default minimum per-sample depth for a site to be reported.
pub const DEFAULT_MIN_DEPTH: u32 = 1;
/// Default minimum variant quality (QUAL) to emit a variant call.
pub const DEFAULT_MIN_QUAL: f32 = 0.0;

/// Call SNPs/indels from a VCF with PL FORMAT fields.
///
/// Reads `input`, applies the consensus Bayesian model, and writes called
/// (variant) sites to `output`.  If `all_sites` is true, every site is
/// emitted with the called GT/GQ appended.
pub fn call(
    input: &mut dyn BufRead,
    output: &mut dyn Write,
    theta: f64,
    min_depth: u32,
    min_qual: f32,
    all_sites: bool,
) -> Result<()> {
    let af = theta / (1.0 + theta);

    for line_result in input.lines() {
        let line = line_result.map_err(RsomicsError::Io)?;

        // Pass through header lines unchanged
        if line.starts_with('#') {
            writeln!(output, "{line}").map_err(RsomicsError::Io)?;
            continue;
        }

        let fields: Vec<&str> = line.splitn(10, '\t').collect();
        if fields.len() < 8 {
            // Malformed record — pass through
            writeln!(output, "{line}").map_err(RsomicsError::Io)?;
            continue;
        }

        let format_col = if fields.len() > 8 { fields[8] } else { "" };
        let sample_col = if fields.len() > 9 { fields[9] } else { "" };

        // Find PL and DP indices in the FORMAT column
        let pl_idx = format_col.split(':').position(|f| f == "PL");
        let dp_idx = format_col.split(':').position(|f| f == "DP");

        // If no PL field, emit as-is (already-called VCF or symbolic allele)
        if pl_idx.is_none() {
            if all_sites {
                writeln!(output, "{line}").map_err(RsomicsError::Io)?;
            }
            continue;
        }

        let pl_idx = pl_idx.unwrap();

        // Check depth filter: skip site if DP < min_depth when filtering is active.
        if min_depth > 1 {
            let dp = dp_idx
                .and_then(|i| sample_col.split(':').nth(i))
                .and_then(|s| s.parse::<u32>().ok())
                .unwrap_or(0);
            if dp < min_depth {
                continue;
            }
        }

        // Parse PL values from the sample column
        let pl_str = sample_col.split(':').nth(pl_idx).unwrap_or(".");
        if pl_str == "." {
            if all_sites {
                writeln!(output, "{line}").map_err(RsomicsError::Io)?;
            }
            continue;
        }

        let pls: Vec<u32> = pl_str.split(',').filter_map(|s| s.parse().ok()).collect();
        if pls.is_empty() {
            if all_sites {
                writeln!(output, "{line}").map_err(RsomicsError::Io)?;
            }
            continue;
        }

        // Count alleles from ALT field
        let alt = fields[4];
        let n_alleles = if alt == "." || alt.is_empty() {
            1usize
        } else {
            alt.split(',').count() + 1 // REF + ALT alleles
        };

        // Build prior probabilities for each diploid genotype (i, j) i ≤ j
        // P(GT=hom_ref) = (1-af)^2
        // P(GT=het)     = 2 * af * (1-af)  (for each het combination)
        // P(GT=hom_alt) = af^2             (for each hom_alt combination)
        let paf = af;
        let p_ref = 1.0 - paf * (n_alleles as f64 - 1.0);
        let prior_rr = p_ref * p_ref;
        let prior_ra = 2.0 * p_ref * paf;
        let prior_aa = paf * paf;

        let n_pl = pls.len();

        // Compute unnormalised posterior for each genotype
        let mut posteriors = Vec::with_capacity(n_pl);
        let mut gt_idx = 0usize;
        let mut best_post = f64::NEG_INFINITY;
        let mut best_gt_i = 0usize;
        let mut best_gt_j = 0usize;

        'outer: for j in 0..n_alleles {
            for i in 0..=j {
                if gt_idx >= n_pl {
                    break 'outer;
                }
                let pl = pls[gt_idx];
                let likelihood = 10.0_f64.powf(-(pl as f64) / 10.0);

                // Prior: RR vs het vs hom_alt
                let prior = if i == 0 && j == 0 {
                    prior_rr
                } else if i == 0 {
                    // ref/alt het
                    prior_ra
                } else if i == j {
                    // hom alt
                    prior_aa
                } else {
                    // alt/alt het (multi-allelic)
                    prior_aa * 0.5
                };

                let post = likelihood * prior;
                posteriors.push(post);

                if post > best_post {
                    best_post = post;
                    best_gt_i = i;
                    best_gt_j = j;
                }

                gt_idx += 1;
            }
        }

        let is_variant = best_gt_i != 0 || best_gt_j != 0;

        // Compute genotype quality (GQ) = -10 log10(1 - P(best_gt) / sum)
        let total_post: f64 = posteriors.iter().sum();
        let gq = if total_post > 0.0 {
            let p_best = best_post / total_post;
            let p_err = 1.0 - p_best;
            if p_err <= 0.0 {
                99.0_f32
            } else {
                (-10.0_f32 * (p_err as f32).log10()).min(99.0)
            }
        } else {
            0.0_f32
        };

        // QUAL = GQ of the best non-ref call; 0 for homref
        let qual = if is_variant { gq } else { 0.0 };

        if !is_variant && !all_sites {
            continue;
        }

        if is_variant && qual < min_qual {
            continue;
        }

        // Build GT string: "0/0", "0/1", "1/1", etc.
        let gt_str = format!("{best_gt_i}/{best_gt_j}");

        // Rebuild sample column with GT prepended (or replacing existing GT)
        let format_fields: Vec<&str> = format_col.split(':').collect();
        let sample_fields: Vec<&str> = sample_col.split(':').collect();

        let has_gt = format_fields.first() == Some(&"GT");
        let (new_format, new_sample) = if has_gt {
            // Replace existing GT and append GQ if not present
            let gq_idx = format_fields.iter().position(|&f| f == "GQ");
            if let Some(gi) = gq_idx {
                let mut sf: Vec<String> = sample_fields.iter().map(|s| s.to_string()).collect();
                sf[0] = gt_str.clone();
                sf[gi] = format!("{gq:.0}");
                (format_col.to_string(), sf.join(":"))
            } else {
                let mut sf: Vec<String> = sample_fields.iter().map(|s| s.to_string()).collect();
                sf[0] = gt_str.clone();
                (format_col.to_string(), sf.join(":"))
            }
        } else {
            // Prepend GT:GQ to FORMAT and sample
            let new_fmt = format!("GT:GQ:{format_col}");
            let new_smp = format!("{gt_str}:{gq:.0}:{sample_col}");
            (new_fmt, new_smp)
        };

        // Reconstruct the QUAL field
        let qual_str = if is_variant {
            format!("{qual:.1}")
        } else {
            ".".to_string()
        };

        writeln!(
            output,
            "{chrom}\t{pos}\t{id}\t{ref_}\t{alt}\t{qual}\t{filter}\t{info}\t{fmt}\t{smp}",
            chrom = fields[0],
            pos = fields[1],
            id = fields[2],
            ref_ = fields[3],
            alt = fields[4],
            qual = qual_str,
            filter = fields[6],
            info = fields[7],
            fmt = new_format,
            smp = new_sample,
        )
        .map_err(RsomicsError::Io)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    const HEADER: &str = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";

    fn call_str(vcf: &str, all_sites: bool) -> String {
        let mut out = Vec::new();
        call(
            &mut BufReader::new(vcf.as_bytes()),
            &mut out,
            DEFAULT_THETA,
            DEFAULT_MIN_DEPTH,
            DEFAULT_MIN_QUAL,
            all_sites,
        )
        .unwrap();
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn header_passthrough() {
        let vcf = HEADER.to_owned();
        let out = call_str(&vcf, false);
        assert!(out.contains("##fileformat=VCFv4.2"));
        assert!(out.contains("#CHROM"));
    }

    #[test]
    fn hom_ref_not_emitted_by_default() {
        // PL = 0,10,100 → REF/REF is most likely
        let vcf = format!("{HEADER}chr1\t100\t.\tA\tG\t.\t.\t.\tDP:PL\t30:0,10,100\n");
        let out = call_str(&vcf, false);
        // Variant line should not appear (hom ref)
        let data_lines: Vec<&str> = out.lines().filter(|l| !l.starts_with('#')).collect();
        assert!(data_lines.is_empty(), "hom-ref should not be emitted");
    }

    #[test]
    fn hom_ref_emitted_all_sites() {
        let vcf = format!("{HEADER}chr1\t100\t.\tA\tG\t.\t.\t.\tDP:PL\t30:0,10,100\n");
        let out = call_str(&vcf, true);
        let data_lines: Vec<&str> = out.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(data_lines.len(), 1);
    }

    #[test]
    fn het_call_emitted() {
        // PL = 100,0,80 → REF/ALT het is most likely
        let vcf = format!("{HEADER}chr1\t200\t.\tA\tG\t.\t.\t.\tDP:PL\t30:100,0,80\n");
        let out = call_str(&vcf, false);
        let data_lines: Vec<&str> = out.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(data_lines.len(), 1);
        assert!(data_lines[0].contains("0/1"), "should call 0/1 het");
    }

    #[test]
    fn hom_alt_call_emitted() {
        // PL = 100,80,0 → ALT/ALT hom_alt is most likely
        let vcf = format!("{HEADER}chr1\t300\t.\tA\tG\t.\t.\t.\tDP:PL\t30:100,80,0\n");
        let out = call_str(&vcf, false);
        let data_lines: Vec<&str> = out.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(data_lines.len(), 1);
        assert!(data_lines[0].contains("1/1"), "should call 1/1 hom_alt");
    }
}
