use clap::Parser;
use rsomics_common::CommonFlags;

/// Bayesian SNP/indel caller from mpileup genotype likelihoods (bcftools call -c port).
///
/// Reads a VCF produced by rsomics-vcf-mpileup or bcftools mpileup, applies the
/// consensus Bayesian diploid model, and outputs variant sites with GT/GQ annotations.
#[derive(Parser, Debug)]
#[command(name = "rsomics-vcf-call", version, author)]
pub struct Cli {
    /// Input VCF file (or bgzipped VCF); use '-' for stdin.
    #[arg(default_value = "-")]
    pub input: String,

    /// Write output to FILE instead of stdout.
    #[arg(short = 'o', long)]
    pub output: Option<String>,

    /// Site-frequency prior theta (controls the prior probability of a variant).
    #[arg(short = 'P', long, default_value_t = rsomics_vcf_call::DEFAULT_THETA)]
    pub theta: f64,

    /// Minimum per-sample depth required to report a site.
    #[arg(short = 'd', long, default_value_t = rsomics_vcf_call::DEFAULT_MIN_DEPTH)]
    pub min_depth: u32,

    /// Minimum variant QUAL score to emit a call.
    #[arg(short = 'Q', long, default_value_t = rsomics_vcf_call::DEFAULT_MIN_QUAL)]
    pub min_qual: f32,

    /// Output all sites, including homozygous-reference calls.
    #[arg(short = 'a', long)]
    pub all_sites: bool,

    #[command(flatten)]
    pub common: CommonFlags,
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
