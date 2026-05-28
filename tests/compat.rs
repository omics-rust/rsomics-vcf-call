use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn ours() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-vcf-call"))
}

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/small.vcf")
}

fn bcftools_version() -> Option<(u32, u32)> {
    let out = Command::new("bcftools").arg("--version").output().ok()?;
    if !out.status.success() {
        return None;
    }
    let text = String::from_utf8_lossy(&out.stdout);
    let first = text.lines().next()?;
    // "bcftools 1.21"
    let ver = first.split_whitespace().nth(1)?;
    let mut parts = ver.split('.');
    let maj: u32 = parts.next()?.parse().ok()?;
    let min: u32 = parts.next()?.parse().ok()?;
    Some((maj, min))
}

/// Data (non-header) lines only — headers carry tool-specific metadata.
fn data_lines(vcf: &[u8]) -> Vec<String> {
    String::from_utf8_lossy(vcf)
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(str::to_owned)
        .collect()
}

/// For each data line keep only CHROM, POS, REF, ALT — bcftools may add
/// INFO fields and differ in QUAL formatting by design (its consensus caller
/// sets QUAL differently). We verify the called sites match, not the exact
/// score encoding.
fn site_keys(lines: &[String]) -> Vec<(String, String, String, String)> {
    lines
        .iter()
        .map(|l| {
            let f: Vec<&str> = l.splitn(10, '\t').collect();
            (
                f.first().unwrap_or(&"").to_string(),
                f.get(1).unwrap_or(&"").to_string(),
                f.get(3).unwrap_or(&"").to_string(),
                f.get(4).unwrap_or(&"").to_string(),
            )
        })
        .collect()
}

#[test]
fn variant_sites_match_bcftools() {
    let (maj, min) = match bcftools_version() {
        Some(v) => v,
        None => {
            eprintln!("SKIP vcf-call compat: bcftools not found");
            return;
        }
    };
    if maj < 1 || (maj == 1 && min < 10) {
        eprintln!("SKIP vcf-call compat: bcftools {maj}.{min} (need >= 1.10)");
        return;
    }

    let fix = fixture();

    // ours: default mode (variant sites only)
    let ours_out = ours().arg(&fix).output().unwrap();
    assert!(
        ours_out.status.success(),
        "rsomics-vcf-call failed: {}",
        String::from_utf8_lossy(&ours_out.stderr)
    );

    // bcftools call -c (consensus model, variant sites only)
    let theirs = Command::new("bcftools")
        .args(["call", "-c"])
        .arg(&fix)
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .output()
        .unwrap();
    if !theirs.status.success() {
        // bcftools call -c can segfault on some platforms (e.g. arm64 macOS
        // with certain VCF formats lacking GL/REF_QUAL tags). Skip gracefully.
        eprintln!(
            "SKIP vcf-call compat: bcftools call returned {}",
            theirs.status
        );
        return;
    }

    let ours_sites = site_keys(&data_lines(&ours_out.stdout));
    let their_sites = site_keys(&data_lines(&theirs.stdout));

    assert_eq!(
        ours_sites, their_sites,
        "called sites differ between rsomics-vcf-call and bcftools call -c"
    );
}

#[test]
fn all_sites_count_matches_input() {
    // --all-sites should emit exactly as many data lines as the input
    let fix = fixture();
    let input_lines = std::fs::read_to_string(&fix)
        .unwrap()
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .count();

    let ours_out = ours().arg("--all-sites").arg(&fix).output().unwrap();
    assert!(
        ours_out.status.success(),
        "rsomics-vcf-call --all-sites failed: {}",
        String::from_utf8_lossy(&ours_out.stderr)
    );

    let out_lines = data_lines(&ours_out.stdout).len();
    assert_eq!(
        out_lines, input_lines,
        "--all-sites emitted {out_lines} data lines but input had {input_lines}"
    );
}
