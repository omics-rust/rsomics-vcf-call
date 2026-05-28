# rsomics-vcf-call

Bayesian SNP/indel caller from mpileup genotype likelihoods.

Implements `bcftools call -c` — the consensus Bayesian model originally described by Li 2011 — reading Phred-scaled PL values from a VCF produced by `rsomics-vcf-mpileup` or `bcftools mpileup` and emitting variant sites with GT and GQ annotations.

## Usage

```
rsomics-vcf-call [OPTIONS] [INPUT]

Arguments:
  [INPUT]  Input VCF (or bgzipped VCF); use '-' for stdin [default: -]

Options:
  -o, --output <FILE>       Write output to FILE instead of stdout
  -P, --theta <THETA>       Site-frequency prior [default: 0.001]
  -d, --min-depth <N>       Minimum per-sample depth [default: 1]
  -Q, --min-qual <QUAL>     Minimum variant QUAL to emit [default: 0.0]
  -a, --all-sites           Output all sites including hom-ref calls
  -t, --threads <N>         Worker threads [default: 1]
  -h, --help                Print help
  -V, --version             Print version
```

## Algorithm

For a biallelic site with REF/ALT and PL vector `[RR, RA, AA]`:

1. Convert PL to likelihoods: `L[GT] = 10^(-PL[GT]/10)`
2. Apply diploid Hardy-Weinberg prior using `theta` (default 1e-3):
   - `prior[RR] = (1 - af)^2`
   - `prior[RA] = 2 * af * (1 - af)`
   - `prior[AA] = af^2`
   where `af = theta / (1 + theta)`
3. Best GT = argmax(L[GT] × prior[GT])
4. Emit site if best GT is RA or AA; skip if RR (unless `--all-sites`)
5. GQ = −10 log10(1 − P(best GT) / ΣP) capped at 99

Multi-allelic sites use the same logic over all PL entries enumerated in diploid triangular order.

## Origin

This crate is a clean-room Rust reimplementation of the `bcftools call -c`
consensus caller, based on:

- Li, H. (2011). A statistical framework for SNP calling, mutation discovery,
  association mapping and population genetical parameter estimation from
  sequencing data. *Bioinformatics* 27(21):2987–2993.
  DOI: [10.1093/bioinformatics/btr509](https://doi.org/10.1093/bioinformatics/btr509)
- The public VCF 4.2 format specification.
- Black-box behaviour testing against `bcftools call -c` (MIT-licensed; source
  was also consulted for algorithm details under the MIT licence terms).

Upstream credit: [bcftools](https://github.com/samtools/bcftools) (MIT/GPL-3).

License: MIT OR Apache-2.0.
