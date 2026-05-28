# rsomics-vcf-call

Bayesian SNP/indel calling from mpileup genotype likelihoods.

Reads a VCF with `PL` FORMAT fields (as produced by `rsomics-vcf-mpileup` or
`bcftools mpileup`) and applies the consensus diploid Bayesian model to call
variant sites.

## Usage

```
rsomics-vcf-call [OPTIONS] INPUT
rsomics-vcf-mpileup in.bam -f ref.fa | rsomics-vcf-call - -o calls.vcf
```

## Origin

This crate is an independent Rust reimplementation of `bcftools call -c`
(consensus model) based on:

- Li, H. (2011). A statistical framework for SNP calling, mutation discovery,
  association mapping and population genetical parameter estimation from
  sequencing data. *Bioinformatics* 27(21):2987-2993.
  DOI: 10.1093/bioinformatics/btr509
- The VCF 4.2 format specification
- Black-box behaviour testing against bcftools call -c

No GPL source code was used as a reference during implementation.
The bcftools source (MIT license) was consulted for algorithm details only.

License: MIT OR Apache-2.0.
Upstream credit: bcftools <https://github.com/samtools/bcftools> (MIT).
