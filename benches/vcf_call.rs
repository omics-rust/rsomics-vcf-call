use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_vcf_call::{DEFAULT_MIN_DEPTH, DEFAULT_MIN_QUAL, DEFAULT_THETA, call};
use std::io::{BufReader, BufWriter};

static FIXTURE: &str = include_str!("../tests/golden/small.vcf");

fn bench_call(c: &mut Criterion) {
    c.bench_function("vcf_call_small", |b| {
        b.iter(|| {
            let mut out = BufWriter::new(Vec::new());
            call(
                &mut BufReader::new(FIXTURE.as_bytes()),
                &mut out,
                DEFAULT_THETA,
                DEFAULT_MIN_DEPTH,
                DEFAULT_MIN_QUAL,
                false,
            )
            .unwrap();
        });
    });
}

criterion_group!(benches, bench_call);
criterion_main!(benches);
