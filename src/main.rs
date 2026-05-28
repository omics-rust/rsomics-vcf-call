mod cli;

use std::io::{BufReader, BufWriter};
use std::process::ExitCode;

use clap::Parser;
use rsomics_vcf_call::call;

use cli::Cli;

fn main() -> ExitCode {
    let args = Cli::parse();

    let mut inp: Box<dyn std::io::Read> = if args.input == "-" {
        Box::new(std::io::stdin().lock())
    } else {
        match std::fs::File::open(&args.input) {
            Ok(f) => Box::new(f),
            Err(e) => {
                eprintln!("rsomics-vcf-call: cannot open '{}': {e}", args.input);
                return ExitCode::FAILURE;
            }
        }
    };

    let mut out: Box<dyn std::io::Write> = match &args.output {
        Some(path) => match std::fs::File::create(path) {
            Ok(f) => Box::new(BufWriter::new(f)),
            Err(e) => {
                eprintln!("rsomics-vcf-call: cannot create '{path}': {e}");
                return ExitCode::FAILURE;
            }
        },
        None => Box::new(BufWriter::new(std::io::stdout().lock())),
    };

    let mut reader = BufReader::new(&mut *inp);
    match call(
        &mut reader,
        &mut *out,
        args.theta,
        args.min_depth,
        args.min_qual,
        args.all_sites,
    ) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("rsomics-vcf-call: {e}");
            ExitCode::FAILURE
        }
    }
}
