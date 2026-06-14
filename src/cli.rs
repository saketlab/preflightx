//! Native command-line front-end. Each subcommand calls one analysis-module
//! function, then branches on `--json`. Kept out of the wasm build (clap +
//! whitelist I/O are native-only).

use anyhow::Result;
use clap::{Parser, Subcommand};
use std::path::PathBuf;

use crate::{assaydetect, chemistry, composition, fastq, report, sex, size, whitelists};

#[derive(Parser)]
#[command(name = "preflightx", version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,

    /// Number of reads to sample from each fastq (default 50000).
    #[arg(long, global = true, default_value_t = 50_000)]
    n_reads: usize,

    /// Emit JSON instead of the human-readable report.
    #[arg(long, global = true)]
    json: bool,

    /// Species for the sex check (a bundled panel; e.g. human, mouse).
    #[arg(long, global = true, default_value = "human")]
    species: String,

    /// Reads to sample for the sex check (Y/XIST are sparse, so deeper than --n-reads).
    #[arg(long, global = true, default_value_t = 200_000)]
    sex_reads: usize,
}

#[derive(Subcommand)]
enum Cmd {
    /// Run every check (chemistry + bulk-vs-sc + size + composition) and report.
    Check { r1: PathBuf, r2: Option<PathBuf> },
    /// 10x Chromium chemistry detection (R1 length + TSO + optional whitelist).
    Chem { r1: PathBuf, r2: PathBuf },
    /// Bulk RNA-seq vs single-cell library classification.
    Bulk { r1: PathBuf },
    /// Read count, length, and quality sanity.
    Size { r1: PathBuf, r2: Option<PathBuf> },
    /// GC content, low-complexity, ambiguous-base composition.
    Composition { r1: PathBuf },
    /// Sample-level genetic sex from cDNA reads (Y-gene + XIST k-mers; alignment-free).
    /// Works for single-cell, single-end bulk, and paired-end bulk: barcode reads are
    /// skipped and every cDNA read is scanned (both mates for paired-end bulk).
    Sex { r1: PathBuf, r2: Option<PathBuf> },
    /// Generic seqspec-driven detection: rank any catalogued single-cell assay
    /// against the reads (geometry + linker anchors, no whitelist needed).
    DetectAssay { r1: PathBuf, r2: PathBuf },
    /// List the embedded seqspec assay catalogue.
    Assays,
    /// Print the library structure of a catalogued assay (by name or spec file).
    Describe { assay: String },
    /// Download 10x barcode whitelists (one-time, ~180 MB) from 10x's public repo
    /// into ~/.cache/preflightx/whitelists/. Required for full v2/v3/v4 disambiguation.
    FetchWhitelists,
    /// Print the cache directory used for whitelists.
    WhitelistDir,
}

pub fn run() -> Result<()> {
    let cli = Cli::parse();

    match cli.cmd {
        Cmd::Check { r1, r2 } => {
            let cache = whitelists::cache_dir()?;
            let chem = if let Some(ref r2p) = r2 {
                Some(chemistry::detect(&r1, r2p, cli.n_reads, Some(&cache))?)
            } else {
                None
            };
            let bulk = size::bulk_or_sc(&r1, cli.n_reads)?;
            let sz = size::size_sanity(&r1, r2.as_deref(), cli.n_reads)?;
            // Composition reads the cDNA (longer read of a pair); sex auto-selects its
            // own reads (cDNA read for single-cell, both mates for paired-end bulk).
            let cdna_read = match r2.as_deref() {
                Some(r2p) => fastq::longer_read(&r1, r2p, 256)?,
                None => &r1,
            };
            let comp = composition::composition(cdna_read, cli.n_reads)?;
            let sx = sex::detect_sex_auto(&r1, r2.as_deref(), &cli.species, cli.sex_reads)?;
            let r = report::FullReport {
                r1: r1.display().to_string(),
                r2: r2.as_ref().map(|p| p.display().to_string()),
                chemistry: chem,
                bulk_vs_sc: bulk,
                size: sz,
                composition: comp,
                sex: Some(sx),
            };
            if cli.json {
                println!("{}", serde_json::to_string_pretty(&r)?);
            } else {
                r.print();
            }
        }
        Cmd::Chem { r1, r2 } => {
            let cache = whitelists::cache_dir()?;
            let c = chemistry::detect(&r1, &r2, cli.n_reads, Some(&cache))?;
            if cli.json {
                println!("{}", serde_json::to_string_pretty(&c)?);
            } else {
                c.print();
            }
        }
        Cmd::Bulk { r1 } => {
            let b = size::bulk_or_sc(&r1, cli.n_reads)?;
            if cli.json {
                println!("{}", serde_json::to_string_pretty(&b)?);
            } else {
                b.print();
            }
        }
        Cmd::Size { r1, r2 } => {
            let s = size::size_sanity(&r1, r2.as_deref(), cli.n_reads)?;
            if cli.json {
                println!("{}", serde_json::to_string_pretty(&s)?);
            } else {
                s.print();
            }
        }
        Cmd::Composition { r1 } => {
            let c = composition::composition(&r1, cli.n_reads)?;
            if cli.json {
                println!("{}", serde_json::to_string_pretty(&c)?);
            } else {
                c.print();
            }
        }
        Cmd::Sex { r1, r2 } => {
            // Auto-selects cDNA reads: the cDNA read for single-cell, both mates for
            // paired-end bulk, or the sole file for single-end.
            let s = sex::detect_sex_auto(&r1, r2.as_deref(), &cli.species, cli.sex_reads)?;
            if cli.json {
                println!("{}", serde_json::to_string_pretty(&s)?);
            } else {
                s.print();
            }
        }
        Cmd::DetectAssay { r1, r2 } => {
            let d = assaydetect::detect(&r1, &r2, cli.n_reads)?;
            if cli.json {
                println!("{}", serde_json::to_string_pretty(&d)?);
            } else {
                d.print();
            }
        }
        Cmd::Assays => {
            let cat = assaydetect::catalogue();
            if cli.json {
                println!("{}", serde_json::to_string_pretty(&cat)?);
            } else {
                println!("{} assays in catalogue:", cat.len());
                for e in &cat {
                    let det = if e.detectable { "detect" } else { "  --  " };
                    let quant = if e.umi_based { "UMI" } else { "plate" };
                    println!(
                        "  [{}] {:<28} {:<24} {:<5} anchors={} ({})",
                        det,
                        e.display_name,
                        e.spec_file.trim_end_matches(".spec.yaml"),
                        quant,
                        e.anchors,
                        e.modalities.join(",")
                    );
                }
            }
        }
        Cmd::Describe { assay } => match assaydetect::find_assay(&assay) {
            Some(a) => print!("{}", assaydetect::describe(&a)),
            None => {
                eprintln!("no catalogued assay matching '{assay}'; try `preflightx assays`.");
                std::process::exit(1);
            }
        },
        Cmd::FetchWhitelists => {
            whitelists::fetch_all()?;
        }
        Cmd::WhitelistDir => {
            println!("{}", whitelists::cache_dir()?.display());
        }
    }
    Ok(())
}
