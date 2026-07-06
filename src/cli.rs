//! Native command-line front-end. Each subcommand calls one analysis-module
//! function, then branches on `--json`. Kept out of the wasm build (clap +
//! whitelist I/O are native-only).

use anyhow::Result;
use clap::{Parser, Subcommand};
use std::path::{Path, PathBuf};

use crate::sex::SexReport;
use crate::{
    assaydetect, chemistry, composition, fastq, nucpanel, report, sex, size, snpanel, species,
    speciespanel, whitelists,
};

/// Sparse-signal detectors (sex, species) under-call on shallow samples. Run `score(n)`,
/// and if the result is still `weak`, resample once at 2x depth before accepting it.
fn deepen_if_weak<T>(
    n: usize,
    weak: impl Fn(&T) -> bool,
    mut score: impl FnMut(usize) -> Result<T>,
) -> Result<T> {
    let first = score(n)?;
    if weak(&first) {
        score(n * 2)
    } else {
        Ok(first)
    }
}

/// Print a cache directory path, as `{"<key>": "<path>"}` under `--json`.
fn print_dir(json: bool, key: &str, dir: &Path) {
    if json {
        println!("{}", serde_json::json!({ key: dir.display().to_string() }));
    } else {
        println!("{}", dir.display());
    }
}

/// Accept `r1,r2` as a single comma-separated positional, or `r1 r2` as two.
/// A lone arg with no comma is single-end (R2 = None).
fn split_pair(r1: String, r2: Option<String>) -> (PathBuf, Option<PathBuf>) {
    match r2 {
        Some(r2) => (r1.into(), Some(r2.into())),
        None => match r1.split_once(',') {
            Some((a, b)) => (a.into(), Some(b.into())),
            None => (r1.into(), None),
        },
    }
}

/// Species + sex from one cDNA read pass (resampled at 2x if either is weak).
/// Shared by `check` and `detect-all`. Degrades to no-species/`--species`-default
/// sex when no panels are cached.
fn species_and_sex(
    r1: &Path,
    r2: Option<&Path>,
    declared: Option<&str>,
    sex_reads: usize,
) -> Result<(Option<species::SpeciesReport>, sex::SexReport)> {
    let catalogue = speciespanel::Catalogue::load()?;
    if catalogue.is_empty() {
        let sp = declared.unwrap_or("human");
        let sx = deepen_if_weak(sex_reads, SexReport::is_uncertain, |n| {
            sex::detect_sex_auto(r1, r2, sp, n)
        })?;
        Ok((None, sx))
    } else {
        // The intron channel rescues single-nucleus species-ID (loaded once; empty if uncached).
        let introns = snpanel::IntronCatalogue::load()
            .unwrap_or_else(|_| snpanel::IntronCatalogue::from_texts(&[]));
        let (sp, sx) = deepen_if_weak(
            sex_reads,
            |(sp, sx): &(species::SpeciesReport, sex::SexReport)| {
                sp.is_unknown() || sx.is_uncertain()
            },
            |n| {
                let reads = sex::sample_cdna(r1, r2, n)?;
                let sp = species::detect_species(&reads, &catalogue, declared);
                let sp = species::apply_intron_channel(sp, &reads, &introns);
                let sx =
                    sex::report_for_reads(catalogue.for_species(&sp.called), &sp.called, &reads);
                Ok((sp, sx))
            },
        )?;
        Ok((Some(sp), sx))
    }
}

/// Species to route the ATAC-vs-RNA panel by: the detected call when confident, else
/// the declared `--species`.
fn atac_species_hint<'a>(
    sp: Option<&'a species::SpeciesReport>,
    declared: Option<&'a str>,
) -> Option<&'a str> {
    sp.filter(|s| !s.is_unknown())
        .map(|s| s.called.as_str())
        .or(declared)
}

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

    /// Also write a self-contained HTML report to this path (`check` only).
    #[arg(long, global = true, value_name = "PATH")]
    html: Option<PathBuf>,

    /// Declared species (e.g. human, mouse, Homo_sapiens).
    #[arg(long, global = true)]
    species: Option<String>,

    /// Reads to sample for the sex check (Y/XIST are sparse, so deeper than --n-reads).
    #[arg(long, global = true, default_value_t = 200_000)]
    sex_reads: usize,

    /// `detect-all`: show the full evidence ledger and notes (default: one line per check).
    #[arg(long, global = true)]
    detailed: bool,
}

#[derive(Subcommand)]
enum Cmd {
    /// Run every check (chemistry + bulk-vs-sc + size + composition) and report.
    Check { r1: PathBuf, r2: Option<PathBuf> },
    /// Detect chemistry + species + sex in one pass. R1/R2 may be two args or a
    /// single comma-separated `r1.fq,r2.fq`. Single-end input skips chemistry.
    DetectAll { r1: String, r2: Option<String> },
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
    /// Closed-set organism ID from cDNA reads (exact k-mer match vs cached panels).
    /// Flags a mislabel when the call disagrees with `--species`.
    Species { r1: PathBuf, r2: Option<PathBuf> },
    /// Single-cell vs single-nucleus from cDNA reads (intron + mitochondrial k-mer
    /// content; alignment-free). Needs the cached nucleus panel (`fetch-panels`).
    Nucleus { r1: PathBuf, r2: Option<PathBuf> },
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
    /// Download the species/sex `.spanel` bundle (~130 MB) into the panel cache.
    FetchPanels,
    /// Print the cache directory used for species/sex panels.
    PanelDir,
}

pub fn run() -> Result<()> {
    let cli = Cli::parse();

    match cli.cmd {
        Cmd::Check { r1, r2 } => {
            let cache = whitelists::cache_dir()?;
            // Species first: routes the sex panel and the chemistry ATAC panel to the
            // called organism; both scan the same cDNA pass.
            let (species_rep, sx) =
                species_and_sex(&r1, r2.as_deref(), cli.species.as_deref(), cli.sex_reads)?;
            let sp_hint = atac_species_hint(species_rep.as_ref(), cli.species.as_deref());
            let (chem, assay_detection) = if let Some(ref r2p) = r2 {
                let (c, a) = report::FullReport::chem_with_fallback(
                    &r1,
                    r2p,
                    cli.n_reads,
                    Some(&cache),
                    sp_hint,
                )?;
                (Some(c), a)
            } else {
                (None, None)
            };
            let bulk = size::bulk_or_sc(&r1, cli.n_reads)?;
            let sz = size::size_sanity(&r1, r2.as_deref(), cli.n_reads)?;
            // Composition reads the cDNA (longer read of a pair).
            let cdna_read = match r2.as_deref() {
                Some(r2p) => fastq::longer_read(&r1, r2p, 256)?,
                None => &r1,
            };
            let comp = composition::composition(cdna_read, cli.n_reads)?;
            // sc-vs-sn + 3'/5' (each None if its panel is uncached; samples cDNA only if a panel is loaded).
            let (nucleus, end) = report::panel_verdicts(&r1, r2.as_deref(), cli.n_reads, sp_hint);
            let r = report::FullReport {
                r1: r1.display().to_string(),
                r2: r2.as_ref().map(|p| p.display().to_string()),
                chemistry: chem,
                assay_detection,
                bulk_vs_sc: bulk,
                size: sz,
                composition: comp,
                species: species_rep,
                sex: Some(sx),
                nucleus,
                end,
            };
            if let Some(path) = &cli.html {
                std::fs::write(path, r.to_html())?;
                eprintln!("wrote HTML report to {}", path.display());
            }
            if cli.json {
                println!("{}", serde_json::to_string_pretty(&r)?);
            } else {
                r.print();
            }
        }
        Cmd::DetectAll { r1, r2 } => {
            let (r1, r2) = split_pair(r1, r2);
            // Species first, so chemistry's ATAC panel can pick the called organism.
            let (species, sex) =
                species_and_sex(&r1, r2.as_deref(), cli.species.as_deref(), cli.sex_reads)?;
            let sp_hint = atac_species_hint(species.as_ref(), cli.species.as_deref());
            let chemistry = r2
                .as_deref()
                .map(|r2p| {
                    let cache = whitelists::cache_dir()?;
                    chemistry::detect(&r1, r2p, cli.n_reads, Some(&cache), sp_hint)
                })
                .transpose()?;
            let r = report::DetectAllReport {
                r1: r1.display().to_string(),
                r2: r2.as_ref().map(|p| p.display().to_string()),
                chemistry,
                species,
                sex,
            };
            match (cli.json, cli.detailed) {
                (true, true) => println!("{}", serde_json::to_string_pretty(&r)?),
                (true, false) => println!("{}", serde_json::to_string_pretty(&r.summary_json())?),
                (false, true) => r.print(),
                (false, false) => r.print_summary(),
            }
        }
        Cmd::Chem { r1, r2 } => {
            let cache = whitelists::cache_dir()?;
            let c = chemistry::detect(&r1, &r2, cli.n_reads, Some(&cache), cli.species.as_deref())?;
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
            let sp = cli.species.as_deref().unwrap_or("human");
            let catalogue = speciespanel::Catalogue::load()?;
            let s = deepen_if_weak(cli.sex_reads, SexReport::is_uncertain, |n| {
                sex::detect_sex_auto_cat(&r1, r2.as_deref(), &catalogue, sp, n)
            })?;
            if cli.json {
                println!("{}", serde_json::to_string_pretty(&s)?);
            } else {
                s.print();
            }
        }
        Cmd::Species { r1, r2 } => {
            let catalogue = speciespanel::Catalogue::load()?;
            let introns = snpanel::IntronCatalogue::load_or_empty();
            let s = deepen_if_weak(cli.sex_reads, species::SpeciesReport::is_unknown, |n| {
                species::detect_species_auto(
                    &r1,
                    r2.as_deref(),
                    &catalogue,
                    &introns,
                    cli.species.as_deref(),
                    n,
                )
            })?;
            if cli.json {
                println!("{}", serde_json::to_string_pretty(&s)?);
            } else {
                s.print();
            }
        }
        Cmd::Nucleus { r1, r2 } => {
            // Nucleus vs cell is a cDNA-read property; sample cDNA as sex/species do.
            let reads = sex::sample_cdna(&r1, r2.as_deref(), cli.n_reads)?;
            let nuc = nucpanel::classify(&reads);
            if cli.json {
                println!("{}", serde_json::to_string_pretty(&nuc)?);
            } else {
                nuc.print();
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
            Some(a) if cli.json => {
                println!(
                    "{}",
                    serde_json::to_string_pretty(&assaydetect::describe_json(&a))?
                )
            }
            Some(a) => print!("{}", assaydetect::describe(&a)),
            None => {
                eprintln!("no catalogued assay matching '{assay}'; try `preflightx assays`.");
                std::process::exit(1);
            }
        },
        Cmd::FetchWhitelists => {
            whitelists::fetch_all()?;
            print_dir(cli.json, "whitelist_dir", &whitelists::cache_dir()?);
        }
        Cmd::WhitelistDir => print_dir(cli.json, "whitelist_dir", &whitelists::cache_dir()?),
        Cmd::FetchPanels => {
            whitelists::fetch_panels()?;
            print_dir(cli.json, "panel_dir", &whitelists::panel_cache_dir()?);
        }
        Cmd::PanelDir => print_dir(cli.json, "panel_dir", &whitelists::panel_cache_dir()?),
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::split_pair;
    use std::path::PathBuf;

    #[test]
    fn split_pair_modes() {
        let p = |s: &str| PathBuf::from(s);
        assert_eq!(
            split_pair("a.fq".into(), Some("b.fq".into())),
            (p("a.fq"), Some(p("b.fq")))
        );
        assert_eq!(
            split_pair("a.fq,b.fq".into(), None),
            (p("a.fq"), Some(p("b.fq")))
        );
        assert_eq!(split_pair("a.fq".into(), None), (p("a.fq"), None));
    }
}
