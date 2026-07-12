//! Sample-level genetic-sex inference from a cDNA (R2) read, alignment-free.
//!
//! Pure analysis: sample the leading reads, count how many hit a Y-linked gene k-mer
//! (male signal) vs the XIST/Xist k-mer set (female signal), per the bundled
//! species panel ([`crate::sexpanel`]). The call uses the **Y/XIST hit ratio**, a
//! depth- and ambient-robust statistic calibrated on 56 human + 17 mouse real samples
//! (see `experiments/sex_panel/`): true females sit at the ambient floor (ratio ≤ ~0.2),
//! true males well above it (≥ ~0.5). False-male is the costly error (X/Y gametolog
//! homology), so a male call additionally needs an absolute Y-hit floor, and marginal
//! cases honestly report `uncertain`/`unknown` rather than guess.

use anyhow::Result;
use colored::Colorize;
use serde::Serialize;
use std::path::Path;

use crate::fastq;
use crate::libstruct;
use crate::sexpanel::{self, Category, SexPanel};

// Decision thresholds (calibrated on real data; see experiments/sex_panel).
const MIN_HITS: u64 = 8; // minimum (Y + XIST) hit-reads to make any call
const MALE_MIN: u64 = 8; // absolute Y-hit-read floor for a male call
const R_MALE: f64 = 0.50; // Y/XIST >= this (and >= MALE_MIN) -> male
const R_FEMALE: f64 = 0.20; // Y/XIST <= this -> female  (between => uncertain)

fn ratio(y: u64, xist: u64) -> f64 {
    y as f64 / xist.max(1) as f64
}

#[derive(Serialize, Debug)]
pub struct SexReport {
    pub species: String,
    pub panel_status: String, // "validated" | "experimental" | "unsupported"
    pub verdict: String,      // "male" | "female" | "uncertain" | "unknown"
    pub confidence: String,   // "high" | "medium" | "low"
    pub reads_scanned: u64,
    pub y_hit_reads: u64,
    pub xist_hit_reads: u64,
    pub y_xist_ratio: f64,
    /// Per-Y-gene hit-read counts, descending. Signal across *many* Y genes ⇒ a real Y
    /// chromosome; concentration in one gene ⇒ likely stray/ambient.
    pub y_per_gene: Vec<(String, u64)>,
    /// True when both signals are simultaneously strong (pooled both-sex sample, or XXY).
    pub mixed_sex_suspected: bool,
    pub notes: Vec<String>,
}

/// The Y/XIST-ratio decision rule.
pub fn classify(y: u64, xist: u64) -> &'static str {
    if y + xist < MIN_HITS {
        return "unknown";
    }
    let r = ratio(y, xist);
    if r >= R_MALE && y >= MALE_MIN {
        "male"
    } else if r <= R_FEMALE {
        "female"
    } else {
        "uncertain"
    }
}

fn confidence(verdict: &str, y: u64, xist: u64) -> &'static str {
    if verdict == "male" || verdict == "female" {
        let r = ratio(y, xist);
        if y + xist >= 20 && (r >= 1.0 || r <= 0.10) {
            return "high";
        }
        return "medium";
    }
    "low"
}

/// Per-read marker tallies: the pure scoring input, separated from the I/O scan so the
/// scoring can be unit-tested without FASTQ fixtures.
struct Tally {
    scanned: u64,
    too_short: u64,
    y_hit_reads: u64,
    xist_hit_reads: u64,
    gene_reads: Vec<u64>, // indexed by SexPanel::genes
}

/// Accumulate marker hits for one read into `t`. `seen_genes` and
/// `read_genes` are caller-owned scratch reused across reads to avoid allocation.
fn count_read(
    panel: &SexPanel,
    seq: &[u8],
    t: &mut Tally,
    seen_genes: &mut [bool],
    read_genes: &mut Vec<u16>,
) {
    t.scanned += 1;
    read_genes.clear();
    let mut hit_male = false;
    let mut hit_female = false;
    let scanned = crate::sexpanel::scan_encoded_windows(seq, panel.k, |key| {
        if let Some((cat, idx)) = panel.get_encoded(key) {
            match cat {
                Category::Female => hit_female = true,
                Category::Male => {
                    hit_male = true;
                    let slot = idx as usize;
                    if let Some(already_seen) = seen_genes.get_mut(slot) {
                        if !*already_seen {
                            *already_seen = true;
                            read_genes.push(idx);
                        }
                    }
                }
            }
        }
    });
    if !scanned {
        t.too_short += 1;
        return;
    }
    t.y_hit_reads += hit_male as u64;
    t.xist_hit_reads += hit_female as u64;
    for &idx in read_genes.iter() {
        t.gene_reads[idx as usize] += 1;
        seen_genes[idx as usize] = false;
    }
}

fn new_tally(panel: &SexPanel) -> Tally {
    Tally {
        scanned: 0,
        too_short: 0,
        y_hit_reads: 0,
        xist_hit_reads: 0,
        gene_reads: vec![0; panel.genes.len()],
    }
}

/// Stream the leading `n_reads` of `cdna` and accumulate Y / XIST marker hits into
/// `t` (the only I/O). Called once per cDNA read, so paired-end bulk sums both mates.
fn scan_into(panel: &SexPanel, cdna: &Path, n_reads: usize, t: &mut Tally) -> Result<()> {
    let mut seen_genes = vec![false; panel.genes.len()];
    let mut read_genes: Vec<u16> = Vec::with_capacity(panel.genes.len());
    fastq::for_each_record(cdna, n_reads, |seq, _| {
        count_read(panel, seq, t, &mut seen_genes, &mut read_genes);
    })?;
    Ok(())
}

/// Score genetic sex from already-sampled in-memory `reads` (shares the read pass
/// with the species scan in `check`). `panel == None` or an empty (no-sex-pool)
/// panel yields a clean "unavailable" no-call.
pub fn report_for_reads(panel: Option<&SexPanel>, species: &str, reads: &[Vec<u8>]) -> SexReport {
    let panel = match panel {
        Some(p) if !p.is_empty() => p,
        _ => {
            return no_panel(
                species,
                format!("sex unavailable: no XIST/Y panel for '{species}'"),
            )
        }
    };
    if panel.mechanism != "XY_XIST" {
        return no_panel(
            species,
            format!(
                "panel mechanism '{}' is not supported by the current decision rule (only XY_XIST)",
                panel.mechanism
            ),
        );
    }
    let mut t = new_tally(panel);
    let mut seen_genes = vec![false; panel.genes.len()];
    let mut read_genes: Vec<u16> = Vec::with_capacity(panel.genes.len());
    for seq in reads {
        count_read(panel, seq, &mut t, &mut seen_genes, &mut read_genes);
    }
    analyze(panel, species, &t)
}

/// Pure scoring: turn a [`Tally`] into a verdict + per-gene evidence ledger + notes.
fn analyze(panel: &SexPanel, species: &str, t: &Tally) -> SexReport {
    let verdict = classify(t.y_hit_reads, t.xist_hit_reads).to_string();
    let conf = confidence(&verdict, t.y_hit_reads, t.xist_hit_reads).to_string();
    let yx_ratio = ratio(t.y_hit_reads, t.xist_hit_reads);

    let mut y_per_gene: Vec<(String, u64)> = t
        .gene_reads
        .iter()
        .enumerate()
        .filter(|(_, &c)| c > 0)
        .map(|(i, &c)| (panel.genes[i].clone(), c))
        .collect();
    y_per_gene.sort_by(|a, b| b.1.cmp(&a.1));

    // Mixed-sex / XXY: both signals genuinely strong and near-balanced. The XIST floor sits
    // well above the ambient level males carry (~8-50 hits), so a normal male with ambient
    // XIST is not flagged.
    let mixed =
        t.y_hit_reads >= MALE_MIN && t.xist_hit_reads >= 25 && (0.5..=2.0).contains(&yx_ratio);

    let mut notes = Vec::new();
    if !panel.has_y() {
        notes.push(format!(
            "the '{species}' panel has no Y-chromosome markers (Y not assembled/annotated); males are not callable — only female via XIST"
        ));
    }
    if panel.status != "validated" {
        notes.push(format!(
            "panel for '{species}' is {}; interpret with care (see experiments/sex_panel)",
            panel.status
        ));
    }
    if t.scanned > 0 && t.too_short * 2 > t.scanned {
        notes.push(format!(
            "{:.0}% of reads are shorter than k={} (barcode/UMI read?); run on the cDNA read (R2)",
            100.0 * t.too_short as f64 / t.scanned as f64,
            panel.k
        ));
    }
    if verdict == "unknown" && t.too_short * 2 <= t.scanned {
        if t.y_hit_reads + t.xist_hit_reads == 0 && t.scanned >= 1000 {
            notes.push(format!(
                "0 of {} reads matched the '{species}' panel; is the species correct?",
                t.scanned
            ));
        } else {
            notes.push(
                "too few marker hits to call; increase --sex-reads (XIST/Y can be lowly expressed)"
                    .into(),
            );
        }
    }
    if mixed {
        notes.push(
            "both Y and XIST signals strong; pooled both-sex sample or XXY (Klinefelter)?".into(),
        );
    }
    if verdict == "male" && y_per_gene.len() <= 1 {
        notes.push(
            "Y signal concentrated in a single gene; possibly stray/ambient rather than a Y chromosome".into(),
        );
    }

    SexReport {
        species: species.to_string(),
        panel_status: panel.status.clone(),
        verdict,
        confidence: conf,
        reads_scanned: t.scanned,
        y_hit_reads: t.y_hit_reads,
        xist_hit_reads: t.xist_hit_reads,
        y_xist_ratio: fastq::round2(yx_ratio),
        y_per_gene,
        mixed_sex_suspected: mixed,
        notes,
    }
}

/// A no-call report for a species/mechanism without a usable panel.
fn no_panel(species: &str, note: String) -> SexReport {
    SexReport {
        species: species.to_string(),
        panel_status: "unsupported".into(),
        verdict: "unknown".into(),
        confidence: "low".into(),
        reads_scanned: 0,
        y_hit_reads: 0,
        xist_hit_reads: 0,
        y_xist_ratio: 0.0,
        y_per_gene: vec![],
        mixed_sex_suspected: false,
        notes: vec![note],
    }
}

/// Infer genetic sex by accumulating Y/XIST hits across one or more cDNA `reads`,
/// sampling up to `n_reads` leading reads from each. Paired-end bulk passes both
/// mates; single-cell passes the single cDNA read.
pub fn detect_sex_reads(reads: &[&Path], species: &str, n_reads: usize) -> Result<SexReport> {
    let panel = match sexpanel::for_species(species) {
        Some(p) => p,
        None => {
            return Ok(no_panel(
                species,
                format!(
                    "no sex-marker panel for species '{species}'; available: {}",
                    sexpanel::species_list().join(", ")
                ),
            ));
        }
    };
    // The decision rule assumes the XY + XIST mechanism; a panel declaring a different
    // dosage-compensation system (e.g. ZW, where the call logic inverts) is not yet handled.
    if panel.mechanism != "XY_XIST" {
        return Ok(no_panel(
            species,
            format!(
                "panel mechanism '{}' is not supported by the current decision rule (only XY_XIST)",
                panel.mechanism
            ),
        ));
    }
    let mut tally = new_tally(panel);
    for path in reads {
        scan_into(panel, path, n_reads, &mut tally)?;
    }
    let mut report = analyze(panel, species, &tally);
    if reads.len() > 1 {
        report.notes.push(format!(
            "scanned {} cDNA reads (paired-end, no barcode read detected)",
            reads.len()
        ));
    }
    Ok(report)
}

/// Genetic sex for `species`, preferring its cached `.spanel` panel and falling
/// back to an embedded one (human/mouse) — so `sex --species <cached>` works.
#[cfg(not(target_arch = "wasm32"))]
pub fn detect_sex_auto_cat(
    r1: &Path,
    r2: Option<&Path>,
    cat: &crate::speciespanel::Catalogue,
    species: &str,
    n_reads: usize,
) -> Result<SexReport> {
    if cat.is_empty() {
        return detect_sex_auto(r1, r2, species, n_reads);
    }
    let reads = sample_cdna(r1, r2, n_reads)?;
    let panel = cat
        .for_species(species)
        .or_else(|| sexpanel::for_species(species));
    Ok(report_for_reads(panel, species, &reads))
}

/// Infer genetic sex from a single cDNA read.
pub fn detect_sex(cdna: &Path, species: &str, n_reads: usize) -> Result<SexReport> {
    detect_sex_reads(&[cdna], species, n_reads)
}

/// Infer genetic sex from `r1` and an optional `r2`, auto-selecting which read(s)
/// carry cDNA: a fixed-length single-cell barcode/UMI read is skipped, and every
/// other read is scanned. So single-cell scans the cDNA read; single-end and
/// paired-end bulk scan their full read(s).
pub fn detect_sex_auto(
    r1: &Path,
    r2: Option<&Path>,
    species: &str,
    n_reads: usize,
) -> Result<SexReport> {
    let reads = cdna_reads(r1, r2)?;
    detect_sex_reads(&reads, species, n_reads)
}

/// Sample up to `n` reads from the cDNA read(s) of `r1`/`r2` into memory, so the
/// species and sex scans can share one read pass (one download for remote inputs).
pub fn sample_cdna(r1: &Path, r2: Option<&Path>, n: usize) -> Result<Vec<Vec<u8>>> {
    let mut reads = Vec::new();
    for p in cdna_reads(r1, r2)? {
        let (seqs, _) = fastq::sample(p, n)?;
        reads.extend(seqs);
    }
    Ok(reads)
}

/// The cDNA reads among `r1`/`r2`: drops a read whose length profile is a known
/// single-cell barcode geometry. Falls back to `r1` if both somehow look like
/// barcode reads (so a call is still attempted rather than scanning nothing).
pub(crate) fn cdna_reads<'a>(r1: &'a Path, r2: Option<&'a Path>) -> Result<Vec<&'a Path>> {
    let mut out = Vec::new();
    for p in std::iter::once(r1).chain(r2) {
        if !is_barcode_read(p)? {
            out.push(p);
        }
    }
    if out.is_empty() {
        out.push(r1);
    }
    Ok(out)
}

/// Sample `path` and report whether it is a fixed-length barcode/UMI read.
fn is_barcode_read(path: &Path) -> Result<bool> {
    let (seqs, _) = fastq::sample(path, 512)?;
    if seqs.is_empty() {
        return Ok(false);
    }
    let n = seqs.len() as f64;
    let mean = seqs.iter().map(|s| s.len()).sum::<usize>() as f64 / n;
    let var = seqs
        .iter()
        .map(|s| {
            let d = s.len() as f64 - mean;
            d * d
        })
        .sum::<f64>()
        / n;
    Ok(libstruct::is_fixed_barcode_read(mean, var.sqrt()))
}

impl SexReport {
    /// A borderline call that a deeper read sample might resolve.
    pub fn is_uncertain(&self) -> bool {
        self.verdict == "uncertain"
    }

    // Weak = uncertain, or unknown with sub-floor hits. Zero-hit unknown is not weak.
    pub fn is_weak(&self) -> bool {
        if self.verdict == "uncertain" {
            return true;
        }
        let hits = self.y_hit_reads + self.xist_hit_reads;
        self.verdict == "unknown" && hits > 0 && hits < MIN_HITS
    }

    pub fn print(&self) {
        println!("{}", "Genetic sex (sample-level)".bold().cyan());
        let v = match self.verdict.as_str() {
            "male" => self.verdict.bright_blue(),
            "female" => self.verdict.magenta(),
            "uncertain" => self.verdict.yellow(),
            _ => self.verdict.dimmed(),
        };
        let c = match self.confidence.as_str() {
            "high" => self.confidence.green(),
            "medium" => self.confidence.yellow(),
            _ => self.confidence.red(),
        };
        println!(
            "  VERDICT: {}   confidence: {}   (species {}, panel {})",
            v.bold(),
            c,
            self.species,
            self.panel_status
        );
        println!(
            "  Y-panel hits: {}   XIST hits: {}   Y/XIST ratio: {:.2}   reads: {}",
            self.y_hit_reads, self.xist_hit_reads, self.y_xist_ratio, self.reads_scanned
        );
        if !self.y_per_gene.is_empty() {
            let ledger: Vec<String> = self
                .y_per_gene
                .iter()
                .map(|(g, c)| format!("{g}:{c}"))
                .collect();
            println!(
                "  Y genes hit ({}): {}",
                self.y_per_gene.len(),
                ledger.join(", ")
            );
        }
        if self.mixed_sex_suspected {
            println!("  {}", "mixed-sex / XXY suspected".yellow());
        }
        for n in &self.notes {
            println!("    {} {}", "note:".dimmed(), n);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sexpanel::for_species;
    use crate::testutil::write_fastq;

    #[test]
    fn sex_unavailable_when_no_xist_or_y_pool() {
        // A panel with only an autosomal pool (e.g. a bird/fish with no XIST/Y) must
        // report sex unavailable, never guess a verdict.
        let cat = crate::speciespanel::Catalogue::from_texts(&[
            "#species\tGallus_gallus\n#k\t10\nAAAAAAAAAA\tautosomal\n",
        ]);
        let panel = cat.for_species("Gallus_gallus");
        let reads: Vec<Vec<u8>> = (0..10).map(|_| vec![b'A'; 50]).collect();
        let r = report_for_reads(panel, "Gallus_gallus", &reads);
        assert_eq!(r.verdict, "unknown");
        assert_eq!(r.panel_status, "unsupported");
    }

    #[test]
    fn analyze_scores_a_tally_without_io() {
        // The pure scoring seam: no FASTQ needed, just hand-built tallies.
        let panel = for_species("human").unwrap();
        let yi = panel
            .genes
            .iter()
            .position(|g| !g.eq_ignore_ascii_case("xist"))
            .unwrap();
        let mut gene_reads = vec![0u64; panel.genes.len()];
        gene_reads[yi] = 30;
        let t = Tally {
            scanned: 100_000,
            too_short: 0,
            y_hit_reads: 40,
            xist_hit_reads: 1,
            gene_reads,
        };
        let r = analyze(panel, "human", &t);
        assert_eq!(r.verdict, "male");
        assert_eq!(r.y_hit_reads, 40);
        assert_eq!(r.y_per_gene.len(), 1);
    }

    #[test]
    fn zero_hits_hints_species_mismatch() {
        // Zero panel matches across many reads means the wrong species/panel, not low
        // expression, so the note should question the species, not suggest more reads.
        let panel = for_species("human").unwrap();
        let t = Tally {
            scanned: 50_000,
            too_short: 0,
            y_hit_reads: 0,
            xist_hit_reads: 0,
            gene_reads: vec![0; panel.genes.len()],
        };
        let r = analyze(panel, "human", &t);
        assert_eq!(r.verdict, "unknown");
        assert!(
            r.notes.iter().any(|n| n.contains("species")),
            "notes={:?}",
            r.notes
        );
        assert!(
            !r.notes.iter().any(|n| n.contains("sex-reads")),
            "notes={:?}",
            r.notes
        );
    }

    #[test]
    fn classify_rule() {
        assert_eq!(classify(0, 0), "unknown");
        assert_eq!(classify(3, 2), "unknown"); // < MIN_HITS
        assert_eq!(classify(40, 2), "male"); // ratio 20, y>=8
        assert_eq!(classify(2, 40), "female"); // ratio 0.05
        assert_eq!(classify(10, 25), "uncertain"); // ratio 0.4, in band
                                                   // a female with low ambient Y must NOT be called male even if ratio is high but
                                                   // Y count is below the absolute floor
        assert_eq!(classify(5, 0), "unknown"); // total 5 < MIN_HITS
        assert_eq!(classify(6, 2), "uncertain"); // total 8, ratio 3 but y=6 < MALE_MIN -> not male
    }

    #[test]
    fn is_weak_covers_subthreshold_unknown_not_zerohit() {
        let panel = for_species("human").unwrap();
        let mk = |y, x| {
            analyze(
                panel,
                "human",
                &Tally {
                    scanned: 200_000,
                    too_short: 0,
                    y_hit_reads: y,
                    xist_hit_reads: x,
                    gene_reads: vec![0; panel.genes.len()],
                },
            )
        };
        // partial evidence under the floor -> weak (depth-limited, deepen it)
        assert!(
            mk(5, 0).is_weak(),
            "sub-threshold unknown y=5 should be weak"
        );
        assert!(
            mk(3, 2).is_weak(),
            "sub-threshold unknown total=5 should be weak"
        );
        // no evidence at all -> NOT weak (species mismatch; more reads won't help)
        assert!(!mk(0, 0).is_weak(), "zero-hit unknown must not be weak");
        // a confident call is not weak
        assert!(!mk(40, 2).is_weak(), "clean male is not weak");
        // uncertain is still weak
        assert!(mk(10, 25).is_weak(), "uncertain should be weak");
    }

    fn read_with_kmer(kmer: &[u8], total_len: usize) -> Vec<u8> {
        let mut r = vec![b'A'; total_len];
        r[..kmer.len()].copy_from_slice(kmer);
        r
    }

    #[test]
    fn detects_male_from_y_kmers() {
        let panel = for_species("human").unwrap();
        let ymer = panel.sample(Category::Male).unwrap().to_vec();
        let reads: Vec<Vec<u8>> = (0..50).map(|_| read_with_kmer(&ymer, 90)).collect();
        let p = write_fastq("sex_male", &reads);
        let r = detect_sex(&p, "human", 1000).unwrap();
        std::fs::remove_file(&p).ok();
        assert!(r.y_hit_reads >= 40, "y_hit_reads={}", r.y_hit_reads);
        assert_eq!(r.xist_hit_reads, 0);
        assert_eq!(r.verdict, "male");
    }

    #[test]
    fn detects_female_from_xist_kmers() {
        let panel = for_species("human").unwrap();
        let xmer = panel.sample(Category::Female).unwrap().to_vec();
        let reads: Vec<Vec<u8>> = (0..50).map(|_| read_with_kmer(&xmer, 90)).collect();
        let p = write_fastq("sex_female", &reads);
        let r = detect_sex(&p, "human", 1000).unwrap();
        std::fs::remove_file(&p).ok();
        assert!(r.xist_hit_reads >= 40);
        assert_eq!(r.y_hit_reads, 0);
        assert_eq!(r.verdict, "female");
    }

    #[test]
    fn matches_reverse_complement_strand() {
        // A read carrying the revcomp of a Y k-mer must still hit (cDNA strand unknown).
        let panel = for_species("human").unwrap();
        let ymer = sexpanel::revcomp(panel.sample(Category::Male).unwrap());
        let reads: Vec<Vec<u8>> = (0..20).map(|_| read_with_kmer(&ymer, 90)).collect();
        let p = write_fastq("sex_rc", &reads);
        let r = detect_sex(&p, "human", 1000).unwrap();
        std::fs::remove_file(&p).ok();
        assert!(r.y_hit_reads >= 15);
    }

    #[test]
    fn accumulates_marker_hits_over_reads() {
        // Two cDNA files (paired-end bulk) accumulate into one tally.
        let panel = for_species("human").unwrap();
        let ymer = panel.sample(Category::Male).unwrap().to_vec();
        let reads: Vec<Vec<u8>> = (0..30).map(|_| read_with_kmer(&ymer, 90)).collect();
        let p1 = write_fastq("sex_acc1", &reads);
        let p2 = write_fastq("sex_acc2", &reads);
        let r = detect_sex_reads(&[&p1, &p2], "human", 1000).unwrap();
        std::fs::remove_file(&p1).ok();
        std::fs::remove_file(&p2).ok();
        assert!(r.y_hit_reads >= 50, "y_hit_reads={}", r.y_hit_reads);
        assert_eq!(r.verdict, "male");
    }

    #[test]
    fn auto_skips_barcode_read_in_sc_pair() {
        // R1 is a fixed 28 bp barcode-geometry read (10x v3): it must be skipped, and
        // only the cDNA R2 scanned.
        let panel = for_species("human").unwrap();
        let ymer = panel.sample(Category::Male).unwrap().to_vec();
        let r1: Vec<Vec<u8>> = (0..40).map(|_| vec![b'A'; 28]).collect();
        let r2: Vec<Vec<u8>> = (0..40).map(|_| read_with_kmer(&ymer, 90)).collect();
        let p1 = write_fastq("sex_sc_r1", &r1);
        let p2 = write_fastq("sex_sc_r2", &r2);
        let r = detect_sex_auto(&p1, Some(&p2), "human", 1000).unwrap();
        std::fs::remove_file(&p1).ok();
        std::fs::remove_file(&p2).ok();
        assert!(
            r.reads_scanned <= 45,
            "scanned R1 too? reads_scanned={}",
            r.reads_scanned
        );
        assert!(r.y_hit_reads >= 30);
        assert_eq!(r.verdict, "male");
    }

    #[test]
    fn auto_scans_both_bulk_mates() {
        // Paired-end bulk: neither mate is a barcode geometry, so both are scanned.
        let panel = for_species("human").unwrap();
        let ymer = panel.sample(Category::Male).unwrap().to_vec();
        let reads: Vec<Vec<u8>> = (0..30)
            .map(|i| read_with_kmer(&ymer, 80 + i % 16))
            .collect();
        let p1 = write_fastq("sex_bulk_r1", &reads);
        let p2 = write_fastq("sex_bulk_r2", &reads);
        let r = detect_sex_auto(&p1, Some(&p2), "human", 1000).unwrap();
        std::fs::remove_file(&p1).ok();
        std::fs::remove_file(&p2).ok();
        assert!(
            r.reads_scanned >= 55,
            "both mates scanned? reads_scanned={}",
            r.reads_scanned
        );
        assert!(r.y_hit_reads >= 50, "y_hit_reads={}", r.y_hit_reads);
    }

    #[test]
    fn unsupported_species_is_unknown_not_error() {
        let reads = vec![vec![b'A'; 90]; 10];
        let p = write_fastq("sex_nospecies", &reads);
        let r = detect_sex(&p, "tardigrade", 1000).unwrap();
        std::fs::remove_file(&p).ok();
        assert_eq!(r.panel_status, "unsupported");
        assert_eq!(r.verdict, "unknown");
        assert!(!r.notes.is_empty());
    }

    const HUMAN_FEMALE_R2: &str =
        "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/028/SRR31131528/SRR31131528_2.fastq.gz";

    #[test]
    #[ignore = "network: streams R2 from ENA"]
    fn real_female_srp541456_reads_female() {
        let reads = crate::fastq::sample(std::path::Path::new(HUMAN_FEMALE_R2), 50_000)
            .expect("stream R2 from ENA")
            .0;
        let r = report_for_reads(for_species("human"), "human", &reads);
        assert_eq!(r.verdict, "female", "report={r:?}");
    }

    #[test]
    #[ignore = "network: streams from ENA"]
    fn real_macaque_male_reads_male() {
        let cat = crate::speciespanel::Catalogue::load().unwrap();
        let base = "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR215/032/SRR21519032/SRR21519032";
        let (r1, r2) = (format!("{base}_1.fastq.gz"), format!("{base}_2.fastq.gz"));
        let r = detect_sex_auto_cat(
            std::path::Path::new(&r1),
            Some(std::path::Path::new(&r2)),
            &cat,
            "Macaca_fascicularis",
            40_000,
        )
        .unwrap();
        assert_eq!(r.verdict, "male", "report={r:?}");
    }
}
