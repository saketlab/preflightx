use serde::Serialize;

pub fn canonical_species(name: &str) -> String {
    let n = name.trim().to_ascii_lowercase().replace(' ', "_");
    match n.as_str() {
        "human" | "homo_sapiens" => "homo_sapiens".to_string(),
        "mouse" | "mus_musculus" => "mus_musculus".to_string(),
        _ => n,
    }
}

pub fn same_species(a: &str, b: &str) -> bool {
    canonical_species(a) == canonical_species(b)
}

#[derive(Serialize, Debug)]
pub struct SpeciesReport {
    pub declared: Option<String>,
    pub called: String,
    /// "high" | "medium" | "low" — how much to trust `called`. Always "low" when
    /// `called == "unknown"`.
    pub confidence: String,
    /// True when the runner-up is within `AMBIG_FRACTION` of the top species; `called`
    /// is still the top but the tie should be treated as suspect.
    pub ambiguous: bool,
    /// Why `called == "unknown"` (below hit floor / no catalogue); None when called.
    pub abstained_reason: Option<String>,
    pub called_hit_reads: u64,
    pub called_hit_fraction: f64,
    pub runner_up: Option<String>,
    pub runner_up_hit_reads: u64,
    pub margin: f64,
    pub reads_scanned: u64,
    pub reads_with_hit: u64,
    pub n_species_in_panel: usize,
    pub mislabel: bool,
    pub notes: Vec<String>,
}

#[cfg(not(target_arch = "wasm32"))]
mod detect {
    use super::*;
    use crate::speciespanel::Catalogue;

    const MIN_HIT_READS: u64 = 10;
    const MIN_HIT_FRACTION: f64 = 0.0005;
    // Hit-fraction floor for "high" (with a clear lead); clearing only MIN is "low".
    const HIGH_HIT_FRACTION: f64 = 0.005; // 10x the calling floor
                                          // Runner-up within this fraction of the top's hit-reads -> `ambiguous`.
    const AMBIG_FRACTION: f64 = 0.80;
    // Top must lead the runner-up by at least this ratio for "high".
    const HIGH_LEAD_RATIO: f64 = 2.0;
    // Lead-ratio path: this many absolute hits with a decisive lead call (and reach "medium")
    // even at a low hit-fraction, for intron-heavy single-nucleus RNA.
    const MIN_LEAD_HITS: u64 = 30;
    const STRONG_LEAD_RATIO: f64 = 3.0;

    fn decisive_lead(top_c: u64, runner_c: u64) -> bool {
        top_c >= MIN_LEAD_HITS
            && (runner_c == 0 || top_c as f64 >= STRONG_LEAD_RATIO * runner_c as f64)
    }

    /// Bucket the call's trustworthiness from its absolute hit-fraction and its lead
    /// over the runner-up. `unknown`/ambiguous calls never reach "high".
    fn classify_confidence(
        top_c: u64,
        runner_c: u64,
        top_frac: f64,
        ambiguous: bool,
    ) -> &'static str {
        if ambiguous {
            return "low";
        }
        let lead_ok = runner_c == 0 || (top_c as f64) >= HIGH_LEAD_RATIO * runner_c as f64;
        if top_frac >= HIGH_HIT_FRACTION && lead_ok {
            "high"
        } else if top_frac >= MIN_HIT_FRACTION * 2.0 || decisive_lead(top_c, runner_c) {
            "medium"
        } else {
            "low"
        }
    }

    /// Per-species hit-read counts over `reads`: each read votes at most once per species.
    /// The shared scan seam for both the exon (`.spanel`) and intron (`.snpanel`) channels.
    pub(crate) fn count_hits(
        reads: &[Vec<u8>],
        k: usize,
        n_species: usize,
        hit_of: impl Fn(u128) -> Option<u16>,
    ) -> (Vec<u64>, u64, u64) {
        let mut counts = vec![0u64; n_species];
        let mut scanned = 0u64;
        let mut reads_with_hit = 0u64;
        let mut hit: Vec<u16> = Vec::new();
        let mut seen = vec![false; n_species];
        for read in reads {
            scanned += 1;
            hit.clear();
            crate::sexpanel::scan_encoded_windows(read, k, |key| {
                if let Some(idx) = hit_of(key) {
                    if let Some(already) = seen.get_mut(idx as usize) {
                        if !*already {
                            *already = true;
                            hit.push(idx);
                        }
                    }
                }
            });
            if !hit.is_empty() {
                reads_with_hit += 1;
                for &i in &hit {
                    counts[i as usize] += 1;
                    seen[i as usize] = false;
                }
            }
        }
        (counts, scanned, reads_with_hit)
    }

    /// Top and runner-up from per-species counts: `(top_idx, top_count, runner_name, runner_count)`.
    fn rank(counts: &[u64], name_of: impl Fn(u16) -> String) -> (u16, u64, Option<String>, u64) {
        let mut ranked: Vec<(u16, u64)> = counts
            .iter()
            .enumerate()
            .map(|(i, &c)| (i as u16, c))
            .filter(|(_, c)| *c > 0)
            .collect();
        ranked.sort_by(|a, b| b.1.cmp(&a.1));
        let (top_idx, top_c) = ranked.first().copied().unwrap_or((0, 0));
        let (runner_up, runner_c) = match ranked.get(1) {
            Some(&(i, c)) => (Some(name_of(i)), c),
            None => (None, 0),
        };
        (top_idx, top_c, runner_up, runner_c)
    }

    fn round4(x: f64) -> f64 {
        crate::fastq::round2(x * 100.0) / 100.0
    }

    pub fn detect_species(
        reads: &[Vec<u8>],
        cat: &Catalogue,
        declared: Option<&str>,
    ) -> SpeciesReport {
        let declared = declared.map(str::to_string);
        if cat.is_empty() {
            return no_catalogue(declared);
        }
        let (counts, scanned, reads_with_hit) = count_hits(reads, cat.k(), cat.len(), |key| {
            cat.autosomal_hit_encoded(key)
        });
        let (top_idx, top_c, runner_up, runner_c) =
            rank(&counts, |i| cat.species_of(i).to_string());

        let scanned_f = scanned.max(1) as f64;
        let top_frac = top_c as f64 / scanned_f;
        let margin = (top_c.saturating_sub(runner_c)) as f64 / scanned_f;

        let decisive = decisive_lead(top_c, runner_c);
        let sparse_but_decisive = decisive && top_frac < MIN_HIT_FRACTION;
        let cleared_floor = (top_c >= MIN_HIT_READS && top_frac >= MIN_HIT_FRACTION) || decisive;
        let called = if cleared_floor {
            cat.species_of(top_idx).to_string()
        } else {
            "unknown".to_string()
        };

        let ambiguous =
            cleared_floor && top_c > 0 && runner_c as f64 >= AMBIG_FRACTION * top_c as f64;
        let confidence = if called == "unknown" {
            "low"
        } else {
            classify_confidence(top_c, runner_c, top_frac, ambiguous)
        }
        .to_string();
        let abstained_reason = if called == "unknown" {
            Some(if top_c == 0 {
                "no catalogued species matched any read".to_string()
            } else {
                format!(
                    "top species hit only {top_c} reads ({:.3}% < {:.3}% floor); off-panel organism or too-shallow/non-cDNA data",
                    top_frac * 100.0,
                    MIN_HIT_FRACTION * 100.0
                )
            })
        } else {
            None
        };

        let mislabel = match &declared {
            Some(d) if called != "unknown" => !same_species(d, &called),
            _ => false,
        };

        let mut notes = Vec::new();
        notes.push(format!(
            "closed-set ID over {} catalogued species; an organism absent from the panel reads as 'unknown'",
            cat.len()
        ));
        if called == "unknown" {
            notes.push(
                "no species cleared the hit threshold; off-panel organism, very shallow data, or non-cDNA reads".into(),
            );
        }
        if sparse_but_decisive {
            notes.push(format!(
                "few exonic hits ({:.3}% of reads) but a decisive {:.0}x lead over the runner-up; typical of intron-heavy single-nucleus RNA",
                top_frac * 100.0,
                top_c as f64 / runner_c.max(1) as f64
            ));
        }
        if ambiguous {
            notes.push(format!(
                "runner-up '{}' is within {:.0}% of the top species — closely-related panel species or a mixture; treat the call as tentative",
                runner_up.as_deref().unwrap_or("?"),
                AMBIG_FRACTION * 100.0
            ));
        }
        if mislabel {
            notes.push(format!(
                "declared '{}' but reads look like '{}' — possible mislabel",
                declared.as_deref().unwrap_or(""),
                called
            ));
        }

        SpeciesReport {
            declared,
            called,
            confidence,
            ambiguous,
            abstained_reason,
            called_hit_reads: top_c,
            called_hit_fraction: round4(top_frac),
            runner_up,
            runner_up_hit_reads: runner_c,
            margin: round4(margin),
            reads_scanned: scanned,
            reads_with_hit,
            n_species_in_panel: cat.len(),
            mislabel,
            notes,
        }
    }

    // Intron-channel thresholds: absolute hit-reads + lead, not the exon hit-fraction floor.
    const SN_MIN_HITS: u64 = 25; // intron hit-reads floor to call from the intron channel
    const SN_LEAD_RATIO: f64 = 3.0; // top must beat the runner-up by this
    const SN_HIGH_HITS: u64 = 75;
    const SN_HIGH_FRACTION: f64 = 0.0015; // 0.15% of reads

    fn sn_decisive(top_c: u64, runner_c: u64) -> bool {
        top_c >= SN_MIN_HITS && (runner_c == 0 || top_c as f64 >= SN_LEAD_RATIO * runner_c as f64)
    }

    /// Rescue a weak exon species call with a decisive intron-channel call. `exon` is the
    /// [`detect_species`] result over the same reads; high-confidence exon calls and an empty
    /// intron catalogue pass through untouched.
    pub fn apply_intron_channel(
        exon: SpeciesReport,
        reads: &[Vec<u8>],
        intron: &crate::snpanel::IntronCatalogue,
    ) -> SpeciesReport {
        if exon.confidence == "high" || intron.is_empty() {
            return exon;
        }
        let (counts, scanned, reads_with_hit) =
            count_hits(reads, intron.k(), intron.len(), |key| {
                intron.hit_encoded(key)
            });
        let (top_idx, top_c, runner_up, runner_c) =
            rank(&counts, |i| intron.species_of(i).to_string());
        if !sn_decisive(top_c, runner_c) {
            return exon;
        }
        let species = intron.species_of(top_idx).to_string();

        // Intron call disagrees with a weak exon call: two organisms, report unknown.
        if exon.called != "unknown" && !same_species(&exon.called, &species) {
            let mut r = exon;
            let ex = r.called.clone();
            r.called = "unknown".into();
            r.confidence = "low".into();
            r.abstained_reason = Some("exon and intron channels disagree".into());
            r.notes.push(format!(
                "exon reads look like '{ex}' but intron reads like '{species}' — possible mixture or contamination"
            ));
            return r;
        }

        let scanned_f = scanned.max(1) as f64;
        let top_frac = top_c as f64 / scanned_f;
        let mut r = exon;
        r.mislabel = matches!(&r.declared, Some(d) if !same_species(d, &species));
        r.confidence = if top_c >= SN_HIGH_HITS && top_frac >= SN_HIGH_FRACTION {
            "high"
        } else {
            "medium"
        }
        .into();
        r.notes = vec![format!(
            "called via the single-nucleus intron channel: {top_c} intron hit-reads ({:.3}% of reads), {:.0}x lead over the runner-up",
            top_frac * 100.0,
            top_c as f64 / runner_c.max(1) as f64
        )];
        if r.mislabel {
            r.notes.push(format!(
                "declared '{}' but reads look like '{species}' — possible mislabel",
                r.declared.as_deref().unwrap_or("")
            ));
        }
        r.called = species;
        r.ambiguous = false;
        r.abstained_reason = None;
        r.called_hit_reads = top_c;
        r.called_hit_fraction = round4(top_frac);
        r.runner_up = runner_up;
        r.runner_up_hit_reads = runner_c;
        r.margin = round4((top_c.saturating_sub(runner_c)) as f64 / scanned_f);
        r.reads_scanned = scanned;
        r.reads_with_hit = reads_with_hit;
        r.n_species_in_panel = intron.len();
        r
    }

    pub fn detect_species_auto(
        r1: &std::path::Path,
        r2: Option<&std::path::Path>,
        cat: &Catalogue,
        intron: &crate::snpanel::IntronCatalogue,
        declared: Option<&str>,
        n: usize,
    ) -> anyhow::Result<SpeciesReport> {
        if cat.is_empty() {
            return Ok(detect_species(&[], cat, declared));
        }
        let reads = crate::sex::sample_cdna(r1, r2, n)?;
        let exon = detect_species(&reads, cat, declared);
        Ok(apply_intron_channel(exon, &reads, intron))
    }

    fn no_catalogue(declared: Option<String>) -> SpeciesReport {
        SpeciesReport {
            declared,
            called: "unknown".into(),
            confidence: "low".into(),
            ambiguous: false,
            abstained_reason: Some("no species panels cached".into()),
            called_hit_reads: 0,
            called_hit_fraction: 0.0,
            runner_up: None,
            runner_up_hit_reads: 0,
            margin: 0.0,
            reads_scanned: 0,
            reads_with_hit: 0,
            n_species_in_panel: 0,
            mislabel: false,
            notes: vec![
                "no species panels cached; run `preflightx fetch-panels` (or copy *.spanel into the panel cache dir)".into(),
            ],
        }
    }
}

#[cfg(not(target_arch = "wasm32"))]
pub use detect::{apply_intron_channel, detect_species, detect_species_auto};

impl SpeciesReport {
    /// No catalogued species cleared the hit threshold — off-panel, or too few reads.
    pub fn is_unknown(&self) -> bool {
        self.called == "unknown"
    }

    pub fn print(&self) {
        use colored::Colorize;
        println!("{}", "Species identification".bold().cyan());
        let called = if self.called == "unknown" {
            self.called.dimmed()
        } else {
            self.called.green().bold()
        };
        print!("  CALLED: {called}");
        let conf = match self.confidence.as_str() {
            "high" => self.confidence.green(),
            "medium" => self.confidence.yellow(),
            _ => self.confidence.red(),
        };
        print!("   confidence: {conf}");
        if self.ambiguous {
            print!("   {}", "[ambiguous]".yellow());
        }
        if let Some(d) = &self.declared {
            print!("   (declared: {d})");
        }
        println!();
        println!(
            "  hit-reads: {} / {} ({:.2}% of reads)   panel: {} species",
            self.called_hit_reads,
            self.reads_scanned,
            self.called_hit_fraction * 100.0,
            self.n_species_in_panel
        );
        if let Some(r) = &self.runner_up {
            println!(
                "  runner-up: {} ({} hit-reads)   margin: {:.2}%",
                r,
                self.runner_up_hit_reads,
                self.margin * 100.0
            );
        }
        if self.mislabel {
            println!(
                "  {}",
                "MISLABEL: declared species disagrees with reads"
                    .red()
                    .bold()
            );
        }
        for n in &self.notes {
            println!("    {} {}", "note:".dimmed(), n);
        }
    }
}

#[cfg(test)]
#[cfg(not(target_arch = "wasm32"))]
mod tests {
    use super::*;

    #[test]
    fn aliases_fold_common_and_latin_names() {
        assert!(same_species("human", "Homo_sapiens"));
        assert!(same_species("Mus musculus", "mouse"));
        assert!(!same_species("human", "Mus_musculus"));
    }

    fn tiny_catalogue() -> crate::speciespanel::Catalogue {
        crate::speciespanel::Catalogue::from_texts(&[
            "#species\tHomo_sapiens\n#k\t10\nAAAAAAAAAA\tautosomal\nACGTACGTAC\tautosomal\n",
            "#species\tMus_musculus\n#k\t10\nGGGAGGGAGG\tautosomal\nTGTGTGTGTG\tautosomal\n",
        ])
    }

    #[test]
    fn calls_top_species_and_flags_mislabel() {
        let cat = tiny_catalogue();
        assert_eq!(cat.len(), 2);
        let mut read = b"AAAAAAAAAA".to_vec();
        read.extend(vec![b'N'; 40]); // pad; Ns never match
        let reads: Vec<Vec<u8>> = (0..50).map(|_| read.clone()).collect();

        let r = detect_species(&reads, &cat, Some("Mus_musculus"));
        assert_eq!(r.called, "Homo_sapiens");
        assert_eq!(r.called_hit_reads, 50);
        assert!(r.mislabel, "expected mislabel, got {r:?}");

        let r2 = detect_species(&reads, &cat, Some("human"));
        assert!(!r2.mislabel);
    }

    #[test]
    fn sparse_but_decisive_lead_calls_at_medium() {
        let cat = tiny_catalogue();
        let pad = |kmer: &[u8]| {
            let mut r = kmer.to_vec();
            r.extend(vec![b'N'; 40]);
            r
        };
        let mut reads: Vec<Vec<u8>> = (0..40).map(|_| pad(b"AAAAAAAAAA")).collect();
        reads.extend((0..5).map(|_| pad(b"GGGAGGGAGG")));
        reads.extend((0..100_000).map(|_| vec![b'N'; 50]));
        let r = detect_species(&reads, &cat, None);
        assert_eq!(r.called, "Homo_sapiens", "{r:?}");
        assert_eq!(r.confidence, "medium");
        assert!(r.notes.iter().any(|n| n.contains("single-nucleus")));
    }

    #[test]
    fn intron_channel_rescues_weak_exon_call() {
        let exon_cat = tiny_catalogue();
        let intron = crate::snpanel::IntronCatalogue::from_texts(&[
            "#species\tFelis_catus\n#k\t10\nATCGATCGAT\tACTB\tdeep_intron\n",
            "#species\tBos_taurus\n#k\t10\nGGGGTTTTCC\tGAPDH\tdeep_intron\n",
        ]);
        let felis_read = {
            let mut r = b"ATCGATCGAT".to_vec();
            r.extend(vec![b'N'; 40]);
            r
        };
        let reads: Vec<Vec<u8>> = (0..90).map(|_| felis_read.clone()).collect();

        let exon = detect_species(&reads, &exon_cat, None);
        assert_eq!(exon.called, "unknown"); // no exon k-mer hit
        let r = apply_intron_channel(exon, &reads, &intron);
        assert_eq!(r.called, "Felis_catus", "{r:?}");
        assert_eq!(r.confidence, "high");
        assert!(r.notes.iter().any(|n| n.contains("intron channel")));
    }

    #[test]
    fn intron_channel_reports_mixture_on_disagreement() {
        let intron = crate::snpanel::IntronCatalogue::from_texts(&[
            "#species\tFelis_catus\n#k\t10\nATCGATCGAT\tACTB\tdeep_intron\n",
        ]);
        let felis_read = {
            let mut r = b"ATCGATCGAT".to_vec();
            r.extend(vec![b'N'; 40]);
            r
        };
        let reads: Vec<Vec<u8>> = (0..60).map(|_| felis_read.clone()).collect();
        let mut exon = super::SpeciesReport {
            declared: None,
            called: "Bos_taurus".into(),
            confidence: "low".into(),
            ambiguous: false,
            abstained_reason: None,
            called_hit_reads: 12,
            called_hit_fraction: 0.0,
            runner_up: None,
            runner_up_hit_reads: 0,
            margin: 0.0,
            reads_scanned: 60,
            reads_with_hit: 12,
            n_species_in_panel: 2,
            mislabel: false,
            notes: vec![],
        };
        exon.called = "Bos_taurus".into();
        let r = apply_intron_channel(exon, &reads, &intron);
        assert_eq!(r.called, "unknown", "{r:?}");
        assert!(r.notes.iter().any(|n| n.contains("mixture")));
    }

    #[test]
    fn off_panel_reads_are_unknown() {
        let cat = tiny_catalogue();
        let reads: Vec<Vec<u8>> = (0..50).map(|_| vec![b'N'; 50]).collect();
        let r = detect_species(&reads, &cat, Some("human"));
        assert_eq!(r.called, "unknown");
        assert!(!r.mislabel); // never flag a mislabel on an unknown call
    }

    #[test]
    fn empty_catalogue_auto_path_does_not_open_fastq() {
        let cat = crate::speciespanel::Catalogue::from_texts(&[]);
        let missing = std::path::Path::new("/tmp/preflightx-missing.fastq");
        let introns = crate::snpanel::IntronCatalogue::from_texts(&[]);
        let r = detect_species_auto(missing, None, &cat, &introns, Some("human"), 1_000).unwrap();
        assert_eq!(r.called, "unknown");
        assert_eq!(r.reads_scanned, 0);
        assert_eq!(r.n_species_in_panel, 0);
    }

    #[test]
    fn strong_call_is_high_confidence_not_ambiguous() {
        let cat = tiny_catalogue();
        let pad = |kmer: &[u8]| {
            let mut r = kmer.to_vec();
            r.extend(vec![b'N'; 40]);
            r
        };
        let reads: Vec<Vec<u8>> = (0..2000).map(|_| pad(b"AAAAAAAAAA")).collect();
        let r = detect_species(&reads, &cat, None);
        assert_eq!(r.called, "Homo_sapiens");
        assert_eq!(r.confidence, "high");
        assert!(!r.ambiguous);
        assert!(r.abstained_reason.is_none());
    }

    #[test]
    fn near_tie_is_ambiguous_and_low_confidence() {
        let cat = tiny_catalogue();
        let pad = |kmer: &[u8]| {
            let mut r = kmer.to_vec();
            r.extend(vec![b'N'; 40]);
            r
        };
        let mut reads: Vec<Vec<u8>> = (0..100).map(|_| pad(b"AAAAAAAAAA")).collect();
        reads.extend((0..90).map(|_| pad(b"GGGAGGGAGG")));
        let r = detect_species(&reads, &cat, None);
        assert_eq!(r.called, "Homo_sapiens"); // still calls the top
        assert!(r.ambiguous, "expected ambiguous, got {r:?}");
        assert_eq!(r.confidence, "low");
        assert!(r.notes.iter().any(|n| n.contains("within")));
    }

    #[test]
    fn off_panel_unknown_sets_abstained_reason() {
        let cat = tiny_catalogue();
        let reads: Vec<Vec<u8>> = (0..50).map(|_| vec![b'N'; 50]).collect();
        let r = detect_species(&reads, &cat, Some("human"));
        assert_eq!(r.called, "unknown");
        assert_eq!(r.confidence, "low");
        assert!(r.abstained_reason.is_some());
        assert!(!r.ambiguous);
    }

    #[test]
    fn margin_reflects_lead_over_runner_up() {
        let cat = tiny_catalogue();
        let pad = |kmer: &[u8]| {
            let mut r = kmer.to_vec();
            r.extend(vec![b'N'; 40]);
            r
        };
        // 40 reads hit human, 10 hit mouse -> margin (40-10)/50 = 0.6.
        let mut reads: Vec<Vec<u8>> = (0..40).map(|_| pad(b"AAAAAAAAAA")).collect();
        reads.extend((0..10).map(|_| pad(b"GGGAGGGAGG")));
        let r = detect_species(&reads, &cat, None);
        assert_eq!(r.called, "Homo_sapiens");
        assert_eq!(r.runner_up.as_deref(), Some("Mus_musculus"));
        assert!((r.margin - 0.6).abs() < 1e-6, "margin={}", r.margin);
        assert!(!r.mislabel); // no declared label -> never a mislabel
    }

    #[test]
    #[ignore = "network: streams from ENA"]
    fn real_mouse_reads_mus_musculus() {
        let cat = crate::speciespanel::Catalogue::load().unwrap();
        let base = "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/091/SRR26554491/SRR26554491";
        let (r1, r2) = (format!("{base}_1.fastq.gz"), format!("{base}_2.fastq.gz"));
        let introns = crate::snpanel::IntronCatalogue::load_or_empty();
        let rep = detect_species_auto(
            std::path::Path::new(&r1),
            Some(std::path::Path::new(&r2)),
            &cat,
            &introns,
            None,
            20_000,
        )
        .unwrap();
        assert_eq!(rep.called, "Mus_musculus", "{rep:?}");
    }

    fn multiome_species(srr: &str, tail: &str) -> SpeciesReport {
        let cat = crate::speciespanel::Catalogue::load().unwrap();
        let base = format!(
            "https://ftp.sra.ebi.ac.uk/vol1/fastq/{}/{tail}/{srr}/{srr}",
            &srr[..6]
        );
        let (r1, r2) = (format!("{base}_1.fastq.gz"), format!("{base}_2.fastq.gz"));
        let introns = crate::snpanel::IntronCatalogue::load_or_empty();
        detect_species_auto(
            std::path::Path::new(&r1),
            Some(std::path::Path::new(&r2)),
            &cat,
            &introns,
            None,
            150_000,
        )
        .unwrap()
    }

    #[test]
    #[ignore = "network: streams from ENA + needs .snpanel cached"]
    fn real_cat_snrna_calls_felis() {
        let r = multiome_species("SRR34653040", "040");
        assert_eq!(r.called, "Felis_catus", "{r:?}");
        assert_eq!(r.confidence, "high");
    }

    #[test]
    #[ignore = "network: streams from ENA + needs .snpanel cached"]
    fn real_horse_snrna_calls_equus() {
        let r = multiome_species("SRR34653030", "030");
        assert_eq!(r.called, "Equus_caballus", "{r:?}");
        assert_eq!(r.confidence, "high");
    }
}
