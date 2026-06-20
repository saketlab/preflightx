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

    pub fn detect_species(
        reads: &[Vec<u8>],
        cat: &Catalogue,
        declared: Option<&str>,
    ) -> SpeciesReport {
        let declared = declared.map(str::to_string);
        if cat.is_empty() {
            return no_catalogue(declared);
        }
        let k = cat.k();
        let mut counts = vec![0u64; cat.len()];
        let mut scanned = 0u64;
        let mut reads_with_hit = 0u64;
        let mut hit: Vec<u16> = Vec::new();
        let mut up: Vec<u8> = Vec::new();

        for read in reads {
            scanned += 1;
            hit.clear();
            crate::sexpanel::scan_windows(read, k, &mut up, |w| {
                if let Some(idx) = cat.autosomal_hit(w) {
                    if !hit.contains(&idx) {
                        hit.push(idx);
                    }
                }
            });
            if !hit.is_empty() {
                reads_with_hit += 1;
                for &i in &hit {
                    counts[i as usize] += 1;
                }
            }
        }

        let mut ranked: Vec<(u16, u64)> = counts
            .iter()
            .enumerate()
            .map(|(i, &c)| (i as u16, c))
            .filter(|(_, c)| *c > 0)
            .collect();
        ranked.sort_by(|a, b| b.1.cmp(&a.1));

        let scanned_f = scanned.max(1) as f64;
        let (top_idx, top_c) = ranked.first().copied().unwrap_or((0, 0));
        let top_frac = top_c as f64 / scanned_f;
        let (runner_up, runner_c) = match ranked.get(1) {
            Some(&(i, c)) => (Some(cat.species_of(i).to_string()), c),
            None => (None, 0),
        };
        let margin = (top_c.saturating_sub(runner_c)) as f64 / scanned_f;

        let called = if top_c >= MIN_HIT_READS && top_frac >= MIN_HIT_FRACTION {
            cat.species_of(top_idx).to_string()
        } else {
            "unknown".to_string()
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
        if mislabel {
            notes.push(format!(
                "declared '{}' but reads look like '{}' — possible mislabel",
                declared.as_deref().unwrap_or(""),
                called
            ));
        }

        let round4 = |x: f64| crate::fastq::round2(x * 100.0) / 100.0;
        SpeciesReport {
            declared,
            called,
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

    pub fn detect_species_auto(
        r1: &std::path::Path,
        r2: Option<&std::path::Path>,
        cat: &Catalogue,
        declared: Option<&str>,
        n: usize,
    ) -> anyhow::Result<SpeciesReport> {
        if cat.is_empty() {
            return Ok(detect_species(&[], cat, declared));
        }
        let reads = crate::sex::sample_cdna(r1, r2, n)?;
        Ok(detect_species(&reads, cat, declared))
    }

    fn no_catalogue(declared: Option<String>) -> SpeciesReport {
        SpeciesReport {
            declared,
            called: "unknown".into(),
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
pub use detect::{detect_species, detect_species_auto};

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
        // 50 reads each carrying a human autosomal k-mer.
        let mut read = b"AAAAAAAAAA".to_vec();
        read.extend(vec![b'N'; 40]); // pad; Ns never match
        let reads: Vec<Vec<u8>> = (0..50).map(|_| read.clone()).collect();

        // Declared as mouse but reads are human -> mislabel.
        let r = detect_species(&reads, &cat, Some("Mus_musculus"));
        assert_eq!(r.called, "Homo_sapiens");
        assert_eq!(r.called_hit_reads, 50);
        assert!(r.mislabel, "expected mislabel, got {r:?}");

        // Declared correctly (alias) -> no mislabel.
        let r2 = detect_species(&reads, &cat, Some("human"));
        assert!(!r2.mislabel);
    }

    #[test]
    fn off_panel_reads_are_unknown() {
        let cat = tiny_catalogue();
        // Reads matching nothing in the panel.
        let reads: Vec<Vec<u8>> = (0..50).map(|_| vec![b'N'; 50]).collect();
        let r = detect_species(&reads, &cat, Some("human"));
        assert_eq!(r.called, "unknown");
        assert!(!r.mislabel); // never flag a mislabel on an unknown call
    }

    #[test]
    fn empty_catalogue_auto_path_does_not_open_fastq() {
        let cat = crate::speciespanel::Catalogue::from_texts(&[]);
        let missing = std::path::Path::new("/tmp/preflightx-missing.fastq");
        let r = detect_species_auto(missing, None, &cat, Some("human"), 1_000).unwrap();
        assert_eq!(r.called, "unknown");
        assert_eq!(r.reads_scanned, 0);
        assert_eq!(r.n_species_in_panel, 0);
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
}
