use colored::Colorize;
use serde::Serialize;

use crate::assaydetect::DetectReport;
use crate::chemistry::ChemReport;
use crate::composition::CompositionReport;
use crate::endpanel::EndReport;
use crate::nucpanel::NucReport;
use crate::sex::SexReport;
use crate::size::{BulkReport, SizeReport};
use crate::species::SpeciesReport;

#[derive(Serialize)]
pub struct FullReport {
    pub r1: String,
    pub r2: Option<String>,
    pub chemistry: Option<ChemReport>,
    /// Fallback seqspec ranking, run only when chemistry finds no 10x/Drop-seq match.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub assay_detection: Option<DetectReport>,
    pub bulk_vs_sc: BulkReport,
    pub size: SizeReport,
    pub composition: CompositionReport,
    /// None when no species panels are cached; serializes to no species key.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub species: Option<SpeciesReport>,
    /// None when sex inference was not requested; serializes to no sex key.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sex: Option<SexReport>,
    /// Single-cell vs single-nucleus (intron/exon + mito k-mer content). None when the
    /// nucleus panel isn't available (e.g. wasm before injection / no cache).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub nucleus: Option<NucReport>,
    /// 3'/5' transcript-end bias. None when the end panel isn't available.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub end: Option<EndReport>,
    pub classification: Classification,
}

/// bulk_vs_sc + nucleus + chemistry, rolled into one call.
#[derive(Serialize, Default, Debug, Clone, PartialEq)]
pub struct Classification {
    /// "bulk" | "single-cell" | "single-nucleus" | "single-undetermined".
    /// Single ATAC is always "single-undetermined"; see [Classification::derive].
    pub library: String,
    /// "RNA" | "ATAC" | "unknown".
    pub modality: String,
    /// e.g. "10x 3' v3". None for bulk.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chemistry: Option<String>,
    /// e.g. "single-nucleus RNA · 10x 3' v3".
    pub label: String,
    pub confidence: String,
}

fn weakest<'a>(a: &'a str, b: &'a str) -> &'a str {
    let rank = |c: &str| match c {
        "high" => 2,
        "medium" => 1,
        _ => 0,
    };
    if rank(a) <= rank(b) {
        a
    } else {
        b
    }
}

impl Classification {
    /// Single if either signal says so; only chemistry catches 10x 5' (bulk-vs-sc calls it
    /// "bulk_or_5prime_broken").
    pub fn derive(
        chemistry: &Option<ChemReport>,
        bulk: &BulkReport,
        nucleus: &Option<NucReport>,
    ) -> Self {
        let chem = chemistry.as_ref();
        let chem_sc = chem.is_some_and(|c| c.is_single_cell);
        let bulk_sc = matches!(bulk.verdict, "single_cell" | "single_cell_uncertain");
        let is_single = chem_sc || bulk_sc;

        // nucleus panel scores transcript k-mers; ATAC reads are gDNA. Most ATAC is run on
        // nuclei, so defaulting to "cell" would usually be wrong.
        if chem.is_some_and(|c| c.modality == "atac") {
            let conf = weakest(
                chem.map_or("low", |c| c.confidence.as_str()),
                bulk.confidence,
            );
            let (library, label) = if bulk_sc {
                ("single-undetermined", "single ATAC")
            } else {
                ("bulk", "bulk ATAC")
            };
            return Classification {
                library: library.into(),
                modality: "ATAC".into(),
                chemistry: Some("ATAC (chromatin accessibility)".into()),
                label: label.into(),
                confidence: conf.into(),
            };
        }
        let modality = if chem_sc {
            "RNA"
        } else {
            match chem.map(|c| c.modality.as_str()) {
                Some("rna") => "RNA",
                _ => "unknown",
            }
        };
        if !is_single {
            let label = if modality == "unknown" {
                "bulk".to_string()
            } else {
                format!("bulk {modality}")
            };
            return Classification {
                library: "bulk".into(),
                modality: modality.into(),
                chemistry: None,
                label,
                confidence: bulk.confidence.to_string(),
            };
        }
        let (library, nuc_conf) = match nucleus {
            Some(n) if n.verdict == "single-nucleus" => ("single-nucleus", n.confidence.as_str()),
            Some(n) if n.verdict == "single-cell" => ("single-cell", n.confidence.as_str()),
            _ => ("single-undetermined", "low"),
        };
        // assay_display is "Not single-cell" on a negative call; assay_id is
        // "unknown" on an abstention.
        let chemistry = chem
            .filter(|c| c.is_single_cell && c.assay_id != "unknown")
            .map(|c| c.assay_display.clone());
        let sc_conf = if chem_sc {
            weakest(
                chem.map_or("low", |c| c.confidence.as_str()),
                bulk.confidence,
            )
        } else {
            bulk.confidence
        };
        let confidence = weakest(nuc_conf, sc_conf).to_string();
        let label = match &chemistry {
            Some(a) => format!("{library} {modality} · {a}"),
            None => format!("{library} {modality}"),
        };
        Classification {
            library: library.to_string(),
            modality: modality.into(),
            chemistry,
            label,
            confidence,
        }
    }
}

/// Trimmed report for the detect-all subcommand: chemistry + species + sex only.
#[derive(Serialize)]
pub struct DetectAllReport {
    pub r1: String,
    pub r2: Option<String>,
    /// None for single-end input (chemistry needs both mates).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chemistry: Option<ChemReport>,
    /// None when no species panels are cached.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub species: Option<SpeciesReport>,
    pub sex: SexReport,
}

impl DetectAllReport {
    pub fn print(&self) {
        println!("{}", "preflightx detect-all".bold().bright_blue());
        println!("  R1: {}", self.r1);
        if let Some(r2) = &self.r2 {
            println!("  R2: {}", r2);
        }
        println!();
        if let Some(chem) = &self.chemistry {
            chem.print();
            println!();
        }
        if let Some(species) = &self.species {
            species.print();
            println!();
        }
        self.sex.print();
        println!();
    }

    /// Default detect-all view; --detailed switches to [print].
    pub fn print_summary(&self) {
        let conf = |c: &str| match c {
            "high" => c.green(),
            "medium" => c.yellow(),
            _ => c.red(),
        };
        println!("{}", "preflightx detect-all".bold().bright_blue());
        if let Some(c) = &self.chemistry {
            let label = if c.modality == "atac" {
                "ATAC (chromatin accessibility)"
            } else {
                &c.assay_display
            };
            println!("  chemistry: {} ({})", label, conf(&c.confidence));
        }
        if let Some(s) = &self.species {
            let mislabel = if s.mislabel {
                "  MISLABEL".red().to_string()
            } else {
                String::new()
            };
            println!(
                "  species:   {} ({}){}",
                s.called,
                conf(&s.confidence),
                mislabel
            );
        }
        println!(
            "  sex:       {} ({})",
            self.sex.verdict,
            conf(&self.sex.confidence)
        );
    }

    /// JSON counterpart of [print_summary]; --detailed serializes the full report.
    pub fn summary_json(&self) -> serde_json::Value {
        use serde_json::json;
        json!({
            "r1": self.r1,
            "r2": self.r2,
            "chemistry": self.chemistry.as_ref().map(|c| json!({
                "assay_id": c.assay_id,
                "assay_display": c.assay_display,
                "confidence": c.confidence,
                "modality": c.modality,
                "is_single_cell": c.is_single_cell,
            })),
            "species": self.species.as_ref().map(|s| json!({
                "called": s.called,
                "confidence": s.confidence,
                "mislabel": s.mislabel,
            })),
            "sex": json!({ "verdict": self.sex.verdict, "confidence": self.sex.confidence }),
        })
    }
}

/// Nucleus and end verdicts, each None when its panel isn't loaded. No FASTQ pass
/// without a panel.
pub fn panel_verdicts(
    r1: &std::path::Path,
    r2: Option<&std::path::Path>,
    n: usize,
    species: Option<&str>,
) -> (Option<NucReport>, Option<EndReport>) {
    if !(crate::nucpanel::have_panels() || crate::endpanel::have_panels()) {
        return (None, None);
    }
    let reads = crate::sex::sample_cdna(r1, r2, n).unwrap_or_default();
    let nucleus = crate::nucpanel::have_panels().then(|| crate::nucpanel::classify(&reads));
    let end = crate::endpanel::have_panels().then(|| crate::endpanel::classify(&reads, species));
    (nucleus, end)
}

impl FullReport {
    /// Run chem; fall back to the generic seqspec detector unless the winning chemistry
    /// has a barcode list. Drop-seq alone has none, so its geometry-only ~20 bp R1 call
    /// would suppress the detector on combinatorial libraries of that length.
    pub fn chem_with_fallback(
        r1: &std::path::Path,
        r2: &std::path::Path,
        n: usize,
        cache: Option<&std::path::Path>,
        species: Option<&str>,
    ) -> anyhow::Result<(ChemReport, Option<DetectReport>)> {
        let chem = crate::chemistry::detect(r1, r2, n, cache, species)?;
        let has_whitelist = chem.is_single_cell && chem.whitelist.is_some();
        let assay = if has_whitelist {
            None
        } else {
            Some(crate::assaydetect::detect(r1, r2, n)?)
        };
        Ok((chem, assay))
    }

    /// Via serde_json::Value, so [crate::reporthtml] stays the only renderer.
    pub fn to_html(&self) -> String {
        crate::reporthtml::render(&serde_json::to_value(self).unwrap_or(serde_json::Value::Null))
    }

    pub fn print(&self) {
        println!("{}", "preflightx full report".bold().bright_blue());
        println!("  R1: {}", self.r1);
        if let Some(r2) = &self.r2 {
            println!("  R2: {}", r2);
        }
        let conf = match self.classification.confidence.as_str() {
            "high" => self.classification.label.green(),
            "medium" => self.classification.label.yellow(),
            _ => self.classification.label.red(),
        };
        println!("  {} {}", "classification:".bold(), conf);
        println!();
        if let Some(chem) = &self.chemistry {
            chem.print();
            println!();
        }
        if let Some(assay) = &self.assay_detection {
            assay.print();
            println!();
        }
        self.bulk_vs_sc.print();
        println!();
        self.size.print();
        println!();
        self.composition.print();
        println!();
        if let Some(species) = &self.species {
            species.print();
            println!();
        }
        if let Some(sex) = &self.sex {
            sex.print();
            println!();
        }
        if let Some(nucleus) = &self.nucleus {
            nucleus.print();
            println!();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::composition::CompositionReport;
    use crate::size::{BulkReport, PerFileSize, SizeReport};

    fn minimal(sex: Option<SexReport>) -> FullReport {
        FullReport {
            r1: "r1".into(),
            r2: None,
            chemistry: None,
            assay_detection: None,
            bulk_vs_sc: BulkReport {
                r1_length_dominant: 28,
                r1_length_stddev: 0.0,
                r1_length_distribution: vec![(28, 100)],
                verdict: "single-cell",
                confidence: "high",
                notes: vec![],
            },
            size: SizeReport {
                r1: PerFileSize {
                    path: "r1".into(),
                    n_reads_sampled: 100,
                    mean_length: 28.0,
                    mean_phred: 35.0,
                    low_qual_pct: 0.0,
                },
                r2: None,
                adequate_for_sc: true,
                notes: vec![],
            },
            species: None,
            composition: CompositionReport {
                gc_pct: 50.0,
                n_pct: 0.0,
                base_pct: [25.0; 5],
                low_complexity_pct: 0.0,
                homopolymer_pct: 0.0,
                guess_organism_family: "unknown",
                is_barcode_read: false,
                notes: vec![],
            },
            sex,
            nucleus: None,
            end: None,
            classification: Classification::default(),
        }
    }

    #[test]
    fn omits_sex_when_not_requested() {
        // A None sex serializes to a report with no sex key.
        let json = serde_json::to_string(&minimal(None)).unwrap();
        assert!(
            !json.contains("\"sex\""),
            "expected no sex key when sex is None, got: {json}"
        );
    }

    fn bulk_rep(verdict: &'static str) -> BulkReport {
        BulkReport {
            r1_length_dominant: 28,
            r1_length_stddev: 0.0,
            r1_length_distribution: vec![],
            verdict,
            confidence: "high",
            notes: vec![],
        }
    }

    fn nuc_rep(verdict: &str) -> crate::nucpanel::NucReport {
        crate::nucpanel::NucReport {
            verdict: verdict.into(),
            confidence: "high".into(),
            n_intron: 0,
            n_exon: 0,
            n_mito: 0,
            reads_scanned: 0,
            intron_frac: 0.0,
            mito_pct: 0.0,
            notes: vec![],
        }
    }

    #[test]
    fn classify_bulk() {
        let c = Classification::derive(&None, &bulk_rep("bulk"), &None);
        assert_eq!(c.library, "bulk");
        assert_eq!(c.chemistry, None);
        assert_eq!(c.label, "bulk");
    }

    #[test]
    fn classify_single_cell_vs_nucleus() {
        let cell = Classification::derive(
            &None,
            &bulk_rep("single_cell"),
            &Some(nuc_rep("single-cell")),
        );
        assert_eq!(cell.library, "single-cell");
        let nuc = Classification::derive(
            &None,
            &bulk_rep("single_cell"),
            &Some(nuc_rep("single-nucleus")),
        );
        assert_eq!(nuc.library, "single-nucleus");
    }

    #[test]
    fn classify_single_undetermined_when_nucleus_missing() {
        let c = Classification::derive(&None, &bulk_rep("single_cell"), &None);
        assert_eq!(c.library, "single-undetermined");
        let amb =
            Classification::derive(&None, &bulk_rep("single_cell"), &Some(nuc_rep("ambiguous")));
        assert_eq!(amb.library, "single-undetermined");
    }

    fn chem_rep(is_single_cell: bool, modality: &str, confidence: &str) -> ChemReport {
        ChemReport {
            r1_length_dominant: 28,
            r1_length_distribution: vec![],
            reads_sampled: 0,
            assay_id: "none".into(),
            // matches what chemistry::report writes on a negative call
            assay_display: if is_single_cell {
                "10x 3' v3".into()
            } else {
                "Not single-cell".into()
            },
            cb_read: "R1".into(),
            cb_offset: 0,
            cb_len: 16,
            umi_offset: 16,
            umi_len: 12,
            n_cb_candidates: 0,
            is_single_cell,
            early_stop: false,
            whitelist: None,
            whitelist_scores: vec![],
            anchor_scores: vec![],
            tool_hints: vec![],
            confidence: confidence.into(),
            modality: modality.into(),
            atac_promoter_frac: 0.0,
            end: "unknown".into(),
            end_confidence: "low".into(),
            n5: 0,
            n3: 0,
            end_ratio: 0.0,
            interpretation: vec![],
            notes: vec![],
        }
    }

    #[test]
    fn negative_chemistry_never_names_an_assay() {
        let c = Classification::derive(
            &Some(chem_rep(false, "rna", "high")),
            &bulk_rep("single_cell"),
            &Some(nuc_rep("single-cell")),
        );
        assert_eq!(c.chemistry, None);
        assert_eq!(c.label, "single-cell RNA");
        assert!(!c.label.contains("Not single-cell"));
    }

    #[test]
    fn uncertain_single_cell_is_not_bulk() {
        let c = Classification::derive(&None, &bulk_rep("single_cell_uncertain"), &None);
        assert_eq!(c.library, "single-undetermined");
        assert_ne!(c.library, "bulk");
    }

    #[test]
    fn confidence_is_the_weakest_signal() {
        let mut bulk = bulk_rep("single_cell");
        bulk.confidence = "low";
        let c = Classification::derive(
            &Some(chem_rep(true, "rna", "low")),
            &bulk,
            &Some(nuc_rep("single-nucleus")),
        );
        assert_eq!(c.confidence, "low");
    }

    #[test]
    fn single_atac_never_claims_cell_or_nucleus() {
        for verdict in ["single_cell", "single_cell_uncertain"] {
            let c = Classification::derive(
                &Some(chem_rep(false, "atac", "high")),
                &bulk_rep(verdict),
                &None,
            );
            assert_eq!(c.library, "single-undetermined");
            assert_eq!(c.label, "single ATAC");
        }
    }

    #[test]
    fn nucleus_panel_cannot_override_atac() {
        // a nucleus verdict on gDNA reads is noise
        let c = Classification::derive(
            &Some(chem_rep(false, "atac", "high")),
            &bulk_rep("single_cell"),
            &Some(nuc_rep("single-cell")),
        );
        assert_eq!(c.library, "single-undetermined");
    }

    #[test]
    fn includes_sex_when_requested() {
        let sx = SexReport {
            species: "human".into(),
            panel_status: "validated".into(),
            verdict: "female".into(),
            confidence: "high".into(),
            reads_scanned: 200_000,
            y_hit_reads: 1,
            xist_hit_reads: 40,
            y_xist_ratio: 0.02,
            y_per_gene: vec![],
            mixed_sex_suspected: false,
            notes: vec![],
        };
        let json = serde_json::to_string(&minimal(Some(sx))).unwrap();
        assert!(json.contains("\"sex\""), "expected sex key, got: {json}");
    }
}
