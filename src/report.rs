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
    /// Fallback seqspec ranking, run only when `chemistry` finds no 10x/Drop-seq match.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub assay_detection: Option<DetectReport>,
    pub bulk_vs_sc: BulkReport,
    pub size: SizeReport,
    pub composition: CompositionReport,
    /// `None` when no species panels are cached; serializes to no `species` key.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub species: Option<SpeciesReport>,
    /// `None` when sex inference was not requested; serializes to no `sex` key.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sex: Option<SexReport>,
    /// Single-cell vs single-nucleus (intron/exon + mito k-mer content). `None` when the
    /// nucleus panel isn't available (e.g. wasm before injection / no cache).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub nucleus: Option<NucReport>,
    /// 3'/5' transcript-end bias. `None` when the end panel isn't available.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub end: Option<EndReport>,
}

/// Trimmed report for the `detect-all` subcommand: chemistry + species + sex only.
#[derive(Serialize)]
pub struct DetectAllReport {
    pub r1: String,
    pub r2: Option<String>,
    /// `None` for single-end input (chemistry needs both mates).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chemistry: Option<ChemReport>,
    /// `None` when no species panels are cached.
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

    /// One line per check: the call and its confidence, no evidence ledger or notes.
    /// The default `detect-all` view; `--detailed` switches to [`print`].
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
            let mislabel = if s.mislabel { "  MISLABEL".red().to_string() } else { String::new() };
            println!("  species:   {} ({}){}", s.called, conf(&s.confidence), mislabel);
        }
        println!("  sex:       {} ({})", self.sex.verdict, conf(&self.sex.confidence));
    }

    /// The concise JSON counterpart of [`print_summary`]: the call + confidence per check.
    /// `--detailed` serializes the full report instead.
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

/// Nucleus (sc-vs-sn) and end (3'/5') verdicts, each `None` when its panel isn't loaded.
/// Shared by the native `check` and the wasm `run_check` so both surface the panels
/// identically. Samples cDNA only when at least one panel is present, so the common
/// no-panel path (no cached panels / wasm before injection) pays no extra FASTQ pass.
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
    /// Run chem; when it finds no 10x/Drop-seq match, fall back to the generic
    /// seqspec detector. Shared by the CLI and web builds.
    pub fn chem_with_fallback(
        r1: &std::path::Path,
        r2: &std::path::Path,
        n: usize,
        cache: Option<&std::path::Path>,
        species: Option<&str>,
    ) -> anyhow::Result<(ChemReport, Option<DetectReport>)> {
        let chem = crate::chemistry::detect(r1, r2, n, cache, species)?;
        let assay = if chem.is_single_cell {
            None
        } else {
            Some(crate::assaydetect::detect(r1, r2, n)?)
        };
        Ok((chem, assay))
    }

    /// Self-contained HTML report. Routes through `serde_json::Value` so
    /// [`crate::reporthtml`] is the only renderer for both CLI and web.
    pub fn to_html(&self) -> String {
        crate::reporthtml::render(&serde_json::to_value(self).unwrap_or(serde_json::Value::Null))
    }

    pub fn print(&self) {
        println!("{}", "preflightx full report".bold().bright_blue());
        println!("  R1: {}", self.r1);
        if let Some(r2) = &self.r2 {
            println!("  R2: {}", r2);
        }
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
        }
    }

    #[test]
    fn omits_sex_when_not_requested() {
        // A `None` sex serializes to a report with no `sex` key.
        let json = serde_json::to_string(&minimal(None)).unwrap();
        assert!(
            !json.contains("\"sex\""),
            "expected no sex key when sex is None, got: {json}"
        );
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
