use colored::Colorize;
use serde::Serialize;

use crate::assaydetect::DetectReport;
use crate::chemistry::ChemReport;
use crate::composition::CompositionReport;
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
}

impl FullReport {
    /// Run chem; when it finds no 10x/Drop-seq match, fall back to the generic
    /// seqspec detector. Shared by the CLI and web builds.
    pub fn chem_with_fallback(
        r1: &std::path::Path,
        r2: &std::path::Path,
        n: usize,
        cache: Option<&std::path::Path>,
    ) -> anyhow::Result<(ChemReport, Option<DetectReport>)> {
        let chem = crate::chemistry::detect(r1, r2, n, cache)?;
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
