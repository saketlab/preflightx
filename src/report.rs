use colored::Colorize;
use serde::Serialize;

use crate::chemistry::ChemReport;
use crate::composition::CompositionReport;
use crate::size::{BulkReport, SizeReport};

#[derive(Serialize)]
pub struct FullReport {
    pub r1: String,
    pub r2: Option<String>,
    pub chemistry: Option<ChemReport>,
    pub bulk_vs_sc: BulkReport,
    pub size: SizeReport,
    pub composition: CompositionReport,
}

impl FullReport {
    pub fn print(&self) {
        println!(
            "{}",
            "============================================================".dimmed()
        );
        println!("{}", "preflightx full report".bold().bright_blue());
        println!("  R1: {}", self.r1);
        if let Some(r2) = &self.r2 {
            println!("  R2: {}", r2);
        }
        println!(
            "{}",
            "============================================================".dimmed()
        );
        println!();
        if let Some(chem) = &self.chemistry {
            chem.print();
            println!();
        }
        self.bulk_vs_sc.print();
        println!();
        self.size.print();
        println!();
        self.composition.print();
        println!();
        println!(
            "{}",
            "============================================================".dimmed()
        );
    }
}
