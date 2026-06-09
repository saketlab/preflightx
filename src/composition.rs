use anyhow::Result;
use colored::Colorize;
use serde::Serialize;
use std::path::Path;

use crate::fastq;
use crate::libstruct;

#[derive(Serialize, Debug)]
pub struct CompositionReport {
    pub gc_pct: f64,
    pub n_pct: f64,
    /// Overall base composition as percentages, ordered [A, C, G, T, N].
    pub base_pct: [f64; 5],
    pub low_complexity_pct: f64, // % reads with single base dominating >= 80%
    pub homopolymer_pct: f64,    // % reads containing a >= 20-base homopolymer
    pub guess_organism_family: &'static str,
    /// True when the reads are a fixed-length cell-barcode/UMI read (synthetic
    /// sequence) rather than cDNA, in which case the organism hint is suppressed.
    pub is_barcode_read: bool,
    pub notes: Vec<String>,
}

fn pct(num: u64, den: u64) -> f64 {
    if den == 0 {
        0.0
    } else {
        100.0 * num as f64 / den as f64
    }
}

pub fn composition(path: &Path, n: usize) -> Result<CompositionReport> {
    let mut total_bases = 0u64;
    let mut gc = 0u64;
    let mut n_bases = 0u64;
    let mut low_complex = 0u64;
    let mut homopoly = 0u64;
    let mut total_reads = 0u64;
    let mut len_sq_sum = 0u64;
    let mut base_totals = [0u64; 5]; // A, C, G, T, N

    fastq::for_each_record(path, n, |seq, _| {
        total_reads += 1;
        total_bases += seq.len() as u64;
        len_sq_sum += (seq.len() as u64) * (seq.len() as u64);
        let mut counts = [0u64; 5];
        let mut run = 1u32;
        let mut max_run = 0u32;
        let mut prev = 0u8;
        for (i, &b) in seq.iter().enumerate() {
            let b_upper = if b.is_ascii_lowercase() { b - 32 } else { b };
            let idx = match b_upper {
                b'A' => 0,
                b'C' => {
                    gc += 1;
                    1
                }
                b'G' => {
                    gc += 1;
                    2
                }
                b'T' => 3,
                _ => {
                    n_bases += 1;
                    4
                }
            };
            counts[idx] += 1;
            if i > 0 && b_upper == prev && b_upper != b'N' {
                run += 1;
                if run > max_run {
                    max_run = run;
                }
            } else {
                run = 1;
            }
            prev = b_upper;
        }
        for k in 0..5 {
            base_totals[k] += counts[k];
        }
        let len = seq.len() as u64;
        if len > 0 {
            let dom: u64 = counts.iter().take(4).copied().max().unwrap_or(0);
            if dom as f64 / len as f64 >= 0.80 {
                low_complex += 1;
            }
        }
        if max_run >= 20 {
            homopoly += 1;
        }
    })?;

    let gc_pct = pct(gc, total_bases);
    let n_pct = pct(n_bases, total_bases);
    let base_pct = [
        pct(base_totals[0], total_bases),
        pct(base_totals[1], total_bases),
        pct(base_totals[2], total_bases),
        pct(base_totals[3], total_bases),
        pct(base_totals[4], total_bases),
    ];
    let lc_pct = pct(low_complex, total_reads);
    let hp_pct = pct(homopoly, total_reads);

    // A fixed-length read matching a known barcode geometry is synthetic CB/UMI, not
    // cDNA -- its GC% is not organism-informative, so suppress the hint below.
    let mean_len = if total_reads > 0 {
        total_bases as f64 / total_reads as f64
    } else {
        0.0
    };
    let len_var = if total_reads > 0 {
        len_sq_sum as f64 / total_reads as f64 - mean_len * mean_len
    } else {
        0.0
    };
    let len_stddev = len_var.max(0.0).sqrt();
    let is_barcode_read = len_stddev < 1.0
        && !libstruct::sc_geometries_matching(mean_len.round() as usize, 1).is_empty();

    // Very rough organism heuristic. This is a SIGNAL, not a CALL.
    let guess = if is_barcode_read {
        "not applicable (R1 barcode/UMI read)"
    } else {
        match gc_pct {
            x if (38.0..=46.0).contains(&x) => "mammalian (human/mouse/rat range)",
            x if (35.0..=38.0).contains(&x) => "lower-GC mammalian / possible mouse",
            x if (46.0..=55.0).contains(&x) => {
                "higher-GC; possible bacterial / viral contamination"
            }
            x if x > 55.0 => "high-GC; bacterial or GC-rich contamination",
            x if x < 35.0 => "low-GC; AT-rich (mitochondrial, parasites, or low-quality reads)",
            _ => "uncertain",
        }
    };

    let mut notes = Vec::new();
    if is_barcode_read {
        notes.push(
            "Organism hint suppressed: these reads are a fixed-length cell-barcode/UMI read \
             (synthetic sequence), not cDNA. Run composition on the cDNA read (R2) for an organism signal."
                .into(),
        );
    }
    if n_pct > 1.0 {
        notes.push(format!(
            "{:.1}% Ns in sequence -- adapter trimming or base-calling issues",
            n_pct
        ));
    }
    if lc_pct > 5.0 {
        notes.push(format!(
            "{:.1}% reads are low-complexity (single base >= 80%) -- adapter-only reads or polyA tail capture",
            lc_pct
        ));
    }
    if hp_pct > 10.0 {
        notes.push(format!(
            "{:.1}% reads contain >= 20-base homopolymer -- typical polyA-tail-rich 10x 3' chemistry",
            hp_pct
        ));
    }

    Ok(CompositionReport {
        gc_pct: fastq::round2(gc_pct),
        n_pct: fastq::round2(n_pct),
        base_pct: base_pct.map(fastq::round2),
        low_complexity_pct: fastq::round2(lc_pct),
        homopolymer_pct: fastq::round2(hp_pct),
        guess_organism_family: guess,
        is_barcode_read,
        notes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::write_fastq;

    #[test]
    fn barcode_read_suppresses_organism_hint() {
        // Fixed 20bp reads = Drop-seq barcode geometry -> hint suppressed.
        let reads: Vec<Vec<u8>> = (0..200)
            .map(|i| (0..20).map(|j| b"ACGT"[(i + j) % 4]).collect())
            .collect();
        let p = write_fastq("comp_barcode", &reads);
        let r = composition(&p, 1000).unwrap();
        std::fs::remove_file(&p).ok();
        assert!(r.is_barcode_read);
        assert_eq!(
            r.guess_organism_family,
            "not applicable (R1 barcode/UMI read)"
        );
    }

    #[test]
    fn cdna_read_gets_organism_hint() {
        // 100bp reads at 40% GC, no barcode geometry -> real organism hint.
        let read: Vec<u8> = b"GCATT".iter().copied().cycle().take(100).collect();
        let p = write_fastq("comp_cdna", &vec![read; 200]);
        let r = composition(&p, 1000).unwrap();
        std::fs::remove_file(&p).ok();
        assert!(!r.is_barcode_read);
        assert_eq!(r.guess_organism_family, "mammalian (human/mouse/rat range)");
    }
}

impl CompositionReport {
    pub fn print(&self) {
        println!("{}", "Sequence composition".bold().cyan());
        println!("  GC content:          {:.2}%", self.gc_pct);
        println!("  N bases:             {:.2}%", self.n_pct);
        println!(
            "  low-complexity:      {:.2}% of reads",
            self.low_complexity_pct
        );
        println!(
            "  long homopolymer:    {:.2}% of reads",
            self.homopolymer_pct
        );
        let hint = if self.is_barcode_read {
            self.guess_organism_family.dimmed()
        } else {
            self.guess_organism_family.yellow()
        };
        println!("  {} {}", "ORGANISM HINT:".bold(), hint);
        for n in &self.notes {
            println!("    {} {}", "note:".dimmed(), n);
        }
    }
}
