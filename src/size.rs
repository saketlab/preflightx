use anyhow::Result;
use colored::Colorize;
use serde::Serialize;
use std::collections::HashMap;
use std::path::Path;

use crate::fastq;
use crate::libstruct;

#[derive(Serialize, Debug)]
pub struct BulkReport {
    pub r1_length_dominant: usize,
    pub r1_length_stddev: f64,
    pub r1_length_distribution: Vec<(usize, usize)>,
    pub verdict: &'static str,
    pub confidence: &'static str,
    pub notes: Vec<String>,
}

#[derive(Serialize, Debug)]
pub struct SizeReport {
    pub r1: PerFileSize,
    pub r2: Option<PerFileSize>,
    pub adequate_for_sc: bool,
    pub notes: Vec<String>,
}

#[derive(Serialize, Debug)]
pub struct PerFileSize {
    pub path: String,
    pub n_reads_sampled: usize,
    pub mean_length: f64,
    pub mean_phred: f64,
    pub low_qual_pct: f64, // % reads with mean PHRED < 20
}

pub fn bulk_or_sc(r1: &Path, n: usize) -> Result<BulkReport> {
    let (seqs, _) = fastq::sample(r1, n)?;
    let mut hist: HashMap<usize, usize> = HashMap::new();
    let mut tot = 0u64;
    let mut tot_sq = 0u64;
    for s in &seqs {
        *hist.entry(s.len()).or_insert(0) += 1;
        tot += s.len() as u64;
        tot_sq += (s.len() as u64) * (s.len() as u64);
    }
    let n_f = seqs.len() as f64;
    let mean = tot as f64 / n_f.max(1.0);
    let var = (tot_sq as f64 / n_f.max(1.0)) - mean * mean;
    let stddev = var.max(0.0).sqrt();
    let mut dist: Vec<(usize, usize)> = hist.into_iter().collect();
    dist.sort_by(|a, b| b.1.cmp(&a.1));
    let dominant = dist.first().map(|(l, _)| *l).unwrap_or(0);

    let mut notes = Vec::new();

    // A fixed-length R1 matching any catalogued seqspec barcode geometry
    // (Drop-seq=20, 10x=26/28, ...) is single-cell.
    let geo_matches: Vec<&'static str> = libstruct::sc_geometries_matching(dominant, 1)
        .iter()
        .map(|g| g.label())
        .collect();

    let (verdict, confidence) = if dominant >= 50 && stddev < 2.0 {
        notes.push("Long R1 with low stddev: bulk RNA-seq, 5'-chemistry broken-paired, or plate-based full-length scRNA-seq (Smart-seq2/3)".into());
        notes.push("Full-length / plate-based protocols (Smart-seq2/3) carry no in-read cell barcode and cannot be distinguished from bulk by read geometry alone".into());
        ("bulk_or_5prime_broken", "medium")
    } else if dominant >= 50 {
        notes.push("Long R1 with high stddev: bulk RNA-seq with adapter trimming, or plate-based full-length scRNA-seq (Smart-seq2/3) — indistinguishable from bulk by geometry (no in-read barcode)".into());
        ("bulk", "high")
    } else if !geo_matches.is_empty() && stddev < 1.0 {
        notes.push(format!(
            "R1 is a fixed-length {}bp barcode read matching a known single-cell geometry: {}",
            dominant,
            geo_matches.join(", "),
        ));
        ("single_cell", "high")
    } else if !geo_matches.is_empty() {
        notes.push(format!(
            "R1 length {}bp matches {} but stddev {:.2} is high; check if reads were trimmed",
            dominant,
            geo_matches.join(", "),
            stddev,
        ));
        ("single_cell_uncertain", "low")
    } else {
        ("uncertain", "low")
    };

    Ok(BulkReport {
        r1_length_dominant: dominant,
        r1_length_stddev: fastq::round2(stddev),
        r1_length_distribution: dist.into_iter().take(5).collect(),
        verdict,
        confidence,
        notes,
    })
}

pub fn size_sanity(r1: &Path, r2: Option<&Path>, n: usize) -> Result<SizeReport> {
    let r1_summary = per_file(r1, n)?;
    let r2_summary = r2.map(|p| per_file(p, n)).transpose()?;
    let mut notes = Vec::new();
    let adequate = r1_summary.n_reads_sampled >= 1000 && r1_summary.mean_phred >= 25.0;
    if r1_summary.n_reads_sampled < n {
        notes.push(format!(
            "Sampled only {} reads from R1 (asked for {}); file may be small.",
            r1_summary.n_reads_sampled, n
        ));
    }
    if r1_summary.mean_phred < 25.0 {
        notes.push(format!(
            "R1 mean PHRED is {:.1} (<25); low base quality reduces downstream accuracy.",
            r1_summary.mean_phred
        ));
    }
    if r1_summary.low_qual_pct > 10.0 {
        notes.push(format!(
            "{:.1}% of R1 reads have mean PHRED < 20.",
            r1_summary.low_qual_pct
        ));
    }

    Ok(SizeReport {
        r1: r1_summary,
        r2: r2_summary,
        adequate_for_sc: adequate,
        notes,
    })
}

fn per_file(path: &Path, n: usize) -> Result<PerFileSize> {
    let mut n_seen = 0;
    let mut total_len = 0u64;
    let mut total_phred = 0.0f64;
    let mut low_q_count = 0u64;
    crate::fastq::for_each_record(path, n, |seq, qual| {
        n_seen += 1;
        total_len += seq.len() as u64;
        if let Some(q) = qual {
            let mp = fastq::mean_phred(q);
            total_phred += mp;
            if mp < 20.0 {
                low_q_count += 1;
            }
        }
    })?;
    let mean_len = total_len as f64 / n_seen.max(1) as f64;
    let mean_phred = total_phred / n_seen.max(1) as f64;
    let low_q_pct = 100.0 * low_q_count as f64 / n_seen.max(1) as f64;
    Ok(PerFileSize {
        path: path.display().to_string(),
        n_reads_sampled: n_seen,
        mean_length: fastq::round2(mean_len),
        mean_phred: fastq::round2(mean_phred),
        low_qual_pct: fastq::round2(low_q_pct),
    })
}

impl BulkReport {
    pub fn print(&self) {
        println!("{}", "Bulk vs single-cell classification".bold().cyan());
        println!("  R1 dominant length:  {} bp", self.r1_length_dominant);
        println!("  R1 length stddev:    {:.2}", self.r1_length_stddev);
        let v = match self.confidence {
            "high" => self.verdict.green(),
            "medium" => self.verdict.yellow(),
            _ => self.verdict.red(),
        };
        println!("  {} {} ({})", "VERDICT:".bold(), v.bold(), self.confidence);
        for n in &self.notes {
            println!("    {} {}", "note:".dimmed(), n);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::{fixed_len_reads, write_fastq};

    #[test]
    fn dropseq_20bp_is_single_cell() {
        let p = write_fastq("size_dropseq", &fixed_len_reads(20, 200));
        let r = bulk_or_sc(&p, 1000).unwrap();
        std::fs::remove_file(&p).ok();
        assert_eq!(r.r1_length_dominant, 20);
        assert_eq!(r.verdict, "single_cell");
        assert_eq!(r.confidence, "high");
    }

    #[test]
    fn long_r1_is_bulk() {
        let p = write_fastq("size_bulk", &fixed_len_reads(90, 200));
        let r = bulk_or_sc(&p, 1000).unwrap();
        std::fs::remove_file(&p).ok();
        assert_eq!(r.verdict, "bulk_or_5prime_broken");
    }
}

impl SizeReport {
    pub fn print(&self) {
        println!("{}", "Size + quality sanity".bold().cyan());
        let p = &self.r1;
        println!(
            "  R1 ({}):  {} reads sampled  mean_len={:.1}  mean_PHRED={:.1}  low_qual={:.1}%",
            p.path, p.n_reads_sampled, p.mean_length, p.mean_phred, p.low_qual_pct
        );
        if let Some(p2) = &self.r2 {
            println!(
                "  R2 ({}):  {} reads sampled  mean_len={:.1}  mean_PHRED={:.1}  low_qual={:.1}%",
                p2.path, p2.n_reads_sampled, p2.mean_length, p2.mean_phred, p2.low_qual_pct
            );
        }
        let label = if self.adequate_for_sc {
            "adequate".green()
        } else {
            "marginal".yellow()
        };
        println!("  {} {}", "QUALITY:".bold(), label.bold());
        for n in &self.notes {
            println!("    {} {}", "note:".dimmed(), n);
        }
    }
}
