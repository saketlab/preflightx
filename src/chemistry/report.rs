//! Chemistry-detection report assembly: build the `ChemReport` from an `Analysis`,
//! render the human-readable interpretation, and `print()`. Split out to keep
//! `chemistry.rs` focused on sampling + scoring.

use super::*;
use crate::assay::{AssaySpec, TSO_PRESENT};
use colored::{ColoredString, Colorize};
use std::collections::HashMap;

/// Turn a scored `Analysis` into the serializable report: winner geometry, tool
/// hints, whitelist/anchor tables, the evidence ledger, and human-facing notes.
/// `n_requested` is the read budget asked for; `early_stop` is true when SPRT
/// decided before reaching it (so a short sample is explained, not mistaken for a
/// truncated input).
pub(crate) fn build_report(
    a: Analysis,
    swapped: bool,
    n_requested: usize,
    early_stop: bool,
) -> ChemReport {
    let Analysis {
        dominant,
        len_dist,
        candidates,
        cb_sets,
        anchor_pcts,
        wl_scores,
        scored,
        confidence,
        n_reads,
        tso_pct,
        max_wl_pct,
        is_single_cell,
    } = a;

    // A non-single-cell library has no real chemistry to report, and its scored
    // winner is a force-fit; never present it as a confident call.
    let confidence = if is_single_cell { confidence } else { "low" };
    let has_wl = !wl_scores.is_empty();
    let cb_sizes: HashMap<(usize, usize), usize> =
        cb_sets.iter().map(|(k, v)| (*k, v.len())).collect();
    let winner = scored.first().map(|(s, _)| *s);
    let runner = scored.get(1).map(|(s, _)| *s);

    // 3'v2 and 5'v1 share both geometry and onlist; once the TSO no longer penalizes
    // absence (one-sided LLR), the top two can tie and differ only in 3'/5' polarity.
    // That polarity isn't recoverable from the reads, so report both rather than guess.
    let gap = scored.first().map(|(_, s)| *s).unwrap_or(0.0)
        - scored.get(1).map(|(_, s)| *s).unwrap_or(0.0);
    let polarity_unresolved = match (winner, runner) {
        (Some(w), Some(r)) => {
            is_single_cell
                && gap.abs() < 1.0
                && same_geometry(w, r)
                && w.id.contains("3prime") != r.id.contains("3prime")
        }
        _ => false,
    };
    // Can't claim "high" for a specific polarity when the polarity is unresolved.
    let confidence = if polarity_unresolved && confidence == "high" {
        "medium"
    } else {
        confidence
    };

    let mut notes: Vec<String> = Vec::new();
    if !is_single_cell {
        notes.push(format!(
            "No single-cell signal detected: R1 is {dominant}bp (single-cell barcode reads are typically ~24-30bp), \
             TSO appears in {tso_pct:.0}% of R2 (background ~1%), and whitelist overlap is {max_wl_pct:.0}%. \
             This looks like bulk or genomic sequencing, not single-cell."
        ));
        if let Some(w) = winner {
            notes.push(format!(
                "If it were single-cell, the closest geometry would be {} ({}bp: CB{}+UMI{}) — but the evidence above does not support a single-cell call.",
                w.display_name, w.r1_expected_len(), w.r1_cb_len, w.r1_umi_len,
            ));
        }
    }
    if polarity_unresolved {
        if let (Some(w), Some(r)) = (winner, runner) {
            notes.insert(
                0,
                format!(
                    "3' vs 5' is unresolved: {} and {} share the same {dominant}bp barcode geometry and the same onlist, and the R2 TSO is at background, so the reads alone can't tell them apart. Check the library-prep records to confirm which.",
                    w.display_name, r.display_name,
                ),
            );
        }
    }
    if swapped {
        notes.push("Input files appeared swapped (R1 was long, R2 was short); auto-oriented to treat the shorter file as R1 (barcodes).".into());
    }
    if early_stop {
        notes.push(format!(
            "Stopped early after {n_reads} reads: the TSO signal was already decisive (sequential test), so more reads were not needed."
        ));
    } else if n_reads < n_requested {
        notes.push(format!(
            "Input ended after {n_reads} reads (fewer than the {n_requested} requested) — a short file or a sliced fetch."
        ));
    }
    // The assay-specific notes below only make sense for a real single-cell call;
    // skip them entirely when the reads are not single-cell.
    if is_single_cell {
        // Chemistries like Drop-seq use random (non-whitelisted) barcodes, so there is no
        // list to confirm them against and `fetch-whitelists` cannot help — say so instead
        // of suggesting a download that won't change the verdict.
        if winner.map(|w| w.whitelist.is_none()).unwrap_or(false) {
            notes.push(format!(
                "{} uses random (non-whitelisted) cell barcodes; no whitelist exists to confirm them, so confidence is capped at MEDIUM (geometry-only).",
                winner.map(|w| w.display_name).unwrap_or("This chemistry"),
            ));
        } else if !has_wl {
            notes.push(
                "No whitelist cache; run `preflightx fetch-whitelists` for full disambiguation."
                    .into(),
            );
        }
        if confidence != "high" && scored.len() > 1 {
            let alts: Vec<&str> = scored.iter().skip(1).take(3).map(|(s, _)| s.id).collect();
            notes.push(format!("Other candidates: {}.", alts.join(", ")));
        }
        if let Some(w) = winner {
            let expected = w.r1_expected_len();
            if expected != dominant && dominant > 0 {
                notes.push(format!(
                    "R1 reads are {}bp but {} expects {}bp; library may have been sequenced at non-standard cycles.",
                    dominant, w.id, expected,
                ));
            }
        }
    }

    let mut assay_id = winner.map(|w| w.id).unwrap_or("unknown").to_string();
    let mut assay_display = winner
        .map(|w| w.display_name)
        .unwrap_or("unknown")
        .to_string();
    if !is_single_cell {
        assay_id = "none".to_string();
        assay_display = "Not single-cell".to_string();
    } else if polarity_unresolved {
        if let (Some(w), Some(r)) = (winner, runner) {
            assay_display = format!("{} or {}", w.display_name, r.display_name);
        }
    }
    let cb_offset = winner.map(|w| w.r1_cb_offset).unwrap_or(0);
    let cb_len = winner.map(|w| w.r1_cb_len).unwrap_or(0);
    let umi_offset = winner.map(|w| w.umi_offset()).unwrap_or(0);
    let umi_len = winner.map(|w| w.r1_umi_len).unwrap_or(0);
    let whitelist_file =
        winner.and_then(|w| w.whitelist.as_ref().map(|wl| wl.filename.to_string()));
    let tool_hints = winner
        .map(|w| {
            let mut hints = vec![ToolHint {
                tool: "starsolo".into(),
                hint: w.starsolo_hint(),
            }];
            hints.extend(w.tool_hints.iter().map(|(t, h)| ToolHint {
                tool: t.to_string(),
                hint: h.to_string(),
            }));
            hints
        })
        .unwrap_or_default();

    let n_cb = winner
        .and_then(|w| cb_sets.get(&(w.r1_cb_offset, w.r1_cb_len)))
        .map(|s| s.len())
        .unwrap_or(0);

    let mut wl_scores_out: Vec<WhitelistScore> = wl_scores
        .iter()
        .map(|(f, (hits, pct))| WhitelistScore {
            filename: f.to_string(),
            hits: *hits,
            pct: fastq::round2(*pct),
        })
        .collect();
    wl_scores_out.sort_by(|a, b| {
        b.pct
            .partial_cmp(&a.pct)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let anchor_scores_out: Vec<AnchorScore> = anchor_pcts
        .iter()
        .map(|(seq, pct)| AnchorScore {
            motif: String::from_utf8_lossy(&seq[..seq.len().min(13)]).into_owned(),
            hit_pct: fastq::round2(*pct),
        })
        .collect();

    let mut interpretation = if is_single_cell {
        build_interpretation(
            winner,
            &scored,
            dominant,
            &candidates,
            &anchor_pcts,
            &wl_scores,
            &cb_sizes,
            n_reads,
            has_wl,
            confidence,
        )
    } else {
        build_not_single_cell_interpretation(dominant, tso_pct, max_wl_pct, n_reads, winner)
    };
    if polarity_unresolved {
        if let (Some(w), Some(r)) = (winner, runner) {
            interpretation.push(Evidence {
                signal: "3' vs 5' polarity".into(),
                observation: format!(
                    "{} and {} tie: same {dominant}bp geometry and onlist, TSO at background",
                    w.id, r.id,
                ),
                llr_contribution: None,
                conclusion:
                    "The reads carry no signal that separates 3' from 5' here, so both are reported. The R2 TSO would normally break this tie when present.".into(),
            });
        }
    }

    ChemReport {
        r1_length_dominant: dominant,
        r1_length_distribution: len_dist.into_iter().take(5).collect(),
        reads_sampled: n_reads,
        assay_id,
        assay_display,
        cb_read: if swapped { "R2" } else { "R1" }.into(),
        cb_offset,
        cb_len,
        umi_offset,
        umi_len,
        n_cb_candidates: n_cb,
        is_single_cell,
        early_stop,
        whitelist: whitelist_file,
        whitelist_scores: wl_scores_out,
        anchor_scores: anchor_scores_out,
        tool_hints,
        confidence: confidence.to_string(),
        interpretation,
        notes,
    }
}

/// Build the evidence ledger explaining *why* `winner` was chosen, signal by signal.
/// Mirrors the decision flow in `detect()`: length gives candidates, TSO gives
/// polarity, whitelist gives version, then the final scored verdict.
#[allow(clippy::too_many_arguments)]
fn build_interpretation(
    winner: Option<&AssaySpec>,
    scored: &[(&AssaySpec, f64)],
    dominant: usize,
    candidates: &[&AssaySpec],
    anchor_pcts: &HashMap<&'static [u8], f64>,
    wl_scores: &HashMap<&str, (usize, f64)>,
    cb_sizes: &HashMap<(usize, usize), usize>,
    n_reads: usize,
    has_wl: bool,
    confidence: &str,
) -> Vec<Evidence> {
    let mut ev = Vec::new();
    let winner = match winner {
        Some(w) => w,
        None => return ev,
    };

    let cand_ids: Vec<&str> = candidates.iter().map(|s| s.id).collect();
    ev.push(Evidence {
        signal: "R1 length".into(),
        observation: format!("{dominant} bp dominant across {n_reads} reads"),
        llr_contribution: None,
        conclusion: format!(
            "The barcode read length fixes the cell-barcode + UMI layout, so only {} \
             chemistries with a ~{dominant} bp R1 (\u{00b1}2 bp) are possible. The next \
             signals pick among them: {}",
            candidates.len(),
            cand_ids.join(", "),
        ),
    });

    let (tso_ll, wl_ll) = llr_components(winner, anchor_pcts, n_reads, wl_scores, cb_sizes);

    // TSO polarity, only for assays that model the motif.
    if !winner.r2_anchors.is_empty() {
        let pct = winner
            .r2_anchors
            .iter()
            .filter_map(|a| anchor_pcts.get(a.seq).copied())
            .fold(0.0_f64, f64::max);
        // Wording follows the observed level, not the winner's model: only call it a
        // TSO signal when it actually contributed evidence (one-sided, so absence = 0).
        let conclusion = if tso_ll > 1.0 {
            format!(
                "TSO present in {pct:.1}% of R2 (background ~{:.0}%): the template-switch oligo is detected, which marks 3' chemistry.",
                P_TSO_BACKGROUND * 100.0,
            )
        } else {
            format!(
                "TSO in only {pct:.1}% of R2, at the ~{:.0}% background: no template-switch signal, so the TSO can't tell 3' from 5' here.",
                P_TSO_BACKGROUND * 100.0,
            )
        };
        let needle = String::from_utf8_lossy(&TSO_PRESENT.seq[..TSO_PRESENT.min_match_len]);
        ev.push(Evidence {
            signal: "TSO motif (R2)".into(),
            observation: format!(
                "{}-mer {needle} in {pct:.1}% of R2[0..{}]",
                TSO_PRESENT.min_match_len, TSO_PRESENT.scan_end,
            ),
            llr_contribution: Some(fastq::round2(tso_ll)),
            conclusion,
        });
    }

    if let Some(wl) = &winner.whitelist {
        if has_wl {
            let win_pct = wl_scores.get(wl.filename).map(|(_, p)| *p).unwrap_or(0.0);
            let alt = wl_scores
                .iter()
                .filter(|(f, _)| **f != wl.filename)
                .max_by(|a, b| {
                    a.1 .1
                        .partial_cmp(&b.1 .1)
                        .unwrap_or(std::cmp::Ordering::Equal)
                });
            let alt_str = alt
                .map(|(f, (_, p))| format!("next pool {} at {:.1}%", f, p))
                .unwrap_or_else(|| "no other pool scored".into());
            let decisive = match alt {
                Some((_, (_, p))) => win_pct >= 50.0 && win_pct >= 3.0 * p.max(0.1),
                None => win_pct >= 50.0,
            };
            let conclusion = if decisive {
                format!(
                    "{:.1}% of barcodes in the {} pool ({}); decisive.",
                    win_pct, wl.filename, alt_str,
                )
            } else {
                format!(
                    "{} matched {:.1}% ({}); not decisive, version from geometry and registry order.",
                    wl.filename, win_pct, alt_str,
                )
            };
            ev.push(Evidence {
                signal: "Whitelist overlap".into(),
                observation: format!("{:.1}% of distinct CBs found in {}", win_pct, wl.filename),
                llr_contribution: Some(fastq::round2(wl_ll)),
                conclusion,
            });
        } else {
            ev.push(Evidence {
                signal: "Whitelist overlap".into(),
                observation: "no whitelist cache".into(),
                llr_contribution: None,
                conclusion: "version within family not checked against barcodes; resolved by length and registry order. Run `preflightx fetch-whitelists` to disambiguate.".into(),
            });
        }
    }

    let best = scored.first().map(|(_, s)| *s).unwrap_or(0.0);
    let second = scored.get(1).map(|(_, s)| *s).unwrap_or(0.0);
    let gap = best - second;
    let runner = scored.get(1).map(|(s, _)| *s);
    let tie = runner.is_some_and(|r| gap.abs() < 1.0 && same_geometry(winner, r));
    let conclusion = if tie {
        format!(
            "{} and {} are indistinguishable by sequence (same geometry and whitelist); \
             chose {} by registry order. Confidence {}.",
            winner.id,
            runner.map(|r| r.id).unwrap_or("?"),
            winner.id,
            confidence,
        )
    } else {
        format!(
            "{}: {:+.0} nats total ({:.3}/read), {:+.0} nats ahead of {}. Confidence {} \
             (high requires a whitelist and a clear lead).",
            winner.id,
            best,
            best / n_reads.max(1) as f64,
            gap,
            runner.map(|r| r.id).unwrap_or("(none)"),
            confidence,
        )
    };
    ev.push(Evidence {
        signal: "Verdict".into(),
        observation: format!("winner {} ({})", winner.id, winner.display_name),
        llr_contribution: Some(fastq::round2(best)),
        conclusion,
    });

    ev
}

/// Evidence ledger for a library with no single-cell signal. Walks the same three
/// signals the positive path uses, but each step says why it *fails* the test, so
/// the reader sees exactly how we concluded "not single-cell" rather than a verdict.
fn build_not_single_cell_interpretation(
    dominant: usize,
    tso_pct: f64,
    max_wl_pct: f64,
    n_reads: usize,
    winner: Option<&AssaySpec>,
) -> Vec<Evidence> {
    let closest = winner
        .map(|w| {
            format!(
                " The closest catalogued geometry is {} ({} bp).",
                w.id,
                w.r1_expected_len()
            )
        })
        .unwrap_or_default();
    vec![
        Evidence {
            signal: "Barcode geometry".into(),
            observation: format!("R1 is {dominant} bp dominant across {n_reads} reads"),
            llr_contribution: None,
            conclusion: format!(
                "Single-cell barcode reads are short (~24-30 bp: a cell barcode + UMI). \
                 {dominant} bp does not match any catalogued single-cell geometry, so no \
                 barcode read is present.{closest}"
            ),
        },
        Evidence {
            signal: "TSO (template-switch oligo)".into(),
            observation: format!("TSO 13-mer in {tso_pct:.1}% of R2 (chance background ~1%)"),
            llr_contribution: None,
            conclusion: "A 10x 3' library carries the TSO in ~20% of R2 reads. At background \
                         level here, there is no 3' single-cell signature."
                .into(),
        },
        Evidence {
            signal: "Barcode whitelist".into(),
            observation: format!("{max_wl_pct:.1}% of R1 barcodes match any 10x onlist"),
            llr_contribution: None,
            conclusion: "Real cell barcodes hit their onlist almost completely. Near-zero \
                         overlap means R1 is not a barcode read."
                .into(),
        },
        Evidence {
            signal: "Verdict".into(),
            observation: "no barcode geometry + no TSO + no whitelist overlap".into(),
            llr_contribution: None,
            conclusion: "All three single-cell signals are absent. This is bulk or genomic \
                         sequencing, not single-cell — no chemistry is reported."
                .into(),
        },
    ]
}

impl ChemReport {
    pub fn print(&self) {
        println!("{}", "Chemistry detection".bold().cyan());
        println!("  R1 dominant length:  {} bp", self.r1_length_dominant);
        if !self.r1_length_distribution.is_empty() {
            print!("  R1 length top-5:    ");
            for (l, c) in &self.r1_length_distribution {
                print!("{}bp={}  ", l, c);
            }
            println!();
        }
        println!("  Reads sampled:       {}", self.reads_sampled);

        let conf = self.confidence.as_str();
        let color = |s: &str| -> ColoredString {
            match conf {
                "high" => s.green(),
                "medium" => s.yellow(),
                _ => s.red(),
            }
        };
        println!(
            "  {} {}  ({})",
            "ASSAY:".bold(),
            color(&self.assay_display).bold(),
            color(conf)
        );
        println!(
            "  Geometry:  CB=R1[{}..{}]  UMI=R1[{}..{}]",
            self.cb_offset,
            self.cb_offset + self.cb_len,
            self.umi_offset,
            self.umi_offset + self.umi_len,
        );
        if let Some(wl) = &self.whitelist {
            println!("  Whitelist: {}", wl);
        }
        println!("  CB candidates (distinct): {}", self.n_cb_candidates);
        if !self.whitelist_scores.is_empty() {
            println!("  Whitelist overlap:");
            for ws in &self.whitelist_scores {
                println!(
                    "    {:<40}  {:>10} hits  {:>6.1}%",
                    ws.filename, ws.hits, ws.pct
                );
            }
        }
        for a in &self.anchor_scores {
            println!("  Anchor [{}...]: {:.1}%", a.motif, a.hit_pct);
        }
        if !self.interpretation.is_empty() {
            println!("  {}", "Why this call:".bold());
            for (i, e) in self.interpretation.iter().enumerate() {
                let llr = e
                    .llr_contribution
                    .map(|v| format!("  [LLR {:+.1}]", v).dimmed().to_string())
                    .unwrap_or_default();
                println!("    {}. {}{}", i + 1, e.signal.bold(), llr);
                println!("       observed: {}", e.observation);
                println!("       {} {}", "implies:".cyan(), e.conclusion);
            }
        }
        if !self.tool_hints.is_empty() {
            println!("  Tool hints:");
            for h in &self.tool_hints {
                println!("    {:12}  {}", h.tool, h.hint);
            }
        }
        for note in &self.notes {
            println!("  {} {}", "note:".dimmed(), note);
        }
    }
}
