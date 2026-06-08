//! Small statistical primitives shared by the two detectors (`chemistry.rs` and
//! `assaydetect.rs`). Both score reads with the same Bernoulli log-likelihood
//! model, so the model lives here once rather than being copied into each.

use std::collections::HashMap;

/// Log-likelihood ratio of observing `k` of `n` Bernoulli trials under the signal
/// rate `p1` versus the background rate `p0`, in nats.
///
/// Returns 0 when there is no information to gain: either no trials (`n == 0`) or
/// the two rates coincide (`p1 == p0`). The combined guard makes this safe for
/// both callers — chemistry's whitelist/TSO scoring (where `p1` may equal `p0`
/// for absent anchors) and assaydetect's anchor scoring (where `n` may be 0).
pub fn bernoulli_llr(k: usize, n: usize, p1: f64, p0: f64) -> f64 {
    if n == 0 || (p1 - p0).abs() < f64::EPSILON {
        return 0.0;
    }
    (k as f64) * (p1 / p0).ln() + ((n - k) as f64) * ((1.0 - p1) / (1.0 - p0)).ln()
}

/// The most common length among `seqs` (0 for an empty slice).
pub fn dominant_len(seqs: &[Vec<u8>]) -> usize {
    let mut hist: HashMap<usize, usize> = HashMap::new();
    for s in seqs {
        *hist.entry(s.len()).or_insert(0) += 1;
    }
    hist.into_iter()
        .max_by_key(|(_, c)| *c)
        .map(|(l, _)| l)
        .unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn llr_is_zero_with_no_information() {
        assert_eq!(bernoulli_llr(0, 0, 0.8, 0.01), 0.0);
        assert_eq!(bernoulli_llr(5, 10, 0.5, 0.5), 0.0);
    }

    #[test]
    fn llr_positive_when_signal_dominates() {
        // All trials hit at a high signal rate vs a tiny background -> strong positive.
        assert!(bernoulli_llr(100, 100, 0.8, 0.01) > 0.0);
        // No hits where the signal expects many -> negative.
        assert!(bernoulli_llr(0, 100, 0.8, 0.01) < 0.0);
    }

    #[test]
    fn dominant_len_picks_the_mode() {
        let seqs = vec![vec![0u8; 16], vec![0u8; 16], vec![0u8; 28]];
        assert_eq!(dominant_len(&seqs), 16);
        assert_eq!(dominant_len(&[]), 0);
    }
}
