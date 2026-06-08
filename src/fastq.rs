use anyhow::{Context, Result};
use needletail::{parse_fastx_file, parse_fastx_reader, FastxReader};
use std::path::Path;
use std::sync::OnceLock;

static HTTP_CLIENT: OnceLock<reqwest::blocking::Client> = OnceLock::new();

fn http_client() -> &'static reqwest::blocking::Client {
    HTTP_CLIENT.get_or_init(reqwest::blocking::Client::new)
}

/// True if `path` is an http(s) URL rather than a local filesystem path.
pub fn is_remote(path: &Path) -> bool {
    let s = path.to_string_lossy();
    s.starts_with("http://") || s.starts_with("https://")
}

/// Open a fastx reader for either a local file or an http(s) URL.
///
/// Remote inputs are streamed: we only pull bytes as records are consumed, so
/// sampling the first N reads of a multi-GB remote fastq downloads just the
/// leading chunk, not the whole file. Compression (gzip/bz2/xz) is detected
/// from the stream's magic bytes by needletail, exactly as for local files.
fn open(path: &Path) -> Result<Box<dyn FastxReader>> {
    if is_remote(path) {
        let url = path.to_string_lossy();
        let resp = http_client()
            .get(url.as_ref())
            .send()
            .with_context(|| format!("fetching fastq: {url}"))?
            .error_for_status()
            .with_context(|| format!("fetching fastq: {url}"))?;
        parse_fastx_reader(resp).with_context(|| format!("parsing remote fastq: {url}"))
    } else {
        parse_fastx_file(path).with_context(|| format!("opening fastq: {}", path.display()))
    }
}

/// Iterate up to `n` records from a fastq(.gz) file or URL, calling `cb` with (seq, qual).
/// `qual` is the raw PHRED+33 byte slice; subtract 33 for the score.
pub fn for_each_record<F>(path: &Path, n: usize, mut cb: F) -> Result<usize>
where
    F: FnMut(&[u8], Option<&[u8]>),
{
    let mut reader = open(path)?;
    let mut seen = 0usize;
    while let Some(rec) = reader.next() {
        let r = rec.with_context(|| format!("parsing record in {}", path.display()))?;
        cb(&r.seq(), r.qual());
        seen += 1;
        if seen >= n {
            break;
        }
    }
    Ok(seen)
}

/// Collect first `n` sequences and (optionally) their qualities as owned Vecs.
pub fn sample(path: &Path, n: usize) -> Result<(Vec<Vec<u8>>, Vec<Vec<u8>>)> {
    let mut seqs = Vec::with_capacity(n);
    let mut quals = Vec::with_capacity(n);
    for_each_record(path, n, |seq, q| {
        seqs.push(seq.to_vec());
        if let Some(qb) = q {
            quals.push(qb.to_vec());
        }
    })?;
    Ok((seqs, quals))
}

/// Stream paired reads from two FASTQ files in batches of `batch_size`.
/// Calls `batch_fn(new_r1_seqs, new_r2_seqs)` after each batch, handing ownership of
/// the batch so callers can move reads out (e.g. `Vec::append`) without cloning.
/// Return `false` from the callback to stop early; returns total pairs read.
pub fn stream_paired_batched<F>(
    r1: &Path,
    r2: &Path,
    n_max: usize,
    batch_size: usize,
    mut batch_fn: F,
) -> Result<usize>
where
    F: FnMut(Vec<Vec<u8>>, Vec<Vec<u8>>) -> bool,
{
    let mut r1_rd = open(r1)?;
    let mut r2_rd = open(r2)?;
    let mut total = 0usize;

    loop {
        let mut b1: Vec<Vec<u8>> = Vec::with_capacity(batch_size);
        let mut b2: Vec<Vec<u8>> = Vec::with_capacity(batch_size);

        while b1.len() < batch_size && total < n_max {
            let s1 = match r1_rd.next() {
                Some(r) => r
                    .with_context(|| format!("parsing {}", r1.display()))?
                    .seq()
                    .to_vec(),
                None => break,
            };
            let s2 = match r2_rd.next() {
                Some(r) => r
                    .with_context(|| format!("parsing {}", r2.display()))?
                    .seq()
                    .to_vec(),
                None => break,
            };
            b1.push(s1);
            b2.push(s2);
            total += 1;
        }

        if b1.is_empty() {
            break;
        }
        let cont = batch_fn(b1, b2);
        if !cont || total >= n_max {
            break;
        }
    }

    Ok(total)
}

/// Mean PHRED of a quality byte slice (subtracting +33 offset).
pub fn mean_phred(qual: &[u8]) -> f64 {
    if qual.is_empty() {
        return 0.0;
    }
    let sum: u64 = qual.iter().map(|&b| b.saturating_sub(33) as u64).sum();
    sum as f64 / qual.len() as f64
}

/// Round to 2 decimal places for report output.
pub(crate) fn round2(x: f64) -> f64 {
    (x * 100.0).round() / 100.0
}
