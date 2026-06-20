use anyhow::{Context, Result};
use needletail::{parse_fastx_reader, FastxReader};
use std::path::Path;

#[cfg(not(target_arch = "wasm32"))]
use needletail::parse_fastx_file;
#[cfg(not(target_arch = "wasm32"))]
use std::sync::OnceLock;

#[cfg(not(target_arch = "wasm32"))]
static HTTP_CLIENT: OnceLock<reqwest::blocking::Client> = OnceLock::new();

#[cfg(not(target_arch = "wasm32"))]
fn http_client() -> &'static reqwest::blocking::Client {
    HTTP_CLIENT.get_or_init(reqwest::blocking::Client::new)
}

/// In-browser virtual file system. The wasm-bindgen entry point ([`crate::wasm`])
/// loads the uploaded FASTQ bytes here keyed by a synthetic name (e.g. `"r1"`),
/// then calls the unchanged `&Path`-based analysis functions with that name as
/// the path. `open()` resolves such a path to an in-memory reader, so the entire
/// analysis layer is shared with the native build with no signature changes.
#[cfg(target_arch = "wasm32")]
pub mod vfs {
    use std::cell::RefCell;
    use std::collections::HashMap;
    use std::sync::Arc;

    thread_local! {
        static FILES: RefCell<HashMap<String, Arc<[u8]>>> = RefCell::new(HashMap::new());
    }

    /// Store `bytes` (raw file contents, possibly gzipped) under `name`.
    pub fn put(name: &str, bytes: Vec<u8>) {
        FILES.with(|f| f.borrow_mut().insert(name.to_string(), Arc::from(bytes)));
    }

    /// Fetch a stored buffer by name, if present.
    pub fn get(name: &str) -> Option<Arc<[u8]>> {
        FILES.with(|f| f.borrow().get(name).cloned())
    }

    /// Drop a single stored buffer.
    pub fn remove(name: &str) {
        FILES.with(|f| f.borrow_mut().remove(name));
    }

    /// True if any stored key starts with `prefix` (used to detect cached onlists).
    pub fn has_prefix(prefix: &str) -> bool {
        FILES.with(|f| f.borrow().keys().any(|k| k.starts_with(prefix)))
    }

    /// Drop all stored buffers (call after a run to free memory).
    pub fn clear() {
        FILES.with(|f| f.borrow_mut().clear());
    }
}

/// True if `path` is an http(s) URL rather than a local filesystem path.
pub fn is_remote(path: &Path) -> bool {
    let s = path.to_string_lossy();
    s.starts_with("http://") || s.starts_with("https://")
}

/// True if `bytes` begin with the gzip magic number.
pub fn is_gzip(bytes: &[u8]) -> bool {
    bytes.len() >= 2 && bytes[0] == 0x1f && bytes[1] == 0x8b
}

/// GET `url`, retrying transient failures (connection drops, timeouts, 5xx) with
/// linear backoff; remote FASTQ hosts (SRA/ENA) flake intermittently. 4xx fails fast.
#[cfg(not(target_arch = "wasm32"))]
fn fetch_with_retry(url: &str) -> Result<reqwest::blocking::Response> {
    const ATTEMPTS: u32 = 4;
    for attempt in 1..=ATTEMPTS {
        let err = match http_client()
            .get(url)
            .send()
            .and_then(|r| r.error_for_status())
        {
            Ok(ok) => return Ok(ok),
            // 4xx is permanent; send errors (no status) fall through to retry.
            Err(e) if e.status().map(|s| s.is_client_error()).unwrap_or(false) => {
                return Err(e).with_context(|| format!("fetching fastq: {url}"));
            }
            Err(e) => anyhow::Error::new(e),
        };
        if attempt == ATTEMPTS {
            return Err(err)
                .with_context(|| format!("fetching fastq after {ATTEMPTS} attempts: {url}"));
        }
        eprintln!(
            "warning: fetch of {url} failed (attempt {attempt}/{ATTEMPTS}): {err:#}; retrying…"
        );
        std::thread::sleep(std::time::Duration::from_secs(2 * attempt as u64));
    }
    unreachable!()
}

/// Open a fastx reader for a local file, an http(s) URL, or (on wasm) an
/// in-memory buffer from [`vfs`].
///
/// Remote inputs are streamed: we only pull bytes as records are consumed, so
/// sampling the first N reads of a multi-GB remote fastq downloads just the
/// leading chunk, not the whole file. Compression (gzip/bz2/xz) is detected
/// from the stream's magic bytes by needletail, exactly as for local files.
#[cfg(not(target_arch = "wasm32"))]
fn open(path: &Path) -> Result<Box<dyn FastxReader>> {
    if is_remote(path) {
        let url = path.to_string_lossy();
        let resp = fetch_with_retry(url.as_ref())?;
        parse_fastx_reader(resp).with_context(|| format!("parsing remote fastq: {url}"))
    } else {
        parse_fastx_file(path).with_context(|| format!("opening fastq: {}", path.display()))
    }
}

/// Wasm build: there is no filesystem and no blocking network. Input must have
/// been registered in [`vfs`] by the bindgen entry point.
///
/// needletail's built-in compression sniffing is gated behind its `compression`
/// feature, which also pulls in the C-backed bz2/xz decoders that don't build
/// for wasm. We therefore link only `flate2` and detect+decompress gzip here by
/// hand (the common case for FASTQ); the decoder is streaming, so sampling the
/// first N reads never inflates the whole buffer. Plain (uncompressed) FASTQ is
/// passed straight through; bz2/xz inputs are unsupported in the browser.
#[cfg(target_arch = "wasm32")]
fn open(path: &Path) -> Result<Box<dyn FastxReader>> {
    use std::io::Cursor;
    let key = path.to_string_lossy();
    let bytes = vfs::get(&key).with_context(|| {
        format!("no in-memory fastq registered as {key:?} (wasm build has no filesystem)")
    })?;
    if is_gzip(&bytes) {
        let dec = flate2::read::MultiGzDecoder::new(Cursor::new(bytes));
        parse_fastx_reader(dec).with_context(|| format!("parsing in-memory gzip fastq: {key}"))
    } else {
        parse_fastx_reader(Cursor::new(bytes))
            .with_context(|| format!("parsing in-memory fastq: {key}"))
    }
}

/// Next record, or `None` at end of stream. A truncated final record reads as clean EOF
/// (inputs are sliced/streamed, so a cut-off tail is expected); real malformation errors.
fn next_record<'a>(
    reader: &'a mut dyn FastxReader,
    path: &Path,
) -> Result<Option<needletail::parser::SequenceRecord<'a>>> {
    match reader.next() {
        Some(Ok(r)) => Ok(Some(r)),
        None => Ok(None),
        Some(Err(e)) if is_truncation(&e) => Ok(None),
        Some(Err(e)) => {
            Err(anyhow::Error::new(e)).with_context(|| format!("parsing {}", path.display()))
        }
    }
}

/// Truncation: a truncated FASTQ record (`UnexpectedEnd`) or gzip EOF (`Io`, which
/// needletail 0.5.1 only exposes via the message). `truncated_gzip_error_shape_is_pinned`
/// guards the message match against a version bump.
fn is_truncation(e: &needletail::errors::ParseError) -> bool {
    use needletail::errors::ParseErrorKind::{Io, UnexpectedEnd};
    matches!(e.kind, UnexpectedEnd) || (matches!(e.kind, Io) && e.msg.contains("unexpected end of file"))
}

/// Iterate up to `n` records from a fastq(.gz) file or URL, calling `cb` with (seq, qual).
/// `qual` is the raw PHRED+33 byte slice; subtract 33 for the score.
pub fn for_each_record<F>(path: &Path, n: usize, mut cb: F) -> Result<usize>
where
    F: FnMut(&[u8], Option<&[u8]>),
{
    let mut reader = open(path)?;
    let mut seen = 0usize;
    while let Some(r) = next_record(reader.as_mut(), path)? {
        cb(&r.seq(), r.qual());
        seen += 1;
        if seen >= n {
            break;
        }
    }
    Ok(seen)
}

/// Owned sampled reads: `(sequences, qualities)`.
pub type SampledReads = (Vec<Vec<u8>>, Vec<Vec<u8>>);

/// Collect first `n` sequences and (optionally) their qualities as owned Vecs.
pub fn sample(path: &Path, n: usize) -> Result<SampledReads> {
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

/// Of two read files, the one whose first `probe` reads are longer on average — i.e. the
/// cDNA insert rather than the short barcode/UMI read. Assay-agnostic (works for 10x,
/// Parse, Drop-seq, …); ties (e.g. paired full-length reads) return `a`.
pub fn longer_read<'a>(a: &'a Path, b: &'a Path, probe: usize) -> Result<&'a Path> {
    fn mean_len(p: &Path, probe: usize) -> Result<f64> {
        let mut total = 0u64;
        let mut n = 0u64;
        for_each_record(p, probe, |seq, _| {
            total += seq.len() as u64;
            n += 1;
        })?;
        Ok(if n == 0 { 0.0 } else { total as f64 / n as f64 })
    }
    Ok(if mean_len(a, probe)? >= mean_len(b, probe)? {
        a
    } else {
        b
    })
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
            let s1 = next_record(r1_rd.as_mut(), r1)?;
            let s2 = next_record(r2_rd.as_mut(), r2)?;
            match (s1, s2) {
                (Some(s1), Some(s2)) => {
                    b1.push(s1.seq().to_vec());
                    b2.push(s2.seq().to_vec());
                }
                (None, None) => break,
                _ => {
                    anyhow::bail!(
                        "paired FASTQ inputs ended at different lengths: {} and {} are out of sync",
                        r1.display(),
                        r2.display()
                    );
                }
            }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::write_fastq;

    /// A gzipped FASTQ cut mid-stream. Large enough to span several decode buffers, else
    /// the cut fixture decodes to nothing.
    fn truncated_gzip_fixture() -> Vec<u8> {
        use flate2::{write::GzEncoder, Compression};
        use std::io::Write;
        let fq = crate::testutil::fastq_bytes(&crate::testutil::fixed_len_reads(30, 50_000));
        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        enc.write_all(&fq).unwrap();
        let mut gz = enc.finish().unwrap();
        gz.truncate(gz.len() * 7 / 10);
        gz
    }

    #[test]
    fn truncated_gzip_samples_complete_records_then_stops() {
        let mut p = std::env::temp_dir();
        p.push(format!("preflightx_{}_trunc.fastq.gz", std::process::id()));
        std::fs::write(&p, truncated_gzip_fixture()).unwrap();
        let result = sample(&p, 10_000);
        std::fs::remove_file(&p).ok();

        let (seqs, _) = result.expect("truncated gzip must sample cleanly, not error");
        assert!(!seqs.is_empty(), "expected the complete leading records");
        assert!(seqs.len() < 50_000, "truncation should drop the tail records");
    }

    /// Pins the error shape `is_truncation` matches on, so a needletail bump fails here.
    #[test]
    fn truncated_gzip_error_shape_is_pinned() {
        use needletail::errors::ParseErrorKind;
        let mut p = std::env::temp_dir();
        p.push(format!("preflightx_{}_shape.fastq.gz", std::process::id()));
        std::fs::write(&p, truncated_gzip_fixture()).unwrap();
        let mut reader = needletail::parse_fastx_file(&p).unwrap();
        let err = loop {
            match reader.next() {
                Some(Ok(_)) => continue,
                Some(Err(e)) => break e,
                None => panic!("expected a truncation error from the cut gzip"),
            }
        };
        std::fs::remove_file(&p).ok();
        assert!(is_truncation(&err));
        assert!(
            matches!(err.kind, ParseErrorKind::Io) && err.msg.contains("unexpected end of file"),
            "needletail truncation error shape changed: {:?} / {}",
            err.kind,
            err.msg
        );
    }

    #[test]
    fn paired_stream_rejects_mismatched_mate_lengths() {
        let r1 = write_fastq("fastq_pair_r1", &[vec![b'A'; 20], vec![b'C'; 20]]);
        let r2 = write_fastq("fastq_pair_r2", &[vec![b'G'; 20]]);
        let err = stream_paired_batched(&r1, &r2, 10, 10, |_, _| true).unwrap_err();
        std::fs::remove_file(&r1).ok();
        std::fs::remove_file(&r2).ok();
        assert!(
            err.to_string().contains("out of sync"),
            "unexpected error: {err:#}"
        );
    }
}
