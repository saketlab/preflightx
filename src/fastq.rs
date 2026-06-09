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

    thread_local! {
        static FILES: RefCell<HashMap<String, Vec<u8>>> = RefCell::new(HashMap::new());
    }

    /// Store `bytes` (raw file contents, possibly gzipped) under `name`.
    pub fn put(name: &str, bytes: Vec<u8>) {
        FILES.with(|f| f.borrow_mut().insert(name.to_string(), bytes));
    }

    /// Fetch a stored buffer by name, if present.
    pub fn get(name: &str) -> Option<Vec<u8>> {
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

/// Pull the next record, or `None` at end of stream. A parse error becomes `None`
/// on wasm — a Range-sliced remote gzip can truncate its final record, so we keep
/// what was already read — but stays a hard error on native.
#[cfg_attr(target_arch = "wasm32", allow(unused_variables))]
fn next_record<'a>(
    reader: &'a mut dyn FastxReader,
    path: &Path,
) -> Result<Option<needletail::parser::SequenceRecord<'a>>> {
    match reader.next() {
        Some(Ok(r)) => Ok(Some(r)),
        None => Ok(None),
        Some(Err(_e)) => {
            #[cfg(target_arch = "wasm32")]
            return Ok(None);
            #[cfg(not(target_arch = "wasm32"))]
            return Err(anyhow::Error::new(_e))
                .with_context(|| format!("parsing {}", path.display()));
        }
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
    while let Some(r) = next_record(reader.as_mut(), path)? {
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
            let s1 = match next_record(r1_rd.as_mut(), r1)? {
                Some(r) => r.seq().to_vec(),
                None => break,
            };
            let s2 = match next_record(r2_rd.as_mut(), r2)? {
                Some(r) => r.seq().to_vec(),
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
