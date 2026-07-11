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

/// Logs per-pull setup vs throughput when `PREFLIGHTX_FETCH_TIMING` is set.
#[cfg(not(target_arch = "wasm32"))]
struct TimedReader<R> {
    inner: R,
    url: String,
    setup: std::time::Duration,
    start: std::time::Instant,
    bytes: u64,
}

#[cfg(not(target_arch = "wasm32"))]
impl<R: std::io::Read> std::io::Read for TimedReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.bytes += n as u64;
        Ok(n)
    }
}

#[cfg(not(target_arch = "wasm32"))]
impl<R> Drop for TimedReader<R> {
    fn drop(&mut self) {
        let body = self.start.elapsed().as_secs_f64();
        let mb = self.bytes as f64 / 1e6;
        eprintln!(
            "fetch-timing {}: setup(connect+headers)={:.0}ms  pulled={:.2}MB in {:.0}ms ({:.1} MB/s)",
            self.url,
            self.setup.as_secs_f64() * 1000.0,
            mb,
            body * 1000.0,
            if body > 0.0 { mb / body } else { 0.0 },
        );
    }
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

/// True if `path` is an http(s)/ftp URL rather than a local filesystem path.
pub fn is_remote(path: &Path) -> bool {
    let s = path.to_string_lossy();
    s.starts_with("http://") || s.starts_with("https://") || s.starts_with("ftp://")
}

/// Hosts to reach over FTP instead of HTTP (see `prefer_ftp`).
#[cfg(not(target_arch = "wasm32"))]
const FTP_PREFERRED_HOSTS: &[&str] = &["ftp.sra.ebi.ac.uk"];

/// FTP hosts that failed this run; later opens skip straight to HTTP.
#[cfg(not(target_arch = "wasm32"))]
static FTP_DEAD_HOSTS: OnceLock<std::sync::Mutex<std::collections::HashSet<String>>> =
    OnceLock::new();

#[cfg(not(target_arch = "wasm32"))]
fn ftp_dead_hosts() -> &'static std::sync::Mutex<std::collections::HashSet<String>> {
    FTP_DEAD_HOSTS.get_or_init(Default::default)
}

/// Rewrite http(s) to `ftp://` for FTP-preferred hosts. `PREFLIGHTX_NO_FTP` disables it.
#[cfg(not(target_arch = "wasm32"))]
fn prefer_ftp(url: &str) -> String {
    if std::env::var_os("PREFLIGHTX_NO_FTP").is_some() {
        return url.to_string();
    }
    for scheme in ["https://", "http://"] {
        if let Some(rest) = url.strip_prefix(scheme) {
            let host = rest.split('/').next().unwrap_or("");
            if FTP_PREFERRED_HOSTS.contains(&host) {
                return format!("ftp://{rest}");
            }
        }
    }
    url.to_string()
}

#[cfg(not(target_arch = "wasm32"))]
fn is_ftp(url: &str) -> bool {
    url.starts_with("ftp://")
}

/// Split `ftp://host/dir/file` into `("host", "/dir/file")`.
#[cfg(not(target_arch = "wasm32"))]
fn split_ftp(url: &str) -> Result<(String, String)> {
    let rest = url.strip_prefix("ftp://").context("not an ftp url")?;
    let slash = rest.find('/').context("ftp url has no path")?;
    Ok((rest[..slash].to_string(), rest[slash..].to_string()))
}

/// True if `bytes` begin with the gzip magic number.
pub fn is_gzip(bytes: &[u8]) -> bool {
    bytes.len() >= 2 && bytes[0] == 0x1f && bytes[1] == 0x8b
}

/// GET `url` with linear-backoff retries on transient failures. 4xx fails fast.
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
            "warning: fetch of {url} failed (attempt {attempt}/{ATTEMPTS}): {err:#}; retrying"
        );
        std::thread::sleep(std::time::Duration::from_secs(2 * attempt as u64));
    }
    unreachable!()
}

/// Holds the control connection open so the RETR stream isn't torn down.
#[cfg(not(target_arch = "wasm32"))]
struct FtpReader<R> {
    data: R,
    _ctrl: suppaftp::FtpStream,
}

#[cfg(not(target_arch = "wasm32"))]
impl<R: std::io::Read> std::io::Read for FtpReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        self.data.read(buf)
    }
}

/// Anonymous FTP stream. Binary mode is required or gzip bytes get mangled.
#[cfg(not(target_arch = "wasm32"))]
fn ftp_reader(url: &str) -> Result<Box<dyn std::io::Read + Send>> {
    use std::net::ToSocketAddrs;
    use suppaftp::{types::FileType, FtpStream};
    let (host, path) = split_ftp(url)?;
    // Bounded connect so a firewalled port fails fast.
    let addr = (host.as_str(), 21)
        .to_socket_addrs()
        .with_context(|| format!("resolving FTP host {host}"))?
        .next()
        .with_context(|| format!("no address for FTP host {host}"))?;
    let mut ftp = FtpStream::connect_timeout(addr, std::time::Duration::from_secs(10))
        .with_context(|| format!("FTP connect {host}"))?;
    ftp.login("anonymous", "anonymous")
        .context("FTP anonymous login")?;
    ftp.transfer_type(FileType::Binary)
        .context("FTP set binary mode")?;
    let data = ftp
        .retr_as_stream(&path)
        .with_context(|| format!("FTP RETR {path}"))?;
    // DataStream is unbuffered; buffer it for the gzip decoder's small reads.
    let data = std::io::BufReader::with_capacity(64 * 1024, data);
    Ok(Box::new(FtpReader { data, _ctrl: ftp }))
}

#[cfg(not(target_arch = "wasm32"))]
fn http_reader(url: &str) -> Result<Box<dyn std::io::Read + Send>> {
    Ok(Box::new(fetch_with_retry(url)?))
}

/// Byte reader for a remote URL: try FTP where preferred, fall back to HTTP and
/// memoize the dead host.
#[cfg(not(target_arch = "wasm32"))]
fn remote_reader(orig: &str) -> Result<Box<dyn std::io::Read + Send>> {
    let url = prefer_ftp(orig);
    if !is_ftp(&url) {
        return http_reader(&url);
    }
    let host = split_ftp(&url).map(|(h, _)| h).unwrap_or_default();
    if !is_ftp(orig) && ftp_dead_hosts().lock().unwrap().contains(&host) {
        return http_reader(orig);
    }
    match ftp_reader(&url) {
        Ok(r) => Ok(r),
        Err(e) if !is_ftp(orig) => {
            eprintln!("warning: FTP fetch of {url} failed ({e:#}); falling back to HTTP");
            ftp_dead_hosts().lock().unwrap().insert(host);
            http_reader(orig)
        }
        Err(e) => Err(e),
    }
}

/// Per-run cache of remote bytes, keyed by URL, so a mate opened by many checks
/// downloads once.
#[cfg(not(target_arch = "wasm32"))]
mod bodycache {
    use super::Result;
    use std::collections::HashMap;
    use std::io::Read;
    use std::sync::{Arc, Mutex, OnceLock};

    struct Body {
        buf: Vec<u8>,
        upstream: Option<Box<dyn Read + Send>>,
    }

    static CACHE: OnceLock<Mutex<HashMap<String, Arc<Mutex<Body>>>>> = OnceLock::new();

    fn cache() -> &'static Mutex<HashMap<String, Arc<Mutex<Body>>>> {
        CACHE.get_or_init(Default::default)
    }

    /// Own cursor over the shared buffer; pulls from the network only past what's cached.
    pub struct CachedReader {
        body: Arc<Mutex<Body>>,
        pos: usize,
    }

    impl Read for CachedReader {
        fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
            let mut body = self.body.lock().unwrap();
            while body.buf.len() <= self.pos && body.upstream.is_some() {
                let mut chunk = [0u8; 64 * 1024];
                let n = body.upstream.as_mut().unwrap().read(&mut chunk)?;
                if n == 0 {
                    body.upstream = None;
                } else {
                    body.buf.extend_from_slice(&chunk[..n]);
                }
            }
            let avail = body.buf.len().saturating_sub(self.pos);
            let k = avail.min(out.len());
            out[..k].copy_from_slice(&body.buf[self.pos..self.pos + k]);
            self.pos += k;
            Ok(k)
        }
    }

    /// `fetch` runs only on first open of `url`, with the map lock released so
    /// concurrent opens of other URLs don't serialize.
    pub fn reader(
        url: &str,
        fetch: impl FnOnce() -> Result<Box<dyn Read + Send>>,
    ) -> Result<CachedReader> {
        if let Some(b) = cache().lock().unwrap().get(url).cloned() {
            return Ok(CachedReader { body: b, pos: 0 });
        }
        let up = fetch()?;
        let b = Arc::new(Mutex::new(Body {
            buf: Vec::new(),
            upstream: Some(up),
        }));
        cache().lock().unwrap().insert(url.to_string(), b.clone());
        Ok(CachedReader { body: b, pos: 0 })
    }
}

/// Open a fastx reader for a local file or http(s)/ftp URL. Remote inputs stream,
/// cached per-URL via [`bodycache`].
#[cfg(not(target_arch = "wasm32"))]
fn open(path: &Path) -> Result<Box<dyn FastxReader>> {
    if !is_remote(path) {
        return parse_fastx_file(path).with_context(|| format!("opening fastq: {}", path.display()));
    }
    let orig = path.to_string_lossy();
    let t0 = std::time::Instant::now();
    let raw: Box<dyn std::io::Read + Send> =
        Box::new(bodycache::reader(&orig, || remote_reader(&orig))?);
    let setup = t0.elapsed();
    let reader: Box<dyn std::io::Read + Send> =
        if std::env::var_os("PREFLIGHTX_FETCH_TIMING").is_some() {
            Box::new(TimedReader {
                inner: raw,
                url: orig.to_string(),
                setup,
                start: std::time::Instant::now(),
                bytes: 0,
            })
        } else {
            raw
        };
    parse_fastx_reader(reader).with_context(|| format!("parsing remote fastq: {orig}"))
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
    matches!(e.kind, UnexpectedEnd)
        || (matches!(e.kind, Io) && e.msg.contains("unexpected end of file"))
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

/// Run two fallible producers concurrently. Overlaps remote R1/R2 fetches.
#[cfg(not(target_arch = "wasm32"))]
pub fn par_pair<A, B, FA, FB>(fa: FA, fb: FB) -> Result<(A, B)>
where
    FA: FnOnce() -> Result<A> + Send,
    FB: FnOnce() -> Result<B> + Send,
    A: Send,
    B: Send,
{
    std::thread::scope(|s| {
        let ha = s.spawn(fa);
        let b = fb()?;
        let a = ha
            .join()
            .map_err(|_| anyhow::anyhow!("paired fetch thread panicked"))??;
        Ok((a, b))
    })
}

/// Wasm has no threads; run sequentially.
#[cfg(target_arch = "wasm32")]
pub fn par_pair<A, B, FA, FB>(fa: FA, fb: FB) -> Result<(A, B)>
where
    FA: FnOnce() -> Result<A>,
    FB: FnOnce() -> Result<B>,
{
    Ok((fa()?, fb()?))
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

/// First `n` sequences only (no quality copy). The common detector input.
pub fn sample_seqs(path: &Path, n: usize) -> Result<Vec<Vec<u8>>> {
    let mut seqs = Vec::with_capacity(n);
    for_each_record(path, n, |seq, _| seqs.push(seq.to_vec()))?;
    Ok(seqs)
}

/// Of two read files, the one whose first `probe` reads are longer on average: the
/// cDNA insert rather than the short barcode/UMI read. Assay-agnostic (works for 10x,
/// Parse, Drop-seq); ties (e.g. paired full-length reads) return `a`.
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
    let (la, lb) = par_pair(|| mean_len(a, probe), || mean_len(b, probe))?;
    Ok(if la >= lb { a } else { b })
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
        assert!(
            seqs.len() < 50_000,
            "truncation should drop the tail records"
        );
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
    fn prefer_ftp_rewrites_only_ena_hosts() {
        // ENA host over http(s) → ftp, path preserved.
        assert_eq!(
            prefer_ftp("https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/x.fastq.gz"),
            "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/x.fastq.gz"
        );
        // Non-ENA host untouched.
        assert_eq!(
            prefer_ftp("https://example.com/a.fastq.gz"),
            "https://example.com/a.fastq.gz"
        );
        // Already ftp / local path untouched.
        assert_eq!(prefer_ftp("ftp://ftp.sra.ebi.ac.uk/x"), "ftp://ftp.sra.ebi.ac.uk/x");
        assert_eq!(prefer_ftp("/local/x.fastq.gz"), "/local/x.fastq.gz");
    }

    #[test]
    fn split_ftp_separates_host_and_path() {
        let (h, p) = split_ftp("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/x.gz").unwrap();
        assert_eq!(h, "ftp.sra.ebi.ac.uk");
        assert_eq!(p, "/vol1/fastq/x.gz");
        assert!(split_ftp("ftp://nohost").is_err());
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
