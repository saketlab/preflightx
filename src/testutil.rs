//! Shared test fixtures, compiled only under cfg(test).

use std::io::Write;
use std::path::PathBuf;

/// reads as minimal FASTQ bytes (uniform 'I' quality).
pub fn fastq_bytes(reads: &[Vec<u8>]) -> Vec<u8> {
    let mut buf = Vec::new();
    for (i, r) in reads.iter().enumerate() {
        writeln!(buf, "@r{i}").unwrap();
        buf.extend_from_slice(r);
        buf.extend_from_slice(b"\n+\n");
        buf.extend_from_slice(&vec![b'I'; r.len()]);
        buf.push(b'\n');
    }
    buf
}

/// Write reads as a minimal FASTQ to a temp file and return its path. The caller is
/// responsible for removing the file.
pub fn write_fastq(name: &str, reads: &[Vec<u8>]) -> PathBuf {
    let mut p = std::env::temp_dir();
    p.push(format!("preflightx_{}_{}.fastq", std::process::id(), name));
    std::fs::write(&p, fastq_bytes(reads)).unwrap();
    p
}

/// count fixed-length reads cycling through ACGT (only the length matters).
pub fn fixed_len_reads(len: usize, count: usize) -> Vec<Vec<u8>> {
    (0..count)
        .map(|i| (0..len).map(|j| b"ACGT"[(i + j) % 4]).collect())
        .collect()
}
