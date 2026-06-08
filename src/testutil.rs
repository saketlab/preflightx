//! Shared test fixtures, compiled only under `cfg(test)`.

use std::io::Write;
use std::path::PathBuf;

/// Write `reads` as a minimal FASTQ (uniform `'I'` quality) to a temp file and
/// return its path. The caller is responsible for removing the file.
pub fn write_fastq(name: &str, reads: &[Vec<u8>]) -> PathBuf {
    let mut p = std::env::temp_dir();
    p.push(format!("preflightx_{}_{}.fastq", std::process::id(), name));
    let mut f = std::fs::File::create(&p).unwrap();
    for (i, r) in reads.iter().enumerate() {
        writeln!(f, "@r{i}").unwrap();
        f.write_all(r).unwrap();
        writeln!(f).unwrap();
        writeln!(f, "+").unwrap();
        f.write_all(&vec![b'I'; r.len()]).unwrap();
        writeln!(f).unwrap();
    }
    p
}

/// `count` fixed-length reads cycling through ACGT (only the length matters).
pub fn fixed_len_reads(len: usize, count: usize) -> Vec<Vec<u8>> {
    (0..count)
        .map(|i| (0..len).map(|j| b"ACGT"[(i + j) % 4]).collect())
        .collect()
}
