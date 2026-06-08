# preflightx

Pre-flight sanity checks for single-cell RNA-seq FASTQ pairs. Samples the first
~50k reads (never the whole file) to answer: what chemistry, bulk or single-cell,
is the quality sane, and which catalogued assay does the geometry match.

## Install

```bash
cargo install --git https://github.com/saketlab/preflightx --locked
```

## Usage

```bash
preflightx check  R1.fq.gz R2.fq.gz     # full report
preflightx chem   R1.fq.gz R2.fq.gz     # 10x chemistry only
preflightx bulk   R1.fq.gz              # bulk vs single-cell (R1 only)
preflightx detect-assay R1.fq.gz R2.fq.gz   # rank seqspec catalogue
preflightx assays                       # list catalogue
preflightx --json check R1.fq.gz R2.fq.gz   # machine-readable

preflightx --n-reads 200000 check R1.fq.gz R2.fq.gz   # sample more (default 50k)
```

An R1/R2 argument may be a local path or an `http(s)` URL — only the leading reads
are streamed, not the whole file.

## Subcommands

| Command | Purpose | Whitelists? |
|---|---|---|
| `check` | All checks, one report | optional (chem) |
| `chem` | 10x chemistry (R1 length + R2 TSO + whitelist) | optional |
| `bulk` | Bulk vs single-cell | no |
| `size` | Read count, length, PHRED | no |
| `composition` | GC, N, low-complexity, homopolymers | no |
| `detect-assay` | Rank seqspec catalogue by geometry + anchors | no |
| `assays` / `describe <name>` | List / inspect catalogued assays | no |
| `fetch-whitelists` / `whitelist-dir` | Manage the optional 10x cache | n/a |

## Example

```
$ preflightx check R1.fq.gz R2.fq.gz

Chemistry detection
  R1 dominant length:  28 bp
  Best match:          3M-february-2018.txt.gz (99.3%)
  R2 TSO %:            1.2% -> 3' (no TSO)
  VERDICT: 3v3 (high)

Bulk vs single-cell
  R1 length 28 bp, stddev 0.05
  VERDICT: single_cell_10x (high)

Size + quality
  R1: 50000 reads  mean_len=28.0  mean_PHRED=36.2
  R2: 50000 reads  mean_len=91.0  mean_PHRED=35.5
  QUALITY: adequate

Composition
  GC 41.3%  N 0.05%  low-complexity 0.40%
  ORGANISM HINT: mammalian
```


## Whitelists (optional)

```bash
preflightx fetch-whitelists   # ~180 MB into ~/.cache/preflightx/whitelists/
```

Only `chem` uses these, and only to disambiguate v2/v3/v4. The 10x whitelists are
**not** redistributed here — they are fetched on demand from the
[scg_lib_structs](https://teichlab.github.io/scg_lib_structs/) mirror, which serves
the canonical 10x files unmodified. Their use is governed by 10x Genomics'
[Software License Agreement](https://www.10xgenomics.com/legal/software-license).
Override the cache path with `PREFLIGHTX_CACHE=<dir>`.

## Acknowledgements

preflightx depends on [seqspec](https://github.com/pachterlab/seqspec) (Booeshaghi
et al., *Bioinformatics* 2024) for its assay-agnostic understanding of
single-cell library structures. 

## License

MIT (see `LICENSE`). 

