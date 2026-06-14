# preflightx

Pre-flight sanity checks for single-cell RNA-seq FASTQ pairs. Samples the first
~50k reads (never the whole file) to answer: what chemistry, bulk or single-cell,
is the quality sane, which catalogued assay does the geometry match, and
(optionally) the genetic sex of the sample.

## Install

```bash
cargo install --git https://github.com/saketlab/preflightx --locked
```

## Usage

```bash
preflightx check  R1.fq.gz R2.fq.gz     # full report
preflightx chem   R1.fq.gz R2.fq.gz     # 10x chemistry only
preflightx bulk   R1.fq.gz              # bulk vs single-cell (R1 only)
preflightx sex    R1.fq.gz R2.fq.gz     # genetic sex (human/mouse)
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
| `sex` | Genetic sex from Y-gene + XIST k-mers (human/mouse) | no |
| `detect-assay` | Rank seqspec catalogue by geometry + anchors (+ onlist overlap when cached) | optional (non-10x) |
| `assays` / `describe <name>` | List / inspect catalogued assays | no |
| `fetch-whitelists` / `whitelist-dir` | Manage the optional barcode cache (10x + combinatorial onlists) | n/a |

## Example

Streamed straight from ENA (only the leading reads are pulled), here ERR4901831,
a mouse 10x run:

```
$ R1=https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR490/001/ERR4901831/ERR4901831_1.fastq.gz
$ R2=https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR490/001/ERR4901831/ERR4901831_2.fastq.gz

$ preflightx chem $R1 $R2
Chemistry detection
  R1 dominant length:  28 bp
  R1 length top-5:    28bp=50000
  Reads sampled:       50000
  ASSAY: 10x Chromium 3' v3  (high)
  Geometry:  CB=R1[0..16]  UMI=R1[16..28]
  Whitelist: 3M-february-2018.txt.gz
  CB candidates (distinct): 18113
  Whitelist overlap:
    3M-february-2018.txt.gz                        16413 hits    90.6%
    737K-august-2016.txt.gz                          265 hits     1.5%
    3M-3pgex-may-2023.txt.gz                         163 hits     0.9%
    3M-5pgex-jan-2023.txt.gz                          47 hits     0.3%
  Anchor [AAGCAGTGGTATC...]: 12.5%
  (... per-signal evidence ledger omitted ...)
  Tool hints:
    starsolo      --soloType CB_UMI_Simple --soloCBwhitelist 3M-february-2018.txt.gz --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12
    kb            -x 10xv3
    alevin        --chromiumV3
    cellranger    --chemistry SC3Pv3

$ preflightx sex $R1 $R2 --species mouse
Genetic sex (sample-level)
  VERDICT: male   confidence: medium   (species mouse, panel experimental)
  Y-panel hits: 17   XIST hits: 0   Y/XIST ratio: 17.00   reads: 200000
  Y genes hit (4): Ddx3y:11, Eif2s3y:3, Uty:2, Kdm5d:1
    note: panel for 'mouse' is experimental; interpret with care (see experiments/sex_panel)
```


## Genetic sex (experimental)

```bash
preflightx sex R1.fq.gz R2.fq.gz            # single-cell, or paired-end bulk
preflightx sex reads.fq.gz --species mouse  # single-end bulk
```

This subcommand counts if the reads either 'map' to the XIST (X silencers) and 
what is the number of reads 'mapping' to Y-genes. There is no real mapping, but
a k-mer based match (alignment free).  Y genes are often lowly expressed, so it samples 
deeper by default (`--sex-reads`, 200k) and returns `uncertain`/`unknown` when the 
signal is weak.


## Whitelists (optional)

```bash
preflightx fetch-whitelists   # ~180 MB into ~/.cache/preflightx/whitelists/
```

`chem` uses the 10x whitelists to disambiguate v2/v3/v4; `detect-assay` uses non-10x
barcode onlists (inDrop, BD Rhapsody, PIP-seq v2/v3, Quartz-seq2) to confirm a
combinatorial assay by barcode overlap on top of geometry. Both are optional: without a
cache, detection falls back to geometry, TSO, and linker anchors. Nothing is redistributed here — every list is fetched
on demand from the [scg_lib_structs](https://teichlab.github.io/scg_lib_structs/) mirror,
which serves the canonical files unmodified. The 10x lists are governed by 10x Genomics'
[Software License Agreement](https://www.10xgenomics.com/legal/software-license).
Override the cache path with `PREFLIGHTX_CACHE=<dir>`.

## Acknowledgements

preflightx depends on [seqspec](https://github.com/pachterlab/seqspec) (Booeshaghi
et al., *Bioinformatics* 2024) for its assay-agnostic understanding of
single-cell library structures. 

## License

MIT (see `LICENSE`). 

