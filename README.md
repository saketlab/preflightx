# preflightx

Pre-flight sanity checks for single-cell RNA-seq FASTQ files. `preflightx` samples
only the leading reads to identify library chemistry, bulk vs single-cell geometry,
read quality, sequence composition, seqspec assay structure, species, and optional
genetic sex. It ships as a native Rust CLI and as an in-browser WebAssembly build.

## Install

```bash
cargo install --git https://github.com/saketlab/preflightx --locked
```

## Web app

The browser build runs locally on uploaded or accession-derived FASTQ slices and does
not upload reads to a server:

```text
https://preflightx.saketlab.org/
```


## Usage

```bash
preflightx check R1.fq.gz R2.fq.gz          # full paired-end report
preflightx check reads.fq.gz                # single-end / bulk-style report
preflightx chem R1.fq.gz R2.fq.gz           # 10x chemistry only
preflightx bulk R1.fq.gz                    # bulk vs single-cell geometry
preflightx size R1.fq.gz R2.fq.gz           # read count, length, PHRED summary
preflightx composition R2.fq.gz             # GC/N/complexity on cDNA reads
preflightx detect-assay R1.fq.gz R2.fq.gz   # rank seqspec catalogue matches
preflightx species R1.fq.gz R2.fq.gz        # closed-set species ID from panels
preflightx sex R1.fq.gz R2.fq.gz --species mouse
preflightx assays                           # list embedded seqspec assays
preflightx describe 10x_rna_v3              # inspect one assay structure
preflightx --json check R1.fq.gz R2.fq.gz   # machine-readable output

preflightx --n-reads 200000 check R1.fq.gz R2.fq.gz   # default: 50k
preflightx --sex-reads 500000 sex R1.fq.gz R2.fq.gz    # default: 200k
```

An input may be a local path or an `http(s)` URL. Native remote inputs are streamed;
analysis stops after the requested sample instead of downloading the whole file.

## Subcommands

| Command | Purpose | Optional cache |
|---|---|---|
| `check` | Chemistry, generic assay, bulk/sc, size, composition, species, sex | whitelists + panels |
| `chem` | 10x chemistry from R1 geometry, R2 TSO, and barcode overlap | 10x whitelists |
| `bulk` | Bulk vs in-read cell-barcode geometry | no |
| `size` | Read count, length, and mean PHRED sanity | no |
| `composition` | GC, N, low-complexity, and homopolymer profile | no |
| `detect-assay` | Rank embedded seqspec assays by geometry, linker anchors, and onlists | non-10x onlists |
| `species` | Closed-set organism ID from autosomal k-mers | `.spanel` panels |
| `sex` | Genetic sex from Y-marker and XIST/Xist k-mers | embedded or cached panels |
| `assays` / `describe <name>` | List or inspect catalogued assay structures | no |
| `fetch-whitelists` / `whitelist-dir` | Manage 10x and non-10x barcode onlists | n/a |
| `fetch-panels` / `panel-dir` | Manage species/sex `.spanel` panels | n/a |

## Barcode Whitelists And Panels

```bash
preflightx fetch-whitelists   # barcode onlists into the cache
preflightx fetch-panels       # species/sex panels into the cache
preflightx whitelist-dir
preflightx panel-dir
```

`chem` uses 10x whitelists to disambiguate v2/v3/v4 chemistries. `detect-assay`
uses hosted non-10x onlists for assays such as inDrop, BD Rhapsody, PIP-seq, and
Quartz-seq2. Without a barcode cache, detection falls back to geometry, TSO, and
linker anchors.

`species` and cached-panel `sex` use `.spanel` files: autosomal k-mers for organism  
along with XIST/Y pools for compatible XY/XIST sex calls. Panels are stored in
`PREFLIGHTX_CACHE=<dir>`.

## Development

```bash
cargo fmt --check
cargo clippy --all-targets -- -D warnings
cargo test
./scripts/build-wasm.sh
```

## Acknowledgements

preflightx embeds seqspec assay definitions from
[pachterlab/seqspec](https://github.com/pachterlab/seqspec) and uses the public
[scg_lib_structs](https://teichlab.github.io/scg_lib_structs/) mirror for barcode
onlists. Vendor barcode lists remain governed by their upstream licenses.

## License

MIT (see `LICENSE`).
