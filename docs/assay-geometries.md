# preflightx: resolved assay read geometries

_Generated from the vendored seqspec catalogue (`vendor/seqspec`) by
`cargo test gen_assay_docs -- --ignored`. Do not edit by hand._

preflightx derives each assay's **observed barcode-read geometry** straight from
its seqspec description: it picks the barcode-bearing read (preferring the RNA
cell-barcode read for multimodal assays), strips leading library primers/adapters,
and spans the barcode / UMI / linker structure up to the cDNA insert. `[ CB n ]` is
a cell barcode, `[ UMI n ]` a UMI, `[ link n ]` a fixed inter-barcode linker, and
`─▶ cdna` marks where the biological insert begins.

## Read structure (representative assays)

```
10x-RNA-v3   (rna-R1.fastq.gz, 28 bp)
  [ CB 16 ]─[ UMI 12 ]

10x-RNA-5prime   (R1.fastq.gz, 26 bp)
  [ CB 16 ]─[ UMI 10 ]

Drop-seq/Illumina   (R1.fastq.gz, 20 bp)
  [ CB 12 ]─[ UMI 8 ]  ─▶ poly_t

SPLiT-seq   (R2.fastq.gz, 94 bp)
  [ UMI 10 ]─[ CB 8 ]─[ link 30 ]─[ CB 8 ]─[ link 30 ]─[ CB 8 ]

inDropv2/Illumina   (R2.fastq.gz, 44 bp)
  [ CB 8 ]─[ probe 22 ]─[ CB 8 ]─[ UMI 6 ]  ─▶ poly_a

BD-Rhapsody-v1/Illumina   (R1.fastq.gz, 60 bp)
  [ CB 9 ]─[ link 12 ]─[ CB 9 ]─[ link 13 ]─[ CB 9 ]─[ UMI 8 ]  ─▶ poly_t

scifi-RNA-Seq   (R1.fastq.gz, 84 bp)
  [ CB 16 ]─[ s7 14 ]─[ truseq_r1 33 ]─[ UMI 8 ]─[ CB 13 ]  ─▶ polyt

10x-ATAC-RNA-MULTI/Illumina   (rna-R1.fastq.gz, 28 bp)
  [ CB 16 ]─[ UMI 12 ]

```

## All 52 catalogued assays (45 with barcode geometry, 7 plate-based)

| Assay | spec | modalities | barcode read |
|-------|------|------------|--------------|
| 10x-ATAC | `10x_atac` | atac | I2.fastq.gz, 16 bp (CB 16) |
| 10x-RNA-5prime | `10x_rna_5prime` | rna | R1.fastq.gz, 26 bp (CB 16 + UMI 10) |
| 10x-ATAC-RNA | `10x_rna_atac` | rna, atac | rna-R1.fastq.gz, 28 bp (CB 16 + UMI 12) |
| 10x-ATAC-RNA-MULTI/Illumina | `10x_rna_atac_multi` | rna, atac, tag | rna-R1.fastq.gz, 28 bp (CB 16 + UMI 12) |
| 10x-RNA-v1 | `10x_rna_v1` | rna | I1.fastq.gz, 14 bp (CB 14) |
| 10x-RNA-v2 | `10x_rna_v2` | rna | R1.fastq.gz, 26 bp (CB 16 + UMI 10) |
| 10x-RNA-v3 | `10x_rna_v3` | rna | rna-R1.fastq.gz, 28 bp (CB 16 + UMI 12) |
| 10xCRISPR | `10xcrispr` | crispr, rna | rna-R1.fastq.gz, 28 bp (CB 16 + UMI 12) |
| 10xFB-3prime | `10xfb_3prime` | rna, protein | R1.fastq.gz, 28 bp (CB 16 + UMI 12) |
| 10xFB protein and RNA assay 5' | `10xfb_5prime` | rna, protein | R1.fastq.gz, 26 bp (CB 16 + UMI 10) |
| 10xFB-VDJ-5prime | `10xfb_vdj_5prime` | rna, vdj | R1.fastq.gz, 26 bp (CB 16 + UMI 10) |
| Example assay | `10xv3_scrnaseq_element_adept_truseq_dual` | rna | _plate-based, no in-read cell barcode_ |
| Example assay | `10xv3_scrnaseq_illumina_truseq_dual` | rna | _plate-based, no in-read cell barcode_ |
| BD-Rhapsody-EB | `bd_rhapsody_eb` | rna | R1.fastq.gz, 43 bp (CB 27 + UMI 8) |
| BD-Rhapsody-v1/Illumina | `bd_rhapsody_v1` | rna | R1.fastq.gz, 60 bp (CB 27 + UMI 8) |
| CEL-Seq | `cel_seq` | rna | R1.fastq.gz, 8 bp (CB 8) |
| CEL-Seq2/Illumina | `cel_seq2` | rna | R1.fastq.gz, 12 bp (CB 6 + UMI 6) |
| DOGMAseq-DIG/Illumina | `dogmaseq_dig` | protein, tag, rna, atac | R1.fastq.gz, 28 bp (CB 16 + UMI 12) |
| DOGMAseq-LLL/Illumina | `dogmaseq_lll` | protein, rna, atac | R1.fastq.gz, 28 bp (CB 16 + UMI 12) |
| Drop-seq/Illumina | `drop_seq` | rna | R1.fastq.gz, 20 bp (CB 12 + UMI 8) |
| inDropv2/Illumina | `indropv2` | rna | R2.fastq.gz, 44 bp (CB 16 + UMI 6) |
| ISSAAC-seq | `issaac_seq` | RNA, ATAC | RNA-I2.fastq.gz, 16 bp (CB 16) |
| MARS-seq | `mars_seq` | RNA | R2.fastq.gz, 10 bp (CB 6 + UMI 4) |
| mcSCRB-seq | `mcscrb_seq` | RNA | R1.fastq.gz, 16 bp (CB 6 + UMI 10) |
| Microwell-seq | `microwell_seq` | RNA | R1.fastq.gz, 54 bp (CB 18 + UMI 6) |
| WT Mega v2 | `parse_wt_v2` | rna | R2.fastq.gz, 86 bp (CB 24 + UMI 10) |
| Pi-ATAC-seq | `pi_atac_seq` | ATAC | _plate-based, no in-read cell barcode_ |
| PIPseq V2 | `pipseqv2` | RNA | R1.fastq.gz, 50 bp (CB 24 + UMI 12) |
| PIPseq V3 | `pipseqv3` | RNA | R1.fastq.gz, 51 bp (CB 28 + UMI 12) |
| Quartz-seq | `quartz_seq` | RNA | R2.fastq.gz, 6 bp (CB 6) |
| Quartz-seq2 | `quartz_seq2` | RNA | R1.fastq.gz, 23 bp (CB 15 + UMI 8) |
| CRISPRa_screening | `sccrispra` | rna, crispr | R2.fastq.gz, 28 bp (CB 16 + UMI 12) |
| sci-RNA-seq | `sci_rna_seq` | RNA | R1.fastq.gz, 18 bp (CB 10 + UMI 8) |
| sci-RNA-seq3 | `sci_rna_seq3` | RNA | R1.fastq.gz, 33 bp (CB 19 + UMI 8) |
| scifi-RNA-Seq | `scifi_rna_seq` | RNA | R1.fastq.gz, 84 bp (CB 29 + UMI 8) |
| Seq-well | `seq_well` | RNA | R1.fastq.gz, 20 bp (CB 12 + UMI 8) |
| Seq-Well-S3 | `seq_well_s3` | RNA | R1.fastq.gz, 20 bp (CB 12 + UMI 8) |
| SHARE-Seq | `share_seq` | atac, rna | R2.fastq.gz, 84 bp (CB 24) |
| Smart-seq2 | `smart_seq2` | RNA | _plate-based, no in-read cell barcode_ |
| Smart-seq3 | `smart_seq3` | RNA_end, RNA_internal | _plate-based, no in-read cell barcode_ |
| sn-m3C-seq | `sn_m3c_seq` | methyl, hic | R1.fastq.gz, 8 bp (CB 8) |
| snmCT-seq | `snmctseq` | methyl, rna | R1.fastq.gz, 8 bp (CB 8) |
| SPLiT-seq | `split_seq` | RNA | R2.fastq.gz, 94 bp (CB 24 + UMI 10) |
| STRT-seq | `strt_seq` | RNA | R1.fastq.gz, 6 bp (CB 6) |
| STRT-seq-2i | `strt_seq_2i` | RNA | R1.fastq.gz, 33 bp (CB 5 + UMI 6) |
| STRT-seq-C1 | `strt_seq_c1` | RNA | R2.fastq.gz, 8 bp (CB 8) |
| SugarSeq/Illumina | `sugarseq` | rna, tag | R2.fastq.gz, 28 bp (CB 16 + UMI 12) |
| SureCell | `surecell` | RNA | R1.fastq.gz, 59 bp (CB 18 + UMI 8) |
| 10x Genomics | `tang2009` | RNA | _plate-based, no in-read cell barcode_ |
| template | `template` | mode | _plate-based, no in-read cell barcode_ |
| VASA-seq-drop | `vasa_seq_drop` | RNA | R2.fastq.gz, 14 bp (CB 8 + UMI 6) |
| VASA-seq | `vasa_seq_plate` | RNA | R2.fastq.gz, 6 bp (CB 6) |

