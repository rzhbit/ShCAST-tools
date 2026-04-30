# ShCAST

Targeted Genetic Screening in Bacteria with a Cas12k-guided Transposase

ShCAST is an analysis pipeline for sequencing data generated from Cas12k-guided transposase experiments in bacteria. It filters paired-end reads by library-specific handle sequences, extracts UMI sequences, trims adapters, maps reads to a reference genome, and reports insertion-site counts.

This repository contains the analysis code associated with the study **"Targeted genetic screening in bacteria with a Cas12k-guided transposase"**.

Contact: renzh6@mail2.sysu.edu.cn

## Features

- Filter paired-end FASTQ reads by LE/RE Cas12k handle sequences.
- Extract UMI sequences and append them to FASTQ read names.
- Retain failed reads for quality control and troubleshooting.
- Trim handle/adaptor sequences with `cutadapt`.
- Map filtered reads to a reference genome with `bowtie2`.
- Parse SAM alignments and count insertion sites by genomic position and UMI support.

## Repository Structure

| File | Description |
| --- | --- |
| `ShCAST_run.sh` | Recommended end-to-end pipeline wrapper. Runs quality filtering, read filtering, adapter trimming, mapping, and insertion-site counting. |
| `ShCAST_FilterReads.py` | Filters paired-end FASTQ files by seed sequences, extracts UMI sequences, and writes passed/failed read pairs. |
| `ShCAST_CountSite.py` | Parses mapped SAM files, infers insertion sites, merges nearby positions, and counts supporting UMIs/reads. |
| `LICENSE` | Apache License 2.0. |

## Requirements

Install the following tools and make sure they are available in your `PATH`:

- Python 3
- fastp
- cutadapt
- bowtie2
- samtools
- Standard Unix utilities: `grep`, `sort`, `cat`, `rm`

Python scripts use only Python standard-library modules.

## Input Data

The pipeline expects paired-end gzipped FASTQ files. The wrapper script assumes the following naming pattern:

```text
sample.R1.fq.gz
sample.R2.fq.gz
```

You also need a pre-built Bowtie2 index for the target reference genome.

Build an index if needed:

```bash
bowtie2-build reference.fa reference_index
```

## Quick Start

Run the full pipeline:

```bash
bash ShCAST_run.sh sample.R1.fq.gz sample.R2.fq.gz LE reference_index
```

Arguments:

| Argument | Description |
| --- | --- |
| `sample.R1.fq.gz` | Read 1 FASTQ file. |
| `sample.R2.fq.gz` | Read 2 FASTQ file. |
| `LE` / `RE` | Library mode. Use `LE` for left-end libraries and `RE` for right-end libraries. |
| `reference_index` | Bowtie2 index prefix. |

## Workflow

The complete workflow contains the following steps:

1. Quality filter raw paired-end reads with `fastp`.
2. Identify reads containing the expected Cas12k handle and READ2 seed.
3. Extract the UMI from the beginning of the seed-containing read.
4. Write passed and failed paired reads into separate FASTQ files.
5. Trim remaining handle/adaptor sequences with `cutadapt`.
6. Map filtered reads to the reference genome with `bowtie2`.
7. Filter mapped reads with `samtools`.
8. Infer insertion sites from SAM flags, read orientation, mate position, and alignment position.
9. Merge nearby insertion positions within a small window and count supporting UMIs/reads.

## Script Usage

### 1. End-to-end pipeline

```bash
bash ShCAST_run.sh <reads.R1.fq.gz> <reads.R2.fq.gz> <mode> <bowtie2_index>
```

Example:

```bash
bash ShCAST_run.sh sample.R1.fq.gz sample.R2.fq.gz LE reference_index
```

### 2. Filter reads and extract UMIs

```bash
python ShCAST_FilterReads.py \
  --fq1 clean_sample.R1.fq.gz \
  --fq2 clean_sample.R2.fq.gz \
  --mode LE \
  --biotin_primer_base 70 \
  --umi_length 8
```

Key parameters:

| Parameter | Default | Description |
| --- | --- | --- |
| `--fq1` | Required | Read 1 FASTQ file. Must be gzip-compressed. |
| `--fq2` | Required | Read 2 FASTQ file. Must be gzip-compressed. |
| `--mode` | Required | Library mode: `LE` or `RE`. |
| `--biotin_primer_base` | `70` | Search window length for the READ1 Cas12k handle. |
| `--umi_length` | `8` | Length of the UMI extracted from the beginning of the seed-containing read. |

Library-specific seed sequences:

| Mode | READ1 seed |
| --- | --- |
| `LE` | `TAATCTGTCACCGACGACAGATAATTTGTCACTGT` |
| `RE` | `CAATAATTTGTCACAACGACATATAATTAGTCACT` |

Shared READ2 seed:

```text
AATATCAGCTTCCAGCTGCCTT
```

Read-filtering outputs:

| Output | Description |
| --- | --- |
| `FliterPass_<input_R1>` / `FliterPass_<input_R2>` | Passed paired reads with UMI appended to the read name. |
| `FliterFail_<input_R1>` / `FliterFail_<input_R2>` | Failed paired reads that do not contain the expected seed combination. |

Note: The output prefix is intentionally documented as `Fliter`, matching the current script output names.

### 3. Count insertion sites

```bash
python ShCAST_CountSite.py \
  --samfilename sample_Cas12k_bowtie2.sam \
  --outfile sample_Cas12k_insertion_sites.txt
```

The counting script first generates an intermediate insertion-site file, then writes a merged count table.

Intermediate files:

| Output | Description |
| --- | --- |
| `insertion_site_<sample>_Cas12k_bowtie2.txt` | Raw insertion-site records parsed from the SAM file. |
| `sort.insertion_site_<sample>_Cas12k_bowtie2.txt` | Sorted insertion-site records used for counting. |

Final output:

```text
#Chr    Pos    strand    umi_count    reads_count
```

| Column | Description |
| --- | --- |
| `Chr` | Reference sequence/chromosome name. |
| `Pos` | Representative insertion-site position. |
| `strand` | Insertion-site strand. |
| `umi_count` | Number of unique UMIs supporting the site. |
| `reads_count` | Number of reads supporting the site. |

## Main Pipeline Outputs

For an input sample named `sample.R1.fq.gz` / `sample.R2.fq.gz`, the full pipeline generates files such as:

| Output | Description |
| --- | --- |
| `clean_sample.R1.fq.gz` / `clean_sample.R2.fq.gz` | Quality-filtered reads generated by `fastp`. |
| `sample.report.html` | `fastp` quality-control report. |
| `FliterPass_clean_sample.R1.fq.gz` / `FliterPass_clean_sample.R2.fq.gz` | Reads passing seed/UMI filtering. |
| `FliterFail_clean_sample.R1.fq.gz` / `FliterFail_clean_sample.R2.fq.gz` | Reads failing seed/UMI filtering. |
| `Cutadapt_FliterPass_clean_sample.R1.fq.gz` / `Cutadapt_FliterPass_clean_sample.R2.fq.gz` | Adapter-trimmed reads used for mapping. |
| `sample_Cas12k_bowtie2.sam` | Mapped and filtered SAM file. |
| `sample_Cas12k_insertion_sites.txt` | Final insertion-site count table. |

## Notes

- Run the pipeline in a clean working directory when possible. Intermediate files are written to the current directory and repeated runs with the same sample name may overwrite existing outputs.
- Make sure the library mode is correct. `LE` and `RE` use different READ1 seed sequences.
- The wrapper script derives the sample prefix by removing `.R1.fq.gz` from the R1 filename.
- The counting step merges positions by strand and groups nearby positions within a 6 bp window.
- Large datasets may require adjusting thread counts in `ShCAST_run.sh` for `fastp`, `bowtie2`, and `samtools`.

## License

This project is licensed under the Apache License 2.0. See [LICENSE](LICENSE) for details.

This software was independently developed by Ze-Hui Ren and integrates multiple open-source tools into a unified analysis pipeline. Third-party software remains the property of its original authors.
