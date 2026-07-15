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

- Python 3.10 or newer
- fastp
- cutadapt
- bowtie2
- samtools

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
bash ShCAST_run.sh sample.R1.fq.gz sample.R2.fq.gz LE reference_index 70
```

Arguments:

| Argument | Description |
| --- | --- |
| `sample.R1.fq.gz` | Read 1 FASTQ file. |
| `sample.R2.fq.gz` | Read 2 FASTQ file. |
| `LE` / `RE` | Library mode. Use `LE` for left-end libraries and `RE` for right-end libraries. |
| `reference_index` | Bowtie2 index prefix. |
| `biotin_primer_base` | Optional experiment-specific handle length/search position. Defaults to `70`. The wrapper passes it to read filtering and `cutadapt -u`. |

## Workflow

The complete workflow contains the following steps:

1. Quality filter raw paired-end reads with `fastp`.
2. Identify reads containing the expected Cas12k handle and READ2 seed.
3. Extract the UMI from the beginning of the seed-containing read.
4. Validate paired FASTQ records and write passed and failed read pairs into separate files.
5. Trim remaining handle/adaptor sequences with `cutadapt`.
6. Map filtered reads to the reference genome with `bowtie2`.
7. Filter mapped reads with `samtools`.
8. Infer insertion sites from SAM flags, read orientation, mate position, alignment position, and CIGAR reference span.
9. Sort sites by reference, strand, and position, then merge nearby positions and count supporting UMIs/reads.

## Script Usage

### 1. End-to-end pipeline

```bash
bash ShCAST_run.sh <reads.R1.fq.gz> <reads.R2.fq.gz> <mode> <bowtie2_index> [biotin_primer_base]
```

Example:

```bash
bash ShCAST_run.sh sample.R1.fq.gz sample.R2.fq.gz LE reference_index 70
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

Input validation:

- R1 and R2 must contain the same number of FASTQ records.
- Corresponding R1 and R2 records must have matching read identifiers; conventional `/1` and `/2` suffixes are accepted.
- Every FASTQ record must contain exactly four lines, with headers beginning with `@` and separator lines beginning with `+`.
- Sequence and quality strings must have equal lengths.
- A seed-containing read must be at least as long as the requested UMI length.
- Incomplete, malformed, or unpaired input causes the script to stop with an explanatory error instead of silently dropping reads.

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

Output files are created in the same directory as their corresponding input files. For example, `/data/sample.R1.fq.gz` produces `/data/FliterPass_sample.R1.fq.gz` and `/data/FliterFail_sample.R1.fq.gz`.

### 3. Count insertion sites

```bash
python ShCAST_CountSite.py \
  --samfilename sample_Cas12k_bowtie2.sam \
  --outfile sample_Cas12k_insertion_sites.txt
```

The counting script generates a raw insertion-site file, creates a sorted intermediate file, and then writes a merged count table. Sorting is performed in Python and does not require external Unix sorting utilities.

Intermediate files:

| Output | Description |
| --- | --- |
| `insertion_site_<sample>_Cas12k_bowtie2.txt` | Raw insertion-site records parsed from the SAM file. |
| `sort.insertion_site_<sample>_Cas12k_bowtie2.txt` | Records sorted by reference sequence, strand, and genomic position for counting. |

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

For forward alignments, the insertion position is based on the SAM alignment start. For reverse alignments, the endpoint is calculated from the CIGAR operations that consume reference bases (`M`, `D`, `N`, `=`, and `X`), rather than from read-sequence length alone. This correctly handles insertions, deletions, clipping, and other non-trivial alignments.

Nearby positions are merged only when they belong to the same reference sequence and strand. An input with no qualifying alignments produces a valid output file containing only the header.

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
- The wrapper accepts FASTQ paths from other directories. Intermediate R1/R2 files are written beside their corresponding inputs, and the report and final results are written beside R1. Input filenames must share a sample prefix and end with `.R1.fq.gz` and `.R2.fq.gz`.
- Make sure the library mode is correct. `LE` and `RE` use different READ1 seed sequences.
- The wrapper script derives the sample prefix by removing `.R1.fq.gz` from the R1 filename.
- The counting step keeps reference sequences and strands separate and groups consecutive positions within a 6 bp window.
- SAM records with unmapped/invalid mate coordinates or mate distances greater than 10,000 bp are excluded from insertion-site counting.
- Insertion-site sorting uses a temporary SQLite database, so large runs require sufficient free space in the system temporary directory but do not keep all site records in memory.
- Large datasets may require adjusting thread counts in `ShCAST_run.sh` for `fastp`, `bowtie2`, and `samtools`.

## License

This project is licensed under the Apache License 2.0. See [LICENSE](LICENSE) for details.

This software was independently developed by Ze-Hui Ren and integrates multiple open-source tools into a unified analysis pipeline. Third-party software remains the property of its original authors.
