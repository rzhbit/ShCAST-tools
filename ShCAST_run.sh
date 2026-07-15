#!/bin/bash
set -euo pipefail

if [ "$#" -lt 4 ] || [ "$#" -gt 5 ]; then
    echo "Usage: $0 <reads.R1.fq.gz> <reads.R2.fq.gz> <mode> <index> [biotin_primer_base]"
    echo "mode: LE or RE"
    echo "biotin_primer_base: non-negative integer (default: 70)"
    exit 1
fi

R1=$1
R2=$2
mode=$3
index=$4
biotin_primer_base=${5:-70}
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

if [ "$mode" != "LE" ] && [ "$mode" != "RE" ]; then
    echo "Error: mode must be LE or RE" >&2
    exit 1
fi

case "$biotin_primer_base" in
    ''|*[!0-9]*)
        echo "Error: biotin_primer_base must be a non-negative integer" >&2
        exit 1
        ;;
esac

if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "Error: one or both input FASTQ files do not exist" >&2
    exit 1
fi

R1_dir=$(cd -- "$(dirname -- "$R1")" && pwd)
R2_dir=$(cd -- "$(dirname -- "$R2")" && pwd)
R1_name=$(basename -- "$R1")
R2_name=$(basename -- "$R2")

if [[ "$R1_name" != *.R1.fq.gz ]] || [[ "$R2_name" != *.R2.fq.gz ]]; then
    echo "Error: input names must end with .R1.fq.gz and .R2.fq.gz" >&2
    exit 1
fi

sample=${R1_name%.R1.fq.gz}
if [ "${R2_name%.R2.fq.gz}" != "$sample" ]; then
    echo "Error: R1 and R2 sample names do not match" >&2
    exit 1
fi

clean_R1="$R1_dir/clean_$R1_name"
clean_R2="$R2_dir/clean_$R2_name"
pass_R1="$R1_dir/FliterPass_clean_$R1_name"
pass_R2="$R2_dir/FliterPass_clean_$R2_name"
trimmed_R1="$R1_dir/Cutadapt_FliterPass_clean_$R1_name"
trimmed_R2="$R2_dir/Cutadapt_FliterPass_clean_$R2_name"
result_prefix="$R1_dir/$sample"

echo "---- $sample ----"

echo "quality filter"
fastp \
  -q 20 -u 40 -n 2 -l 30 -w 8 \
  -i "$R1" -I "$R2" \
  -o "$clean_R1" -O "$clean_R2" \
  -h "$result_prefix.report.html"

echo "filter reads"
python "$script_dir/ShCAST_FilterReads.py" \
  --fq1 "$clean_R1" --fq2 "$clean_R2" --mode "$mode" \
  --biotin_primer_base "$biotin_primer_base"

R1_handle_tail=GTCACTGTACA
R2_handle_tail=TGTACAGTGAC

echo "cut handle adapter"
cutadapt -u "$biotin_primer_base" \
  -g "$R1_handle_tail" -a AAGGCAGCTGGAAGCTGATATT \
  -G ATCAGCTTCCAGCTGCCTT -A "$R2_handle_tail" \
  -n 2 -m 30 -q 20,15 -e 0.2 --trimmed-only \
  -o "$trimmed_R1" -p "$trimmed_R2" \
  "$pass_R1" "$pass_R2"

echo "mapping"
bowtie2 -q -p 16 -x "$index" -1 "$trimmed_R1" -2 "$trimmed_R2" | \
  samtools view -@ 16 -F 4 -f 64 -h -o "${result_prefix}_Cas12k_bowtie2.sam"

echo "count insertion sites"
python "$script_dir/ShCAST_CountSite.py" \
  --samfilename "${result_prefix}_Cas12k_bowtie2.sam" \
  --outfile "${result_prefix}_Cas12k_insertion_sites.txt"
