#!/bin/bash
set -e

R1=$1
R2=$2
mode=$3
index=$4

if [ -z $mode ] || [ $mode == "LE" ] || [ $mode == "RE" ] || [ -z $index ]; then
    echo "Usage: $0 <reads.R1.fq.gz> <reads.R2.fq.gz> <mode> <index>"
    echo "mode: LE or RE"
    exit 1
fi

R1_filetail=.R1.fq.gz
R2_filetail=.R2.fq.gz

# ----
filename=${R1/${R1_filetail}/}
# 判断文件是否存在
if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "File does not exist, skip"
    continue
fi
echo "---- ${filename} ----"


# ---- raw2clean ----
echo "quality filter"
fastp \
  -q 20 -u 40 -n 2 -l 30 -w 8 \
  -i $R1 -I $R2 \
  -o clean_$R1 -O  clean_$R2 \
  -h $filename.report.html
# ---- End raw2clean ----


# ---- Filter Reads ----
# 过滤FASTQ文件，获取UMI序列
echo "filter reads"
python -m ShCAST_FilterReads.py --fq1 clean_$R1 --fq2 clean_$R2 --mode $mode
# ---- End Reads ----


# ---- cutadapt ----
# 去除handle adapter
R1_handle_tail=GTCACTGTACA
R2_handle_tail=TGTACAGTGAC

echo "cut handle adapter"
cutadapt -u $biotin_primer_base \
  -g $R1_handle_tail -a AAGGCAGCTGGAAGCTGATATT \
  -G ATCAGCTTCCAGCTGCCTT -A $R2_handle_tail \
  -n 2 -m 30 -q 20,15 -e 0.2 --trimmed-only \
  -o Cutadapt_FliterPass_clean_$R1 -p Cutadapt_FliterPass_clean_$R2 \
  FliterPass_clean_$R1 FliterPass_clean_$R2
# ---- cutadapt ----


# ---- mapping ----
echo "mapping"
bowtie2 -q -p 16 -x $index -1 Cutadapt_FliterPass_clean_$R1 -2 Cutadapt_FliterPass_clean_$R2 -S | \
  samtools view -@ 16 -F4 -f 64 -S -o ${filename}_Cas12k_bowtie2.sam
# ---- mapping ----


# ----count insertion sites ----
echo "count insertion sites"
python -m ShCAST_CountSite.py --samfilename ${filename}_Cas12k_bowtie2.sam --outfile ${filename}_Cas12k_insertion_sites.txt
# ---- count insertion sites ----
