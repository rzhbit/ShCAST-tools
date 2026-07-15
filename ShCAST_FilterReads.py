#!/usr/bin/env python3
"""Filter paired FASTQ reads by Cas12k handles and extract UMIs."""

import argparse
import gzip
import logging
import re
from itertools import islice, zip_longest
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

SEED1_LE = r"\w?AATCTGTCACCGACGACAGATAATTTGTCACTGT|T\w?ATCTGTCACCGACGACAGATAATTTGTCACTGT|TA\w?TCTGTCACCGACGACAGATAATTTGTCACTGT|TAA\w?CTGTCACCGACGACAGATAATTTGTCACTGT|TAAT\w?TGTCACCGACGACAGATAATTTGTCACTGT|TAATC\w?GTCACCGACGACAGATAATTTGTCACTGT|TAATCT\w?TCACCGACGACAGATAATTTGTCACTGT|TAATCTG\w?CACCGACGACAGATAATTTGTCACTGT|TAATCTGT\w?ACCGACGACAGATAATTTGTCACTGT|TAATCTGTC\w?CCGACGACAGATAATTTGTCACTGT|TAATCTGTCA\w?CGACGACAGATAATTTGTCACTGT|TAATCTGTCAC\w?GACGACAGATAATTTGTCACTGT|TAATCTGTCACC\w?ACGACAGATAATTTGTCACTGT|TAATCTGTCACCG\w?CGACAGATAATTTGTCACTGT|TAATCTGTCACCGA\w?GACAGATAATTTGTCACTGT|TAATCTGTCACCGAC\w?ACAGATAATTTGTCACTGT|TAATCTGTCACCGACG\w?CAGATAATTTGTCACTGT|TAATCTGTCACCGACGA\w?AGATAATTTGTCACTGT|TAATCTGTCACCGACGAC\w?GATAATTTGTCACTGT|TAATCTGTCACCGACGACA\w?ATAATTTGTCACTGT|TAATCTGTCACCGACGACAG\w?TAATTTGTCACTGT|TAATCTGTCACCGACGACAGA\w?AATTTGTCACTGT|TAATCTGTCACCGACGACAGAT\w?ATTTGTCACTGT|TAATCTGTCACCGACGACAGATA\w?TTTGTCACTGT|TAATCTGTCACCGACGACAGATAA\w?TTGTCACTGT|TAATCTGTCACCGACGACAGATAAT\w?TGTCACTGT|TAATCTGTCACCGACGACAGATAATT\w?GTCACTGT|TAATCTGTCACCGACGACAGATAATTT\w?TCACTGT|TAATCTGTCACCGACGACAGATAATTTG\w?CACTGT|TAATCTGTCACCGACGACAGATAATTTGT\w?ACTGT|TAATCTGTCACCGACGACAGATAATTTGTC\w?CTGT|TAATCTGTCACCGACGACAGATAATTTGTCA\w?TGT|TAATCTGTCACCGACGACAGATAATTTGTCAC\w?GT|TAATCTGTCACCGACGACAGATAATTTGTCACT\w?T|TAATCTGTCACCGACGACAGATAATTTGTCACTG\w?"
SEED1_RE = r"\w?AATAATTTGTCACAACGACATATAATTAGTCACT|C\w?ATAATTTGTCACAACGACATATAATTAGTCACT|CA\w?TAATTTGTCACAACGACATATAATTAGTCACT|CAA\w?AATTTGTCACAACGACATATAATTAGTCACT|CAAT\w?ATTTGTCACAACGACATATAATTAGTCACT|CAATA\w?TTTGTCACAACGACATATAATTAGTCACT|CAATAA\w?TTGTCACAACGACATATAATTAGTCACT|CAATAAT\w?TGTCACAACGACATATAATTAGTCACT|CAATAATT\w?GTCACAACGACATATAATTAGTCACT|CAATAATTT\w?TCACAACGACATATAATTAGTCACT|CAATAATTTG\w?CACAACGACATATAATTAGTCACT|CAATAATTTGT\w?ACAACGACATATAATTAGTCACT|CAATAATTTGTC\w?CAACGACATATAATTAGTCACT|CAATAATTTGTCA\w?AACGACATATAATTAGTCACT|CAATAATTTGTCAC\w?ACGACATATAATTAGTCACT|CAATAATTTGTCACA\w?CGACATATAATTAGTCACT|CAATAATTTGTCACAA\w?GACATATAATTAGTCACT|CAATAATTTGTCACAAC\w?ACATATAATTAGTCACT|CAATAATTTGTCACAACG\w?CATATAATTAGTCACT|CAATAATTTGTCACAACGA\w?ATATAATTAGTCACT|CAATAATTTGTCACAACGAC\w?TATAATTAGTCACT|CAATAATTTGTCACAACGACA\w?ATAATTAGTCACT|CAATAATTTGTCACAACGACAT\w?TAATTAGTCACT|CAATAATTTGTCACAACGACATA\w?AATTAGTCACT|CAATAATTTGTCACAACGACATAT\w?ATTAGTCACT|CAATAATTTGTCACAACGACATATA\w?TTAGTCACT|CAATAATTTGTCACAACGACATATAA\w?TAGTCACT|CAATAATTTGTCACAACGACATATAAT\w?AGTCACT|CAATAATTTGTCACAACGACATATAATT\w?GTCACT|CAATAATTTGTCACAACGACATATAATTA\w?TCACT|CAATAATTTGTCACAACGACATATAATTAG\w?CACT|CAATAATTTGTCACAACGACATATAATTAGT\w?ACT|CAATAATTTGTCACAACGACATATAATTAGTC\w?CT|CAATAATTTGTCACAACGACATATAATTAGTCA\w?T|CAATAATTTGTCACAACGACATATAATTAGTCAC\w?"
SEED2 = r"AA\w?ATCAGCTTCCAGCTGCCTT|AA[GCT]TATCAGCTTCCAGCTGCCTT|AAT\w?TCAGCTTCCAGCTGCCTT|AAT\w?ATCAGCTTCCAGCTGCCTT|AATA\w?CAGCTTCCAGCTGCCTT|AATA\w?TCAGCTTCCAGCTGCCTT|AATAT\w?AGCTTCCAGCTGCCTT|AATAT\w?CAGCTTCCAGCTGCCTT|AATATC\w?GCTTCCAGCTGCCTT|AATATC\w?AGCTTCCAGCTGCCTT|AATATCA\w?CTTCCAGCTGCCTT|AATATCA\w?GCTTCCAGCTGCCTT|AATATCAG\w?TTCCAGCTGCCTT|AATATCAG\w?CTTCCAGCTGCCTT|AATATCAGC\w?TCCAGCTGCCTT|AATATCAGC\w?TTCCAGCTGCCTT|AATATCAGCT\w?CCAGCTGCCTT|AATATCAGCT\w?TCCAGCTGCCTT|AATATCAGCTT\w?CAGCTGCCTT|AATATCAGCTT\w?CCAGCTGCCTT|AATATCAGCTTC\w?AGCTGCCTT|AATATCAGCTTC\w?CAGCTGCCTT|AATATCAGCTTCC\w?GCTGCCTT|AATATCAGCTTCC\w?AGCTGCCTT|AATATCAGCTTCCA\w?CTGCCTT|AATATCAGCTTCCA\w?GCTGCCTT|AATATCAGCTTCCAG\w?TGCCTT|AATATCAGCTTCCAG\w?CTGCCTT|AATATCAGCTTCCAGC\w?GCCTT|AATATCAGCTTCCAGC\w?TGCCTT|AATATCAGCTTCCAGCT\w?CCTT|AATATCAGCTTCCAGCT\w?GCCTT|AATATCAGCTTCCAGCTG\w?CTT|AATATCAGCTTCCAGCTG\w?CCTT|AATATCAGCTTCCAGCTGC\w?TT|AATATCAGCTTCCAGCTGC\w?CTT|AATATCAGCTTCCAGCTGCC\w?T|AATATCAGCTTCCAGCTGCC[AGC]TT|AATATCAGCTTCCAGCTGCCT\w?"


def read_fastq_records(handle, source_name):
    """Yield validated four-line FASTQ records."""
    iterator = iter(handle)
    record_number = 0
    while True:
        lines = list(islice(iterator, 4))
        if not lines:
            return
        record_number += 1
        if len(lines) != 4:
            raise ValueError(f"Incomplete FASTQ record {record_number} in {source_name}")
        lines = [line.rstrip("\r\n") for line in lines]
        if not lines[0].startswith("@") or not lines[2].startswith("+"):
            raise ValueError(f"Malformed FASTQ record {record_number} in {source_name}")
        if len(lines[1]) != len(lines[3]):
            raise ValueError(f"Sequence/quality length mismatch in record {record_number} of {source_name}")
        yield lines


def prefixed_path(filename, prefix):
    path = Path(filename)
    return path.with_name(prefix + path.name)


def read_identifier(header):
    """Return a paired-read identifier without a conventional /1 or /2 suffix."""
    identifier = header.split(maxsplit=1)[0]
    if identifier.endswith(("/1", "/2")):
        return identifier[:-2]
    return identifier


def rename_record(record, umi, identifier, metadata_header):
    metadata = metadata_header.split(maxsplit=1)
    record[0] = f"{identifier}_{umi}" + (f" {metadata[1]}" if len(metadata) == 2 else "")


def write_record(handle, record):
    handle.write("\n".join(record) + "\n")


def filter_fq(fq1, fq2, biotin_primer_base, umi_length, seed1, seed2):
    """Filter paired FASTQ records and append the extracted UMI to read names."""
    if biotin_primer_base < 0 or umi_length < 1:
        raise ValueError("biotin_primer_base must be non-negative and umi_length must be positive")

    fq1_path, fq2_path = Path(fq1), Path(fq2)
    output_paths = (
        prefixed_path(fq1_path, "FliterPass_"),
        prefixed_path(fq2_path, "FliterPass_"),
        prefixed_path(fq1_path, "FliterFail_"),
        prefixed_path(fq2_path, "FliterFail_"),
    )
    pass_count = fail_count = 0
    seed1_pattern = re.compile(seed1)
    seed2_pattern = re.compile(seed2)

    logging.info("Filtering FASTQ files: %s, %s", fq1_path, fq2_path)
    with gzip.open(fq1_path, "rt") as input1, gzip.open(fq2_path, "rt") as input2, \
            gzip.open(output_paths[0], "wt") as pass1, gzip.open(output_paths[1], "wt") as pass2, \
            gzip.open(output_paths[2], "wt") as fail1, gzip.open(output_paths[3], "wt") as fail2:
        pairs = zip_longest(
            read_fastq_records(input1, str(fq1_path)),
            read_fastq_records(input2, str(fq2_path)),
        )
        for pair_number, (record1, record2) in enumerate(pairs, 1):
            if record1 is None or record2 is None:
                raise ValueError(f"Paired FASTQ files contain different numbers of records (at pair {pair_number})")
            if read_identifier(record1[0]) != read_identifier(record2[0]):
                raise ValueError(
                    f"Paired FASTQ read names differ at pair {pair_number}: "
                    f"{record1[0]!r} != {record2[0]!r}"
                )

            header1 = record1[0]
            header2 = record2[0]
            identifier = header1.split(maxsplit=1)[0]

            r1_has_handle = seed1_pattern.search(record1[1][: biotin_primer_base + 15])
            r2_has_handle = seed1_pattern.search(record2[1][: biotin_primer_base + 15])
            r1_has_seed = seed2_pattern.search(record1[1][:35])
            r2_has_seed = seed2_pattern.search(record2[1][:35])

            if r1_has_handle and r2_has_seed:
                if len(record2[1]) < umi_length:
                    raise ValueError(f"Read sequence is shorter than UMI length at pair {pair_number}")
                umi = record2[1][:umi_length]
                rename_record(record1, umi, identifier, header1)
                rename_record(record2, umi, identifier, header2)
                write_record(pass1, record1)
                write_record(pass2, record2)
                pass_count += 1
            elif r2_has_handle and r1_has_seed:
                if len(record1[1]) < umi_length:
                    raise ValueError(f"Read sequence is shorter than UMI length at pair {pair_number}")
                umi = record1[1][:umi_length]
                rename_record(record2, umi, identifier, header1)
                rename_record(record1, umi, identifier, header2)
                write_record(pass1, record2)
                write_record(pass2, record1)
                pass_count += 1
            else:
                write_record(fail1, record1)
                write_record(fail2, record2)
                fail_count += 1

    logging.info("Filter result: %d passed pairs, %d failed pairs", pass_count, fail_count)
    return tuple(str(path) for path in output_paths)


def main():
    parser = argparse.ArgumentParser(description="Filter FASTQ files based on seed patterns")
    parser.add_argument("--fq1", required=True, help="First gzip-compressed FASTQ file")
    parser.add_argument("--fq2", required=True, help="Second gzip-compressed FASTQ file")
    parser.add_argument("--mode", required=True, choices=("LE", "RE"), help="Library mode")
    parser.add_argument(
        "--biotin_primer_base",
        type=int,
        default=70,
        help="Experiment-specific biotin-primer length (default: 70)",
    )
    parser.add_argument("--umi_length", type=int, default=8, help="UMI length (default: 8)")
    args = parser.parse_args()

    seed1 = SEED1_LE if args.mode == "LE" else SEED1_RE
    pass_fq1, pass_fq2, fail_fq1, fail_fq2 = filter_fq(
        args.fq1,
        args.fq2,
        args.biotin_primer_base,
        args.umi_length,
        seed1,
        SEED2,
    )
    logging.info("Output files:")
    logging.info("Passed reads: %s, %s", pass_fq1, pass_fq2)
    logging.info("Failed reads: %s, %s", fail_fq1, fail_fq2)


if __name__ == "__main__":
    main()
