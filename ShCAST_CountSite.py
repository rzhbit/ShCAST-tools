#!/usr/bin/env python3
"""Extract and count Cas12k insertion sites from a SAM alignment file."""

import argparse
import logging
import sqlite3
import tempfile
from collections import Counter
from contextlib import closing
from pathlib import Path
from time import localtime, strftime

TIME_FORMAT = "%Y-%m-%d %H:%M:%S"
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")


def reference_span(cigar: str, sequence: str) -> int:
    """Return the number of reference bases consumed by a SAM alignment."""
    if cigar == "*":
        return len(sequence)
    number = ""
    span = 0
    for char in cigar:
        if char.isdigit():
            number += char
        else:
            if not number:
                raise ValueError(f"Invalid CIGAR string: {cigar}")
            if char in "MDN=X":
                span += int(number)
            number = ""
    if number:
        raise ValueError(f"Invalid CIGAR string: {cigar}")
    return span


def analyse_sam(sam_filename: str) -> Path:
    sam_path = Path(sam_filename)
    if not sam_path.is_file():
        raise FileNotFoundError(f"SAM file does not exist: {sam_path}")

    output_path = sam_path.with_name(f"insertion_site_{sam_path.stem}.txt")
    sorted_path = output_path.with_name(f"sort.{output_path.name}")
    alignment_count = proper_count = 0

    logging.info("Analysing SAM file %s", sam_path)
    with tempfile.TemporaryDirectory(prefix="shcast_count_") as temp_dir:
        database_path = Path(temp_dir) / "sites.sqlite3"
        with closing(sqlite3.connect(database_path)) as database, output_path.open(
            "w", encoding="utf-8", newline=""
        ) as raw_output, sam_path.open(encoding="utf-8") as source:
            database.execute(
                "CREATE TABLE sites (chrom TEXT, position INTEGER, strand TEXT, umi TEXT)"
            )
            for line_number, line in enumerate(source, 1):
                if line.startswith("@"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 11:
                    raise ValueError(f"Malformed SAM record at line {line_number}")

                alignment_count += 1
                flag = int(fields[1])
                position = int(fields[3])
                mate_position = int(fields[7])
                if position <= 0 or mate_position <= 0 or abs(mate_position - position) > 10_000:
                    continue

                if flag & 0x20:  # mate is reverse-complemented
                    if mate_position - position < -10:
                        continue
                    site_position = position
                    strand = "+"
                elif flag & 0x10:  # this read is reverse-complemented
                    if mate_position - position > 10:
                        continue
                    site_position = position + reference_span(fields[5], fields[9])
                    strand = "-"
                else:
                    continue

                umi = fields[0].rsplit("_", 1)[-1]
                record = (fields[2], site_position, strand, umi)
                raw_output.write("\t".join(map(str, record)) + "\n")
                database.execute("INSERT INTO sites VALUES (?, ?, ?, ?)", record)
                proper_count += 1

            database.commit()
            with sorted_path.open("w", encoding="utf-8", newline="") as sorted_output:
                rows = database.execute(
                    "SELECT chrom, position, strand, umi FROM sites "
                    "ORDER BY chrom, strand, position"
                )
                for record in rows:
                    sorted_output.write("\t".join(map(str, record)) + "\n")

    logging.info(
        "SAM file contains %d alignments; %d support insertion sites",
        alignment_count,
        proper_count,
    )
    return sorted_path


def write_group(output, group: list[tuple[str, int, str, str]]) -> None:
    representative = Counter((chrom, pos, strand) for chrom, pos, strand, _ in group).most_common(1)[0][0]
    umi_count = len({umi for *_, umi in group})
    output.write(f"{representative[0]}\t{representative[1]}\t{representative[2]}\t{umi_count}\t{len(group)}\n")


def count_sites(result_file: str | Path, output_file: str) -> None:
    result_path = Path(result_file)
    group = []
    previous_key = None
    previous_position = None

    with result_path.open(encoding="utf-8") as source, Path(output_file).open(
        "w", encoding="utf-8", newline=""
    ) as output:
        output.write("#Chr\tPos\tstrand\tumi_count\treads_count\n")
        for line_number, line in enumerate(source, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 4:
                raise ValueError(f"Malformed insertion-site record at line {line_number}")
            record = (fields[0], int(fields[1]), fields[2], fields[3])
            key = (record[0], record[2])
            if group and (key != previous_key or record[1] - previous_position > 6):
                write_group(output, group)
                group = []
            group.append(record)
            previous_key = key
            previous_position = record[1]
        if group:
            write_group(output, group)


def main() -> None:
    parser = argparse.ArgumentParser(description="Count insertion sites")
    parser.add_argument("--samfilename", required=True, help="Input SAM file")
    parser.add_argument("--outfile", required=True, help="Output count table")
    args = parser.parse_args()

    start_time = strftime(TIME_FORMAT, localtime())
    sorted_file = analyse_sam(args.samfilename)
    count_sites(sorted_file, args.outfile)
    logging.info("Finished (started %s)", start_time)


if __name__ == "__main__":
    main()
