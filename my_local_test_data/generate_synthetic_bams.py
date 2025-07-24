#!/usr/bin/env python3
"""Generate synthetic BAM files with one read per peak.

Usage:
  python generate_synthetic_bams.py <peaks_csv> <output_bam>

The script reads a peaks file (CSV with columns chr,start,end,...) and emits
SAM records with 100M alignments covering the start coordinate of each peak.
It then converts the SAM to BAM using samtools (must be on PATH) and indexes
it. The resulting BAM has minimal yet valid alignments so that tools like
DiffBind can count reads overlapping the supplied peaks.

We purposely keep the read depth minimal (1 read per peak) to reduce file size.
"""
import csv
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

def build_sam(peaks_csv: Path) -> str:
    """Return SAM text (str)."""
    chrom_max_pos: dict[str, int] = defaultdict(int)
    sam_body: list[str] = []
    seq = "A" * 100
    qual = "I" * 100  # high quality
    with peaks_csv.open() as handle:
        reader = csv.DictReader(handle)
        fieldnames = set(reader.fieldnames or [])
        if not {"chr", "start", "end"}.issubset(fieldnames):
            raise ValueError("Peaks file must contain chr,start,end columns")
        for idx, row in enumerate(reader, start=1):
            chrom = row["chr"].strip()
            start = int(row["start"]) + 1  # SAM is 1-based
            chrom_max_pos[chrom] = max(chrom_max_pos[chrom], start + 100)
            sam_body.append(
                f"read{idx}\t0\t{chrom}\t{start}\t255\t100M\t*\t0\t0\t{seq}\t{qual}\n"
            )
    header = ["@HD\tVN:1.6\tSO:coordinate\n"]
    for chrom, ln in chrom_max_pos.items():
        header.append(f"@SQ\tSN:{chrom}\tLN:{ln}\n")
    return "".join(header + sam_body)

def main():
    if len(sys.argv) != 3:
        print("Usage: python generate_synthetic_bams.py <peaks_csv> <output_bam>", file=sys.stderr)
        sys.exit(1)

    peaks_csv = Path(sys.argv[1]).resolve()
    output_bam = Path(sys.argv[2]).resolve()

    if not peaks_csv.exists():
        sys.exit(f"Peaks file {peaks_csv} not found")

    if output_bam.exists():
        print(f"{output_bam} already exists – skipping generation", file=sys.stderr)
        sys.exit(0)

    sam_text = build_sam(peaks_csv)

    # Pipe SAM to samtools view -> BAM
    samtools_view = subprocess.Popen([
        "samtools", "view", "-bS", "-o", str(output_bam), "-"  # read SAM from stdin
    ], stdin=subprocess.PIPE)
    samtools_view.communicate(input=sam_text.encode())
    samtools_view.wait()
    if samtools_view.returncode != 0:
        sys.exit("samtools view failed")

    # Create BAI index
    subprocess.check_call(["samtools", "index", str(output_bam)])
    print(f"Generated {output_bam} (+ .bai)")

if __name__ == "__main__":
    main() 