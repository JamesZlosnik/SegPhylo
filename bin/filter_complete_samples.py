#!/usr/bin/env python3
"""
Filter multifasta files for a tri-segmented virus.

The complete-trio outputs retain only samples present in all three segments
(L, M, S) and meeting all segment-specific length thresholds. The per-segment
outputs retain every sample that passes the threshold for that segment.

Headers must follow the format:  >SampleName|SEGMENT
Outputs:
  - filtered_L.fasta, filtered_M.fasta, filtered_S.fasta
  - segment_filtered_L.fasta, segment_filtered_M.fasta, segment_filtered_S.fasta
  - filtering_report.txt
"""

import argparse
import sys
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--l_fasta",    required=True, help="Multifasta for segment L")
    parser.add_argument("--m_fasta",    required=True, help="Multifasta for segment M")
    parser.add_argument("--s_fasta",    required=True, help="Multifasta for segment S")
    parser.add_argument("--min_len_l",  type=int, default=6000,
                        help="Minimum length for segment L [default: 6000]")
    parser.add_argument("--min_len_m",  type=int, default=3000,
                        help="Minimum length for segment M [default: 3000]")
    parser.add_argument("--min_len_s",  type=int, default=1000,
                        help="Minimum length for segment S [default: 1000]")
    parser.add_argument(
        "--dropped_strains",
        help="Optional file of sample names to exclude, one per line. "
             "Blank lines and lines beginning with '#' are ignored.",
    )
    return parser.parse_args()


def load_segment(fpath):
    recs = {}
    for rec in SeqIO.parse(str(fpath), "fasta"):
        parts = rec.id.split("|")
        if len(parts) < 2:
            print(f"WARNING: unexpected header format (no '|'): {rec.id}", file=sys.stderr)
        sample = parts[0]
        recs[sample] = rec
    return recs


def apply_length_filter(recs, min_len, seg):
    """Remove sequences below min_len. Returns (passing, failed_dict)."""
    passing = {}
    failed  = {}
    for sample, rec in recs.items():
        seq_len = len(str(rec.seq).replace("-", ""))
        if seq_len >= min_len:
            passing[sample] = rec
        else:
            failed[sample] = seq_len
    return passing, failed


def load_dropped_strains(fpath):
    if not fpath:
        return set()

    dropped = set()
    with open(fpath) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            dropped.add(line)
    return dropped


def main():
    args = parse_args()

    min_lengths = {"L": args.min_len_l, "M": args.min_len_m, "S": args.min_len_s}
    dropped_strains = load_dropped_strains(args.dropped_strains)

    raw = {
        "L": load_segment(args.l_fasta),
        "M": load_segment(args.m_fasta),
        "S": load_segment(args.s_fasta),
    }

    if dropped_strains:
        for seg in ["L", "M", "S"]:
            for sample in dropped_strains:
                raw[seg].pop(sample, None)

    # Apply per-segment length filters
    segs          = {}
    length_failed = {}
    for seg in ["L", "M", "S"]:
        passing, failed      = apply_length_filter(raw[seg], min_lengths[seg], seg)
        segs[seg]            = passing
        length_failed[seg]   = failed

    complete  = set(segs["L"]) & set(segs["M"]) & set(segs["S"])
    all_samps = set(raw["L"])  | set(raw["M"])  | set(raw["S"])
    discarded = all_samps - complete

    with open("filtering_report.txt", "w") as fh:
        fh.write("=== Segment Filtering Report ===\n\n")

        fh.write("Minimum length thresholds:\n")
        for seg in ["L", "M", "S"]:
            fh.write(f"  {seg}: {min_lengths[seg]} bp\n")

        fh.write("\nSamples per segment (raw):\n")
        for seg in ["L", "M", "S"]:
            fh.write(f"  {seg}: {len(raw[seg])}\n")

        if dropped_strains:
            fh.write(f"\nDropped by explicit exclusion list: {len(dropped_strains)}\n")
            for sample in sorted(dropped_strains):
                fh.write(f"  {sample}\n")

        fh.write("\nFailed minimum length filter:\n")
        any_failed = False
        for seg in ["L", "M", "S"]:
            if length_failed[seg]:
                any_failed = True
                for sample, length in sorted(length_failed[seg].items()):
                    fh.write(f"  {sample} [{seg}]: {length} bp (min {min_lengths[seg]} bp)\n")
        if not any_failed:
            fh.write("  None\n")

        fh.write(f"\nRetained (complete trios, length-passing): {len(complete)}\n")
        fh.write(f"Discarded (incomplete or too short):        {len(discarded)}\n")

        fh.write("\nRetained per segment (length-passing):\n")
        for seg in ["L", "M", "S"]:
            fh.write(f"  {seg}: {len(segs[seg])}\n")

        if discarded:
            fh.write("\nDiscarded samples:\n")
            for s in sorted(discarded):
                reasons = []
                for seg in ["L", "M", "S"]:
                    if s in length_failed[seg]:
                        reasons.append(f"{seg} too short ({length_failed[seg][s]} bp)")
                    elif s not in segs[seg]:
                        reasons.append(f"{seg} missing")
                fh.write(f"  {s}: {'; '.join(reasons)}\n")

        fh.write("\nRetained samples:\n")
        for s in sorted(complete):
            fh.write(f"  {s}\n")

    for seg, out_fname in [
        ("L", "filtered_L.fasta"),
        ("M", "filtered_M.fasta"),
        ("S", "filtered_S.fasta"),
    ]:
        records = [segs[seg][s] for s in sorted(complete)]
        SeqIO.write(records, out_fname, "fasta")

    for seg, out_fname in [
        ("L", "segment_filtered_L.fasta"),
        ("M", "segment_filtered_M.fasta"),
        ("S", "segment_filtered_S.fasta"),
    ]:
        records = [segs[seg][s] for s in sorted(segs[seg])]
        SeqIO.write(records, out_fname, "fasta")

    n_length = sum(len(v) for v in length_failed.values())
    print(
        f"Filtering complete: {len(complete)} retained, "
        f"{len(discarded)} discarded "
        f"({n_length} length failures, "
        f"{len(discarded) - len(set(s for v in length_failed.values() for s in v))} incomplete trios)."
    )


if __name__ == "__main__":
    main()
