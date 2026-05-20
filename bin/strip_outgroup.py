#!/usr/bin/env python3
"""
Remove the OUTGROUP sequence from an aligned fasta file.
Used when --remove_reference is set, mirroring augur align --remove-reference.

Output: the input alignment minus any sequence with ID 'OUTGROUP'
"""

import argparse
import sys
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment", required=True, help="Input aligned fasta")
    parser.add_argument("--output",    required=True, help="Output fasta (without OUTGROUP)")
    return parser.parse_args()


def main():
    args = parse_args()

    records = [r for r in SeqIO.parse(args.alignment, "fasta") if r.id != "OUTGROUP"]

    if not records:
        print("ERROR: no sequences remain after removing OUTGROUP.", file=sys.stderr)
        sys.exit(1)

    SeqIO.write(records, args.output, "fasta")
    print(f"Removed OUTGROUP — wrote {len(records)} sequences to {args.output}")


if __name__ == "__main__":
    main()
