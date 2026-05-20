#!/usr/bin/env python3
"""
Mask terminal positions in an aligned fasta by replacing nucleotides with 'N'.
Gap characters (-) at masked positions are left as-is.

Output: {segment}_masked.fasta
"""

import argparse
from Bio import SeqIO
from Bio.Seq import Seq


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment",       required=True,
                        help="Aligned fasta file")
    parser.add_argument("--segment",         required=True,
                        help="Segment name (L, M, or S) — used for output filename")
    parser.add_argument("--mask_from_start", type=int, default=30,
                        help="Positions to mask from the start of the alignment [default: 30]")
    parser.add_argument("--mask_from_end",   type=int, default=50,
                        help="Positions to mask from the end of the alignment  [default: 50]")
    return parser.parse_args()


def mask_sequence(seq_str, mask_start, mask_end):
    seq = list(seq_str)
    aln_len = len(seq)
    for i in range(min(mask_start, aln_len)):
        if seq[i] != "-":
            seq[i] = "N"
    for i in range(max(0, aln_len - mask_end), aln_len):
        if seq[i] != "-":
            seq[i] = "N"
    return "".join(seq)


def main():
    args = parse_args()

    records = list(SeqIO.parse(args.alignment, "fasta"))
    if not records:
        raise ValueError(f"No sequences found in {args.alignment}")

    for rec in records:
        rec.seq = Seq(mask_sequence(str(rec.seq), args.mask_from_start, args.mask_from_end))

    out_fname = f"{args.segment}_masked.fasta"
    SeqIO.write(records, out_fname, "fasta")
    print(
        f"Masked {args.mask_from_start} bp from start and "
        f"{args.mask_from_end} bp from end → {out_fname} "
        f"({len(records)} sequences)"
    )


if __name__ == "__main__":
    main()
