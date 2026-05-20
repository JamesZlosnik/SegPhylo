#!/usr/bin/env python3
"""
Prepare the root/outgroup sequence for use in the pipeline.

Two modes:
  Default          -- prepend OUTGROUP to sample fasta → {segment}_with_outgroup.fasta
  --outgroup_only  -- write OUTGROUP alone             → outgroup_ref.fasta
                      (used as the MAFFT reference for --addfragments)

Supported root formats: FASTA (.fasta .fa .fna .fas) or GenBank (.gb .gbk .genbank)
"""

import argparse
import sys
from pathlib import Path
from Bio import SeqIO


FASTA_EXTS   = {'.fasta', '.fa', '.fna', '.fas'}
GENBANK_EXTS = {'.gb', '.gbk', '.genbank'}


def detect_format(fpath):
    suffix = Path(fpath).suffix.lower()
    if suffix in FASTA_EXTS:
        return 'fasta'
    if suffix in GENBANK_EXTS:
        return 'genbank'
    print(
        f"WARNING: unrecognised extension '{suffix}' for {fpath}, assuming FASTA.",
        file=sys.stderr,
    )
    return 'fasta'


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--segment",      required=True,
                        help="Segment name (L, M, or S)")
    parser.add_argument("--fasta",        required=True,
                        help="Filtered segment multifasta (samples)")
    parser.add_argument("--root_fasta",   required=True,
                        help="Root/outgroup sequence (FASTA or GenBank)")
    parser.add_argument("--outgroup_only", action="store_true",
                        help="Write only the outgroup to outgroup_ref.fasta "
                             "(for use as MAFFT --addfragments reference)")
    return parser.parse_args()


def main():
    args = parse_args()

    fmt  = detect_format(args.root_fasta)
    recs = list(SeqIO.parse(args.root_fasta, fmt))

    if not recs:
        print(f"ERROR: no sequences found in {args.root_fasta}", file=sys.stderr)
        sys.exit(1)
    if len(recs) > 1:
        print(
            f"WARNING: {len(recs)} records found in root file — using the first one.",
            file=sys.stderr,
        )

    root             = recs[0]
    root.id          = "OUTGROUP"
    root.name        = ""
    root.description = ""

    if args.outgroup_only:
        SeqIO.write([root], "outgroup_ref.fasta", "fasta")
        print(f"Written outgroup reference to outgroup_ref.fasta [{fmt}]")
    else:
        samples  = list(SeqIO.parse(args.fasta, "fasta"))
        out_path = f"{args.segment}_with_outgroup.fasta"
        SeqIO.write([root] + samples, out_path, "fasta")
        print(
            f"Written {len(samples) + 1} sequences "
            f"(1 outgroup [{fmt}] + {len(samples)} samples) to {out_path}"
        )


if __name__ == "__main__":
    main()
