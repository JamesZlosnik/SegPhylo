#!/usr/bin/env python3
"""
Concatenate three per-segment MAFFT alignments (L, M, S) into a single
alignment and write an IQ-TREE2-compatible RAxML-style partition file.

Outputs:
  - concatenated_aligned.fasta
  - partitions.txt
"""

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--l_aln", required=True, help="Aligned segment L fasta")
    parser.add_argument("--m_aln", required=True, help="Aligned segment M fasta")
    parser.add_argument("--s_aln", required=True, help="Aligned segment S fasta")
    return parser.parse_args()


def normalize_id(seq_id):
    """
    Normalize sequence IDs so segments match across files.
    Example:
        PP_006VP2A.1|L -> PP_006VP2A.1
    """
    return seq_id.split("|")[0]


def load_alignment(fname):
    """
    Load FASTA alignment into dict: {normalized_id: sequence}
    """
    seqs = {}
    for rec in SeqIO.parse(fname, "fasta"):
        nid = normalize_id(rec.id)
        seqs[nid] = str(rec.seq)
    return seqs


def main():
    args = parse_args()

    # Load alignments
    alns = {
        "L": load_alignment(args.l_aln),
        "M": load_alignment(args.m_aln),
        "S": load_alignment(args.s_aln),
    }

    # Check lengths
    lengths = {}
    for seg, seqs in alns.items():
        if not seqs:
            raise ValueError(f"No sequences found in {seg} alignment")
        lengths[seg] = len(next(iter(seqs.values())))

    # Check alignment consistency
    for seg, seqs in alns.items():
        bad = [sid for sid, seq in seqs.items() if len(seq) != lengths[seg]]
        if bad:
            raise ValueError(f"Inconsistent alignment length in {seg}: {bad[:5]}")

    # Union of taxa across all segments
    all_ids = set().union(*[set(alns[seg]) for seg in ["L", "M", "S"]])

    # Debug info
    print("=== Input summary ===")
    for seg in ["L", "M", "S"]:
        print(f"{seg}: {len(alns[seg])} sequences (length={lengths[seg]})")
    print(f"Total unique taxa: {len(all_ids)}")

    # Missing data report
    missing = {
        seg: len(all_ids - set(alns[seg]))
        for seg in ["L", "M", "S"]
    }
    print("Missing taxa per segment:", missing)

    # Write partition file
    pos = 1
    with open("partitions.txt", "w") as fh:
        for seg in ["L", "M", "S"]:
            end = pos + lengths[seg] - 1
            fh.write(f"GTR+G, {seg} = {pos}-{end}\n")
            pos = end + 1

    # Concatenate sequences
    records = []
    for sid in sorted(all_ids):
        concat_seq = "".join(
            alns[seg].get(sid, "-" * lengths[seg]) for seg in ["L", "M", "S"]
        )
        records.append(SeqRecord(Seq(concat_seq), id=sid, name="", description=""))

    # Write concatenated alignment
    SeqIO.write(records, "concatenated_aligned.fasta", "fasta")

    total_len = sum(lengths.values())
    print(
        f"\n✅ Concatenated alignment written: "
        f"{len(records)} taxa, {total_len} sites "
        f"(L={lengths['L']}, M={lengths['M']}, S={lengths['S']})"
    )


if __name__ == "__main__":
    main()
