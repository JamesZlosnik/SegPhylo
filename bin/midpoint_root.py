#!/usr/bin/env python3
"""
Apply midpoint rooting to a Newick tree file.
Mirrors augur refine --root mid_point behaviour.
"""

import argparse
from Bio import Phylo


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree",   required=True, help="Input Newick tree file")
    parser.add_argument("--output", required=True, help="Output midpoint-rooted Newick file")
    return parser.parse_args()


def main():
    args = parse_args()

    tree = Phylo.read(args.tree, "newick")
    tree.root_at_midpoint()
    Phylo.write(tree, args.output, "newick")
    print(f"Midpoint-rooted tree written to {args.output}")


if __name__ == "__main__":
    main()
