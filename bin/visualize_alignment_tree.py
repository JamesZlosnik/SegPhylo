#!/usr/bin/env python3
"""
Phylogenetic tree + barcode-style visualization of the concatenated alignment.

Left panel  : phylogenetic tree drawn with branch lengths, tips top-to-bottom
Right panel : barcode alignment, rows matched exactly to the tree leaf order

- OUTGROUP row: colored by its own nucleotide sequence (reference view)
- Sample rows:  white where identical to OUTGROUP, colored where different
- Dashed vertical lines mark L / M / S segment boundaries
- Shaded bands mark masked terminal regions within each segment
- Columns are auto-binned when the alignment exceeds 5000 positions
- Optional --highlight_samples file colours those labels in red and marks
  their tree tips with a filled circle

Works with trees that include or exclude the OUTGROUP sequence.

Outputs: alignment_tree.png, alignment_tree.svg
"""

import argparse
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap
from Bio import AlignIO, Phylo

HIGHLIGHT_COLOR = "#E53935"   # red for highlighted sample labels / tips


# ─────────────────────────────────────────────────────────────────────────────
# Argument parsing
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment",          required=True,
                        help="Concatenated aligned FASTA (always includes OUTGROUP)")
    parser.add_argument("--partition",          required=True,
                        help="RAxML-style partition file")
    parser.add_argument("--tree",               required=True,
                        help="Rooted concatenated tree in Newick format")
    parser.add_argument("--mask_from_start",    type=int, default=30,
                        help="Positions masked from start of each segment [default: 30]")
    parser.add_argument("--mask_from_end",      type=int, default=50,
                        help="Positions masked from end of each segment   [default: 50]")
    parser.add_argument("--label_taxa_max",     type=int, default=250,
                        help="Hide sample labels when n_taxa exceeds this [default: 250]")
    parser.add_argument("--highlight_samples",  default=None,
                        help="Plain-text file, one sample name per line; those labels "
                             "are highlighted in red on the alignment and marked in the "
                             "tree with a filled circle.")
    return parser.parse_args()


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def load_partitions(partition_file):
    """Parse a RAxML-style partition file into {segment: (start_0idx, end_0idx)}."""
    partitions = {}
    with open(partition_file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            seg_part       = line.split(",")[1]
            seg            = seg_part.split("=")[0].strip()
            start_s, end_s = seg_part.split("=")[1].strip().split("-")
            partitions[seg] = (int(start_s) - 1, int(end_s) - 1)
    return partitions


def encode_seq(seq):
    """Encode a nucleotide string to uint8 array. 0 reserved for 'same as outgroup'."""
    nuc_map = {"A": 1, "T": 2, "G": 3, "C": 4, "-": 5}
    return np.array([nuc_map.get(c.upper(), 6) for c in seq], dtype=np.uint8)


def load_highlights(path):
    """Return a set of sample names to highlight (skips blanks and #-comments)."""
    names = set()
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                names.add(line)
    return names


# ─────────────────────────────────────────────────────────────────────────────
# Tree geometry
# ─────────────────────────────────────────────────────────────────────────────

def compute_node_positions(tree, row_map):
    """
    Compute {clade: (x, y)} for every node.

    x = cumulative branch length from the root (root = 0).
    y = for terminals: their row index in row_map;
        for internal nodes: mean y of their direct children.

    Both x and y are in data coordinates that map directly to the
    alignment display matrix.
    """
    pos = {}

    def _set_x(clade, x_parent):
        x = x_parent + (clade.branch_length or 0.0)
        pos[id(clade)] = {"x": x, "y": None, "clade": clade}
        for child in clade.clades:
            _set_x(child, x)

    _set_x(tree.root, 0.0)

    def _set_y(clade):
        if clade.is_terminal():
            y = row_map.get(clade.name)
            if y is None:
                # Leaf present in tree but absent from alignment — skip gracefully
                y = 0
        else:
            child_ys = [_set_y(c) for c in clade.clades]
            y = sum(child_ys) / len(child_ys)
        pos[id(clade)]["y"] = y
        return y

    _set_y(tree.root)
    return pos


def draw_tree(ax, tree, pos, highlight, n_rows):
    """
    Draw the phylogenetic tree on ax.

    Branch lengths are normalised so the deepest tip sits at x = 1.
    Root is at x = 0 (left); tips at x ≈ 1 (right, bordering the alignment).
    Y axis: row 0 at visual top, row n_rows-1 at visual bottom.
    """
    max_x = max(p["x"] for p in pos.values()) or 1.0

    for clade in tree.find_clades(order="level"):
        p  = pos[id(clade)]
        xc = p["x"] / max_x
        yc = p["y"]

        if clade.clades:
            child_ys = [pos[id(ch)]["y"] for ch in clade.clades]
            # Vertical clade bar connecting children
            ax.plot([xc, xc], [min(child_ys), max(child_ys)],
                    color="#555555", linewidth=0.7, solid_capstyle="butt", zorder=2)
            # Horizontal branch to each child
            for child in clade.clades:
                pc = pos[id(child)]
                xch = pc["x"] / max_x
                ych = pc["y"]
                ax.plot([xc, xch], [ych, ych],
                        color="#555555", linewidth=0.7, zorder=2)
        else:
            # Terminal: draw a marker
            is_outgroup = clade.name == "OUTGROUP"
            is_highlight = clade.name in highlight
            color  = HIGHLIGHT_COLOR if is_highlight else ("#1565C0" if is_outgroup else "#777777")
            marker = "D" if is_outgroup else ("*" if is_highlight else "o")
            size   = 20  if is_highlight else (14 if is_outgroup else 6)
            ax.scatter([xc], [yc], s=size, c=color, marker=marker,
                       zorder=4, linewidths=0)

    # Axes styling
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(n_rows - 0.5, -0.5)   # top-to-bottom: row 0 at top
    ax.axis("off")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    # ── Highlights ────────────────────────────────────────────────────────────
    highlight = load_highlights(args.highlight_samples) if args.highlight_samples else set()

    # ── Alignment ─────────────────────────────────────────────────────────────
    aln_records = list(AlignIO.read(args.alignment, "fasta"))
    seq_dict    = {rec.id: rec for rec in aln_records}

    if "OUTGROUP" not in seq_dict:
        print("ERROR: OUTGROUP sequence not found in alignment.", file=sys.stderr)
        sys.exit(1)

    # ── Tree ──────────────────────────────────────────────────────────────────
    tree          = Phylo.read(args.tree, "newick")
    leaf_order    = [c.name for c in tree.get_terminals()]
    # Handle trees that may or may not include OUTGROUP
    in_tree       = set(leaf_order)
    outgroup_in_tree = "OUTGROUP" in in_tree

    if outgroup_in_tree:
        # Use the tree's natural leaf order for all rows (OUTGROUP wherever it lands)
        all_rows = leaf_order
    else:
        # OUTGROUP is the alignment reference but not in the tree
        # Prepend it as a fixed top row with no tree branch
        all_rows = ["OUTGROUP"] + leaf_order

    # Warn about any alignment samples absent from the tree
    aln_samples  = set(seq_dict.keys()) - {"OUTGROUP"}
    tree_samples = in_tree - {"OUTGROUP"}
    missing = aln_samples - tree_samples
    if missing:
        print(
            f"WARNING: {len(missing)} sample(s) present in the alignment but absent "
            f"from the tree will appear at the bottom of the plot: "
            f"{', '.join(sorted(missing)[:5])}{'...' if len(missing) > 5 else ''}",
            file=sys.stderr,
        )
        all_rows = all_rows + [s for s in sorted(missing)]

    # Samples in tree but absent from alignment (can't draw their alignment row)
    extra = tree_samples - aln_samples
    if extra:
        print(
            f"WARNING: {len(extra)} sample(s) in tree but absent from alignment "
            f"(ignored): {', '.join(sorted(extra)[:5])}",
            file=sys.stderr,
        )

    n_rows  = len(all_rows)
    aln_len = len(str(seq_dict["OUTGROUP"].seq))

    # row_map maps each leaf name to its display row index
    row_map = {name: i for i, name in enumerate(all_rows)}

    # ── Display matrix ─────────────────────────────────────────────────────────
    outgroup_enc = encode_seq(str(seq_dict["OUTGROUP"].seq))
    display = np.zeros((n_rows, aln_len), dtype=np.uint8)

    for i, name in enumerate(all_rows):
        rec = seq_dict.get(name)
        if rec is None:
            continue  # sample not in alignment — leave as zeros (white)
        if name == "OUTGROUP":
            display[i] = outgroup_enc
        else:
            enc  = encode_seq(str(rec.seq))
            diff = enc != outgroup_enc
            display[i][diff]  = enc[diff]
            display[i][~diff] = 0

    # ── Bin wide alignments ───────────────────────────────────────────────────
    max_cols = 5000
    if aln_len > max_cols:
        bin_size = aln_len // max_cols
        new_len  = aln_len // bin_size
        binned   = np.zeros((n_rows, new_len), dtype=np.uint8)
        for j in range(new_len):
            window = display[:, j * bin_size:(j + 1) * bin_size]
            for i in range(n_rows):
                nonzero = window[i][window[i] > 0]
                if len(nonzero):
                    vals, counts = np.unique(nonzero, return_counts=True)
                    binned[i, j] = vals[np.argmax(counts)]
        display  = binned
        scale    = bin_size
        bin_note = f" (1 px = {bin_size} bp)"
    else:
        scale    = 1
        bin_note = ""

    # ── Partitions ────────────────────────────────────────────────────────────
    partitions = load_partitions(args.partition)

    # ── Colormap ──────────────────────────────────────────────────────────────
    colors = ["#FFFFFF", "#4CAF50", "#F44336", "#FF9800", "#2196F3", "#424242", "#BDBDBD"]
    cmap   = ListedColormap(colors)

    # ── Adaptive figure sizing ────────────────────────────────────────────────
    if n_rows <= 100:
        row_h       = 0.20
        show_labels = True
        label_fs    = max(6, min(9, 120 // n_rows))
    elif n_rows <= 250:
        row_h       = 0.14
        show_labels = True
        label_fs    = max(4, min(6, 120 // n_rows))
    else:
        row_h       = 0.08
        show_labels = n_rows <= args.label_taxa_max
        label_fs    = 3

    fig_height = max(5, n_rows * row_h + 3.0)

    # ── Figure layout: tree (left) + alignment (right) ────────────────────────
    fig = plt.figure(figsize=(30, fig_height))
    gs  = gridspec.GridSpec(
        1, 2,
        width_ratios=[1, 3],
        wspace=0.01,
        left=0.01,
        right=0.99,
        top=0.92,
        bottom=0.08,
    )
    ax_tree = fig.add_subplot(gs[0])
    ax_aln  = fig.add_subplot(gs[1])

    # ── Draw tree ─────────────────────────────────────────────────────────────
    pos = compute_node_positions(tree, row_map)
    draw_tree(ax_tree, tree, pos, highlight, n_rows)

    # Label the tree panel
    ax_tree.set_title("Concatenated tree", fontsize=9, pad=4)

    # ── Draw alignment ────────────────────────────────────────────────────────
    ax_aln.imshow(display, aspect="auto", cmap=cmap, vmin=0, vmax=6,
                  interpolation="none", origin="upper")

    # Separator after OUTGROUP row (wherever it is)
    og_row = row_map["OUTGROUP"]
    ax_aln.axhline(y=og_row + 0.5, color="black", linewidth=1.5, zorder=5)

    # Masked region overlays
    for seg, (seg_start, seg_end) in sorted(partitions.items(), key=lambda x: x[1][0]):
        seg_len        = seg_end - seg_start + 1
        mask_s         = min(args.mask_from_start, seg_len)
        x0 = seg_start / scale - 0.5
        x1 = (seg_start + mask_s) / scale - 0.5
        ax_aln.axvspan(x0, x1, color="black", alpha=0.20, zorder=3, linewidth=0)
        mask_e = max(seg_start, seg_end - args.mask_from_end + 1)
        x0 = mask_e / scale - 0.5
        x1 = (seg_end + 1) / scale - 0.5
        ax_aln.axvspan(x0, x1, color="black", alpha=0.20, zorder=3, linewidth=0)

    # Segment boundary dashed lines and labels
    for seg, (start, end) in sorted(partitions.items(), key=lambda x: x[1][0]):
        x_start = start / scale
        x_end   = end   / scale
        ax_aln.axvline(x=x_start, color="black", linewidth=1.2,
                       linestyle="--", alpha=0.8, zorder=4)
        ax_aln.text(
            x_start + (x_end - x_start) / 2, -0.9, seg,
            ha="center", va="top", fontsize=10, fontweight="bold",
            transform=ax_aln.get_xaxis_transform(),
        )

    # Y-axis labels
    if show_labels:
        ax_aln.set_yticks(range(n_rows))
        ax_aln.set_yticklabels(all_rows, fontsize=label_fs)
        for tick, name in zip(ax_aln.get_yticklabels(), all_rows):
            if name == "OUTGROUP":
                tick.set_fontweight("bold")
                tick.set_color("#1565C0")
            elif name in highlight:
                tick.set_color(HIGHLIGHT_COLOR)
                tick.set_fontweight("bold")
    else:
        ax_aln.set_yticks([og_row])
        ax_aln.set_yticklabels(["OUTGROUP"], fontsize=6, fontweight="bold")
        ax_aln.tick_params(axis="y", length=0)
        n_samples = n_rows - 1
        ax_aln.text(
            -0.01, 0.5,
            f"{n_samples} samples (labels hidden — use alignment FASTA for IDs)",
            transform=ax_aln.transAxes, fontsize=7, va="center",
            ha="right", rotation=90, color="grey",
        )

    # X-axis ticks
    tick_pos    = np.linspace(0, display.shape[1] - 1, 10, dtype=int)
    tick_labels = [str(int(tp * scale)) for tp in tick_pos]
    ax_aln.set_xticks(tick_pos)
    ax_aln.set_xticklabels(tick_labels, fontsize=7, rotation=45, ha="right")
    ax_aln.set_xlabel(f"Alignment position{bin_note}", fontsize=9)

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor="#FFFFFF", edgecolor="grey", label="Same as OUTGROUP"),
        mpatches.Patch(facecolor="#4CAF50", label="A"),
        mpatches.Patch(facecolor="#F44336", label="T"),
        mpatches.Patch(facecolor="#FF9800", label="G"),
        mpatches.Patch(facecolor="#2196F3", label="C"),
        mpatches.Patch(facecolor="#424242", label="Gap (-)"),
        mpatches.Patch(facecolor="#BDBDBD", label="N / ambiguous"),
        mpatches.Patch(facecolor="black",   alpha=0.20,
                       label=f"Masked (±{args.mask_from_start}/{args.mask_from_end} bp)"),
    ]
    if highlight:
        legend_elements.append(
            mpatches.Patch(facecolor=HIGHLIGHT_COLOR, label="Highlighted sample")
        )
    ax_aln.legend(
        handles=legend_elements, loc="upper right",
        bbox_to_anchor=(1.0, 1.02), fontsize=7,
        framealpha=0.9, ncol=len(legend_elements),
    )

    fig.suptitle(
        "Concatenated alignment — differences relative to OUTGROUP\n"
        "Tree tip order matches alignment rows  |  dashed lines = segment boundaries  |  "
        "shaded bands = masked regions",
        fontsize=10,
    )

    plt.savefig("alignment_tree.png", dpi=150, bbox_inches="tight")
    plt.savefig("alignment_tree.svg", bbox_inches="tight")
    plt.close()

    n_samples = sum(1 for r in all_rows if r != "OUTGROUP")
    hl_note   = f"  ({len(highlight)} highlighted)" if highlight else ""
    print(
        f"Saved alignment_tree.png and alignment_tree.svg  "
        f"({n_samples} samples, {aln_len} sites{bin_note}{hl_note})"
    )


if __name__ == "__main__":
    main()
