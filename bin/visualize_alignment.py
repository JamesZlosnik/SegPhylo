#!/usr/bin/env python3
"""
Barcode-style visualization of the concatenated alignment.

- OUTGROUP row: colored by its own nucleotide sequence (reference view)
- Sample rows:  white where identical to OUTGROUP, colored where different
- Dashed vertical lines mark L / M / S segment boundaries
- Grey shaded bands mark the masked terminal regions within each segment
- Columns are auto-binned when the alignment exceeds 5000 positions

Outputs: alignment_barcode.png, alignment_barcode.svg
"""

import argparse
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
from Bio import AlignIO


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment",       required=True,
                        help="Concatenated aligned fasta")
    parser.add_argument("--partition",       required=True,
                        help="RAxML-style partition file")
    parser.add_argument("--mask_from_start", type=int, default=30,
                        help="Positions masked from start of each segment [default: 30]")
    parser.add_argument("--mask_from_end",   type=int, default=50,
                        help="Positions masked from end of each segment   [default: 50]")
    parser.add_argument("--label_taxa_max",  type=int, default=250,
                        help="Hide sample labels when n_taxa exceeds this (default: 250)")
    return parser.parse_args()


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


def main():
    args = parse_args()

    # ── Load alignment ────────────────────────────────────────────────────────
    aln_records = list(AlignIO.read(args.alignment, "fasta"))
    outgroup    = None
    samples     = []
    for rec in aln_records:
        if rec.id == "OUTGROUP":
            outgroup = rec
        else:
            samples.append(rec)

    if outgroup is None:
        print("ERROR: OUTGROUP sequence not found in alignment.", file=sys.stderr)
        sys.exit(1)

    aln_len    = len(str(outgroup.seq))
    n_samples  = len(samples)
    sample_ids = [rec.id for rec in samples]

    # ── Load partitions ───────────────────────────────────────────────────────
    partitions = load_partitions(args.partition)

    # ── Build display matrix ──────────────────────────────────────────────────
    outgroup_enc = encode_seq(str(outgroup.seq))
    n_rows       = n_samples + 1
    display      = np.zeros((n_rows, aln_len), dtype=np.uint8)
    display[0]   = outgroup_enc

    for i, rec in enumerate(samples):
        enc              = encode_seq(str(rec.seq))
        diff             = enc != outgroup_enc
        display[i + 1][diff]  = enc[diff]
        display[i + 1][~diff] = 0

    # ── Bin columns for very wide alignments ──────────────────────────────────
    max_cols = 5000
    if aln_len > max_cols:
        bin_size = aln_len // max_cols
        new_len  = aln_len // bin_size
        binned   = np.zeros((n_rows, new_len), dtype=np.uint8)
        for j in range(new_len):
            window = display[:, j * bin_size:(j + 1) * bin_size]
            for i in range(n_rows):
                nonzero = window[i][window[i] > 0]
                if len(nonzero) > 0:
                    vals, counts = np.unique(nonzero, return_counts=True)
                    binned[i, j] = vals[np.argmax(counts)]
        display  = binned
        scale    = bin_size
        bin_note = f" (1 px = {bin_size} bp)"
    else:
        scale    = 1
        bin_note = ""

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
    fig, ax    = plt.subplots(figsize=(22, fig_height))

    ax.imshow(display, aspect="auto", cmap=cmap, vmin=0, vmax=6,
              interpolation="none", origin="upper")

    # Separator between OUTGROUP row and samples
    ax.axhline(y=0.5, color="black", linewidth=1.5, zorder=5)

    # ── Masked region overlays ────────────────────────────────────────────────
    # Draw a semi-transparent grey band over the masked terminal positions
    # within each segment. Positions are converted to display (binned) coords.
    for seg, (seg_start, seg_end) in sorted(partitions.items(), key=lambda x: x[1][0]):
        seg_len = seg_end - seg_start + 1

        # Start mask band
        mask_start_end = min(args.mask_from_start, seg_len)
        x0 = (seg_start) / scale - 0.5
        x1 = (seg_start + mask_start_end) / scale - 0.5
        ax.axvspan(x0, x1, color="black", alpha=0.20, zorder=3, linewidth=0)

        # End mask band
        mask_end_start = max(seg_start, seg_end - args.mask_from_end + 1)
        x0 = (mask_end_start) / scale - 0.5
        x1 = (seg_end + 1) / scale - 0.5
        ax.axvspan(x0, x1, color="black", alpha=0.20, zorder=3, linewidth=0)

    # ── Segment boundary lines and labels ─────────────────────────────────────
    for seg, (start, end) in sorted(partitions.items(), key=lambda x: x[1][0]):
        x_start = start / scale
        x_end   = end   / scale
        ax.axvline(x=x_start, color="black", linewidth=1.2,
                   linestyle="--", alpha=0.8, zorder=4)
        ax.text(
            x_start + (x_end - x_start) / 2, -0.9,
            seg, ha="center", va="top", fontsize=10, fontweight="bold",
            transform=ax.get_xaxis_transform(),
        )

    # ── Y-axis labels ─────────────────────────────────────────────────────────
    if show_labels:
        ax.set_yticks(range(n_rows))
        ax.set_yticklabels(["OUTGROUP"] + sample_ids, fontsize=label_fs)
        ax.get_yticklabels()[0].set_fontweight("bold")
    else:
        ax.set_yticks([0])
        ax.set_yticklabels(["OUTGROUP"], fontsize=6, fontweight="bold")
        ax.tick_params(axis="y", length=0)
        ax.text(-0.01, 0.5,
                f"{n_samples} samples (labels hidden — use alignment FASTA for IDs)",
                transform=ax.transAxes, fontsize=7, va="center",
                ha="right", rotation=90, color="grey")

    # ── X-axis ────────────────────────────────────────────────────────────────
    tick_pos    = np.linspace(0, display.shape[1] - 1, 10, dtype=int)
    tick_labels = [str(int(tp * scale)) for tp in tick_pos]
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_labels, fontsize=7, rotation=45, ha="right")
    ax.set_xlabel(f"Alignment position{bin_note}", fontsize=9)

    # ── Legend ────────────────────────────────────────────────────────────────
    legend_elements = [
        mpatches.Patch(facecolor="#FFFFFF", edgecolor="grey", label="Same as OUTGROUP"),
        mpatches.Patch(facecolor="#4CAF50", label="A"),
        mpatches.Patch(facecolor="#F44336", label="T"),
        mpatches.Patch(facecolor="#FF9800", label="G"),
        mpatches.Patch(facecolor="#2196F3", label="C"),
        mpatches.Patch(facecolor="#424242", label="Gap (-)"),
        mpatches.Patch(facecolor="#BDBDBD", label="N / ambiguous"),
        mpatches.Patch(facecolor="black",   alpha=0.20, label=f"Masked (±{args.mask_from_start}/{args.mask_from_end} bp)"),
    ]
    ax.legend(
        handles=legend_elements, loc="upper right",
        bbox_to_anchor=(1.0, 1.02), fontsize=7,
        framealpha=0.9, ncol=len(legend_elements),
    )

    ax.set_title(
        "Concatenated alignment — differences relative to OUTGROUP\n"
        "OUTGROUP row shows its own sequence; dashed lines = segment boundaries; "
        "shaded bands = masked regions",
        fontsize=10, pad=12,
    )

    plt.tight_layout()
    fig.savefig("alignment_barcode.png", dpi=150, bbox_inches="tight")
    fig.savefig("alignment_barcode.svg", bbox_inches="tight")
    plt.close()
    print(
        f"Saved alignment_barcode.png and alignment_barcode.svg  "
        f"({n_samples} samples, {aln_len} sites{bin_note})"
    )


if __name__ == "__main__":
    main()
