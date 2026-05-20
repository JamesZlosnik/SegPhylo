# Output Files

All outputs are published to `--outdir` (default: `results/`). The directory is created automatically if it does not exist.

---

## Directory structure

```
results/
├── 01_filtered/
│   ├── filtered_L.fasta
│   ├── filtered_M.fasta
│   ├── filtered_S.fasta
│   ├── segment_filtered_L.fasta
│   ├── segment_filtered_M.fasta
│   ├── segment_filtered_S.fasta
│   └── filtering_report.txt
│
├── 02_alignments/
│   ├── L_aligned.fasta
│   ├── M_aligned.fasta
│   ├── S_aligned.fasta
│   ├── L_masked.fasta
│   ├── M_masked.fasta
│   ├── S_masked.fasta
│   ├── concatenated_aligned.fasta
│   └── partitions.txt
│
├── 03_trees/
│   ├── concatenated.nwk          ← primary output
│   ├── concatenated.log
│   ├── concatenated.iqtree
│   ├── concatenated.contree
│   └── per_segment/
│       ├── L.nwk
│       ├── L.log
│       ├── L.iqtree
│       ├── L.contree
│       ├── M.nwk
│       ├── M.log
│       ├── M.iqtree
│       ├── M.contree
│       ├── S.nwk
│       ├── S.log
│       ├── S.iqtree
│       └── S.contree
│
└── 04_visualization/
    ├── alignment_barcode.png
    └── alignment_barcode.svg
```

> **Note on `--remove_reference`:** When `--remove_reference true` is set, the OUTGROUP sequence is absent from the `.nwk` tree files. The alignment files in `02_alignments/` always retain the OUTGROUP row, as they are also used for the barcode visualisation.

> **Note on `--root midpoint`:** The `.nwk` files are post-processed by midpoint rooting after IQ-TREE2 finishes. The `.iqtree`, `.log`, and `.contree` files reflect the raw unrooted IQ-TREE2 output.

---

## 01_filtered/

Output of `FILTER_COMPLETE_SAMPLES`.

### `filtered_L.fasta`, `filtered_M.fasta`, `filtered_S.fasta`

Filtered multifasta files containing only samples that passed all three filtering stages:

1. Not in the `--dropped_strains` exclusion list
2. Meet the minimum length threshold for their segment
3. Present in all three segments (complete trio)

Sequences are sorted alphabetically by sample name. Headers retain the original `SampleName|SEGMENT` format.

### `segment_filtered_L.fasta`, `segment_filtered_M.fasta`, `segment_filtered_S.fasta`

Per-segment filtered FASTAs containing every sample that:

1. Is not in the `--dropped_strains` exclusion list
2. Meets the minimum length threshold for that segment

These files may contain samples that are missing one or both of the other segments. They are used for per-segment tree building only when `--segment_trees_all_passing true` is set. The concatenated tree never uses incomplete samples.

### `filtering_report.txt`

Plain-text summary of the filtering step. Contains:

- Number of explicitly dropped strains (if `--dropped_strains` was provided)
- Minimum length thresholds applied
- Raw sample counts per segment after explicit exclusions
- Samples that failed the length threshold (with observed and minimum lengths)
- Samples discarded for incomplete trios (with which segments were present)
- Per-segment retained counts
- Final count of complete-trio retained samples and full list

Example:
```
=== Segment Filtering Report ===

Minimum length thresholds:
  L: 6000 bp
  M: 3000 bp
  S: 1000 bp

Samples per segment (raw, after explicit exclusions):
  L: 210
  M: 216
  S: 218

Dropped by explicit exclusion list: 2
  SampleDropped1
  SampleDropped2

Failed minimum length filter:
  SampleX [L]: 4821 bp (min 6000 bp)

Retained (complete trios, length-passing): 207
Discarded (incomplete or too short):        11

Retained per segment (length-passing):
  L: 209
  M: 213
  S: 215

Discarded samples:
  SampleX: L too short (4821 bp)
  SampleY: M missing

Retained samples:
  SampleA
  SampleB
  ...
```

---

## 02_alignments/

Outputs of `MAFFT_ALIGN`, `MASK_ALIGNMENT`, and `CONCATENATE_ALIGNMENTS`.

### `L_aligned.fasta`, `M_aligned.fasta`, `S_aligned.fasta`

Reference-guided MAFFT alignments for each segment. The `OUTGROUP` sequence appears first. Alignment width is fixed to the reference sequence length (`--keeplength`). Sequences that do not cover the full reference length have their uncovered ends filled with gap characters (`-`).

### `L_masked.fasta`, `M_masked.fasta`, `S_masked.fasta`

Per-segment alignments after terminal masking. The first `--mask_from_start` bp and last `--mask_from_end` bp of each segment alignment have nucleotides replaced with `N`. Gap characters at masked positions are left unchanged.

These are the alignments used for tree building (unless `--remove_reference true`, in which case the OUTGROUP row is additionally stripped before IQ-TREE2 runs — but the published files here always include OUTGROUP).

### `concatenated_aligned.fasta`

The three masked segment alignments concatenated in **L → M → S** order, always including the OUTGROUP row. Sequence IDs have the `|SEGMENT` suffix stripped so they match across segments. Only taxa present in all three segment alignments are retained in the concatenated alignment. Used for both the concatenated tree and the barcode visualisation.

### `partitions.txt`

RAxML-style partition file defining the start and end positions of each segment within the concatenated alignment. Used by IQ-TREE2 to fit an independent substitution model per segment. The model string comes from `--iqtree_model`.

Example:
```
GTR+G, L = 1-6420
GTR+G, M = 6421-9561
GTR+G, S = 9562-11470
```

---

## 03_trees/

Outputs of `IQTREE_SEGMENT`, `IQTREE_CONCATENATED`, and (if applicable) `MIDPOINT_ROOT` / `MIDPOINT_ROOT_CONCAT`.

### `concatenated.nwk` ← primary output

The full-genome phylogenetic tree in Newick format, built from the concatenated alignment using a partitioned model (one model per segment). This is the main deliverable of the pipeline.

- With `--root outgroup` (default): rooted on the `OUTGROUP` sequence.
- With `--root midpoint`: midpoint-rooted using Biopython after IQ-TREE2 finishes.
- With `--remove_reference true`: `OUTGROUP` is absent from the tree.

### `L.nwk`, `M.nwk`, `S.nwk`

Per-segment phylogenetic trees in Newick format, each built independently from its masked alignment. Rooting follows the same logic as `concatenated.nwk`.

By default, these trees use the complete-trio sample set. With `--segment_trees_all_passing true`, each segment tree uses every sample that passes that segment's length threshold, even if other segments are missing.

### `*.iqtree`

IQ-TREE2 summary report for each tree. Contains the substitution model parameters, log-likelihood, AIC/BIC scores, tree topology, and bootstrap support values. Reflects the raw IQ-TREE2 output before any midpoint rooting is applied.

### `*.log`

Full IQ-TREE2 run log. Useful for diagnosing tree-building issues or verifying which model was selected.

### `*.contree`

Consensus tree with bootstrap support values annotated on branches as percentages. Reflects the raw IQ-TREE2 output before any midpoint rooting post-processing.

---

## 04_visualization/

Output of `VISUALIZE_ALIGNMENT`.

### `alignment_barcode.png` / `alignment_barcode.svg`

Barcode-style visualisation of the concatenated alignment. Always includes the OUTGROUP row as a reference regardless of `--remove_reference`. Rows are samples; columns are alignment positions.

| Element | Description |
|---------|-------------|
| OUTGROUP row (top, bold) | Colored by its own nucleotide at each position |
| White cells | Position identical to OUTGROUP |
| Green cells | A — differs from OUTGROUP |
| Red cells | T — differs from OUTGROUP |
| Orange cells | G — differs from OUTGROUP |
| Blue cells | C — differs from OUTGROUP |
| Dark grey cells | Gap (`-`) |
| Light grey cells | N / ambiguous base |
| Dashed vertical lines | Segment boundaries (L / M / S) with segment labels |
| Shaded grey bands | Masked terminal regions (width reflects `--mask_from_start` / `--mask_from_end`) |

**Adaptive sizing:** row height and label visibility scale with the number of taxa:

| Taxa count | Row height | Labels shown |
|------------|------------|--------------|
| ≤ 100 | 0.20 in | Yes (6–9 pt) |
| 101–250 | 0.14 in | Yes (4–6 pt) |
| > 250 | 0.08 in | Only if ≤ `--label_taxa_max` |

For large datasets where labels are hidden, use the row index into `concatenated_aligned.fasta` to identify samples by position.

The PNG is rendered at 150 dpi. The SVG is resolution-independent and recommended for publication or further editing in Inkscape or Illustrator.
