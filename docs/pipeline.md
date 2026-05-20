# Pipeline Steps

## Overview

```
Input multifastas (L, M, S)  +  root sequences (L, M, S)
                │
                ▼
   ┌─────────────────────────────┐
   │  FILTER_COMPLETE_SAMPLES    │  Remove dropped strains, apply min length
   │  + dropped strains filter   │  filter, emit complete + per-segment sets
   └─────────────────────────────┘
                │
        ┌───────┴───────┐
        │  (×3 parallel)│
        ▼               │
   ┌──────────────┐     │
   │  MAFFT_ALIGN │     │  Reference-guided alignment per segment
   └──────────────┘     │
        │               │
        ▼               │
   ┌──────────────────┐ │
   │  MASK_ALIGNMENT  │ │  Mask noisy terminal positions
   └──────────────────┘ │
        │               │
        ├───────────────┘
        │
        │  [if --remove_reference]
        ├──────────────────────────────┐
        │                             │
        ▼                             ▼
   ┌──────────────────┐    ┌──────────────────────┐
   │  STRIP_OUTGROUP  │    │  (keep full alignment │
   │  (×3 parallel)   │    │   for visualization)  │
   └──────────────────┘    └──────────────────────┘
        │                             │
        ▼                             ▼
   ┌──────────────────┐    ┌──────────────────────────┐
   │  IQTREE_SEGMENT  │    │  CONCATENATE_ALIGNMENTS  │
   │  (×3 parallel)   │    │  (always with OUTGROUP)  │
   └──────────────────┘    └──────────────────────────┘
        │                             │
        │  [if --root midpoint]       ├───────────────────────┐
        ▼                             │                       │
   ┌──────────────────┐               │  [if --remove_ref]    │
   │  MIDPOINT_ROOT   │    ┌──────────▼──────────┐           │
   │  (×3 parallel)   │    │ STRIP_OUTGROUP_CONCAT│           │
   └──────────────────┘    └──────────────────────┘           │
        │                             │                       │
        ▼                             ▼                       ▼
   Per-segment .nwk       ┌─────────────────────┐  ┌──────────────────────┐
                          │ IQTREE_CONCATENATED │  │  VISUALIZE_ALIGNMENT │
                          └─────────────────────┘  └──────────────────────┘
                                     │                       │
                          [if --root midpoint]               ▼
                                     ▼             alignment_barcode.png/svg
                          ┌──────────────────────┐
                          │  MIDPOINT_ROOT_CONCAT │
                          └──────────────────────┘
                                     │
                                     ▼
                             concatenated.nwk
```

---

## Step 1 — FILTER_COMPLETE_SAMPLES

**Script:** `bin/filter_complete_samples.py`

Parses headers from each segment multifasta in the format `>SampleName|SEGMENT`. Performs filtering in this order:

1. **Dropped strains** — if `--dropped_strains` is provided, named samples are removed from all segments before any other filtering. Lines beginning with `#` in the file are treated as comments.
2. **Minimum length filter** — sequences shorter than the per-segment threshold are removed (length is measured on the raw sequence with gaps stripped).
3. **Complete trio filter** — only samples present in all three segments after the above filters are retained.

Produces:

- `filtered_L.fasta`, `filtered_M.fasta`, `filtered_S.fasta`: complete-trio samples only
- `segment_filtered_L.fasta`, `segment_filtered_M.fasta`, `segment_filtered_S.fasta`: all samples passing that segment's threshold
- `filtering_report.txt`: counts at each stage, per-sample failure reasons, retained complete-trio samples, and per-segment retained counts

By default, the complete-trio FASTAs feed both per-segment and concatenated trees. When `--segment_trees_all_passing true` is set, the per-segment FASTAs feed the per-segment tree path, but concatenation still retains only samples present in all three segments.

**Runs once.** All three segment fastas are processed together so the intersection can be computed.

---

## Step 2 — MAFFT_ALIGN

**Script:** `bin/add_outgroup.py` (via `--outgroup_only`), then `mafft`

Performs **reference-guided alignment** for each segment independently (three parallel tasks).

The root sequence is first extracted and labelled `OUTGROUP` using `add_outgroup.py --outgroup_only`, producing a single-sequence `outgroup_ref.fasta`. MAFFT then aligns all sample sequences against this reference using `--addfragments --keeplength`:

- `--addfragments` — each sample is aligned pairwise to the reference rather than all sequences being pooled together, avoiding alignment errors propagating across divergent sequences.
- `--keeplength` — the alignment width is fixed to the reference sequence length. Sequences shorter than the reference have their ends filled with gaps (`-`), equivalent to `augur align --fill-gaps`.

The `OUTGROUP` sequence appears first in the output alignment.

> **Why reference-guided?** De novo MAFFT alignment pools all sequences including the outgroup, which can introduce gaps into the reference at positions where any sample has an insertion. Reference-guided alignment keeps the outgroup coordinates stable and is consistent with the approach used by `augur align --reference-sequence`.

---

## Step 3 — MASK_ALIGNMENT

**Script:** `bin/mask_alignment.py`

Masks terminal positions of each per-segment alignment by replacing nucleotides with `N`. Gap characters (`-`) at masked positions are left unchanged.

Default masking (matching [hodcroftlab/andv](https://github.com/hodcroftlab/andv)):
- **30 bp** from the start of each segment alignment
- **50 bp** from the end of each segment alignment

These regions correspond to primer binding sites and areas of lower sequencing coverage common in viral amplicon sequencing, which can distort tree topology if left unmasked.

Runs in parallel for L, M, and S.

---

## Step 4 — STRIP_OUTGROUP *(conditional)*

**Script:** `bin/strip_outgroup.py`

*Only runs when `--remove_reference true`.*

Removes the `OUTGROUP` sequence from each masked alignment before tree building. Mirrors the behaviour of `augur align --remove-reference`. Runs in parallel for L, M, and S.

When this step is active, the OUTGROUP is used solely as a coordinate reference for alignment and does not appear in the output trees.

> **Note:** The concatenated alignment passed to `VISUALIZE_ALIGNMENT` always retains the OUTGROUP row regardless of this setting, so the barcode plot always has a reference to compare against.

---

## Step 5 — IQTREE_SEGMENT

**Tool:** `iqtree2`

Builds a maximum-likelihood phylogenetic tree for each segment independently. Runs in parallel for L, M, and S.

Key flags used:

| Flag | Effect |
|------|--------|
| `-o OUTGROUP` | Root on the outgroup sequence (only set when `--root outgroup` and `--remove_reference false`) |
| `-m GTR+G` | Substitution model (configurable via `--iqtree_model`) |
| `-B 1000` | Ultrafast bootstrap replicates (configurable via `--iqtree_boot`) |
| `-czb` | Collapse zero-length branches into polytomies |
| `-st DNA` | Explicitly declare sequence type as DNA |
| `--redo` | Overwrite any existing IQ-TREE2 output |

Output trees are renamed from `.treefile` to `.nwk`.

---

## Step 6 — MIDPOINT_ROOT *(conditional)*

**Script:** `bin/midpoint_root.py`

*Only runs when `--root midpoint`.*

Applies midpoint rooting to each per-segment tree using Biopython's `Phylo.root_at_midpoint()`. Mirrors the behaviour of `augur refine --root mid_point`. Runs in parallel for L, M, and S.

The midpoint-rooted `.nwk` files overwrite the IQ-TREE2 output at the same paths.

---

## Step 7 — CONCATENATE_ALIGNMENTS

**Script:** `bin/concatenate_alignments.py`

Concatenates the three masked, aligned segment fastas (always the full alignments **with** OUTGROUP) in **L → M → S** order. Also produces an IQ-TREE2-compatible RAxML-style partition file (`partitions.txt`) recording the site boundaries of each segment.

The script:
- Normalises sequence IDs by stripping the `|SEGMENT` suffix so samples match across files
- Retains only taxa present in all three segment alignments for concatenation
- Reports taxa missing from individual segments and excludes them from the concatenated alignment
- Validates that all sequences within each alignment are the same length
- Prints an input summary and missing-data report to the log

The partition file format:
```
GTR+G, L = 1-6420
GTR+G, M = 6421-9561
GTR+G, S = 9562-11470
```

The model string comes from `--iqtree_model`; for example, `--iqtree_model GTR` writes `GTR, L = ...` for each partition.

---

## Step 8 — STRIP_OUTGROUP_CONCAT *(conditional)*

**Script:** `bin/strip_outgroup.py`

*Only runs when `--remove_reference true`.*

Removes the `OUTGROUP` sequence from the concatenated alignment before building the full-genome tree. The visualization step still receives the full concatenated alignment (with OUTGROUP).

---

## Step 9 — IQTREE_CONCATENATED

**Tool:** `iqtree2`

Builds the primary full-genome phylogenetic tree from the concatenated alignment using a **partitioned substitution model** — one independent model is fitted per segment rather than a single model across the whole alignment. The model for each partition is set by `--iqtree_model` and defaults to `GTR+G`.

Uses the same `-czb -st DNA` flags and conditional `-o OUTGROUP` logic as the per-segment trees.

Output tree is written as `concatenated.nwk`.

---

## Step 10 — MIDPOINT_ROOT_CONCAT *(conditional)*

**Script:** `bin/midpoint_root.py`

*Only runs when `--root midpoint`.*

Applies midpoint rooting to the concatenated tree, overwriting `concatenated.nwk`.

---

## Step 11 — VISUALIZE_ALIGNMENT

**Script:** `bin/visualize_alignment.py`

Produces a barcode-style visualisation of the concatenated alignment (always with OUTGROUP as the reference row) for manual inspection.

Coloring scheme:

| Color | Meaning |
|-------|---------|
| White | Identical to OUTGROUP at this position |
| Green | A (differs from OUTGROUP) |
| Red | T (differs from OUTGROUP) |
| Orange | G (differs from OUTGROUP) |
| Blue | C (differs from OUTGROUP) |
| Dark grey | Gap (`-`) |
| Light grey | N / ambiguous |
| Semi-transparent black band | Masked terminal region |

Additional features:
- The **OUTGROUP row** (top) is colored by its own nucleotides so the reference can be read directly
- **Dashed vertical lines** mark L / M / S segment boundaries
- **Shaded bands** indicate masked terminal regions within each segment
- Row height and label visibility adapt automatically to the number of taxa (labels suppressed above `--label_taxa_max`)
- Columns are auto-binned when the alignment exceeds 5000 positions

Output is written as both PNG (150 dpi) and SVG.

---

## Comparison with hodcroftlab/andv

The hodcroftlab/andv Nextstrain build uses Snakemake + augur and targets Nextstrain visualisation. The table below shows how this pipeline maps onto their approach.

| Aspect | hodcroftlab/andv | This pipeline (default) | Closest equivalent |
|--------|-----------------|------------------------|--------------------|
| Alignment | `augur align --reference-sequence --fill-gaps` | MAFFT `--addfragments --keeplength` | Equivalent |
| Reference in tree | Removed (`--remove-reference`) | Kept as OUTGROUP | `--remove_reference true` |
| Tree rooting | Midpoint (`augur refine --root mid_point`) | Outgroup | `--root midpoint` |
| Terminal masking | 30 bp / 50 bp (`augur mask`) | 30 bp / 50 bp | Equivalent (defaults match) |
| Substitution model | GTR (augur/IQ-TREE default) | GTR+G for per-segment and concatenated partitions | `--iqtree_model GTR` |
| Zero-branch collapse | `-czb` | `-czb` | Equivalent |
| Sequence type flag | `-st DNA` | `-st DNA` | Equivalent |
| Dropped strains | `config/dropped_strains.txt` | `--dropped_strains` | Equivalent |
| Min length filter | S=1000, M=3000, L=6000 | S=1000, M=3000, L=6000 | Equivalent (defaults match) |
| Timetree / refine | Yes (augur refine) | Not applicable | — |
| Nextstrain export | Yes | Not applicable | — |

**Command to most closely replicate hodcroftlab/andv tree topology:**
```bash
nextflow run main.nf \
    --remove_reference true \
    --root             midpoint \
    --iqtree_model     GTR \
    --root_l           config/outgroup_L.gb \
    --root_m           config/outgroup_M.gb \
    --root_s           config/outgroup_S.gb \
    [...]
```
