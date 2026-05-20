# Segmented Virus Phylogenetics Pipeline

A Nextflow (DSL2) pipeline for constructing phylogenetic trees from a
tri-segmented virus (L, M, S segments). Builds both per-segment trees and a
single concatenated tree using a partitioned substitution model.

## Pipeline overview

```
Input multifastas (L, M, S)  +  root sequences (L, M, S)
              │
              ▼
   FILTER_COMPLETE_SAMPLES
   Keep only samples present in all three segments.
   Writes filtering_report.txt.
              │
              ▼
      ADD_OUTGROUP (×3)
   Prepend root sequence (labelled "OUTGROUP") to each segment.
              │
              ▼
      MAFFT_ALIGN (×3)
   Align each segment independently.
              │
         ┌────┴────┐
         ▼         ▼
  IQTREE_SEGMENT  CONCATENATE_ALIGNMENTS
  Per-segment     Concatenate L+M+S alignments;
  trees (×3)      generate RAxML partition file.
                        │
                        ▼
               IQTREE_CONCATENATED
               Full-genome tree with
               per-segment GTR+G models.
```

## Requirements

| Tool       | Minimum version | Notes                        |
|------------|-----------------|------------------------------|
| Nextflow   | ≥ 23.04         | DSL2 required                |
| Python     | ≥ 3.8           | + BioPython                  |
| MAFFT      | ≥ 7.0           |                              |
| IQ-TREE2   | ≥ 2.2           | `-B` ultrafast bootstrap     |

All tools can be managed automatically via the `conda` profile (see below).

## Input format

**Segment multifastas** — headers must follow `>SampleName|SEGMENT`:
```
>StrainA|L
ATGCATGC...
>StrainB|L
ATGCATGC...
```

**Root sequences** — one file per segment, FASTA or GenBank format (any header/locus name; relabelled `OUTGROUP` internally):
```
>ReferenceStrain
ATGCATGC...
```
Accepted extensions: `.fasta` `.fa` `.fna` `.fas` `.gb` `.gbk` `.genbank`
If a GenBank file contains multiple records, the first is used.

## Usage

```bash
nextflow run main.nf \
    --l_fasta  sequences_L.fasta \
    --m_fasta  sequences_M.fasta \
    --s_fasta  sequences_S.fasta \
    --root_l   root_L.fasta \
    --root_m   root_M.fasta \
    --root_s   root_S.fasta \
    --outdir   results
```

### With profiles

```bash
# Auto-install tools via conda
nextflow run main.nf -profile conda [...]

# Run on a SLURM cluster
nextflow run main.nf -profile slurm [...]

# Combine profiles
nextflow run main.nf -profile slurm,conda [...]
```

### Override key parameters

```bash
# Use a stricter MAFFT algorithm
--mafft_args "--localpair --maxiterate 1000"

# Change substitution model for per-segment and concatenated trees
--iqtree_model "TVM+F+G4"

# Build per-segment trees from every sample that passes that segment's threshold
--segment_trees_all_passing true

# More bootstrap replicates
--iqtree_boot 5000
```

> **Note:** The *concatenated* tree uses a partition model with one partition per
> segment. `--iqtree_model` sets the substitution model used for each partition.

## Output structure

```
results/
├── 01_filtered/
│   ├── filtered_L.fasta          # Complete-trio, sorted sample sequences
│   ├── filtered_M.fasta
│   ├── filtered_S.fasta
│   ├── segment_filtered_L.fasta  # All L samples passing the L length threshold
│   ├── segment_filtered_M.fasta
│   ├── segment_filtered_S.fasta
│   └── filtering_report.txt      # Retained / discarded sample summary
│
├── 02_alignments/
│   ├── L_aligned.fasta           # Per-segment MAFFT alignments (incl. OUTGROUP)
│   ├── M_aligned.fasta
│   ├── S_aligned.fasta
│   ├── concatenated_aligned.fasta
│   └── partitions.txt            # RAxML-style partition file (L, M, S boundaries)
│
├── 03_trees/
│   ├── concatenated.treefile     ← Primary output: full-genome rooted newick
│   ├── concatenated.*            # All IQ-TREE2 output files (log, iqtree, etc.)
│   └── per_segment/
│       ├── L.treefile
│       ├── M.treefile
│       ├── S.treefile
│       └── *.{log,iqtree,...}
│
├── 04_visualization/
│   ├── alignment_barcode.png     ← Barcode plot (differences vs OUTGROUP)
│   └── alignment_barcode.svg
│
├── pipeline_report.html
├── pipeline_timeline.html
└── pipeline_trace.txt
```

## Notes

- Samples with any missing segment are excluded from the concatenated tree and
  listed in `filtering_report.txt`.
- By default, per-segment trees use the same complete-trio sample set as the
  concatenated tree. Use `--segment_trees_all_passing true` to build each
  per-segment tree from every sample that passes that segment's length threshold.
- The `OUTGROUP` sequence appears at the root in all output trees.
- The partitioned model in the concatenated tree (`partitions.txt`) assigns an
  independent model to each segment using `--iqtree_model`.
- Resume a failed run with `nextflow run main.nf -resume [...]` — Nextflow
  caches completed processes automatically.
