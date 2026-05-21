# Usage Guide

This guide is for running the pipeline end to end. For every available option, see [parameters.md](parameters.md).

## 1. Check Your Inputs

You need one FASTA per segment and one outgroup/reference sequence per segment:

```text
sequences_L.fasta
sequences_M.fasta
sequences_S.fasta
config/outgroup_L.gb
config/outgroup_M.gb
config/outgroup_S.gb
```

Segment FASTA headers must use this pattern:

```fasta
>SampleName|L
ATGC...
```

Use `|M` for M-segment records and `|S` for S-segment records. The sample name before `|` is how records are matched across segments.

Quick path check:

```bash
ls sequences_L.fasta sequences_M.fasta sequences_S.fasta
ls config/outgroup_L.gb config/outgroup_M.gb config/outgroup_S.gb
```

## 2. Run The Default Pipeline

```bash
nextflow run main.nf \
    --l_fasta  sequences_L.fasta \
    --m_fasta  sequences_M.fasta \
    --s_fasta  sequences_S.fasta \
    --root_l   config/outgroup_L.gb \
    --root_m   config/outgroup_M.gb \
    --root_s   config/outgroup_S.gb \
    --outdir   results \
    -profile   conda
```

The main tree is written to:

```text
results/03_trees/concatenated.nwk
```

The main visual outputs are:

```text
results/04_visualization/alignment_tree.png
results/04_visualization/alignment_tree.svg
results/04_visualization/alignment_barcode.png
results/04_visualization/alignment_barcode.svg
```

## 3. Highlight Isolates Of Interest

Create a plain-text file with one sample name per line. Use the sample name without the `|L`, `|M`, or `|S` suffix.

```text
# config/highlight_samples.txt
SampleA
SampleB
SampleC
```

Then pass it with `--highlight_samples`:

```bash
nextflow run main.nf \
    --l_fasta            sequences_L.fasta \
    --m_fasta            sequences_M.fasta \
    --s_fasta            sequences_S.fasta \
    --root_l             config/outgroup_L.gb \
    --root_m             config/outgroup_M.gb \
    --root_s             config/outgroup_S.gb \
    --highlight_samples  config/highlight_samples.txt \
    --outdir             results_highlighted \
    -profile             conda
```

Highlighted samples appear in `alignment_tree.png/svg` with red labels and marked tree tips.

## 4. Match The hodcroftlab/andv-Style Run

This mode removes the reference from trees, midpoint-roots the outputs, and uses `GTR`:

```bash
nextflow run main.nf \
    --l_fasta          sequences_L.fasta \
    --m_fasta          sequences_M.fasta \
    --s_fasta          sequences_S.fasta \
    --root_l           config/outgroup_L.gb \
    --root_m           config/outgroup_M.gb \
    --root_s           config/outgroup_S.gb \
    --remove_reference true \
    --root             midpoint \
    --iqtree_model     GTR \
    --outdir           results_andv \
    -profile           conda
```

Do not combine `--remove_reference true` with `--root outgroup`; the pipeline exits because the OUTGROUP cannot be used for rooting after it has been removed.

## 5. Exclude Dropped Strains

Create a file with one sample name per line:

```text
# config/dropped_strains.txt
ProblemSample1
ProblemSample2
```

Run with:

```bash
nextflow run main.nf \
    --l_fasta          sequences_L.fasta \
    --m_fasta          sequences_M.fasta \
    --s_fasta          sequences_S.fasta \
    --root_l           config/outgroup_L.gb \
    --root_m           config/outgroup_M.gb \
    --root_s           config/outgroup_S.gb \
    --dropped_strains  config/dropped_strains.txt \
    --outdir           results_filtered \
    -profile           conda
```

Dropped strains are removed before length filtering and complete-trio filtering.

## 6. Use All Passing Samples For Per-Segment Trees

By default, per-segment trees and the concatenated tree use the same complete-trio sample set. To let each per-segment tree include every sample that passes that segment's length threshold:

```bash
nextflow run main.nf \
    --l_fasta                    sequences_L.fasta \
    --m_fasta                    sequences_M.fasta \
    --s_fasta                    sequences_S.fasta \
    --root_l                     config/outgroup_L.gb \
    --root_m                     config/outgroup_M.gb \
    --root_s                     config/outgroup_S.gb \
    --segment_trees_all_passing  true \
    --outdir                     results_segment_all \
    -profile                     conda
```

The concatenated tree remains complete-trio only.

## 7. Resume A Run

If a run stops partway through, rerun the same command with `-resume`:

```bash
nextflow run main.nf \
    --l_fasta  sequences_L.fasta \
    --m_fasta  sequences_M.fasta \
    --s_fasta  sequences_S.fasta \
    --root_l   config/outgroup_L.gb \
    --root_m   config/outgroup_M.gb \
    --root_s   config/outgroup_S.gb \
    --outdir   results \
    -profile   conda \
    -resume
```

If report files already exist in the output directory, use a fresh `--outdir` or see [troubleshooting.md](troubleshooting.md).

## 8. Run On SLURM

```bash
nextflow run main.nf \
    --l_fasta  sequences_L.fasta \
    --m_fasta  sequences_M.fasta \
    --s_fasta  sequences_S.fasta \
    --root_l   config/outgroup_L.gb \
    --root_m   config/outgroup_M.gb \
    --root_s   config/outgroup_S.gb \
    --outdir   results_slurm \
    -profile   slurm,conda
```

Adjust the `slurm` profile in `nextflow.config` if your site requires account, partition, or queue settings.
