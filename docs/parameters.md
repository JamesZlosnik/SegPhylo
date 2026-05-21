# Parameter Reference

All parameters are set on the command line with `--param value`. Defaults are defined in `nextflow.config` and can be overridden per-run.

---

## Required parameters

These must be provided on every run. Omitting all of them shows the help message; omitting some produces an error listing every missing parameter.

| Parameter | Description |
|-----------|-------------|
| `--l_fasta` | Multifasta file for segment L. Headers must follow `>SampleName\|L`. |
| `--m_fasta` | Multifasta file for segment M. Headers must follow `>SampleName\|M`. |
| `--s_fasta` | Multifasta file for segment S. Headers must follow `>SampleName\|S`. |
| `--root_l` | Root/outgroup sequence for segment L. Accepts FASTA or GenBank (see [formats](#root-sequence-formats)). |
| `--root_m` | Root/outgroup sequence for segment M. Accepts FASTA or GenBank. |
| `--root_s` | Root/outgroup sequence for segment S. Accepts FASTA or GenBank. |

---

## Filtering parameters

Applied in `FILTER_COMPLETE_SAMPLES` in this order: dropped strains → length filter → complete trio filter.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--dropped_strains` | `null` | Path to a plain-text file of sample names to exclude before any other filtering, one name per line. Lines beginning with `#` are treated as comments. Mirrors the `config/dropped_strains.txt` approach used in hodcroftlab/andv. |
| `--min_len_l` | `6000` | Minimum sequence length (bp) for segment L. Sequences shorter than this are discarded. Length is measured after stripping gap characters. |
| `--min_len_m` | `3000` | Minimum sequence length (bp) for segment M. |
| `--min_len_s` | `1000` | Minimum sequence length (bp) for segment S. |

The filter emits both complete-trio FASTAs and per-segment FASTAs. Complete-trio FASTAs contain only samples passing all three segment filters and present in all three segments. Per-segment FASTAs contain every sample passing the threshold for that segment.

---

## Alignment parameters

Applied in `MAFFT_ALIGN`.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--mafft_args` | `--auto` | Arguments passed to MAFFT in addition to `--keeplength` and `--addfragments`, which are always set by the pipeline. For higher accuracy at the cost of runtime, use `"--localpair --maxiterate 1000"` (L-INS-i, recommended for up to ~200 sequences per segment). |

> **Note:** Do not include `--keeplength` or `--addfragments` in `--mafft_args`; they are always appended automatically.

---

## Masking parameters

Applied in `MASK_ALIGNMENT` after alignment. Masked positions are set to `N`. Defaults match [hodcroftlab/andv](https://github.com/hodcroftlab/andv).

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--mask_from_start` | `30` | Number of positions to mask from the **start** of each segment alignment. |
| `--mask_from_end` | `50` | Number of positions to mask from the **end** of each segment alignment. |

Set either to `0` to disable masking for that end. The same values are used to draw shaded bands on the barcode visualisation.

---

## Reference / rooting parameters

These two parameters work together and control how the root sequence is used in tree building.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--remove_reference` | `false` | Remove the OUTGROUP sequence from the alignment before tree building. Mirrors `augur align --remove-reference`. When `true`, the OUTGROUP is used only as a coordinate reference for alignment and does not appear in the output trees. The barcode visualisation always retains the OUTGROUP row regardless of this setting. |
| `--root` | `outgroup` | Tree rooting method. Accepted values: `outgroup` or `midpoint`. `outgroup` roots the tree on the OUTGROUP sequence using IQ-TREE2's `-o OUTGROUP` flag. `midpoint` applies Biopython midpoint rooting after tree building, mirroring `augur refine --root mid_point`. |

**Valid combinations:**

| `--remove_reference` | `--root` | Behaviour |
|---------------------|----------|-----------|
| `false` | `outgroup` | Default. OUTGROUP in alignment and tree; outgroup-rooted. Most rigorous option when a genuine outgroup is available. |
| `false` | `midpoint` | OUTGROUP in alignment and tree; midpoint-rooted. |
| `true` | `midpoint` | OUTGROUP used for alignment only, then removed. Tree is midpoint-rooted. Closest to hodcroftlab/andv. |
| `true` | `outgroup` | ❌ **Not allowed.** Exits with an error — the OUTGROUP cannot be used for rooting if it has been removed from the alignment. |

---

## Tree-building parameters

Applied in `IQTREE_SEGMENT` and `IQTREE_CONCATENATED`.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--segment_trees_all_passing` | `false` | When `false`, per-segment trees use the same complete-trio sample set as the concatenated tree. When `true`, each per-segment tree uses all samples passing that segment's length threshold. The concatenated tree remains complete-trio only. |
| `--iqtree_model` | `GTR+G` | Substitution model for the per-segment trees and for each concatenated partition. Any model string accepted by IQ-TREE2 is valid (e.g. `GTR`, `GTR+G+I`, `TVM+F+G4`). Use `GTR` to match the hodcroftlab/andv default. |
| `--iqtree_boot` | `1000` | Number of ultrafast bootstrap replicates (`-B`). Applied to both per-segment and concatenated trees. |

> **Note:** The **concatenated tree** always uses a partitioned model, with one partition per segment in `partitions.txt`. `--iqtree_model` controls the model string written for each partition.
>
> The flags `-czb` (collapse zero-length branches), `-st DNA` (explicit DNA sequence type), and `--redo` are always set and cannot be overridden.

---

## Visualisation parameters

Applied in `VISUALIZE_ALIGNMENT`.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--label_taxa_max` | `250` | Maximum number of taxa before sample labels are hidden on the barcode plot. Above this threshold only the OUTGROUP label is shown, with a note on the y-axis. Override to force labels: `--label_taxa_max 600`. |
| `--highlight_samples` | `null` | Optional plain-text file of sample names to highlight in the tree-plus-alignment plot. Names should match the sample ID before `|SEGMENT`. Blank lines and `#` comments are ignored. Highlighted samples use red labels and tree-tip markers. |

---

## General parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | `results` | Directory for all published output files. Created if it does not exist. |
| `--help` | `false` | Print the help message and exit. Also triggered automatically when no required parameters are provided. |

---

## Root sequence formats

The `--root_l`, `--root_m`, and `--root_s` parameters accept files with any of the following extensions:

| Extension | Format parsed as |
|-----------|-----------------|
| `.fasta` `.fa` `.fna` `.fas` | FASTA |
| `.gb` `.gbk` `.genbank` | GenBank |
| Any other | FASTA (with a warning) |

If a root file contains multiple records, the first record is used and a warning is printed to the log.

---

## Example commands

**Show help:**
```bash
nextflow run main.nf --help
```

**Default run (outgroup rooting, GTR+G):**
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

**Closest to hodcroftlab/andv (remove reference, midpoint rooting, GTR):**
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

**With dropped strains exclusion list:**
```bash
nextflow run main.nf \
    --l_fasta         sequences_L.fasta \
    --m_fasta         sequences_M.fasta \
    --s_fasta         sequences_S.fasta \
    --root_l          config/outgroup_L.gb \
    --root_m          config/outgroup_M.gb \
    --root_s          config/outgroup_S.gb \
    --dropped_strains config/dropped_strains.txt \
    --outdir          results \
    -profile          conda
```

**Highlight isolates of interest in the tree-plus-alignment plot:**
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

**Build per-segment trees from all length-passing segment samples:**
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

**High-accuracy alignment (L-INS-i, recommended for ≤200 sequences per segment):**
```bash
nextflow run main.nf \
    --l_fasta      sequences_L.fasta \
    --m_fasta      sequences_M.fasta \
    --s_fasta      sequences_S.fasta \
    --root_l       config/outgroup_L.gb \
    --root_m       config/outgroup_M.gb \
    --root_s       config/outgroup_S.gb \
    --mafft_args   "--localpair --maxiterate 1000" \
    --outdir       results_hq \
    -profile       conda
```

**Midpoint rooting, keep reference in alignment:**
```bash
nextflow run main.nf \
    --l_fasta  sequences_L.fasta \
    --m_fasta  sequences_M.fasta \
    --s_fasta  sequences_S.fasta \
    --root_l   config/outgroup_L.gb \
    --root_m   config/outgroup_M.gb \
    --root_s   config/outgroup_S.gb \
    --root     midpoint \
    --outdir   results_midpoint \
    -profile   conda
```

**Disable terminal masking:**
```bash
nextflow run main.nf \
    --l_fasta         sequences_L.fasta \
    --m_fasta         sequences_M.fasta \
    --s_fasta         sequences_S.fasta \
    --root_l          config/outgroup_L.gb \
    --root_m          config/outgroup_M.gb \
    --root_s          config/outgroup_S.gb \
    --mask_from_start 0 \
    --mask_from_end   0 \
    --outdir          results_unmasked \
    -profile          conda
```

**SLURM cluster:**
```bash
nextflow run main.nf \
    --l_fasta  sequences_L.fasta \
    --m_fasta  sequences_M.fasta \
    --s_fasta  sequences_S.fasta \
    --root_l   config/outgroup_L.gb \
    --root_m   config/outgroup_M.gb \
    --root_s   config/outgroup_S.gb \
    --outdir   results \
    -profile   slurm,conda
```

**Resume a failed or interrupted run:**
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
