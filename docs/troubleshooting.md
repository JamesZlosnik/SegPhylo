# Troubleshooting

Common issues and quick checks for local and cluster runs.

---

## `ModuleNotFoundError: No module named 'Bio'`

This means a task is running with a Python environment that does not include Biopython.

Use the conda profile:

```bash
nextflow run main.nf ... -profile conda
```

The pipeline calls helper scripts with `python ${projectDir}/bin/...` so they use the Python from the active Nextflow conda environment. If this error appears, check the failed task wrapper:

```bash
sed -n '/conda environment/,+5p' work/*/*/.command.run
```

The wrapper should activate a `work/conda/env-*` environment before running `.command.sh`.

---

## Existing report or timeline files block reruns

Nextflow can fail at the end of a rerun if these files already exist:

```text
pipeline_report.html
pipeline_timeline.html
pipeline_trace.txt
```

Use a fresh `--outdir`, remove or rename the old report files, or configure report overwriting in `nextflow.config` if that fits your workflow.

---

## `--root outgroup` with `--remove_reference true`

This combination is intentionally rejected. If the OUTGROUP is removed before tree building, IQ-TREE2 cannot root on it.

Use:

```bash
--remove_reference true --root midpoint
```

or keep the OUTGROUP in the tree:

```bash
--remove_reference false --root outgroup
```

---

## Concatenated tree has fewer samples than per-segment trees

This is expected when:

```bash
--segment_trees_all_passing true
```

Per-segment trees can include samples that pass only that segment's threshold. The concatenated tree is always restricted to samples present in all three segments after filtering.

---

## Input files are not found

Paths are interpreted relative to the directory where `nextflow run main.nf` is executed. Use absolute paths or run from the project directory if in doubt.

Quick check:

```bash
ls sequences_L.fasta sequences_M.fasta sequences_S.fasta
ls config/outgroup_L.gb config/outgroup_M.gb config/outgroup_S.gb
```

---

## Resume behaviour looks stale

Nextflow caches completed tasks in `work/`. After changing code or parameters, Nextflow usually detects what must be rerun. If results still look stale, rerun without `-resume` or use a fresh `--outdir`.
