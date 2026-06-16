# Data Policy

This repository is for code, workflow documentation, configuration, and selected reviewable outputs. It intentionally excludes raw data and sensitive or large intermediate analysis products.

## Excluded From Git

The following file classes are excluded:

- Raw and processed single-cell objects: `.rds`, `.RData`, `.h5ad`, `.loom`
- Matrix files: `.mtx`, `.mtx.gz`
- Large exported tables, logs, and intermediate objects
- Patient-level clinical metadata and merged analysis datasets
- Cell-level observation tables, score matrices, and barcodes
- Local presentation and note files unrelated to the NSCLC GitHub repository

## Server Data Location

The original full working directory was:

```text
/data/podo/Projects/project_jhl/20260215_nsclc
```

The heaviest server-side directory was:

```text
/data/podo/Projects/project_jhl/20260215_nsclc/data
```

That directory was not copied into Git.

## Included Outputs

Included outputs are limited to small files useful for review:

- `run_manifest.json`, `input_manifest.json`
- model metrics and comparison summaries
- selected feature and pruning summaries
- representative figures
- non-sensitive workflow notes

When in doubt, patient-level tables are excluded even if the files are small.
