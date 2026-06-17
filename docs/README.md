# File Structure

1. Keep notebooks lightweight.
   - Load and override configuration values in notebooks only when needed.
2. Keep `src` focused on reusable analysis logic.
3. Use config files to externalize parameters and preserve runnable snapshots.
   - Keep `RunConfig` defaults as fallbacks.
   - Read run-specific values from `configs/*.json`.
   - Append the final configuration to `logs/config_history.jsonl` for each run.
   - Store the final configuration in `adata.uns["run_config"]` when saving AnnData objects.

# Parameter Handling

1. `cfg`: default operating values from the config file
   - Definition: `RunConfig`
   - Input file: `configs/nsclc_v1.json`
   - Included fields:
     - `min_counts`, `max_counts`, `min_genes`, `max_genes`
     - `max_pct_mt`, `max_pct_ribo`, `max_pct_hb`
     - `expected_doublet_rate`
     - `use_*_filter` flags for enabling or disabling each filter

2. `params`: ad hoc notebook overrides
   - Definition: `QCParams` in `src/qc.py`
   - Included fields:
     - The same threshold fields listed above.
     - Optional overrides for the `use_*_filter` flags.

# PBMC QC Notes

1. Basic QC
   - `n_counts`
   - `n_genes`
   - `percentage_mt`
   - doublet score
2. Contamination checks
   - RBC contamination: hemoglobin genes and `percentage_hb`
   - Platelet signal: possible ambient RNA signal using platelet markers
   - Ambient RNA
3. Interpretation
   - Separate cell-type proportion signals from cell-state signals during interpretation.
