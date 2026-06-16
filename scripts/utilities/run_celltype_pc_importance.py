#!/usr/bin/env python

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd
import yaml


DEFAULT_CELLTYPES = [
    "effector CD8+T",
    "CD14 Mono",
    "NK",
    "naive/TCM CD4+T",
]


def log_step(message: str) -> None:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}", flush=True)


def sanitize_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("_")


def load_yaml(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def write_yaml(payload: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        yaml.safe_dump(payload, f, sort_keys=False, allow_unicode=True)


def load_json(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        payload = json.load(f)
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object at {path}, got {type(payload).__name__}")
    return payload


def parse_celltype_npcs_overrides(items: list[str]) -> dict[str, int]:
    overrides: dict[str, int] = {}
    for item in items:
        if "=" not in item:
            raise ValueError(
                f"Invalid --celltype-npcs entry '{item}'. Use the format 'celltype=npcs'."
            )
        celltype, npcs = item.split("=", 1)
        celltype = celltype.strip()
        if not celltype:
            raise ValueError(f"Invalid --celltype-npcs entry '{item}': empty celltype.")
        try:
            overrides[celltype] = int(npcs)
        except ValueError as e:
            raise ValueError(
                f"Invalid --celltype-npcs entry '{item}': '{npcs}' is not an integer."
            ) from e
    return overrides


def resolve_celltype_npcs(
    celltypes: list[str],
    default_npcs: int,
    selection_file: Path | None,
    inline_overrides: list[str],
) -> dict[str, int]:
    resolved = {celltype: int(default_npcs) for celltype in celltypes}

    if selection_file is not None:
        payload = load_json(selection_file)
        for key, value in payload.items():
            if key in resolved:
                resolved[key] = int(value)

    for key, value in parse_celltype_npcs_overrides(inline_overrides).items():
        if key not in resolved:
            raise ValueError(
                f"--celltype-npcs specified '{key}', but it is not present in --celltypes."
            )
        resolved[key] = int(value)

    return resolved


def build_celltype_config(
    base_cfg: dict,
    celltype: str,
    npcs: int,
    outdir: Path,
) -> dict:
    cfg = json.loads(json.dumps(base_cfg))

    cfg.setdefault("celltypes", {})
    cfg["celltypes"]["proportion"] = []
    cfg["celltypes"]["pseudobulk_pca"] = [celltype]
    cfg["celltypes"]["pseudobulk_nmf"] = []

    cfg.setdefault("pca_features", {})
    cfg["pca_features"]["npcs_per_celltype"] = int(npcs)
    cfg["pca_features"]["max_npcs_to_compute"] = max(int(cfg["pca_features"].get("max_npcs_to_compute", npcs)), int(npcs))
    cfg["pca_features"]["pc_selection_file"] = None
    cfg["pca_features"]["pc_selection"] = {}

    cfg.setdefault("nmf_features", {})
    cfg["nmf_features"]["enabled"] = False
    cfg["nmf_features"]["selected_modules"] = []
    cfg["nmf_features"]["selected_modules_file"] = None

    cfg.setdefault("feature_selection", {})
    cfg["feature_selection"]["always_keep_proportions"] = False

    cfg.setdefault("paths", {})
    cfg["paths"]["outdir"] = str(outdir)

    return cfg


def run_command(cmd: list[str], cwd: Path) -> None:
    log_step("Running: " + " ".join(cmd))
    subprocess.run(cmd, cwd=str(cwd), check=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run PCA-only celltype-specific subtype discovery and permutation importance in batch."
    )
    parser.add_argument(
        "--base-config",
        default="configs/version1_subtype_base.yaml",
        help="Base YAML config to clone and modify",
    )
    parser.add_argument(
        "--run-root",
        default="data/20260309_pilot/results/version1/subtype_discovery/pc_only_celltype_sweep",
        help="Directory where per-celltype runs and generated configs will be written",
    )
    parser.add_argument(
        "--celltypes",
        nargs="+",
        default=DEFAULT_CELLTYPES,
        help="Celltypes to run one-by-one",
    )
    parser.add_argument(
        "--npcs",
        type=int,
        default=10,
        help="Default number of PCs to keep for each single-celltype run",
    )
    parser.add_argument(
        "--pc-selection-file",
        default=None,
        help="Optional JSON file mapping celltype to PC count. Missing celltypes fall back to --npcs.",
    )
    parser.add_argument(
        "--celltype-npcs",
        nargs="*",
        default=[],
        help="Optional inline overrides in the format 'celltype=npcs'. These override --pc-selection-file.",
    )
    parser.add_argument(
        "--n-permutations",
        type=int,
        default=50,
        help="Number of permutations per feature for importance scoring",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=123,
        help="Random seed passed to permutation importance",
    )
    parser.add_argument(
        "--python",
        default=sys.executable,
        help="Python executable to use for downstream runs. Defaults to the current interpreter.",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    base_config_path = (repo_root / args.base_config).resolve() if not Path(args.base_config).is_absolute() else Path(args.base_config)
    run_root = (repo_root / args.run_root).resolve() if not Path(args.run_root).is_absolute() else Path(args.run_root)
    config_dir = run_root / "generated_configs"
    run_root.mkdir(parents=True, exist_ok=True)
    config_dir.mkdir(parents=True, exist_ok=True)

    base_cfg = load_yaml(base_config_path)
    pc_selection_file = None
    if args.pc_selection_file:
        pc_selection_file = (
            (repo_root / args.pc_selection_file).resolve()
            if not Path(args.pc_selection_file).is_absolute()
            else Path(args.pc_selection_file)
        )
    celltype_npcs = resolve_celltype_npcs(
        celltypes=list(args.celltypes),
        default_npcs=int(args.npcs),
        selection_file=pc_selection_file,
        inline_overrides=list(args.celltype_npcs),
    )

    summary_tables = []

    for celltype in args.celltypes:
        npcs = int(celltype_npcs[celltype])
        slug = sanitize_name(celltype)
        run_dir = run_root / f"{slug}_pc1to{npcs}"
        config_path = config_dir / f"{slug}_pc1to{npcs}.yaml"

        cfg = build_celltype_config(
            base_cfg=base_cfg,
            celltype=celltype,
            npcs=npcs,
            outdir=run_dir,
        )
        write_yaml(cfg, config_path)

        log_step(f"Prepared config for {celltype} with PC1-{npcs}: {config_path}")
        run_command(
            [args.python, "src/version1/run_subtype_pipeline.py", "--config", str(config_path)],
            cwd=repo_root,
        )
        run_command(
            [
                args.python,
                "src/version1/permutation_feature_importance.py",
                "--run-dir",
                str(run_dir),
                "--n-permutations",
                str(int(args.n_permutations)),
                "--seed",
                str(int(args.seed)),
            ],
            cwd=repo_root,
        )

        summary_path = run_dir / "permutation_feature_importance" / "permutation_importance_summary.csv"
        if summary_path.exists():
            df = pd.read_csv(summary_path)
            df.insert(0, "celltype", celltype)
            df.insert(1, "run_dir", str(run_dir))
            pc_index = df["feature"].str.extract(r"_PC(\d+)$", expand=False)
            df["pc_index"] = pd.to_numeric(pc_index, errors="coerce")
            summary_tables.append(df)

    if summary_tables:
        combined = pd.concat(summary_tables, ignore_index=True)
        sort_col = "importance_fixed_k" if "importance_fixed_k" in combined.columns else "importance"
        combined = combined.sort_values(["celltype", sort_col], ascending=[True, False]).reset_index(drop=True)
        combined.to_csv(run_root / "permutation_importance_all_celltypes.csv", index=False)
        combined.groupby("celltype", as_index=False).head(5).to_csv(
            run_root / "permutation_importance_top5_per_celltype.csv",
            index=False,
        )
        log_step(f"Wrote combined importance tables to {run_root}")
    else:
        log_step("No permutation importance summaries were found to aggregate.")


if __name__ == "__main__":
    main()
