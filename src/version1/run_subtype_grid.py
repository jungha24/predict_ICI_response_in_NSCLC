#!/usr/bin/env python

from __future__ import annotations

import argparse
import copy
import itertools
import json
import subprocess
from pathlib import Path

import pandas as pd
import yaml


def read_yaml(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def write_yaml(obj: dict, path: str | Path) -> None:
    with open(path, "w", encoding="utf-8") as f:
        yaml.safe_dump(obj, f, sort_keys=False, allow_unicode=True)


def set_nested(cfg: dict, dotted_key: str, value):
    parts = dotted_key.split(".")
    cur = cfg
    for p in parts[:-1]:
        if p not in cur or not isinstance(cur[p], dict):
            cur[p] = {}
        cur = cur[p]
    cur[parts[-1]] = value


def main(sweep_config_path: str):
    sweep_cfg = read_yaml(sweep_config_path)

    base_config_path = sweep_cfg["base_config"]
    run_root = Path(sweep_cfg["run_root"])
    grid = sweep_cfg["grid"]

    run_root.mkdir(parents=True, exist_ok=True)

    base_cfg = read_yaml(base_config_path)

    keys = list(grid.keys())
    values = [grid[k] for k in keys]

    summary_rows = []

    for idx, combo in enumerate(itertools.product(*values), start=1):
        cfg = copy.deepcopy(base_cfg)

        combo_dict = dict(zip(keys, combo))
        for dotted_key, value in combo_dict.items():
            set_nested(cfg, dotted_key, value)

        run_name = f"run_{idx:03d}"
        outdir = run_root / run_name
        outdir.mkdir(parents=True, exist_ok=True)

        cfg["paths"]["outdir"] = str(outdir)

        config_out = outdir / "config.yaml"
        write_yaml(cfg, config_out)

        cmd = ["python", "src/version1/run_subtype_pipeline.py", "--config", str(config_out)]
        result = subprocess.run(cmd, capture_output=True, text=True)

        record = {"run_name": run_name, **combo_dict, "returncode": result.returncode}

        metrics_path = outdir / "metrics.json"
        if metrics_path.exists():
            with open(metrics_path, "r", encoding="utf-8") as f:
                metrics = json.load(f)
            record.update(metrics)

        summary_rows.append(record)

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(run_root / "grid_summary.csv", index=False)
    print(f"Saved: {run_root / 'grid_summary.csv'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="Path to sweep YAML config")
    args = parser.parse_args()
    main(args.config)