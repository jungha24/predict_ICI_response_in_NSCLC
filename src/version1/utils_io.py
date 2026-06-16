from __future__ import annotations

import time
import json
import logging
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, Literal, Any

try:
    import scanpy as sc
    import anndata as ad
except Exception as e:
    raise RuntimeError("Please install scanpy/anndata in this environment.") from e


# -------------------------
# Logging / timing utilities
# -------------------------

def setup_logger(
    log_dir: str | Path = "logs",
    name: str = "sc_pipeline",
    level: int = logging.INFO,
) -> logging.Logger:
    Path(log_dir).mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.handlers.clear()

    fmt = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    fh = logging.FileHandler(Path(log_dir) / f"{name}.log")
    fh.setLevel(level)
    fh.setFormatter(fmt)

    sh = logging.StreamHandler()
    sh.setLevel(level)
    sh.setFormatter(fmt)

    logger.addHandler(fh)
    logger.addHandler(sh)
    logger.propagate = False
    return logger


class Timer:
    def __init__(self, logger: logging.Logger, label: str):
        self.logger = logger
        self.label = label
        self.t0 = None

    def __enter__(self):
        self.t0 = time.time()
        self.logger.info(f"[START] {self.label}")
        return self

    def __exit__(self, exc_type, exc, tb):
        dt = time.time() - (self.t0 or time.time())
        self.logger.info(f"[END] {self.label} | {dt:.2f}s")


# -------------------------
# Config / dirs
# -------------------------

@dataclass
class RunConfig:
    project_dir: str = "."
    fig_dir: str = "figures"
    data_dir: str = "data"

    sample_col: str = "sample"
    batch_col: str = "batch"
    patient_col: str = "patient"

    seed: int = 0
    dataset_type: Literal["pbmc", "tumor", "generic"] = "generic"

    # QC thresholds (None -> auto-suggest in QC step)
    min_counts: Optional[float] = None
    use_min_counts_filter: bool = True
    max_counts: Optional[float] = None
    use_max_counts_filter: bool = True
    min_genes: Optional[float] = None
    use_min_genes_filter: bool = True
    max_genes: Optional[float] = None
    max_pct_mt: Optional[float] = None
    use_max_pct_mt_filter: bool = True
    use_max_genes_filter: bool = True

    max_pct_ribo: Optional[float] = None
    use_max_pct_ribo_filter: bool = False

    max_pct_hb: Optional[float] = None
    # auto: pbmc일 때만 고려, on: 항상 고려, off: 미사용
    use_max_pct_hb_filter: Literal["auto", "on", "off"] = "auto"

    # Doublet settings
    run_doublet: bool = True
    expected_doublet_rate: float = 0.06

    # If True, reading huge .h5ad in backed mode (read-only)
    use_backed: bool = False
    n_jobs: int = 1

    def to_json(self, path: str | Path) -> None:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            json.dump(asdict(self), f, indent=2)


def update_config(cfg: RunConfig, overrides: dict[str, Any]) -> RunConfig:
    valid_fields = set(RunConfig.__dataclass_fields__.keys())
    unknown = set(overrides.keys()) - valid_fields
    if unknown:
        unknown_str = ", ".join(sorted(unknown))
        raise KeyError(f"Unknown config keys: {unknown_str}")

    dataset_type = overrides.get("dataset_type")
    if dataset_type is not None and dataset_type not in {"pbmc", "tumor", "generic"}:
        raise ValueError("dataset_type must be one of: pbmc, tumor, generic")
    hb_mode = overrides.get("use_max_pct_hb_filter")
    if hb_mode is not None and hb_mode not in {"auto", "on", "off"}:
        raise ValueError("use_max_pct_hb_filter must be one of: auto, on, off")
    expected_doublet_rate = overrides.get("expected_doublet_rate")
    if expected_doublet_rate is not None:
        if not (0.0 < float(expected_doublet_rate) < 1.0):
            raise ValueError("expected_doublet_rate must be in (0, 1)")

    for key, value in overrides.items():
        setattr(cfg, key, value)
    return cfg


def load_config_json(path: str | Path, base_cfg: Optional[RunConfig] = None) -> RunConfig:
    """
    Load config overrides from JSON and merge into RunConfig defaults.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    with open(path) as f:
        overrides = json.load(f)

    if not isinstance(overrides, dict):
        raise ValueError("Config JSON must contain a JSON object.")

    cfg = base_cfg or RunConfig()
    return update_config(cfg, overrides)


def append_config_history(
    cfg: RunConfig,
    log_path: str | Path = "logs/config_history.jsonl",
    run_name: Optional[str] = None,
    source_config: Optional[str | Path] = None,
    note: Optional[str] = None,
) -> Path:
    """
    Append one config snapshot per run as JSONL for auditability.
    """
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    payload = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "run_name": run_name,
        "source_config": str(source_config) if source_config is not None else None,
        "note": note,
        "config": asdict(cfg),
    }
    with open(log_path, "a") as f:
        f.write(json.dumps(payload) + "\n")
    return log_path


def attach_run_config_to_adata(
    adata,
    cfg: RunConfig,
    source_config: Optional[str | Path] = None,
    note: Optional[str] = None,
) -> None:
    """
    Store the exact run config inside adata.uns for reproducibility.
    """
    adata.uns["run_config"] = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "source_config": str(source_config) if source_config is not None else None,
        "note": note,
        "config": asdict(cfg),
    }


def ensure_dirs(cfg: RunConfig) -> None:
    Path(cfg.project_dir, cfg.fig_dir).mkdir(parents=True, exist_ok=True)
    Path(cfg.project_dir, cfg.data_dir).mkdir(parents=True, exist_ok=True)
    Path(cfg.project_dir, "logs").mkdir(parents=True, exist_ok=True)


# -------------------------
# Data loading / saving  (I/O는 여기로 고정)
# -------------------------

def load_h5ad(path: str | Path, cfg: Optional[RunConfig] = None):
    """
    Load AnnData from .h5ad.
    If cfg.use_backed=True, loads in backed mode (read-only, avoids huge RAM).
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    if cfg and cfg.use_backed:
        return ad.read_h5ad(path, backed="r")
    return ad.read_h5ad(path)


def load_10x_h5(path: str | Path):
    """
    Load 10x Genomics .h5 (filtered_feature_bc_matrix.h5 / raw_feature_bc_matrix.h5).
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    return sc.read_10x_h5(path)


def save_h5ad(adata, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path)


def checkpoint_path(cfg: RunConfig, name: str) -> Path:
    return Path(cfg.project_dir, cfg.data_dir) / f"{name}.h5ad"


def save_checkpoint(adata, cfg: RunConfig, name: str, logger: Optional[logging.Logger] = None) -> Path:
    p = checkpoint_path(cfg, name)
    if logger:
        logger.info(f"Saving checkpoint: {p}")
    save_h5ad(adata, p)
    return p


def load_checkpoint(cfg: RunConfig, name: str):
    p = checkpoint_path(cfg, name)
    return load_h5ad(p, cfg=None)


# -------------------------
# Sanity checks
# -------------------------

def assert_basic_sanity(adata, cfg: RunConfig, logger: logging.Logger) -> None:
    if adata.var_names.has_duplicates:
        logger.warning("var_names has duplicates. Consider .var_names_make_unique().")

    for col in [cfg.sample_col, cfg.batch_col]:
        if col not in adata.obs.columns:
            logger.warning(f"Missing obs column '{col}'. Recommend adding it for multi-sample/batch analyses.")

    X = adata.X
    is_sparse = hasattr(X, "tocsr")
    logger.info(f"X type: {type(X)} | sparse={is_sparse}")
