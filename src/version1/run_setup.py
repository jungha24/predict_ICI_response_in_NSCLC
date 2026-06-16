from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple

from .utils_io import (
    RunConfig,
    load_config_json,
    setup_logger,
    ensure_dirs,
    append_config_history,
)

import logging


def init_run(
    config_path: str | Path,
    logger_name: str,
    run_name: Optional[str] = None,
    note: Optional[str] = None,
    log_dir: str | Path = "logs",
    config_history_path: str | Path = "logs/config_history.jsonl",
) -> Tuple[RunConfig, logging.Logger]:
    """
    One-time run initialization:
    - Load cfg from configs/*.json (override RunConfig defaults)
    - Ensure dirs
    - Setup logger
    - Append config snapshot to JSONL history
    """
    cfg = load_config_json(config_path)

    # ensure project dirs based on cfg.project_dir
    ensure_dirs(cfg)

    logger = setup_logger(log_dir=Path(cfg.project_dir) / log_dir, name=logger_name)

    append_config_history(
        cfg,
        log_path=Path(cfg.project_dir) / config_history_path,
        run_name=run_name,
        source_config=config_path,
        note=note,
    )

    logger.info(f"Loaded config: {config_path}")
    logger.info(f"Run name: {run_name} | note: {note}")
    return cfg, logger