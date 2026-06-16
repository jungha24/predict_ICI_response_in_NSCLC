from __future__ import annotations

import argparse
import copy
import itertools
import os
import resource
import socket
import sys
import threading
import time
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from typing import Iterable, Optional
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
try:
    from lifelines.utils import concordance_index
except ImportError:  # pragma: no cover - optional when running binary-only configs
    concordance_index = None
from sklearn.metrics import average_precision_score, brier_score_loss
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold

from .data import build_search_dataset, drop_bad_columns
from .design import BinaryEndpointSpec, SearchDataset, SurvivalEndpointSpec, resolve_endpoint_specs
from .io_utils import ensure_outdir, log_step, read_yaml, write_json, write_yaml
from .models import (
    bootstrap_stability_cox,
    bootstrap_stability_logistic,
    evaluate_cox_nested_cv,
    evaluate_logistic_nested_cv,
    fit_cox_model,
    fit_logistic_model,
    resolve_stratified_splits,
    safe_auc,
    tune_cox_elastic_net,
    tune_logistic_elastic_net,
)

try:
    import psutil
except ImportError:  # pragma: no cover - optional dependency
    psutil = None


_ACTIVE_RUN_TRACKER = None


def _require_lifelines() -> None:
    if concordance_index is None:
        raise ImportError(
            "lifelines is required for survival/Cox modeling. "
            "Install lifelines or use a binary-only endpoint config."
        )


def _logging_cfg(cfg: dict) -> dict:
    logging_cfg = cfg.get("logging", {}) or {}
    return {
        "verbose": bool(logging_cfg.get("verbose", True)),
        "stage_progress_every": max(1, int(logging_cfg.get("stage_progress_every", 25))),
        "resource_sample_interval_seconds": max(0.1, float(logging_cfg.get("resource_sample_interval_seconds", 0.5))),
    }


def _resolve_baseline_clinical_feature_set(cfg: dict) -> str:
    search_cfg = cfg.get("search", {}) or {}
    baseline_cfg = search_cfg.get("baseline", {}) or {}

    explicit_name = baseline_cfg.get("clinical_feature_set") or baseline_cfg.get("feature_set")
    if explicit_name:
        return str(explicit_name)

    mode = baseline_cfg.get("mode")
    if mode is not None:
        normalized = str(mode).strip().lower()
        mode_map = {
            "clinical_only": "base",
            "base": "base",
            "clinical_with_pd_l1": "with_pd_l1",
            "with_pd_l1": "with_pd_l1",
            "pd_l1": "with_pd_l1",
        }
        if normalized not in mode_map:
            raise ValueError(
                "search.baseline.mode must be one of: "
                "clinical_only, base, clinical_with_pd_l1, with_pd_l1, pd_l1"
            )
        return mode_map[normalized]

    include_pd_l1 = baseline_cfg.get("include_pd_l1")
    if include_pd_l1 is not None:
        return "with_pd_l1" if bool(include_pd_l1) else "base"

    return str(search_cfg.get("baseline_clinical_feature_set", "base"))


def _bytes_to_mb(value: float) -> float:
    return float(value) / (1024.0 * 1024.0)


def _ru_maxrss_bytes(usage: resource.struct_rusage) -> int:
    value = int(getattr(usage, "ru_maxrss", 0) or 0)
    if sys.platform == "darwin":
        return value
    return value * 1024


class RunTracker:
    def __init__(self, outdir: Path, cfg: dict) -> None:
        self.outdir = ensure_outdir(outdir / "runtime")
        self.cfg = cfg
        self.logging_cfg = _logging_cfg(cfg)
        self.parallel_cfg = _parallel_cfg(cfg)
        self.start_wall = time.time()
        self.start_perf = time.perf_counter()
        self.stage_events: list[dict[str, object]] = []
        self.log_path = self.outdir / "execution.log"
        self.timing_json_path = self.outdir / "run_timing.json"
        self.timing_tsv_path = self.outdir / "stage_timing.tsv"
        self.resource_json_path = self.outdir / "run_resource_summary.json"
        self.metadata = {
            "pid": int(os.getpid()),
            "hostname": socket.gethostname(),
            "cpu_count": int(os.cpu_count() or 1),
            "parallel": self.parallel_cfg,
        }
        self._log_handle = self.log_path.open("w", encoding="utf-8")
        self._sampler_stop = threading.Event()
        self._sampler_thread: Optional[threading.Thread] = None
        self._peak_tree_rss_bytes = 0
        self._peak_tree_processes = 1
        self._peak_self_rss_bytes = 0
        self._sampler_method = "resource"
        if psutil is not None:
            self._sampler_method = "psutil_process_tree"
            self._sampler_thread = threading.Thread(target=self._sample_process_tree_rss, daemon=True)
            self._sampler_thread.start()

    def log(self, message: str, *, verbose_only: bool = False) -> None:
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{timestamp}] {message}"
        if (not verbose_only) or self.logging_cfg["verbose"]:
            log_step(message)
        print(line, file=self._log_handle, flush=True)

    @contextmanager
    def stage(self, name: str, metadata: Optional[dict[str, object]] = None):
        start_iso = datetime.now().isoformat(timespec="seconds")
        start_perf = time.perf_counter()
        self.log(f"START {name}")
        try:
            yield
        except Exception as exc:
            end_iso = datetime.now().isoformat(timespec="seconds")
            elapsed = time.perf_counter() - start_perf
            self.stage_events.append(
                {
                    "stage": name,
                    "status": "failed",
                    "start_time": start_iso,
                    "end_time": end_iso,
                    "elapsed_seconds": float(elapsed),
                    "error": str(exc),
                    **(metadata or {}),
                }
            )
            self.log(f"FAIL {name} ({elapsed:.2f}s): {exc}")
            raise
        else:
            end_iso = datetime.now().isoformat(timespec="seconds")
            elapsed = time.perf_counter() - start_perf
            self.stage_events.append(
                {
                    "stage": name,
                    "status": "completed",
                    "start_time": start_iso,
                    "end_time": end_iso,
                    "elapsed_seconds": float(elapsed),
                    **(metadata or {}),
                }
            )
            self.log(f"END {name} ({elapsed:.2f}s)")

    def progress(self, message: str) -> None:
        self.log(message, verbose_only=True)

    def _sample_process_tree_rss(self) -> None:
        assert psutil is not None
        try:
            root = psutil.Process(self.metadata["pid"])
        except Exception:
            return
        interval = self.logging_cfg["resource_sample_interval_seconds"]
        while not self._sampler_stop.is_set():
            rss_total = 0
            process_count = 0
            try:
                processes = [root] + root.children(recursive=True)
            except Exception:
                processes = [root]
            for process in processes:
                try:
                    rss_total += int(process.memory_info().rss)
                    process_count += 1
                except Exception:
                    continue
            self._peak_tree_rss_bytes = max(self._peak_tree_rss_bytes, rss_total)
            self._peak_tree_processes = max(self._peak_tree_processes, process_count)
            try:
                self._peak_self_rss_bytes = max(self._peak_self_rss_bytes, int(root.memory_info().rss))
            except Exception:
                pass
            self._sampler_stop.wait(interval)

    def finalize(self) -> None:
        self._sampler_stop.set()
        if self._sampler_thread is not None:
            self._sampler_thread.join(timeout=2.0)

        end_wall = time.time()
        end_perf = time.perf_counter()
        usage_self = resource.getrusage(resource.RUSAGE_SELF)
        usage_children = resource.getrusage(resource.RUSAGE_CHILDREN)
        peak_self_bytes = max(self._peak_self_rss_bytes, _ru_maxrss_bytes(usage_self))
        peak_children_bytes = _ru_maxrss_bytes(usage_children)
        peak_tree_bytes = max(self._peak_tree_rss_bytes, peak_self_bytes, peak_children_bytes)

        timing_payload = {
            "started_at": datetime.fromtimestamp(self.start_wall).isoformat(timespec="seconds"),
            "finished_at": datetime.fromtimestamp(end_wall).isoformat(timespec="seconds"),
            "wall_time_seconds": float(end_perf - self.start_perf),
            "stage_events": self.stage_events,
        }
        write_json(timing_payload, self.timing_json_path)

        stage_df = pd.DataFrame(self.stage_events)
        stage_df.to_csv(self.timing_tsv_path, sep="\t", index=False)

        resource_payload = {
            **self.metadata,
            "sampler_method": self._sampler_method,
            "rss_peak_self_mb": _bytes_to_mb(peak_self_bytes),
            "rss_peak_children_mb": _bytes_to_mb(peak_children_bytes),
            "rss_peak_process_tree_mb": _bytes_to_mb(peak_tree_bytes),
            "peak_process_count_observed": int(self._peak_tree_processes),
            "wall_time_seconds": float(end_perf - self.start_perf),
        }
        write_json(resource_payload, self.resource_json_path)
        self.log(
            "Runtime summary: "
            f"wall={resource_payload['wall_time_seconds']:.2f}s, "
            f"peak_tree_rss={resource_payload['rss_peak_process_tree_mb']:.1f} MB, "
            f"pid={resource_payload['pid']}, cpu_count={resource_payload['cpu_count']}, "
            f"max_workers={self.parallel_cfg['max_workers']}"
        )
        self._log_handle.close()


def _get_tracker() -> Optional[RunTracker]:
    return _ACTIVE_RUN_TRACKER


@contextmanager
def _tracked_stage(name: str, metadata: Optional[dict[str, object]] = None):
    tracker = _get_tracker()
    if tracker is None:
        yield
        return
    with tracker.stage(name, metadata):
        yield


def _progress_log(cfg: dict, message: str) -> None:
    tracker = _get_tracker()
    if tracker is not None:
        tracker.progress(message)
    elif _logging_cfg(cfg)["verbose"]:
        log_step(message)


def _append_detail_frame(
    bucket: dict[str, list[pd.DataFrame]],
    endpoint_slug: str,
    fold_df: Optional[pd.DataFrame],
    extra_cols: dict[str, object],
) -> None:
    if fold_df is None:
        return
    detail_df = fold_df.copy()
    for column, value in extra_cols.items():
        detail_df.insert(len(detail_df.columns), column, value)
    bucket.setdefault(endpoint_slug, []).append(detail_df)


def _write_detail_bundles(
    bucket: dict[str, list[pd.DataFrame]],
    outdir: Path,
    prefix: str,
) -> None:
    for endpoint_slug, frames in bucket.items():
        if not frames:
            continue
        pd.concat(frames, ignore_index=True).to_csv(outdir / f"{prefix}__{endpoint_slug}.csv", index=False)


def _metric_key_for_endpoint(endpoint_type: str, cfg: dict) -> str:
    stage2_cfg = cfg.get("search", {}).get("stage2", {})
    if endpoint_type == "binary":
        return str(stage2_cfg.get("ranking_metric_binary", "delta_roc_auc_mean"))
    return str(stage2_cfg.get("ranking_metric_survival", "delta_cindex_mean"))


def _single_feature_feature_set_name(feature: str) -> str:
    return f"clinical_plus__{feature}"


def _subset_name(features: Iterable[str]) -> str:
    return "plus__".join(features)


def _raw_candidate_feature_names(dataset: SearchDataset) -> list[str]:
    if "feature_name" not in dataset.feature_catalog.columns:
        return list(dataset.candidate_features)
    return [str(col) for col in dataset.feature_catalog["feature_name"].tolist() if str(col) in dataset.analysis_df.columns]


def _build_dataset_view(
    base_dataset: SearchDataset,
    analysis_df: pd.DataFrame,
    cfg: dict,
    subset_label: str,
) -> SearchDataset:
    analysis_df = analysis_df.copy().reset_index(drop=True)
    raw_candidate_features = _raw_candidate_feature_names(base_dataset)
    qc_cfg = cfg.get("modeling", {}).get("data_qc", {})
    min_non_missing_fraction = float(qc_cfg.get("min_non_missing_fraction", 0.6))
    min_unique_values = int(qc_cfg.get("min_unique_values", 2))
    candidate_features = drop_bad_columns(
        analysis_df,
        columns=raw_candidate_features,
        min_non_missing_fraction=min_non_missing_fraction,
        min_unique_values=min_unique_values,
    )

    feature_catalog = base_dataset.feature_catalog.copy()
    feature_catalog["non_missing_fraction"] = feature_catalog["feature_name"].map(
        lambda col: float(analysis_df[col].notna().mean()) if col in analysis_df.columns else np.nan
    )
    feature_catalog["n_unique"] = feature_catalog["feature_name"].map(
        lambda col: int(analysis_df[col].dropna().nunique()) if col in analysis_df.columns else 0
    )
    feature_catalog["kept_for_modeling"] = feature_catalog["feature_name"].isin(candidate_features)
    feature_catalog = feature_catalog.sort_values(
        ["kept_for_modeling", "feature_level", "feature_name"],
        ascending=[False, True, True],
    ).reset_index(drop=True)

    manifest = {
        **base_dataset.manifest,
        "subset_label": subset_label,
        "n_matched_patients": int(analysis_df.shape[0]),
        "n_candidate_features_before_qc": int(len(raw_candidate_features)),
        "n_candidate_features_after_qc": int(len(candidate_features)),
        "data_qc": {
            "min_non_missing_fraction": min_non_missing_fraction,
            "min_unique_values": min_unique_values,
        },
    }

    return SearchDataset(
        analysis_df=analysis_df,
        clinical_feature_sets={key: [col for col in value if col in analysis_df.columns] for key, value in base_dataset.clinical_feature_sets.items()},
        candidate_features=candidate_features,
        carry_columns=[col for col in base_dataset.carry_columns if col in analysis_df.columns],
        feature_catalog=feature_catalog,
        manifest=manifest,
    )


def _single_endpoint_cfg(cfg: dict, endpoint_name: str, endpoint_type: str) -> dict:
    cfg_copy = copy.deepcopy(cfg)
    endpoints_cfg = cfg_copy.setdefault("modeling", {}).setdefault("endpoints", {})
    if endpoint_type == "binary":
        endpoints_cfg["binary"] = [spec for spec in endpoints_cfg.get("binary", []) if str(spec.get("name", spec.get("outcome_col", ""))) == endpoint_name]
        endpoints_cfg["survival"] = []
    else:
        endpoints_cfg["survival"] = [spec for spec in endpoints_cfg.get("survival", []) if str(spec.get("name", spec.get("time_col", ""))) == endpoint_name]
        endpoints_cfg["binary"] = []
    return cfg_copy


def _outer_validation_cfg(cfg: dict) -> dict:
    outer_cfg = cfg.get("validation", {}).get("outer_search_cv", {}) or {}
    return {
        "enabled": bool(outer_cfg.get("enabled", False)),
        "n_splits": int(outer_cfg.get("n_splits", 3)),
        "n_repeats": int(outer_cfg.get("n_repeats", 1)),
        "candidate_priority": str(outer_cfg.get("candidate_priority", "stage4_then_stage2_then_baseline")),
        "write_detailed_cv_files": bool(outer_cfg.get("write_detailed_cv_files", False)),
    }


def _parallel_cfg(cfg: dict) -> dict:
    parallel_cfg = cfg.get("parallel", {}) or {}
    cpu_count = os.cpu_count() or 1
    max_workers = int(parallel_cfg.get("max_workers", 1))
    if max_workers <= 0:
        max_workers = cpu_count
    stage2_workers = int(parallel_cfg.get("stage2_workers", max_workers))
    if stage2_workers <= 0:
        stage2_workers = max_workers
    stage4_workers = int(parallel_cfg.get("stage4_workers", max_workers))
    if stage4_workers <= 0:
        stage4_workers = max_workers
    return {
        "enabled": bool(parallel_cfg.get("enabled", False)) and max_workers > 1,
        "max_workers": max(1, min(max_workers, cpu_count)),
        "stage2_workers": max(1, min(stage2_workers, cpu_count)),
        "stage4_workers": max(1, min(stage4_workers, cpu_count)),
    }


def _iter_outer_splits(df: pd.DataFrame, endpoint, endpoint_type: str, cfg: dict) -> list[tuple[int, np.ndarray, np.ndarray]]:
    outer_cfg = _outer_validation_cfg(cfg)
    modeling_cfg = cfg.get("modeling", {})
    random_state = int(modeling_cfg.get("random_state", 42))

    if endpoint_type == "binary":
        y = df[endpoint.outcome_col].astype(int).to_numpy()
        resolved_splits = resolve_stratified_splits(y, outer_cfg["n_splits"])
        if resolved_splits < 2:
            return []
        splitter = RepeatedStratifiedKFold(
            n_splits=resolved_splits,
            n_repeats=max(outer_cfg["n_repeats"], 1),
            random_state=random_state,
        )
        return [(fold_idx, train_idx, test_idx) for fold_idx, (train_idx, test_idx) in enumerate(splitter.split(df, y), start=1)]

    y_event = df[endpoint.event_col].astype(int).to_numpy()
    resolved_splits = resolve_stratified_splits(y_event, outer_cfg["n_splits"])
    if resolved_splits < 2:
        return []

    splits: list[tuple[int, np.ndarray, np.ndarray]] = []
    fold_id = 0
    for repeat_idx in range(max(outer_cfg["n_repeats"], 1)):
        splitter = StratifiedKFold(
            n_splits=resolved_splits,
            shuffle=True,
            random_state=random_state + repeat_idx,
        )
        for train_idx, test_idx in splitter.split(df, y_event):
            fold_id += 1
            splits.append((fold_id, train_idx, test_idx))
    return splits


def _select_outer_candidate(
    baseline_df: pd.DataFrame,
    stage2_df: pd.DataFrame,
    subset_results_df: pd.DataFrame,
    endpoint_name: str,
    endpoint_type: str,
    cfg: dict,
) -> dict[str, object]:
    metric_key = _metric_key_for_endpoint(endpoint_type, cfg)
    priority = _outer_validation_cfg(cfg)["candidate_priority"]

    baseline_row = baseline_df[
        (baseline_df["endpoint_name"] == endpoint_name)
        & (baseline_df["endpoint_type"] == endpoint_type)
    ].iloc[0]

    def best_subset() -> Optional[pd.Series]:
        if subset_results_df.empty:
            return None
        sub = subset_results_df[
            (subset_results_df["endpoint_name"] == endpoint_name)
            & (subset_results_df["endpoint_type"] == endpoint_type)
        ].copy()
        if sub.empty:
            return None
        sub = sub.sort_values([metric_key, "subset_size"], ascending=[False, False])
        return sub.iloc[0]

    def best_single() -> Optional[pd.Series]:
        if stage2_df.empty:
            return None
        sub = stage2_df[
            (stage2_df["endpoint_name"] == endpoint_name)
            & (stage2_df["endpoint_type"] == endpoint_type)
            & (stage2_df["status"] == "ok")
        ].copy()
        if sub.empty:
            return None
        sub = sub.sort_values([metric_key, "full_inner_score"], ascending=[False, False])
        return sub.iloc[0]

    ordered_choices = []
    if priority == "stage2_then_stage4_then_baseline":
        ordered_choices = [("stage2_single", best_single()), ("stage4_subset", best_subset())]
    else:
        ordered_choices = [("stage4_subset", best_subset()), ("stage2_single", best_single())]

    for source, row in ordered_choices:
        if row is None:
            continue
        if source == "stage4_subset":
            feature_names = tuple(str(row["feature_names"]).split("|"))
            return {
                "candidate_source": source,
                "feature_names": feature_names,
                "feature_names_joined": "|".join(feature_names),
                "subset_name": row["subset_name"],
                "subset_size": int(row["subset_size"]),
                "train_ranking_metric": row.get(metric_key, np.nan),
            }
        feature_name = str(row["feature_name"])
        return {
            "candidate_source": source,
            "feature_names": (feature_name,),
            "feature_names_joined": feature_name,
            "subset_name": _single_feature_feature_set_name(feature_name),
            "subset_size": 1,
            "train_ranking_metric": row.get(metric_key, np.nan),
        }

    return {
        "candidate_source": "baseline_only",
        "feature_names": tuple(),
        "feature_names_joined": "",
        "subset_name": "baseline_only",
        "subset_size": 0,
        "train_ranking_metric": baseline_row.get(metric_key, np.nan),
    }


def _score_binary_outer_test(
    train_df: pd.DataFrame,
    test_df: pd.DataFrame,
    feature_cols: list[str],
    clinical_cols: list[str],
    endpoint: BinaryEndpointSpec,
    cfg: dict,
    fold_seed: int,
) -> dict[str, object]:
    modeling = cfg["modeling"]
    logistic_cfg = modeling["logistic"]
    resampling = modeling["resampling"]
    alpha, l1_ratio, inner_auc = tune_logistic_elastic_net(
        df_train=train_df,
        feature_cols=feature_cols,
        clinical_cols=clinical_cols,
        outcome_col=endpoint.outcome_col,
        alpha_grid=list(logistic_cfg["alpha_grid"]),
        l1_ratio_grid=list(logistic_cfg["l1_ratio_grid"]),
        n_inner_splits=int(resampling["n_inner_splits"]),
        random_state=fold_seed,
        max_iter=int(logistic_cfg.get("max_iter", 5000)),
        selection_rule=str(logistic_cfg.get("selection_rule", "best")),
        collinearity_cfg=modeling.get("multicollinearity", {}),
    )
    pred, coef_df = fit_logistic_model(
        train_df=train_df,
        test_df=test_df,
        feature_cols=feature_cols,
        clinical_cols=clinical_cols,
        outcome_col=endpoint.outcome_col,
        alpha=alpha,
        l1_ratio=l1_ratio,
        max_iter=int(logistic_cfg.get("max_iter", 5000)),
        random_state=fold_seed,
        collinearity_cfg=modeling.get("multicollinearity", {}),
    )
    y_train = train_df[endpoint.outcome_col].astype(int).to_numpy()
    y_test = test_df[endpoint.outcome_col].astype(int).to_numpy()
    return {
        "inner_score": inner_auc,
        "best_alpha": alpha,
        "best_l1_ratio": l1_ratio,
        "roc_auc": safe_auc(y_test, pred),
        "auprc": average_precision_score(y_test, pred) if len(np.unique(y_test)) > 1 else np.nan,
        "brier": brier_score_loss(y_test, pred),
        "n_train": int(len(train_df)),
        "n_test": int(len(test_df)),
        "events_train": int(y_train.sum()),
        "events_test": int(y_test.sum()),
        "selected_terms": int(coef_df.shape[0]),
    }


def _score_survival_outer_test(
    train_df: pd.DataFrame,
    test_df: pd.DataFrame,
    feature_cols: list[str],
    clinical_cols: list[str],
    endpoint: SurvivalEndpointSpec,
    cfg: dict,
    fold_seed: int,
) -> dict[str, object]:
    _require_lifelines()
    modeling = cfg["modeling"]
    cox_cfg = modeling["cox"]
    resampling = modeling["resampling"]
    penalizer, l1_ratio, inner_cindex = tune_cox_elastic_net(
        df_train=train_df,
        feature_cols=feature_cols,
        clinical_cols=clinical_cols,
        time_col=endpoint.time_col,
        event_col=endpoint.event_col,
        penalizer_grid=list(cox_cfg["penalizer_grid"]),
        l1_ratio_grid=list(cox_cfg["l1_ratio_grid"]),
        n_inner_splits=int(resampling["n_inner_splits"]),
        random_state=fold_seed,
        clinical_penalty_factor=float(cox_cfg.get("clinical_penalty_factor", 1.0)),
        selection_rule=str(cox_cfg.get("selection_rule", "best")),
        collinearity_cfg=modeling.get("multicollinearity", {}),
    )
    risk, coef_df = fit_cox_model(
        train_df=train_df,
        test_df=test_df,
        feature_cols=feature_cols,
        clinical_cols=clinical_cols,
        time_col=endpoint.time_col,
        event_col=endpoint.event_col,
        penalizer=penalizer,
        l1_ratio=l1_ratio,
        clinical_penalty_factor=float(cox_cfg.get("clinical_penalty_factor", 1.0)),
        collinearity_cfg=modeling.get("multicollinearity", {}),
    )
    return {
        "inner_score": inner_cindex,
        "best_penalizer": penalizer,
        "best_l1_ratio": l1_ratio,
        "cindex": concordance_index(
            test_df[endpoint.time_col].to_numpy(),
            -risk,
            test_df[endpoint.event_col].to_numpy(),
        ),
        "n_train": int(len(train_df)),
        "n_test": int(len(test_df)),
        "events_train": int(train_df[endpoint.event_col].sum()),
        "events_test": int(test_df[endpoint.event_col].sum()),
        "selected_terms": int(coef_df.shape[0]),
    }


def _stage2_binary_feature_task(args: tuple[pd.DataFrame, list[str], BinaryEndpointSpec, dict, dict[str, object], str]) -> tuple[dict[str, object], Optional[pd.DataFrame]]:
    endpoint_df, clinical_cols, endpoint, cfg, baseline_metrics, feature = args
    feature_cols = list(clinical_cols) + [feature]
    metric_key = _metric_key_for_endpoint("binary", cfg)
    try:
        fold_df, summary = _evaluate_binary_analysis(endpoint_df, feature_cols, clinical_cols, endpoint, cfg)
    except Exception as exc:
        return ({
            "endpoint_name": endpoint.name,
            "endpoint_type": "binary",
            "feature_name": feature,
            "feature_set": _single_feature_feature_set_name(feature),
            "n_patients": int(endpoint_df.shape[0]),
            "status": "failed",
            "error": str(exc),
        }, None)

    row = {
        "endpoint_name": endpoint.name,
        "endpoint_type": "binary",
        "feature_name": feature,
        "feature_set": _single_feature_feature_set_name(feature),
        "n_patients": int(endpoint_df.shape[0]),
        "status": "ok",
        **summary,
        "baseline_roc_auc_mean": baseline_metrics.get("roc_auc_mean", np.nan),
        "baseline_auprc_mean": baseline_metrics.get("auprc_mean", np.nan),
        "baseline_brier_mean": baseline_metrics.get("brier_mean", np.nan),
    }
    row["delta_roc_auc_mean"] = row["roc_auc_mean"] - row["baseline_roc_auc_mean"] if pd.notna(row["roc_auc_mean"]) and pd.notna(row["baseline_roc_auc_mean"]) else np.nan
    row["delta_auprc_mean"] = row["auprc_mean"] - row["baseline_auprc_mean"] if pd.notna(row["auprc_mean"]) and pd.notna(row["baseline_auprc_mean"]) else np.nan
    row["delta_brier_mean"] = row["brier_mean"] - row["baseline_brier_mean"] if pd.notna(row["brier_mean"]) and pd.notna(row["baseline_brier_mean"]) else np.nan
    row["ranking_metric"] = row.get(metric_key, np.nan)
    try:
        logistic_cfg = cfg["modeling"]["logistic"]
        resampling = cfg["modeling"]["resampling"]
        alpha, l1_ratio, inner_auc = tune_logistic_elastic_net(
            df_train=endpoint_df,
            feature_cols=feature_cols,
            clinical_cols=clinical_cols,
            outcome_col=endpoint.outcome_col,
            alpha_grid=list(logistic_cfg["alpha_grid"]),
            l1_ratio_grid=list(logistic_cfg["l1_ratio_grid"]),
            n_inner_splits=int(resampling["n_inner_splits"]),
            random_state=int(cfg["modeling"].get("random_state", 42)),
            max_iter=int(logistic_cfg.get("max_iter", 5000)),
            selection_rule=str(logistic_cfg.get("selection_rule", "best")),
            collinearity_cfg=cfg["modeling"].get("multicollinearity", {}),
        )
        _, coef_df = fit_logistic_model(
            train_df=endpoint_df,
            test_df=endpoint_df,
            feature_cols=feature_cols,
            clinical_cols=clinical_cols,
            outcome_col=endpoint.outcome_col,
            alpha=alpha,
            l1_ratio=l1_ratio,
            max_iter=int(logistic_cfg.get("max_iter", 5000)),
            random_state=int(cfg["modeling"].get("random_state", 42)),
            collinearity_cfg=cfg["modeling"].get("multicollinearity", {}),
        )
        feature_coef = coef_df.loc[coef_df["raw_feature"] == feature, "coef"]
        row["full_data_coef"] = float(feature_coef.iloc[0]) if len(feature_coef) else np.nan
        row["full_data_coef_sign"] = np.sign(row["full_data_coef"]) if pd.notna(row["full_data_coef"]) else np.nan
        row["best_alpha"] = alpha
        row["best_l1_ratio"] = l1_ratio
        row["full_inner_score"] = inner_auc
    except Exception:
        row["full_data_coef"] = np.nan
        row["full_data_coef_sign"] = np.nan
        row["best_alpha"] = np.nan
        row["best_l1_ratio"] = np.nan
        row["full_inner_score"] = np.nan
    return row, fold_df


def _stage2_survival_feature_task(args: tuple[pd.DataFrame, list[str], SurvivalEndpointSpec, dict, dict[str, object], str]) -> tuple[dict[str, object], Optional[pd.DataFrame]]:
    endpoint_df, clinical_cols, endpoint, cfg, baseline_metrics, feature = args
    feature_cols = list(clinical_cols) + [feature]
    metric_key = _metric_key_for_endpoint("survival", cfg)
    try:
        fold_df, summary = _evaluate_survival_analysis(endpoint_df, feature_cols, clinical_cols, endpoint, cfg)
    except Exception as exc:
        return ({
            "endpoint_name": endpoint.name,
            "endpoint_type": "survival",
            "feature_name": feature,
            "feature_set": _single_feature_feature_set_name(feature),
            "n_patients": int(endpoint_df.shape[0]),
            "status": "failed",
            "error": str(exc),
        }, None)

    row = {
        "endpoint_name": endpoint.name,
        "endpoint_type": "survival",
        "feature_name": feature,
        "feature_set": _single_feature_feature_set_name(feature),
        "n_patients": int(endpoint_df.shape[0]),
        "status": "ok",
        **summary,
        "baseline_cindex_mean": baseline_metrics.get("cindex_mean", np.nan),
    }
    row["delta_cindex_mean"] = row["cindex_mean"] - row["baseline_cindex_mean"] if pd.notna(row["cindex_mean"]) and pd.notna(row["baseline_cindex_mean"]) else np.nan
    row["ranking_metric"] = row.get(metric_key, np.nan)
    try:
        cox_cfg = cfg["modeling"]["cox"]
        resampling = cfg["modeling"]["resampling"]
        penalizer, l1_ratio, inner_cindex = tune_cox_elastic_net(
            df_train=endpoint_df,
            feature_cols=feature_cols,
            clinical_cols=clinical_cols,
            time_col=endpoint.time_col,
            event_col=endpoint.event_col,
            penalizer_grid=list(cox_cfg["penalizer_grid"]),
            l1_ratio_grid=list(cox_cfg["l1_ratio_grid"]),
            n_inner_splits=int(resampling["n_inner_splits"]),
            random_state=int(cfg["modeling"].get("random_state", 42)),
            clinical_penalty_factor=float(cox_cfg.get("clinical_penalty_factor", 1.0)),
            selection_rule=str(cox_cfg.get("selection_rule", "best")),
            collinearity_cfg=cfg["modeling"].get("multicollinearity", {}),
        )
        _, coef_df = fit_cox_model(
            train_df=endpoint_df,
            test_df=endpoint_df,
            feature_cols=feature_cols,
            clinical_cols=clinical_cols,
            time_col=endpoint.time_col,
            event_col=endpoint.event_col,
            penalizer=penalizer,
            l1_ratio=l1_ratio,
            clinical_penalty_factor=float(cox_cfg.get("clinical_penalty_factor", 1.0)),
            collinearity_cfg=cfg["modeling"].get("multicollinearity", {}),
        )
        feature_coef = coef_df.loc[coef_df["raw_feature"] == feature, "coef"]
        row["full_data_coef"] = float(feature_coef.iloc[0]) if len(feature_coef) else np.nan
        row["full_data_coef_sign"] = np.sign(row["full_data_coef"]) if pd.notna(row["full_data_coef"]) else np.nan
        row["best_penalizer"] = penalizer
        row["best_l1_ratio"] = l1_ratio
        row["full_inner_score"] = inner_cindex
    except Exception:
        row["full_data_coef"] = np.nan
        row["full_data_coef_sign"] = np.nan
        row["best_penalizer"] = np.nan
        row["best_l1_ratio"] = np.nan
        row["full_inner_score"] = np.nan
    return row, fold_df


def _stage4_subset_task(
    args: tuple[pd.DataFrame, list[str], object, str, tuple[str, ...], dict, dict[str, object]],
) -> tuple[dict[str, object], Optional[pd.DataFrame]]:
    endpoint_df, clinical_cols, endpoint, endpoint_type, subset, cfg, baseline_metrics = args
    feature_cols = list(clinical_cols) + list(subset)
    metric_key = _metric_key_for_endpoint(endpoint_type, cfg)

    try:
        if endpoint_type == "binary":
            fold_df, summary = _evaluate_binary_analysis(endpoint_df, feature_cols, clinical_cols, endpoint, cfg)
        else:
            fold_df, summary = _evaluate_survival_analysis(endpoint_df, feature_cols, clinical_cols, endpoint, cfg)
    except Exception as exc:
        return ({
            "endpoint_name": endpoint.name,
            "endpoint_type": endpoint_type,
            "feature_names": "|".join(subset),
            "subset_name": _subset_name(subset),
            "subset_size": len(subset),
            "status": "failed",
            "error": str(exc),
        }, None)

    row = {
        "endpoint_name": endpoint.name,
        "endpoint_type": endpoint_type,
        "feature_names": "|".join(subset),
        "subset_name": _subset_name(subset),
        "subset_size": len(subset),
        "n_patients": int(endpoint_df.shape[0]),
        **summary,
    }
    if endpoint_type == "binary":
        row["delta_roc_auc_mean"] = row["roc_auc_mean"] - baseline_metrics.get("roc_auc_mean", np.nan)
        row["delta_auprc_mean"] = row["auprc_mean"] - baseline_metrics.get("auprc_mean", np.nan)
        row["delta_brier_mean"] = row["brier_mean"] - baseline_metrics.get("brier_mean", np.nan)
    else:
        row["delta_cindex_mean"] = row["cindex_mean"] - baseline_metrics.get("cindex_mean", np.nan)
    row["ranking_metric"] = row.get(metric_key, np.nan)
    return row, fold_df


def _evaluate_binary_analysis(
    df: pd.DataFrame,
    feature_cols: list[str],
    clinical_cols: list[str],
    endpoint: BinaryEndpointSpec,
    cfg: dict,
) -> tuple[pd.DataFrame, dict[str, float]]:
    modeling = cfg["modeling"]
    resampling = modeling["resampling"]
    logistic_cfg = modeling["logistic"]
    return evaluate_logistic_nested_cv(
        df=df,
        feature_cols=feature_cols,
        clinical_cols=clinical_cols,
        outcome_col=endpoint.outcome_col,
        alpha_grid=list(logistic_cfg["alpha_grid"]),
        l1_ratio_grid=list(logistic_cfg["l1_ratio_grid"]),
        n_outer_splits=int(resampling["n_outer_splits"]),
        n_outer_repeats=int(resampling["n_outer_repeats"]),
        n_inner_splits=int(resampling["n_inner_splits"]),
        random_state=int(modeling.get("random_state", 42)),
        max_iter=int(logistic_cfg.get("max_iter", 5000)),
        selection_rule=str(logistic_cfg.get("selection_rule", "best")),
        collinearity_cfg=modeling.get("multicollinearity", {}),
    )


def _evaluate_survival_analysis(
    df: pd.DataFrame,
    feature_cols: list[str],
    clinical_cols: list[str],
    endpoint: SurvivalEndpointSpec,
    cfg: dict,
) -> tuple[pd.DataFrame, dict[str, float]]:
    modeling = cfg["modeling"]
    resampling = modeling["resampling"]
    cox_cfg = modeling["cox"]
    return evaluate_cox_nested_cv(
        df=df,
        feature_cols=feature_cols,
        clinical_cols=clinical_cols,
        time_col=endpoint.time_col,
        event_col=endpoint.event_col,
        penalizer_grid=list(cox_cfg["penalizer_grid"]),
        l1_ratio_grid=list(cox_cfg["l1_ratio_grid"]),
        n_outer_splits=int(resampling["n_outer_splits"]),
        n_outer_repeats=int(resampling["n_outer_repeats"]),
        n_inner_splits=int(resampling["n_inner_splits"]),
        random_state=int(modeling.get("random_state", 42)),
        clinical_penalty_factor=float(cox_cfg.get("clinical_penalty_factor", 1.0)),
        selection_rule=str(cox_cfg.get("selection_rule", "best")),
        collinearity_cfg=modeling.get("multicollinearity", {}),
    )


def _run_baseline_metrics(dataset: SearchDataset, endpoint, endpoint_type: str, baseline_clinical_set: str, cfg: dict) -> dict[str, object]:
    clinical_cols = list(dataset.clinical_feature_sets[baseline_clinical_set])
    if endpoint_type == "binary":
        endpoint_df = dataset.analysis_df.dropna(subset=[endpoint.outcome_col]).copy()
        fold_df, summary = _evaluate_binary_analysis(endpoint_df, clinical_cols, clinical_cols, endpoint, cfg)
    else:
        endpoint_df = dataset.analysis_df.dropna(subset=[endpoint.time_col, endpoint.event_col]).copy()
        fold_df, summary = _evaluate_survival_analysis(endpoint_df, clinical_cols, clinical_cols, endpoint, cfg)
    return {
        "endpoint_name": endpoint.name,
        "endpoint_type": endpoint_type,
        "n_patients": int(endpoint_df.shape[0]),
        "clinical_feature_set": baseline_clinical_set,
        **summary,
        "cv_folds": fold_df,
    }


def _map_tasks(task_fn, tasks: list[tuple], max_workers: int) -> Iterable:
    if max_workers <= 1 or len(tasks) <= 1:
        for task in tasks:
            yield task_fn(task)
        return
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        for result in executor.map(task_fn, tasks):
            yield result


def run_stage2_single_feature_scan(
    dataset: SearchDataset,
    cfg: dict,
    outdir: Path,
    write_fold_detail_files: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    baseline_set = _resolve_baseline_clinical_feature_set(cfg)
    if baseline_set not in dataset.clinical_feature_sets:
        raise ValueError(f"Unknown baseline clinical feature set: {baseline_set}")

    binary_endpoints, survival_endpoints = resolve_endpoint_specs(cfg, dataset.analysis_df.columns.tolist())
    parallel_cfg = _parallel_cfg(cfg)
    logging_cfg = _logging_cfg(cfg)
    baseline_rows: list[dict[str, object]] = []
    scan_rows: list[dict[str, object]] = []
    detail_frames_by_endpoint: dict[str, list[pd.DataFrame]] = {}

    # Baseline per endpoint
    for endpoint in binary_endpoints:
        baseline = _run_baseline_metrics(dataset, endpoint, "binary", baseline_set, cfg)
        baseline_fold_df = baseline.pop("cv_folds")
        if write_fold_detail_files:
            _append_detail_frame(
                detail_frames_by_endpoint,
                endpoint.slug,
                baseline_fold_df,
                {
                    "endpoint_name": endpoint.name,
                    "endpoint_type": "binary",
                    "detail_scope": "baseline",
                    "feature_name": pd.NA,
                    "feature_set": baseline_set,
                    "is_baseline": True,
                },
            )
        baseline_rows.append(baseline)
    for endpoint in survival_endpoints:
        baseline = _run_baseline_metrics(dataset, endpoint, "survival", baseline_set, cfg)
        baseline_fold_df = baseline.pop("cv_folds")
        if write_fold_detail_files:
            _append_detail_frame(
                detail_frames_by_endpoint,
                endpoint.slug,
                baseline_fold_df,
                {
                    "endpoint_name": endpoint.name,
                    "endpoint_type": "survival",
                    "detail_scope": "baseline",
                    "feature_name": pd.NA,
                    "feature_set": baseline_set,
                    "is_baseline": True,
                },
            )
        baseline_rows.append(baseline)

    baseline_df = pd.DataFrame(baseline_rows)
    baseline_df.to_csv(outdir / "stage2_baseline_metrics.csv", index=False)

    # Singles
    for endpoint in binary_endpoints:
        endpoint_df = dataset.analysis_df.dropna(subset=[endpoint.outcome_col]).copy()
        clinical_cols = list(dataset.clinical_feature_sets[baseline_set])
        baseline_row = baseline_df[(baseline_df["endpoint_name"] == endpoint.name) & (baseline_df["endpoint_type"] == "binary")].iloc[0]
        tasks = [
            (endpoint_df, clinical_cols, endpoint, cfg, baseline_row.to_dict(), feature)
            for feature in dataset.candidate_features
        ]
        total_tasks = len(tasks)
        _progress_log(
            cfg,
            f"Stage 2 [{endpoint.slug}] starting single-feature scan for {total_tasks} features "
            f"with stage2_workers={parallel_cfg['stage2_workers']}",
        )
        completed = 0
        for row, fold_df in _map_tasks(_stage2_binary_feature_task, tasks, parallel_cfg["stage2_workers"]):
            completed += 1
            scan_rows.append(row)
            if write_fold_detail_files and fold_df is not None and row.get("status") == "ok":
                _append_detail_frame(
                    detail_frames_by_endpoint,
                    endpoint.slug,
                    fold_df,
                    {
                        "endpoint_name": endpoint.name,
                        "endpoint_type": "binary",
                        "detail_scope": "single_feature",
                        "feature_name": row["feature_name"],
                        "feature_set": row["feature_set"],
                        "is_baseline": False,
                    },
                )
            if completed % logging_cfg["stage_progress_every"] == 0 or completed == total_tasks:
                _progress_log(cfg, f"Stage 2 [{endpoint.slug}] progress: {completed}/{total_tasks} features scored")

    for endpoint in survival_endpoints:
        endpoint_df = dataset.analysis_df.dropna(subset=[endpoint.time_col, endpoint.event_col]).copy()
        clinical_cols = list(dataset.clinical_feature_sets[baseline_set])
        baseline_row = baseline_df[(baseline_df["endpoint_name"] == endpoint.name) & (baseline_df["endpoint_type"] == "survival")].iloc[0]
        tasks = [
            (endpoint_df, clinical_cols, endpoint, cfg, baseline_row.to_dict(), feature)
            for feature in dataset.candidate_features
        ]
        total_tasks = len(tasks)
        _progress_log(
            cfg,
            f"Stage 2 [{endpoint.slug}] starting single-feature scan for {total_tasks} features "
            f"with stage2_workers={parallel_cfg['stage2_workers']}",
        )
        completed = 0
        for row, fold_df in _map_tasks(_stage2_survival_feature_task, tasks, parallel_cfg["stage2_workers"]):
            completed += 1
            scan_rows.append(row)
            if write_fold_detail_files and fold_df is not None and row.get("status") == "ok":
                _append_detail_frame(
                    detail_frames_by_endpoint,
                    endpoint.slug,
                    fold_df,
                    {
                        "endpoint_name": endpoint.name,
                        "endpoint_type": "survival",
                        "detail_scope": "single_feature",
                        "feature_name": row["feature_name"],
                        "feature_set": row["feature_set"],
                        "is_baseline": False,
                    },
                )
            if completed % logging_cfg["stage_progress_every"] == 0 or completed == total_tasks:
                _progress_log(cfg, f"Stage 2 [{endpoint.slug}] progress: {completed}/{total_tasks} features scored")

    scan_df = pd.DataFrame(scan_rows)
    scan_df.to_csv(outdir / "stage2_single_feature_scan.csv", index=False)
    if write_fold_detail_files:
        _write_detail_bundles(detail_frames_by_endpoint, outdir, "stage2_cv_details")
    return baseline_df, scan_df


def _rank_stage2_candidates(stage2_df: pd.DataFrame, endpoint_name: str, endpoint_type: str, top_k: int) -> pd.DataFrame:
    sub = stage2_df[(stage2_df["endpoint_name"] == endpoint_name) & (stage2_df["endpoint_type"] == endpoint_type) & (stage2_df["status"] == "ok")].copy()
    if sub.empty:
        return sub
    sub = sub.sort_values(["ranking_metric", "full_inner_score"], ascending=[False, False])
    return sub.head(top_k).reset_index(drop=True)


def _pairwise_corr_prune(candidate_df: pd.DataFrame, dataset_df: pd.DataFrame, threshold: float) -> tuple[list[str], pd.DataFrame]:
    kept: list[str] = []
    dropped_rows: list[dict[str, object]] = []
    features = candidate_df["feature_name"].tolist()
    rank_map = {feat: idx for idx, feat in enumerate(features)}
    values = dataset_df[features].copy()
    values = values.apply(pd.to_numeric, errors="coerce")
    corr = values.corr().abs()
    dropped: set[str] = set()
    for i, left in enumerate(features):
        if left in dropped:
            continue
        kept.append(left)
        for right in features[i + 1:]:
            if right in dropped:
                continue
            corr_val = corr.loc[left, right]
            if pd.notna(corr_val) and corr_val >= threshold:
                dropped.add(right)
                dropped_rows.append({
                    "kept_feature": left,
                    "dropped_feature": right,
                    "abs_correlation": float(corr_val),
                    "reason": "high_correlation",
                })
    return kept, pd.DataFrame(dropped_rows)


def _compute_vif_trace(dataset_df: pd.DataFrame, features: list[str], vif_max: float) -> tuple[list[str], pd.DataFrame]:
    from sklearn.linear_model import LinearRegression

    if len(features) <= 1:
        return features, pd.DataFrame(columns=["step", "feature", "vif", "action"])

    current = dataset_df[features].apply(pd.to_numeric, errors="coerce").copy()
    for col in current.columns:
        current[col] = current[col].fillna(current[col].median())

    def vif_series(df: pd.DataFrame) -> pd.Series:
        scores = {}
        x = df.to_numpy(dtype=float)
        for idx, col in enumerate(df.columns):
            y = x[:, idx]
            x_other = np.delete(x, idx, axis=1)
            if x_other.shape[1] == 0 or np.std(y) <= 1e-12:
                scores[col] = 1.0
                continue
            model = LinearRegression().fit(x_other, y)
            r2 = model.score(x_other, y)
            r2 = min(max(float(r2), 0.0), 0.999999)
            scores[col] = 1.0 / max(1.0 - r2, 1e-6)
        return pd.Series(scores)

    trace = []
    step = 0
    while current.shape[1] > 1:
        step += 1
        vif = vif_series(current)
        max_feature = str(vif.idxmax())
        max_vif = float(vif.max())
        trace.append({"step": step, "feature": max_feature, "vif": max_vif, "action": "keep" if max_vif <= vif_max else "drop"})
        if max_vif <= vif_max:
            break
        current = current.drop(columns=[max_feature])
    return list(current.columns), pd.DataFrame(trace)


def run_stage3_redundancy_pruning(dataset: SearchDataset, stage2_df: pd.DataFrame, cfg: dict, outdir: Path) -> pd.DataFrame:
    stage3_cfg = cfg["search"]["stage3"]
    top_k = int(stage3_cfg.get("top_k_per_endpoint", 12))
    max_per_family = int(stage3_cfg.get("max_per_family", 2))
    corr_threshold = float(stage3_cfg.get("correlation_threshold", 0.8))
    vif_max = float(stage3_cfg.get("vif_max", 10.0))

    feature_catalog = dataset.feature_catalog.copy()
    pruned_rows: list[dict[str, object]] = []
    redundancy_rows: list[pd.DataFrame] = []
    vif_rows: list[pd.DataFrame] = []
    family_cap_rows: list[dict[str, object]] = []

    for (endpoint_name, endpoint_type), sub in stage2_df[stage2_df["status"] == "ok"].groupby(["endpoint_name", "endpoint_type"]):
        ranked = _rank_stage2_candidates(stage2_df, endpoint_name, endpoint_type, max(top_k * 3, top_k))
        if ranked.empty:
            continue
        ranked = ranked.merge(feature_catalog, left_on="feature_name", right_on="feature_name", how="left")
        ranked["family"] = ranked["family"].fillna(ranked["feature_name"])
        ranked["family_cap_group"] = ranked.get("family_cap_group", pd.Series(pd.NA, index=ranked.index, dtype="object"))
        ranked["family_cap_applies"] = ranked.get("family_cap_applies", pd.Series(False, index=ranked.index, dtype=bool))
        ranked["family_cap_applies"] = ranked["family_cap_applies"].fillna(False).astype(bool)
        ranked["feature_level"] = ranked["feature_level"].fillna("patient_feature")

        # family cap
        family_counts: dict[str, int] = {}
        family_kept: list[str] = []
        for rank_idx, row in enumerate(ranked.itertuples(index=False), start=1):
            fam = str(getattr(row, "family"))
            cap_group = getattr(row, "family_cap_group")
            cap_applies = bool(getattr(row, "family_cap_applies"))
            if not cap_applies or pd.isna(cap_group):
                family_kept.append(row.feature_name)
                family_cap_rows.append({
                    "endpoint_name": endpoint_name,
                    "endpoint_type": endpoint_type,
                    "rank_in_stage2_window": rank_idx,
                    "feature_name": row.feature_name,
                    "family": fam,
                    "family_cap_group": pd.NA,
                    "family_cap_applies": False,
                    "family_count_before": pd.NA,
                    "decision": "keep_no_family_cap",
                    "ranking_metric": row.ranking_metric,
                })
                continue

            count = family_counts.get(str(cap_group), 0)
            decision = "keep_within_family_cap" if count < max_per_family else "drop_by_family_cap"
            family_cap_rows.append({
                "endpoint_name": endpoint_name,
                "endpoint_type": endpoint_type,
                "rank_in_stage2_window": rank_idx,
                "feature_name": row.feature_name,
                "family": fam,
                "family_cap_group": str(cap_group),
                "family_cap_applies": True,
                "family_count_before": count,
                "decision": decision,
                "ranking_metric": row.ranking_metric,
            })
            if count < max_per_family:
                family_counts[str(cap_group)] = count + 1
                family_kept.append(row.feature_name)
        family_df = ranked[ranked["feature_name"].isin(family_kept)].copy()
        if family_df.empty:
            continue

        # corr prune
        corr_kept, corr_diag = _pairwise_corr_prune(family_df, dataset.analysis_df, corr_threshold)
        corr_df = family_df[family_df["feature_name"].isin(corr_kept)].copy()
        if not corr_diag.empty:
            corr_diag.insert(0, "endpoint_name", endpoint_name)
            corr_diag.insert(1, "endpoint_type", endpoint_type)
            redundancy_rows.append(corr_diag)

        # vif prune
        vif_kept, vif_trace = _compute_vif_trace(dataset.analysis_df, corr_df["feature_name"].tolist(), vif_max)
        if not vif_trace.empty:
            vif_trace.insert(0, "endpoint_name", endpoint_name)
            vif_trace.insert(1, "endpoint_type", endpoint_type)
            vif_rows.append(vif_trace)
        final_df = corr_df[corr_df["feature_name"].isin(vif_kept)].copy().head(top_k)
        for row in final_df.itertuples(index=False):
            pruned_rows.append({
                "endpoint_name": endpoint_name,
                "endpoint_type": endpoint_type,
                "feature_name": row.feature_name,
                "family": getattr(row, "family", row.feature_name),
                "family_cap_group": getattr(row, "family_cap_group", pd.NA),
                "family_cap_applies": bool(getattr(row, "family_cap_applies", False)),
                "feature_level": getattr(row, "feature_level", "patient_feature"),
                "ranking_metric": row.ranking_metric,
                "full_data_coef": row.full_data_coef,
            })

    pruned_df = pd.DataFrame(pruned_rows)
    family_cap_df = pd.DataFrame(family_cap_rows)
    pruned_df.to_csv(outdir / "stage3_pruned_candidates.csv", index=False)
    pd.concat(redundancy_rows, ignore_index=True).to_csv(outdir / "stage3_redundancy_pairs.csv", index=False) if redundancy_rows else pd.DataFrame(columns=["endpoint_name","endpoint_type","kept_feature","dropped_feature","abs_correlation","reason"]).to_csv(outdir / "stage3_redundancy_pairs.csv", index=False)
    pd.concat(vif_rows, ignore_index=True).to_csv(outdir / "stage3_vif_trace.csv", index=False) if vif_rows else pd.DataFrame(columns=["endpoint_name","endpoint_type","step","feature","vif","action"]).to_csv(outdir / "stage3_vif_trace.csv", index=False)
    family_cap_df.to_csv(outdir / "stage3_family_cap_trace.csv", index=False) if not family_cap_df.empty else pd.DataFrame(
        columns=[
            "endpoint_name",
            "endpoint_type",
            "rank_in_stage2_window",
            "feature_name",
            "family",
            "family_cap_group",
            "family_cap_applies",
            "family_count_before",
            "decision",
            "ranking_metric",
        ]
    ).to_csv(outdir / "stage3_family_cap_trace.csv", index=False)
    return pruned_df


def _evaluate_subset(dataset: SearchDataset, endpoint, endpoint_type: str, baseline_set: str, subset: tuple[str, ...], cfg: dict) -> dict[str, object]:
    clinical_cols = list(dataset.clinical_feature_sets[baseline_set])
    feature_cols = clinical_cols + list(subset)
    if endpoint_type == "binary":
        endpoint_df = dataset.analysis_df.dropna(subset=[endpoint.outcome_col]).copy()
        fold_df, summary = _evaluate_binary_analysis(endpoint_df, feature_cols, clinical_cols, endpoint, cfg)
    else:
        endpoint_df = dataset.analysis_df.dropna(subset=[endpoint.time_col, endpoint.event_col]).copy()
        fold_df, summary = _evaluate_survival_analysis(endpoint_df, feature_cols, clinical_cols, endpoint, cfg)
    return {
        "subset": subset,
        "subset_name": _subset_name(subset),
        "subset_size": len(subset),
        "n_patients": int(endpoint_df.shape[0]),
        **summary,
        "fold_df": fold_df,
    }


def run_stage4_subset_search(
    dataset: SearchDataset,
    pruned_df: pd.DataFrame,
    baseline_df: pd.DataFrame,
    cfg: dict,
    outdir: Path,
    write_fold_detail_files: bool = True,
) -> pd.DataFrame:
    stage4_cfg = cfg["search"]["stage4"]
    parallel_cfg = _parallel_cfg(cfg)
    logging_cfg = _logging_cfg(cfg)
    baseline_set = _resolve_baseline_clinical_feature_set(cfg)
    exhaustive_sizes = sorted({int(x) for x in stage4_cfg.get("exhaustive_subset_sizes", [2, 3]) if int(x) >= 2})
    heuristic_sizes = sorted({int(x) for x in stage4_cfg.get("heuristic_subset_sizes", [4]) if int(x) >= 4})
    beam_width = int(stage4_cfg.get("beam_width", 20))
    max_expansions_per_parent = int(stage4_cfg.get("max_expansions_per_parent", 20))
    max_total_evaluations = int(stage4_cfg.get("max_total_evaluations", 200))

    binary_endpoints, survival_endpoints = resolve_endpoint_specs(cfg, dataset.analysis_df.columns.tolist())
    endpoint_lookup = {(ep.name, "binary"): ep for ep in binary_endpoints} | {(ep.name, "survival"): ep for ep in survival_endpoints}
    result_rows: list[dict[str, object]] = []
    detail_frames_by_endpoint: dict[str, list[pd.DataFrame]] = {}

    for (endpoint_name, endpoint_type), sub in pruned_df.groupby(["endpoint_name", "endpoint_type"]):
        endpoint = endpoint_lookup[(endpoint_name, endpoint_type)]
        candidates = sub.sort_values("ranking_metric", ascending=False)["feature_name"].tolist()
        if len(candidates) < 2:
            continue
        metric_key = _metric_key_for_endpoint(endpoint_type, cfg)
        baseline_row = baseline_df[(baseline_df["endpoint_name"] == endpoint_name) & (baseline_df["endpoint_type"] == endpoint_type)].iloc[0]
        clinical_cols = list(dataset.clinical_feature_sets[baseline_set])
        if endpoint_type == "binary":
            endpoint_df = dataset.analysis_df.dropna(subset=[endpoint.outcome_col]).copy()
        else:
            endpoint_df = dataset.analysis_df.dropna(subset=[endpoint.time_col, endpoint.event_col]).copy()
        seen: set[tuple[str, ...]] = set()
        cache: dict[tuple[str, ...], dict[str, object]] = {}

        def score_subsets(subsets: list[tuple[str, ...]], search_method: str) -> list[dict[str, object]]:
            fresh_subsets = [subset for subset in subsets if subset not in cache]
            if fresh_subsets:
                _progress_log(
                    cfg,
                    f"Stage 4 [{endpoint.slug}] scoring {len(fresh_subsets)} {search_method} subsets "
                    f"with stage4_workers={parallel_cfg['stage4_workers']}",
                )
                tasks = [
                    (endpoint_df, clinical_cols, endpoint, endpoint_type, subset, cfg, baseline_row.to_dict())
                    for subset in fresh_subsets
                ]
                completed = 0
                total = len(fresh_subsets)
                for row, fold_df in _map_tasks(_stage4_subset_task, tasks, parallel_cfg["stage4_workers"]):
                    completed += 1
                    subset = tuple(str(row["feature_names"]).split("|"))
                    row["search_method"] = search_method
                    cache[subset] = row | {"fold_df": fold_df}
                    if write_fold_detail_files and fold_df is not None and row.get("status", "ok") == "ok":
                        _append_detail_frame(
                            detail_frames_by_endpoint,
                            endpoint.slug,
                            fold_df,
                            {
                                "endpoint_name": endpoint.name,
                                "endpoint_type": endpoint_type,
                                "feature_names": row["feature_names"],
                                "subset_name": row["subset_name"],
                                "subset_size": row["subset_size"],
                                "search_method": search_method,
                            },
                        )
                    if completed % logging_cfg["stage_progress_every"] == 0 or completed == total:
                        _progress_log(cfg, f"Stage 4 [{endpoint.slug}] {search_method} progress: {completed}/{total} subsets scored")
            return [{k: v for k, v in cache[subset].items() if k != "fold_df"} for subset in subsets if subset in cache]

        # exhaustive
        for size in exhaustive_sizes:
            if len(candidates) < size:
                continue
            subsets = []
            for subset in itertools.combinations(candidates, size):
                seen.add(subset)
                subsets.append(subset)
            result_rows.extend(score_subsets(subsets, search_method="exhaustive"))

        # heuristic beam for 4+
        previous_level = [tuple([feat]) for feat in candidates[:beam_width]]
        # seed from best exhaustive 3 if available
        if exhaustive_sizes:
            max_seed_size = max(exhaustive_sizes)
            seed_rows = [r for r in result_rows if r["endpoint_name"] == endpoint_name and r["endpoint_type"] == endpoint_type and r["subset_size"] == max_seed_size]
            if seed_rows:
                seed_rows = sorted(seed_rows, key=lambda x: (pd.notna(x.get(metric_key, np.nan)), x.get(metric_key, -np.inf)), reverse=True)[:beam_width]
                previous_level = [tuple(r["feature_names"].split("|")) for r in seed_rows]
        total_heuristic_eval = 0
        for size in heuristic_sizes:
            expanded: list[dict[str, object]] = []
            for parent in previous_level[:beam_width]:
                remaining = [feat for feat in candidates if feat not in parent][:max_expansions_per_parent]
                for feat in remaining:
                    subset = tuple(sorted(parent + (feat,)))
                    if len(subset) != size or subset in seen:
                        continue
                    seen.add(subset)
                    expanded.append(subset)
                    total_heuristic_eval += 1
                    if total_heuristic_eval >= max_total_evaluations:
                        break
                if total_heuristic_eval >= max_total_evaluations:
                    break
            if not expanded:
                break
            expanded_rows = score_subsets(expanded, search_method="heuristic")
            expanded_rows = sorted(expanded_rows, key=lambda x: x.get(metric_key, -np.inf), reverse=True)
            result_rows.extend(expanded_rows)
            previous_level = [tuple(r["feature_names"].split("|")) for r in expanded_rows[:beam_width]]
            if total_heuristic_eval >= max_total_evaluations:
                break

        # save fold files for best subsets of this endpoint
        endpoint_results = pd.DataFrame([r for r in result_rows if r["endpoint_name"] == endpoint_name and r["endpoint_type"] == endpoint_type])
        if not endpoint_results.empty:
            endpoint_results.sort_values(metric_key, ascending=False).head(10).to_csv(outdir / f"stage4_top10__{endpoint_type}__{endpoint_name}.csv", index=False)

    results_df = pd.DataFrame(result_rows)
    if not results_df.empty:
        results_df = results_df.sort_values(["endpoint_name", "endpoint_type", "ranking_metric"], ascending=[True, True, False])
    results_df.to_csv(outdir / "stage4_subset_results.csv", index=False)
    if write_fold_detail_files:
        _write_detail_bundles(detail_frames_by_endpoint, outdir, "stage4_cv_details")
    return results_df


def run_feature_search_inner(
    dataset: SearchDataset,
    cfg: dict,
    outdir: Path,
    *,
    run_stage5: bool = True,
    write_fold_detail_files: bool = True,
    stage_prefix: str = "feature_search",
) -> dict[str, pd.DataFrame]:
    outdir = ensure_outdir(outdir)
    write_yaml(cfg, outdir / "config.yaml")
    dataset.analysis_df.to_csv(outdir / "merged_feature_search_dataset.tsv", sep="\t", index=False)
    dataset.feature_catalog.to_csv(outdir / "feature_catalog_resolved.csv", index=False)
    write_json(dataset.manifest, outdir / "input_manifest.json")

    with _tracked_stage(
        f"{stage_prefix}.stage2_single_feature_scan",
        {"outdir": str(outdir), "n_patients": int(dataset.analysis_df.shape[0]), "n_candidate_features": int(len(dataset.candidate_features))},
    ):
        baseline_df, stage2_df = run_stage2_single_feature_scan(
            dataset,
            cfg,
            outdir,
            write_fold_detail_files=write_fold_detail_files,
        )
    with _tracked_stage(
        f"{stage_prefix}.stage3_redundancy_pruning",
        {"outdir": str(outdir), "n_stage2_rows": int(stage2_df.shape[0])},
    ):
        pruned_df = run_stage3_redundancy_pruning(dataset, stage2_df, cfg, outdir)
    with _tracked_stage(
        f"{stage_prefix}.stage4_subset_search",
        {"outdir": str(outdir), "n_stage3_candidates": int(pruned_df.shape[0])},
    ):
        subset_results_df = run_stage4_subset_search(
            dataset,
            pruned_df,
            baseline_df,
            cfg,
            outdir,
            write_fold_detail_files=write_fold_detail_files,
        )

    best_df = pd.DataFrame()
    coef_df = pd.DataFrame()
    coef_stability_df = pd.DataFrame()
    if run_stage5:
        with _tracked_stage(
            f"{stage_prefix}.stage5_best_subset_refits",
            {"outdir": str(outdir), "n_stage4_subsets": int(subset_results_df.shape[0])},
        ):
            best_df, coef_df, coef_stability_df = run_stage5_best_subset_refits(dataset, baseline_df, subset_results_df, cfg, outdir)

    manifest = {
        **dataset.manifest,
        "n_candidate_features": int(len(dataset.candidate_features)),
        "n_stage2_rows": int(stage2_df.shape[0]),
        "n_stage3_candidates": int(pruned_df.shape[0]),
        "n_stage4_subsets": int(subset_results_df.shape[0]),
        "n_stage5_best_rows": int(best_df.shape[0]),
    }
    write_json(manifest, outdir / "run_manifest.json")
    return {
        "baseline_df": baseline_df,
        "stage2_df": stage2_df,
        "pruned_df": pruned_df,
        "subset_results_df": subset_results_df,
        "best_df": best_df,
        "coef_df": coef_df,
        "coef_stability_df": coef_stability_df,
    }


def _write_search_inputs(dataset: SearchDataset, cfg: dict, outdir: Path) -> None:
    outdir = ensure_outdir(outdir)
    write_yaml(cfg, outdir / "config.yaml")
    dataset.analysis_df.to_csv(outdir / "merged_feature_search_dataset.tsv", sep="\t", index=False)
    dataset.feature_catalog.to_csv(outdir / "feature_catalog_resolved.csv", index=False)
    write_json(dataset.manifest, outdir / "input_manifest.json")


def _filter_to_endpoint_keys(df: pd.DataFrame, endpoint_keys: set[tuple[str, str]]) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    keep_mask = [(str(row.endpoint_name), str(row.endpoint_type)) in endpoint_keys for row in df.itertuples(index=False)]
    return df.loc[keep_mask].copy()


def _load_reused_stage2_results(
    dataset: SearchDataset,
    cfg: dict,
    stage2_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    baseline_path = stage2_dir / "stage2_baseline_metrics.csv"
    stage2_path = stage2_dir / "stage2_single_feature_scan.csv"
    if not baseline_path.exists():
        raise FileNotFoundError(f"Missing reusable stage2 baseline file: {baseline_path}")
    if not stage2_path.exists():
        raise FileNotFoundError(f"Missing reusable stage2 scan file: {stage2_path}")

    baseline_df = pd.read_csv(baseline_path)
    stage2_df = pd.read_csv(stage2_path)
    required_baseline_cols = {"endpoint_name", "endpoint_type", "clinical_feature_set"}
    required_stage2_cols = {"endpoint_name", "endpoint_type", "feature_name", "status", "ranking_metric"}
    missing_baseline = required_baseline_cols.difference(baseline_df.columns)
    missing_stage2 = required_stage2_cols.difference(stage2_df.columns)
    if missing_baseline:
        raise ValueError(f"Reusable stage2 baseline file is missing required columns: {sorted(missing_baseline)}")
    if missing_stage2:
        raise ValueError(f"Reusable stage2 scan file is missing required columns: {sorted(missing_stage2)}")

    baseline_set = _resolve_baseline_clinical_feature_set(cfg)
    observed_sets = {str(value) for value in baseline_df["clinical_feature_set"].dropna().astype(str).unique().tolist()}
    if observed_sets and observed_sets != {baseline_set}:
        raise ValueError(
            f"Reusable stage2 results were computed with clinical feature set(s) {sorted(observed_sets)}, "
            f"but current config resolves baseline set '{baseline_set}'. Re-run stage2 for this baseline."
        )

    binary_endpoints, survival_endpoints = resolve_endpoint_specs(cfg, dataset.analysis_df.columns.tolist())
    endpoint_keys = {
        *((ep.name, "binary") for ep in binary_endpoints),
        *((ep.name, "survival") for ep in survival_endpoints),
    }
    baseline_df = _filter_to_endpoint_keys(baseline_df, endpoint_keys)
    stage2_df = _filter_to_endpoint_keys(stage2_df, endpoint_keys)

    if baseline_df.empty:
        raise ValueError("No reusable stage2 baseline rows remain after filtering to the current endpoint config.")
    if stage2_df.empty:
        raise ValueError("No reusable stage2 scan rows remain after filtering to the current endpoint config.")

    candidate_set = set(dataset.candidate_features)
    before_rows = int(stage2_df.shape[0])
    stage2_df = stage2_df[stage2_df["feature_name"].isin(candidate_set)].copy()
    dropped_rows = before_rows - int(stage2_df.shape[0])
    if stage2_df.empty:
        raise ValueError("No reusable stage2 scan rows remain after aligning to current candidate features.")
    if dropped_rows > 0 and _ACTIVE_RUN_TRACKER is not None:
        _ACTIVE_RUN_TRACKER.log(
            f"Reusable stage2 filtering dropped {dropped_rows} rows because the feature is not in the current candidate set.",
            verbose_only=True,
        )

    return baseline_df.reset_index(drop=True), stage2_df.reset_index(drop=True)


def _coefficient_stability_binary(df: pd.DataFrame, feature_cols: list[str], clinical_cols: list[str], endpoint: BinaryEndpointSpec, cfg: dict) -> pd.DataFrame:
    modeling = cfg["modeling"]
    resampling = modeling["resampling"]
    logistic_cfg = modeling["logistic"]
    from sklearn.model_selection import RepeatedStratifiedKFold

    y = df[endpoint.outcome_col].astype(int).to_numpy()
    splitter = RepeatedStratifiedKFold(
        n_splits=int(resampling["n_outer_splits"]),
        n_repeats=int(resampling["n_outer_repeats"]),
        random_state=int(modeling.get("random_state", 42)),
    )
    coef_rows = []
    for fold_idx, (train_idx, test_idx) in enumerate(splitter.split(df, y), start=1):
        train_df = df.iloc[train_idx].copy()
        alpha, l1_ratio, _ = tune_logistic_elastic_net(
            df_train=train_df,
            feature_cols=feature_cols,
            clinical_cols=clinical_cols,
            outcome_col=endpoint.outcome_col,
            alpha_grid=list(logistic_cfg["alpha_grid"]),
            l1_ratio_grid=list(logistic_cfg["l1_ratio_grid"]),
            n_inner_splits=int(resampling["n_inner_splits"]),
            random_state=int(modeling.get("random_state", 42)) + fold_idx,
            max_iter=int(logistic_cfg.get("max_iter", 5000)),
            selection_rule=str(logistic_cfg.get("selection_rule", "best")),
            collinearity_cfg=modeling.get("multicollinearity", {}),
        )
        _, coef_df = fit_logistic_model(
            train_df=train_df,
            test_df=train_df,
            feature_cols=feature_cols,
            clinical_cols=clinical_cols,
            outcome_col=endpoint.outcome_col,
            alpha=alpha,
            l1_ratio=l1_ratio,
            max_iter=int(logistic_cfg.get("max_iter", 5000)),
            random_state=int(modeling.get("random_state", 42)) + fold_idx,
            collinearity_cfg=modeling.get("multicollinearity", {}),
        )
        coef_df = coef_df[["raw_feature", "coef", "feature_scope"]].copy()
        coef_df["fold"] = fold_idx
        coef_rows.append(coef_df)
    coef_all = pd.concat(coef_rows, ignore_index=True)
    summary = coef_all.groupby(["raw_feature", "feature_scope"], as_index=False).agg(
        coef_mean=("coef", "mean"),
        coef_sd=("coef", "std"),
        sign_consistency=("coef", lambda s: float(np.mean(np.sign(s.fillna(0.0)) == np.sign(np.nanmean(s))))),
        n_folds=("coef", "size"),
    )
    return summary


def _coefficient_stability_survival(df: pd.DataFrame, feature_cols: list[str], clinical_cols: list[str], endpoint: SurvivalEndpointSpec, cfg: dict) -> pd.DataFrame:
    modeling = cfg["modeling"]
    resampling = modeling["resampling"]
    cox_cfg = modeling["cox"]
    from sklearn.model_selection import StratifiedKFold

    y = df[endpoint.event_col].astype(int).to_numpy()
    coef_rows = []
    for repeat_idx in range(int(resampling["n_outer_repeats"])):
        splitter = StratifiedKFold(
            n_splits=int(resampling["n_outer_splits"]),
            shuffle=True,
            random_state=int(modeling.get("random_state", 42)) + repeat_idx,
        )
        for fold_idx, (train_idx, test_idx) in enumerate(splitter.split(df, y), start=1):
            global_fold = repeat_idx * int(resampling["n_outer_splits"]) + fold_idx
            train_df = df.iloc[train_idx].copy()
            penalizer, l1_ratio, _ = tune_cox_elastic_net(
                df_train=train_df,
                feature_cols=feature_cols,
                clinical_cols=clinical_cols,
                time_col=endpoint.time_col,
                event_col=endpoint.event_col,
                penalizer_grid=list(cox_cfg["penalizer_grid"]),
                l1_ratio_grid=list(cox_cfg["l1_ratio_grid"]),
                n_inner_splits=int(resampling["n_inner_splits"]),
                random_state=int(modeling.get("random_state", 42)) + global_fold,
                clinical_penalty_factor=float(cox_cfg.get("clinical_penalty_factor", 1.0)),
                selection_rule=str(cox_cfg.get("selection_rule", "best")),
                collinearity_cfg=modeling.get("multicollinearity", {}),
            )
            _, coef_df = fit_cox_model(
                train_df=train_df,
                test_df=train_df,
                feature_cols=feature_cols,
                clinical_cols=clinical_cols,
                time_col=endpoint.time_col,
                event_col=endpoint.event_col,
                penalizer=penalizer,
                l1_ratio=l1_ratio,
                clinical_penalty_factor=float(cox_cfg.get("clinical_penalty_factor", 1.0)),
                collinearity_cfg=modeling.get("multicollinearity", {}),
            )
            coef_df = coef_df[["raw_feature", "coef", "feature_scope"]].copy()
            coef_df["fold"] = global_fold
            coef_rows.append(coef_df)
    coef_all = pd.concat(coef_rows, ignore_index=True)
    summary = coef_all.groupby(["raw_feature", "feature_scope"], as_index=False).agg(
        coef_mean=("coef", "mean"),
        coef_sd=("coef", "std"),
        sign_consistency=("coef", lambda s: float(np.mean(np.sign(s.fillna(0.0)) == np.sign(np.nanmean(s))))),
        n_folds=("coef", "size"),
    )
    return summary


def run_stage5_best_subset_refits(dataset: SearchDataset, baseline_df: pd.DataFrame, subset_results_df: pd.DataFrame, cfg: dict, outdir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if subset_results_df.empty:
        empty = pd.DataFrame()
        empty.to_csv(outdir / "stage5_best_subsets.csv", index=False)
        empty.to_csv(outdir / "stage5_best_subset_coefficients.csv", index=False)
        empty.to_csv(outdir / "stage5_best_subset_coefficient_stability.csv", index=False)
        empty.to_csv(outdir / "stage5_best_subset_bootstrap_stability.csv", index=False)
        return empty, empty, empty

    stage5_cfg = cfg.get("search", {}).get("stage5", {})
    top_n = int(stage5_cfg.get("top_n_best_subsets_per_endpoint", 3))
    baseline_set = _resolve_baseline_clinical_feature_set(cfg)
    binary_endpoints, survival_endpoints = resolve_endpoint_specs(cfg, dataset.analysis_df.columns.tolist())
    endpoint_lookup = {(ep.name, "binary"): ep for ep in binary_endpoints} | {(ep.name, "survival"): ep for ep in survival_endpoints}

    best_rows = []
    coef_frames = []
    coef_stability_frames = []
    bootstrap_frames = []

    for (endpoint_name, endpoint_type), sub in subset_results_df.groupby(["endpoint_name", "endpoint_type"]):
        metric_key = _metric_key_for_endpoint(endpoint_type, cfg)
        top_sub = sub.sort_values(metric_key, ascending=False).head(top_n)
        endpoint = endpoint_lookup[(endpoint_name, endpoint_type)]
        for row in top_sub.itertuples(index=False):
            subset = tuple(str(row.feature_names).split("|"))
            clinical_cols = list(dataset.clinical_feature_sets[baseline_set])
            feature_cols = clinical_cols + list(subset)
            if endpoint_type == "binary":
                endpoint_df = dataset.analysis_df.dropna(subset=[endpoint.outcome_col]).copy()
                logistic_cfg = cfg["modeling"]["logistic"]
                resampling = cfg["modeling"]["resampling"]
                alpha, l1_ratio, inner_auc = tune_logistic_elastic_net(
                    df_train=endpoint_df,
                    feature_cols=feature_cols,
                    clinical_cols=clinical_cols,
                    outcome_col=endpoint.outcome_col,
                    alpha_grid=list(logistic_cfg["alpha_grid"]),
                    l1_ratio_grid=list(logistic_cfg["l1_ratio_grid"]),
                    n_inner_splits=int(resampling["n_inner_splits"]),
                    random_state=int(cfg["modeling"].get("random_state", 42)),
                    max_iter=int(logistic_cfg.get("max_iter", 5000)),
                    selection_rule=str(logistic_cfg.get("selection_rule", "best")),
                    collinearity_cfg=cfg["modeling"].get("multicollinearity", {}),
                )
                _, coef_df = fit_logistic_model(
                    train_df=endpoint_df,
                    test_df=endpoint_df,
                    feature_cols=feature_cols,
                    clinical_cols=clinical_cols,
                    outcome_col=endpoint.outcome_col,
                    alpha=alpha,
                    l1_ratio=l1_ratio,
                    max_iter=int(logistic_cfg.get("max_iter", 5000)),
                    random_state=int(cfg["modeling"].get("random_state", 42)),
                    collinearity_cfg=cfg["modeling"].get("multicollinearity", {}),
                )
                stability_df = bootstrap_stability_logistic(
                    df=endpoint_df,
                    feature_cols=feature_cols,
                    clinical_cols=clinical_cols,
                    outcome_col=endpoint.outcome_col,
                    alpha=alpha,
                    l1_ratio=l1_ratio,
                    n_bootstrap=int(cfg["modeling"]["stability"]["n_bootstrap"]),
                    selection_threshold=float(cfg["modeling"]["stability"]["selection_threshold"]),
                    random_state=int(cfg["modeling"].get("random_state", 42)),
                    max_iter=int(logistic_cfg.get("max_iter", 5000)),
                    selection_scope=str(cfg["modeling"]["stability"].get("selection_scope", "all")),
                    collinearity_cfg=cfg["modeling"].get("multicollinearity", {}),
                )
                coef_stab_df = _coefficient_stability_binary(endpoint_df, feature_cols, clinical_cols, endpoint, cfg)
                extra = {"best_alpha": alpha, "best_l1_ratio": l1_ratio, "full_inner_score": inner_auc}
            else:
                endpoint_df = dataset.analysis_df.dropna(subset=[endpoint.time_col, endpoint.event_col]).copy()
                cox_cfg = cfg["modeling"]["cox"]
                resampling = cfg["modeling"]["resampling"]
                penalizer, l1_ratio, inner_cindex = tune_cox_elastic_net(
                    df_train=endpoint_df,
                    feature_cols=feature_cols,
                    clinical_cols=clinical_cols,
                    time_col=endpoint.time_col,
                    event_col=endpoint.event_col,
                    penalizer_grid=list(cox_cfg["penalizer_grid"]),
                    l1_ratio_grid=list(cox_cfg["l1_ratio_grid"]),
                    n_inner_splits=int(resampling["n_inner_splits"]),
                    random_state=int(cfg["modeling"].get("random_state", 42)),
                    clinical_penalty_factor=float(cox_cfg.get("clinical_penalty_factor", 1.0)),
                    selection_rule=str(cox_cfg.get("selection_rule", "best")),
                    collinearity_cfg=cfg["modeling"].get("multicollinearity", {}),
                )
                _, coef_df = fit_cox_model(
                    train_df=endpoint_df,
                    test_df=endpoint_df,
                    feature_cols=feature_cols,
                    clinical_cols=clinical_cols,
                    time_col=endpoint.time_col,
                    event_col=endpoint.event_col,
                    penalizer=penalizer,
                    l1_ratio=l1_ratio,
                    clinical_penalty_factor=float(cox_cfg.get("clinical_penalty_factor", 1.0)),
                    collinearity_cfg=cfg["modeling"].get("multicollinearity", {}),
                )
                stability_df = bootstrap_stability_cox(
                    df=endpoint_df,
                    feature_cols=feature_cols,
                    clinical_cols=clinical_cols,
                    time_col=endpoint.time_col,
                    event_col=endpoint.event_col,
                    penalizer=penalizer,
                    l1_ratio=l1_ratio,
                    n_bootstrap=int(cfg["modeling"]["stability"]["n_bootstrap"]),
                    selection_threshold=float(cfg["modeling"]["stability"]["selection_threshold"]),
                    random_state=int(cfg["modeling"].get("random_state", 42)),
                    clinical_penalty_factor=float(cox_cfg.get("clinical_penalty_factor", 1.0)),
                    selection_scope=str(cfg["modeling"]["stability"].get("selection_scope", "all")),
                    collinearity_cfg=cfg["modeling"].get("multicollinearity", {}),
                )
                coef_stab_df = _coefficient_stability_survival(endpoint_df, feature_cols, clinical_cols, endpoint, cfg)
                extra = {"best_penalizer": penalizer, "best_l1_ratio": l1_ratio, "full_inner_score": inner_cindex}

            best_rows.append({
                "endpoint_name": endpoint_name,
                "endpoint_type": endpoint_type,
                "feature_names": row.feature_names,
                "subset_name": row.subset_name,
                "subset_size": row.subset_size,
                "ranking_metric": row.ranking_metric,
                **extra,
            })
            coef_df = coef_df.copy()
            coef_df.insert(0, "subset_name", row.subset_name)
            coef_df.insert(0, "feature_names", row.feature_names)
            coef_df.insert(0, "endpoint_type", endpoint_type)
            coef_df.insert(0, "endpoint_name", endpoint_name)
            coef_frames.append(coef_df)

            stability_df = stability_df.copy()
            stability_df.insert(0, "subset_name", row.subset_name)
            stability_df.insert(0, "feature_names", row.feature_names)
            stability_df.insert(0, "endpoint_type", endpoint_type)
            stability_df.insert(0, "endpoint_name", endpoint_name)
            bootstrap_frames.append(stability_df)

            coef_stab_df = coef_stab_df.copy()
            coef_stab_df.insert(0, "subset_name", row.subset_name)
            coef_stab_df.insert(0, "feature_names", row.feature_names)
            coef_stab_df.insert(0, "endpoint_type", endpoint_type)
            coef_stab_df.insert(0, "endpoint_name", endpoint_name)
            coef_stability_frames.append(coef_stab_df)

    best_df = pd.DataFrame(best_rows)
    coef_all = pd.concat(coef_frames, ignore_index=True) if coef_frames else pd.DataFrame()
    coef_stability_all = pd.concat(coef_stability_frames, ignore_index=True) if coef_stability_frames else pd.DataFrame()
    bootstrap_all = pd.concat(bootstrap_frames, ignore_index=True) if bootstrap_frames else pd.DataFrame()

    best_df.to_csv(outdir / "stage5_best_subsets.csv", index=False)
    coef_all.to_csv(outdir / "stage5_best_subset_coefficients.csv", index=False)
    coef_stability_all.to_csv(outdir / "stage5_best_subset_coefficient_stability.csv", index=False)
    bootstrap_all.to_csv(outdir / "stage5_best_subset_bootstrap_stability.csv", index=False)
    return best_df, coef_all, coef_stability_all


def _outer_endpoint_dataset(base_dataset: SearchDataset, endpoint, endpoint_type: str, cfg: dict) -> SearchDataset:
    if endpoint_type == "binary":
        analysis_df = base_dataset.analysis_df.dropna(subset=[endpoint.outcome_col]).copy()
    else:
        analysis_df = base_dataset.analysis_df.dropna(subset=[endpoint.time_col, endpoint.event_col]).copy()
    return _build_dataset_view(base_dataset, analysis_df, cfg, subset_label=f"outer_endpoint__{endpoint.slug}")


def _evaluate_outer_models(
    train_df: pd.DataFrame,
    test_df: pd.DataFrame,
    clinical_cols: list[str],
    selected_features: tuple[str, ...],
    endpoint,
    endpoint_type: str,
    cfg: dict,
    fold_seed: int,
) -> tuple[dict[str, object], dict[str, object]]:
    baseline_feature_cols = list(clinical_cols)
    selected_feature_cols = list(clinical_cols) + list(selected_features)

    if endpoint_type == "binary":
        baseline_metrics = _score_binary_outer_test(
            train_df=train_df,
            test_df=test_df,
            feature_cols=baseline_feature_cols,
            clinical_cols=clinical_cols,
            endpoint=endpoint,
            cfg=cfg,
            fold_seed=fold_seed,
        )
        selected_metrics = _score_binary_outer_test(
            train_df=train_df,
            test_df=test_df,
            feature_cols=selected_feature_cols,
            clinical_cols=clinical_cols,
            endpoint=endpoint,
            cfg=cfg,
            fold_seed=fold_seed,
        )
    else:
        baseline_metrics = _score_survival_outer_test(
            train_df=train_df,
            test_df=test_df,
            feature_cols=baseline_feature_cols,
            clinical_cols=clinical_cols,
            endpoint=endpoint,
            cfg=cfg,
            fold_seed=fold_seed,
        )
        selected_metrics = _score_survival_outer_test(
            train_df=train_df,
            test_df=test_df,
            feature_cols=selected_feature_cols,
            clinical_cols=clinical_cols,
            endpoint=endpoint,
            cfg=cfg,
            fold_seed=fold_seed,
        )
    return baseline_metrics, selected_metrics


def _summarize_outer_metrics(fold_metrics_df: pd.DataFrame, endpoint_type: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    if fold_metrics_df.empty:
        return pd.DataFrame(), pd.DataFrame()

    metric_cols = ["roc_auc", "auprc", "brier"] if endpoint_type == "binary" else ["cindex"]
    available_metric_cols = [col for col in metric_cols if col in fold_metrics_df.columns]

    summary = fold_metrics_df.groupby(["endpoint_name", "endpoint_type", "model_role"], as_index=False).agg(
        n_folds=("fold", "nunique"),
        **{
            f"{metric}_mean": (metric, "mean")
            for metric in available_metric_cols
        },
        **{
            f"{metric}_sd": (metric, lambda s: float(s.std(ddof=1)) if len(s) > 1 else 0.0)
            for metric in available_metric_cols
        },
    )

    baseline = fold_metrics_df[fold_metrics_df["model_role"] == "baseline"].copy()
    selected = fold_metrics_df[fold_metrics_df["model_role"] == "selected"].copy()
    if baseline.empty or selected.empty:
        return summary, pd.DataFrame()

    merge_cols = ["endpoint_name", "endpoint_type", "fold"]
    baseline = baseline[merge_cols + available_metric_cols].rename(columns={metric: f"baseline_{metric}" for metric in available_metric_cols})
    selected = selected[merge_cols + available_metric_cols].rename(columns={metric: f"selected_{metric}" for metric in available_metric_cols})
    delta_df = baseline.merge(selected, on=merge_cols, how="inner")
    for metric in available_metric_cols:
        delta_df[f"delta_{metric}"] = delta_df[f"selected_{metric}"] - delta_df[f"baseline_{metric}"]

    delta_summary = delta_df.groupby(["endpoint_name", "endpoint_type"], as_index=False).agg(
        n_folds=("fold", "nunique"),
        **{
            f"delta_{metric}_mean": (f"delta_{metric}", "mean")
            for metric in available_metric_cols
        },
        **{
            f"delta_{metric}_sd": (f"delta_{metric}", lambda s: float(s.std(ddof=1)) if len(s) > 1 else 0.0)
            for metric in available_metric_cols
        },
    )
    return summary, delta_summary


def run_outer_search_validation(dataset: SearchDataset, cfg: dict, outdir: Path) -> Optional[Path]:
    outer_cfg = _outer_validation_cfg(cfg)
    if not outer_cfg["enabled"]:
        return None

    outer_outdir = ensure_outdir(outdir / "outer_search_validation")
    binary_endpoints, survival_endpoints = resolve_endpoint_specs(cfg, dataset.analysis_df.columns.tolist())

    fold_metric_rows: list[dict[str, object]] = []
    selection_rows: list[dict[str, object]] = []
    skipped_rows: list[dict[str, object]] = []

    all_endpoints = [(endpoint, "binary") for endpoint in binary_endpoints] + [(endpoint, "survival") for endpoint in survival_endpoints]
    for endpoint, endpoint_type in all_endpoints:
        endpoint_cfg = _single_endpoint_cfg(cfg, endpoint.name, endpoint_type)
        endpoint_dataset = _outer_endpoint_dataset(dataset, endpoint, endpoint_type, endpoint_cfg)
        fold_splits = _iter_outer_splits(endpoint_dataset.analysis_df, endpoint, endpoint_type, cfg)
        if not fold_splits:
            skipped_rows.append({
                "endpoint_name": endpoint.name,
                "endpoint_type": endpoint_type,
                "reason": "insufficient_class_or_event_variation_for_outer_cv",
                "n_patients": int(endpoint_dataset.analysis_df.shape[0]),
            })
            continue

        baseline_set = _resolve_baseline_clinical_feature_set(cfg)
        clinical_cols = list(endpoint_dataset.clinical_feature_sets[baseline_set])
        endpoint_outdir = ensure_outdir(outer_outdir / endpoint.slug)
        _progress_log(
            cfg,
            f"Outer validation [{endpoint.slug}] starting with {len(fold_splits)} folds "
            f"and n={endpoint_dataset.analysis_df.shape[0]} patients",
        )

        for fold_id, train_idx, test_idx in fold_splits:
            fold_dir = ensure_outdir(endpoint_outdir / f"fold_{fold_id:02d}")
            train_dataset = _build_dataset_view(
                endpoint_dataset,
                endpoint_dataset.analysis_df.iloc[train_idx].copy(),
                endpoint_cfg,
                subset_label=f"{endpoint.slug}__outer_train_fold_{fold_id:02d}",
            )
            test_dataset = _build_dataset_view(
                endpoint_dataset,
                endpoint_dataset.analysis_df.iloc[test_idx].copy(),
                endpoint_cfg,
                subset_label=f"{endpoint.slug}__outer_test_fold_{fold_id:02d}",
            )

            _progress_log(
                cfg,
                f"Outer validation {endpoint.slug} fold {fold_id}: "
                f"search on train n={train_dataset.analysis_df.shape[0]}, test n={test_dataset.analysis_df.shape[0]}",
            )
            try:
                with _tracked_stage(
                    f"outer_validation.{endpoint.slug}.fold_{fold_id:02d}.train_search",
                    {
                        "endpoint_name": endpoint.name,
                        "endpoint_type": endpoint_type,
                        "fold": int(fold_id),
                        "n_train": int(train_dataset.analysis_df.shape[0]),
                        "n_test": int(test_dataset.analysis_df.shape[0]),
                    },
                ):
                    inner_results = run_feature_search_inner(
                        train_dataset,
                        endpoint_cfg,
                        fold_dir,
                        run_stage5=False,
                        write_fold_detail_files=outer_cfg["write_detailed_cv_files"],
                        stage_prefix=f"outer_validation.{endpoint.slug}.fold_{fold_id:02d}.train_search",
                    )
                    candidate = _select_outer_candidate(
                        inner_results["baseline_df"],
                        inner_results["stage2_df"],
                        inner_results["subset_results_df"],
                        endpoint_name=endpoint.name,
                        endpoint_type=endpoint_type,
                        cfg=endpoint_cfg,
                    )

                fold_seed = int(cfg["modeling"].get("random_state", 42)) + int(fold_id)
                with _tracked_stage(
                    f"outer_validation.{endpoint.slug}.fold_{fold_id:02d}.outer_test_evaluation",
                    {
                        "endpoint_name": endpoint.name,
                        "endpoint_type": endpoint_type,
                        "fold": int(fold_id),
                        "candidate_source": candidate["candidate_source"],
                        "selected_subset_size": int(candidate["subset_size"]),
                    },
                ):
                    baseline_metrics, selected_metrics = _evaluate_outer_models(
                        train_df=train_dataset.analysis_df,
                        test_df=test_dataset.analysis_df,
                        clinical_cols=clinical_cols,
                        selected_features=tuple(candidate["feature_names"]),
                        endpoint=endpoint,
                        endpoint_type=endpoint_type,
                        cfg=endpoint_cfg,
                        fold_seed=fold_seed,
                    )
            except Exception as exc:
                skipped_rows.append({
                    "endpoint_name": endpoint.name,
                    "endpoint_type": endpoint_type,
                    "reason": "fold_failed",
                    "fold": fold_id,
                    "error": str(exc),
                    "n_patients": int(endpoint_dataset.analysis_df.shape[0]),
                })
                continue

            selection_rows.append({
                "endpoint_name": endpoint.name,
                "endpoint_type": endpoint_type,
                "endpoint_slug": endpoint.slug,
                "fold": fold_id,
                "candidate_source": candidate["candidate_source"],
                "feature_names": candidate["feature_names_joined"],
                "subset_name": candidate["subset_name"],
                "subset_size": int(candidate["subset_size"]),
                "train_ranking_metric": candidate["train_ranking_metric"],
                "n_train_candidates_after_qc": int(len(train_dataset.candidate_features)),
            })

            common_row = {
                "endpoint_name": endpoint.name,
                "endpoint_type": endpoint_type,
                "endpoint_slug": endpoint.slug,
                "fold": fold_id,
                "candidate_source": candidate["candidate_source"],
                "selected_feature_names": candidate["feature_names_joined"],
                "selected_subset_name": candidate["subset_name"],
                "selected_subset_size": int(candidate["subset_size"]),
                "train_ranking_metric": candidate["train_ranking_metric"],
            }
            fold_metric_rows.append(common_row | {"model_role": "baseline"} | baseline_metrics)
            fold_metric_rows.append(common_row | {"model_role": "selected"} | selected_metrics)

    fold_metrics_df = pd.DataFrame(fold_metric_rows)
    selection_df = pd.DataFrame(selection_rows)
    skipped_df = pd.DataFrame(skipped_rows)

    if fold_metrics_df.empty:
        summary_df = pd.DataFrame()
        delta_summary_df = pd.DataFrame()
    else:
        binary_summary_df, binary_delta_df = _summarize_outer_metrics(
            fold_metrics_df[fold_metrics_df["endpoint_type"] == "binary"].copy(),
            endpoint_type="binary",
        )
        survival_summary_df, survival_delta_df = _summarize_outer_metrics(
            fold_metrics_df[fold_metrics_df["endpoint_type"] == "survival"].copy(),
            endpoint_type="survival",
        )
        summary_df = pd.concat([binary_summary_df, survival_summary_df], ignore_index=True)
        delta_summary_df = pd.concat([binary_delta_df, survival_delta_df], ignore_index=True)

    fold_metrics_df.to_csv(outer_outdir / "outer_fold_metrics.csv", index=False)
    selection_df.to_csv(outer_outdir / "outer_selected_candidates.csv", index=False)
    summary_df.to_csv(outer_outdir / "outer_metrics_summary.csv", index=False)
    delta_summary_df.to_csv(outer_outdir / "outer_delta_summary.csv", index=False)
    skipped_df.to_csv(outer_outdir / "outer_skipped_endpoints.csv", index=False)
    write_json(
        {
            "enabled": True,
            "n_endpoints_attempted": int(len(all_endpoints)),
            "n_endpoints_skipped": int(skipped_df.shape[0]),
            "n_outer_folds_total": int(selection_df.shape[0]),
            "settings": outer_cfg,
        },
        outer_outdir / "outer_validation_manifest.json",
    )
    return outer_outdir


def run_feature_search(config_path: str) -> Path:
    global _ACTIVE_RUN_TRACKER
    cfg = read_yaml(config_path)
    outdir = ensure_outdir(cfg["paths"]["output_dir"])
    tracker = RunTracker(outdir, cfg)
    _ACTIVE_RUN_TRACKER = tracker
    tracker.log(
        f"Starting feature-library search with config: {config_path} | "
        f"pid={tracker.metadata['pid']} cpu_count={tracker.metadata['cpu_count']} "
        f"max_workers={tracker.parallel_cfg['max_workers']}"
    )

    try:
        with _tracked_stage("build_search_dataset", {"config_path": str(config_path)}):
            dataset = build_search_dataset(cfg)
        _progress_log(
            cfg,
            f"Dataset ready: n_patients={dataset.analysis_df.shape[0]}, "
            f"n_candidate_features={len(dataset.candidate_features)}",
        )

        tracker.log("Stage 2-5: exploratory full-cohort feature search")
        with _tracked_stage("exploratory_full_cohort", {"output_dir": str(outdir)}):
            inner_results = run_feature_search_inner(
                dataset,
                cfg,
                outdir,
                run_stage5=True,
                write_fold_detail_files=True,
                stage_prefix="exploratory_full_cohort",
            )

        outer_outdir = None
        if _outer_validation_cfg(cfg)["enabled"]:
            tracker.log("Outer validation: wrapping the full feature search inside outer CV")
            with _tracked_stage("outer_search_validation", {"output_dir": str(outdir / 'outer_search_validation')}):
                outer_outdir = run_outer_search_validation(dataset, cfg, outdir)

        manifest = {
            **dataset.manifest,
            "n_candidate_features": int(len(dataset.candidate_features)),
            "n_stage2_rows": int(inner_results["stage2_df"].shape[0]),
            "n_stage3_candidates": int(inner_results["pruned_df"].shape[0]),
            "n_stage4_subsets": int(inner_results["subset_results_df"].shape[0]),
            "n_stage5_best_rows": int(inner_results["best_df"].shape[0]),
            "outer_validation_enabled": bool(_outer_validation_cfg(cfg)["enabled"]),
            "outer_validation_dir": str(outer_outdir) if outer_outdir is not None else None,
            "runtime_dir": str(outdir / "runtime"),
        }
        write_json(manifest, outdir / "run_manifest.json")
        tracker.log(f"Feature-library search completed. Output written to: {outdir}")
        return outdir
    finally:
        tracker.finalize()
        _ACTIVE_RUN_TRACKER = None


def run_feature_search_from_existing_stage2(
    config_path: str,
    stage2_dir: str,
    *,
    run_outer_validation: bool = False,
) -> Path:
    global _ACTIVE_RUN_TRACKER
    cfg = read_yaml(config_path)
    outdir = ensure_outdir(cfg["paths"]["output_dir"])
    tracker = RunTracker(outdir, cfg)
    _ACTIVE_RUN_TRACKER = tracker
    tracker.log(
        f"Starting feature-library search from existing stage2 with config: {config_path} | "
        f"stage2_dir={stage2_dir} pid={tracker.metadata['pid']} cpu_count={tracker.metadata['cpu_count']} "
        f"max_workers={tracker.parallel_cfg['max_workers']}"
    )

    try:
        with _tracked_stage("build_search_dataset", {"config_path": str(config_path)}):
            dataset = build_search_dataset(cfg)
        _progress_log(
            cfg,
            f"Dataset ready: n_patients={dataset.analysis_df.shape[0]}, "
            f"n_candidate_features={len(dataset.candidate_features)}",
        )

        stage2_source_dir = Path(stage2_dir)
        with _tracked_stage("reuse_existing_stage2", {"stage2_dir": str(stage2_source_dir)}):
            baseline_df, stage2_df = _load_reused_stage2_results(dataset, cfg, stage2_source_dir)
        tracker.log(
            f"Reusing stage2 results from {stage2_source_dir}: "
            f"n_baseline_rows={baseline_df.shape[0]} n_stage2_rows={stage2_df.shape[0]}"
        )

        _write_search_inputs(dataset, cfg, outdir)
        baseline_df.to_csv(outdir / "stage2_baseline_metrics.csv", index=False)
        stage2_df.to_csv(outdir / "stage2_single_feature_scan.csv", index=False)

        tracker.log("Stage 3-5: exploratory full-cohort feature search using reused stage2 outputs")
        with _tracked_stage("exploratory_full_cohort.stage3_redundancy_pruning", {"output_dir": str(outdir), "stage2_source_dir": str(stage2_source_dir)}):
            pruned_df = run_stage3_redundancy_pruning(dataset, stage2_df, cfg, outdir)
        with _tracked_stage("exploratory_full_cohort.stage4_subset_search", {"output_dir": str(outdir), "n_stage3_candidates": int(pruned_df.shape[0])}):
            subset_results_df = run_stage4_subset_search(
                dataset,
                pruned_df,
                baseline_df,
                cfg,
                outdir,
                write_fold_detail_files=True,
            )
        with _tracked_stage("exploratory_full_cohort.stage5_best_subset_refits", {"output_dir": str(outdir), "n_stage4_subsets": int(subset_results_df.shape[0])}):
            best_df, coef_df, coef_stability_df = run_stage5_best_subset_refits(dataset, baseline_df, subset_results_df, cfg, outdir)

        outer_outdir = None
        outer_cfg = _outer_validation_cfg(cfg)
        if run_outer_validation and outer_cfg["enabled"]:
            tracker.log("Outer validation: running from scratch with the new v2 family-cap logic")
            with _tracked_stage("outer_search_validation", {"output_dir": str(outdir / 'outer_search_validation')}):
                outer_outdir = run_outer_search_validation(dataset, cfg, outdir)
        elif outer_cfg["enabled"]:
            tracker.log(
                "Outer validation is enabled in config but was skipped in stage3-resume mode. "
                "Pass --run-outer-validation to recompute it from scratch.",
                verbose_only=True,
            )

        manifest = {
            **dataset.manifest,
            "stage2_reused": True,
            "stage2_reused_from": str(stage2_source_dir),
            "n_candidate_features": int(len(dataset.candidate_features)),
            "n_stage2_rows": int(stage2_df.shape[0]),
            "n_stage3_candidates": int(pruned_df.shape[0]),
            "n_stage4_subsets": int(subset_results_df.shape[0]),
            "n_stage5_best_rows": int(best_df.shape[0]),
            "outer_validation_enabled": bool(run_outer_validation and outer_cfg["enabled"]),
            "outer_validation_dir": str(outer_outdir) if outer_outdir is not None else None,
            "runtime_dir": str(outdir / "runtime"),
        }
        write_json(manifest, outdir / "run_manifest.json")
        tracker.log(f"Stage3-resume feature-library search completed. Output written to: {outdir}")
        return outdir
    finally:
        tracker.finalize()
        _ACTIVE_RUN_TRACKER = None


def main(config_path: str) -> None:
    run_feature_search(config_path)


def cli_main() -> None:
    parser = argparse.ArgumentParser(description="Feature-library search pipeline for single-feature scans and subset search.")
    parser.add_argument("--config", required=True, help="Path to YAML config file")
    parser.add_argument(
        "--reuse-stage2-dir",
        help="Existing full-cohort output directory containing stage2_baseline_metrics.csv and stage2_single_feature_scan.csv. "
             "When set, the pipeline skips stage2 and reruns stage3-5 with the current code.",
    )
    parser.add_argument(
        "--run-outer-validation",
        action="store_true",
        help="When reusing stage2, also recompute outer validation from scratch with the current code.",
    )
    args = parser.parse_args()
    if args.reuse_stage2_dir:
        run_feature_search_from_existing_stage2(
            args.config,
            args.reuse_stage2_dir,
            run_outer_validation=bool(args.run_outer_validation),
        )
    else:
        main(args.config)


if __name__ == "__main__":
    cli_main()
