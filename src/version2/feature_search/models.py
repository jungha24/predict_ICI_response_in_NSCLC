"""Shared preprocessing plus elastic-net logistic/Cox model helpers."""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


CLINICAL_CATEGORICAL = {
    "sex",
    "histology",
    "smoking",
    "ecog",
    "egfr_status",
    "io_line",
    "previous_chemo",
    "previous_target",
    "drug_class"
}
PREFERRED_NUMERIC_CLINICAL = {"age", "pd_l1_tps"}


def resolve_stratified_splits(y: np.ndarray, requested_splits: int) -> int:
    if y.size == 0:
        return 0
    counts = np.bincount(y.astype(int), minlength=2)
    if np.count_nonzero(counts) < 2:
        return 0
    return int(min(requested_splits, counts.min()))


def infer_clinical_column_types(train_df: pd.DataFrame, clinical_cols: list[str]) -> tuple[list[str], list[str]]:
    categorical_cols = [column for column in clinical_cols if column in CLINICAL_CATEGORICAL]
    numeric_cols: list[str] = []
    for column in clinical_cols:
        if column in categorical_cols:
            continue
        series = train_df[column]
        non_missing = int(series.notna().sum())
        converted = pd.to_numeric(series, errors="coerce")
        numeric_non_missing = int(converted.notna().sum())
        if non_missing == 0:
            if column in PREFERRED_NUMERIC_CLINICAL:
                numeric_cols.append(column)
            else:
                categorical_cols.append(column)
            continue
        if numeric_non_missing == non_missing:
            numeric_cols.append(column)
            continue
        categorical_cols.append(column)
    return numeric_cols, categorical_cols


def build_transformer(numeric_cols: list[str], categorical_cols: list[str]) -> ColumnTransformer:
    # Keep clinical categorical handling in one place so every resampling split uses the same encoding rules.
    num_pipe = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
        ]
    )
    cat_pipe = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="most_frequent")),
            ("onehot", OneHotEncoder(drop="first", handle_unknown="ignore", sparse_output=False)),
        ]
    )

    transformer = ColumnTransformer(
        transformers=[
            ("num", num_pipe, numeric_cols),
            ("cat", cat_pipe, categorical_cols),
        ],
        remainder="drop",
        verbose_feature_names_out=True,
    )
    return transformer


def fit_transform_fold(
    train_df: pd.DataFrame,
    test_df: pd.DataFrame,
    feature_cols: list[str],
    clinical_cols: list[str],
    collinearity_cfg: Optional[dict] = None,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], list[str]]:
    # Fit preprocessing on the training fold only, then apply the identical transform to the held-out fold.
    immune_cols = [column for column in feature_cols if column not in clinical_cols]
    numeric_clinical, categorical_cols = infer_clinical_column_types(train_df, clinical_cols)
    numeric_cols = immune_cols + numeric_clinical
    transformer = build_transformer(numeric_cols, categorical_cols)

    train_features = train_df[feature_cols].copy()
    test_features = test_df[feature_cols].copy()
    for column in numeric_cols:
        train_features[column] = pd.to_numeric(train_features[column], errors="coerce")
        test_features[column] = pd.to_numeric(test_features[column], errors="coerce")

    x_train = transformer.fit_transform(train_features)
    x_test = transformer.transform(test_features)

    feature_names = list(transformer.get_feature_names_out())
    raw_names = list(numeric_cols)
    if categorical_cols:
        encoder = transformer.named_transformers_["cat"].named_steps["onehot"]
        for raw_col, categories in zip(categorical_cols, encoder.categories_):
            retained_categories = list(categories[1:]) if len(categories) > 1 else []
            raw_names.extend([raw_col] * len(retained_categories))

    x_train_df = pd.DataFrame(x_train, index=train_df.index, columns=feature_names).astype(np.float32)
    x_test_df = pd.DataFrame(x_test, index=test_df.index, columns=feature_names).astype(np.float32)

    keep_mask = x_train_df.std(axis=0, ddof=0) > 1e-12
    x_train_df = x_train_df.loc[:, keep_mask]
    x_test_df = x_test_df.loc[:, keep_mask]
    feature_names = list(x_train_df.columns)
    raw_names = [raw_name for raw_name, keep in zip(raw_names, keep_mask.tolist()) if keep]

    x_train_df, x_test_df, feature_names, raw_names = apply_multicollinearity_pruning(
        x_train_df=x_train_df,
        x_test_df=x_test_df,
        raw_names=raw_names,
        clinical_cols=clinical_cols,
        collinearity_cfg=collinearity_cfg,
    )
    return x_train_df, x_test_df, feature_names, raw_names


def build_coef_df(
    feature_names: list[str],
    raw_names: list[str],
    coef_values: np.ndarray,
    clinical_cols: list[str],
) -> pd.DataFrame:
    coef_df = pd.DataFrame(
        {
            "feature": feature_names,
            "raw_feature": raw_names,
            "coef": coef_values.astype(float),
        }
    )
    coef_df["abs_coef"] = coef_df["coef"].abs()
    coef_df["feature_scope"] = np.where(coef_df["raw_feature"].isin(clinical_cols), "clinical", "immune")
    return coef_df


def feature_scope(raw_name: str, clinical_cols: list[str]) -> str:
    return "clinical" if raw_name in clinical_cols else "immune"


def choose_feature_to_drop(
    feature_a: str,
    feature_b: str,
    raw_name_a: str,
    raw_name_b: str,
    clinical_cols: list[str],
) -> str:
    is_clinical_a = raw_name_a in clinical_cols
    is_clinical_b = raw_name_b in clinical_cols
    if is_clinical_a and not is_clinical_b:
        return feature_b
    if is_clinical_b and not is_clinical_a:
        return feature_a
    return feature_b


def prune_correlated_features(
    x_train_df: pd.DataFrame,
    x_test_df: pd.DataFrame,
    raw_names: list[str],
    clinical_cols: list[str],
    correlation_threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], list[str]]:
    if x_train_df.shape[1] < 2 or correlation_threshold <= 0 or correlation_threshold >= 1:
        return x_train_df, x_test_df, list(x_train_df.columns), raw_names

    corr = x_train_df.corr().abs().fillna(0.0)
    columns = list(x_train_df.columns)
    raw_name_map = dict(zip(columns, raw_names))
    drop_columns: set[str] = set()

    for left_idx, left_column in enumerate(columns):
        if left_column in drop_columns:
            continue
        for right_idx in range(left_idx + 1, len(columns)):
            right_column = columns[right_idx]
            if right_column in drop_columns:
                continue
            if float(corr.iloc[left_idx, right_idx]) < correlation_threshold:
                continue
            drop_column = choose_feature_to_drop(
                left_column,
                right_column,
                raw_name_map[left_column],
                raw_name_map[right_column],
                clinical_cols,
            )
            drop_columns.add(drop_column)

    if not drop_columns:
        return x_train_df, x_test_df, list(x_train_df.columns), raw_names

    keep_columns = [column for column in columns if column not in drop_columns]
    kept_raw_names = [raw_name_map[column] for column in keep_columns]
    return x_train_df[keep_columns], x_test_df[keep_columns], keep_columns, kept_raw_names


def compute_vif_scores(x_train_df: pd.DataFrame) -> pd.Series:
    if x_train_df.shape[1] < 2:
        return pd.Series(1.0, index=x_train_df.columns, dtype=float)

    x = x_train_df.to_numpy(dtype=float)
    scores: dict[str, float] = {}
    for idx, column in enumerate(x_train_df.columns):
        y = x[:, idx]
        x_other = np.delete(x, idx, axis=1)
        if x_other.shape[1] == 0 or np.nanstd(y) <= 1e-12:
            scores[column] = 1.0
            continue
        try:
            coef, _, _, _ = np.linalg.lstsq(x_other, y, rcond=None)
            pred = x_other @ coef
            ss_tot = float(np.sum((y - np.mean(y)) ** 2))
            ss_res = float(np.sum((y - pred) ** 2))
            if ss_tot <= 1e-12:
                scores[column] = 1.0
                continue
            r_squared = 1.0 - (ss_res / ss_tot)
            r_squared = min(max(r_squared, 0.0), 0.999999)
            scores[column] = 1.0 / max(1.0 - r_squared, 1e-6)
        except Exception:
            scores[column] = np.inf
    return pd.Series(scores, dtype=float)


def prune_high_vif_features(
    x_train_df: pd.DataFrame,
    x_test_df: pd.DataFrame,
    raw_names: list[str],
    clinical_cols: list[str],
    max_vif: float,
    max_iter: int,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], list[str]]:
    if x_train_df.shape[1] < 2 or max_vif <= 1 or max_iter <= 0:
        return x_train_df, x_test_df, list(x_train_df.columns), raw_names

    current_train = x_train_df.copy()
    current_test = x_test_df.copy()
    current_raw_names = list(raw_names)

    for _ in range(max_iter):
        vif_scores = compute_vif_scores(current_train)
        if vif_scores.empty or float(vif_scores.max()) <= max_vif:
            break

        above_threshold = vif_scores[vif_scores > max_vif].sort_values(ascending=False)
        raw_name_map = dict(zip(current_train.columns, current_raw_names))
        drop_column = str(above_threshold.index[0])
        current_train = current_train.drop(columns=[drop_column])
        current_test = current_test.drop(columns=[drop_column])
        current_raw_names = [raw_name_map[column] for column in current_train.columns]
        if current_train.shape[1] < 2:
            break

    return current_train, current_test, list(current_train.columns), current_raw_names


def apply_multicollinearity_pruning(
    x_train_df: pd.DataFrame,
    x_test_df: pd.DataFrame,
    raw_names: list[str],
    clinical_cols: list[str],
    collinearity_cfg: Optional[dict],
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], list[str]]:
    cfg = collinearity_cfg or {}
    if not bool(cfg.get("enabled", False)):
        return x_train_df, x_test_df, list(x_train_df.columns), raw_names

    train_df, test_df, feature_names, feature_raw_names = prune_correlated_features(
        x_train_df=x_train_df,
        x_test_df=x_test_df,
        raw_names=raw_names,
        clinical_cols=clinical_cols,
        correlation_threshold=float(cfg.get("correlation_threshold", 0.95)),
    )

    vif_cfg = cfg.get("vif", {})
    if bool(vif_cfg.get("enabled", False)):
        train_df, test_df, feature_names, feature_raw_names = prune_high_vif_features(
            x_train_df=train_df,
            x_test_df=test_df,
            raw_names=feature_raw_names,
            clinical_cols=clinical_cols,
            max_vif=float(vif_cfg.get("max_vif", 10.0)),
            max_iter=int(vif_cfg.get("max_iter", 50)),
        )

    return train_df, test_df, feature_names, feature_raw_names


def penalty_weights(raw_names: list[str], clinical_cols: list[str], base_penalty: float, clinical_penalty_factor: float) -> np.ndarray:
    weights = np.full(len(raw_names), float(base_penalty), dtype=float)
    factor = max(float(clinical_penalty_factor), 0.0)
    for idx, raw_name in enumerate(raw_names):
        if raw_name in clinical_cols:
            weights[idx] = float(base_penalty) * factor
    return weights


def safe_auc(y_true: np.ndarray, y_prob: np.ndarray) -> float:
    if len(np.unique(y_true)) < 2:
        return np.nan
    return float(roc_auc_score(y_true, y_prob))


def summarize_resampling_scores(scores: Sequence[float]) -> dict[str, float]:
    valid = np.asarray([score for score in scores if np.isfinite(score)], dtype=float)
    if valid.size == 0:
        return {
            "mean_score": np.nan,
            "score_sd": np.nan,
            "score_se": np.nan,
            "n_valid_scores": 0,
        }

    mean_score = float(valid.mean())
    score_sd = float(valid.std(ddof=1)) if valid.size > 1 else 0.0
    score_se = float(score_sd / np.sqrt(valid.size))
    return {
        "mean_score": mean_score,
        "score_sd": score_sd,
        "score_se": score_se,
        "n_valid_scores": int(valid.size),
    }


def select_candidate_by_rule(
    candidate_rows: list[dict[str, float]],
    strength_key: str,
    selection_rule: str,
) -> Optional[dict[str, float]]:
    valid_rows = [row for row in candidate_rows if np.isfinite(row.get("mean_score", np.nan))]
    if not valid_rows:
        return None

    best_row = max(valid_rows, key=lambda row: (row["mean_score"], row[strength_key], row["l1_ratio"]))
    if str(selection_rule or "best").lower() != "one_se":
        return best_row

    best_score_se = float(best_row.get("score_se", 0.0)) if np.isfinite(best_row.get("score_se", np.nan)) else 0.0
    threshold = float(best_row["mean_score"]) - best_score_se
    eligible = [row for row in valid_rows if row["mean_score"] >= threshold]
    if not eligible:
        return best_row
    return max(eligible, key=lambda row: (row[strength_key], row["l1_ratio"], row["mean_score"]))


def tune_logistic_elastic_net(
    df_train: pd.DataFrame,
    feature_cols: list[str],
    clinical_cols: list[str],
    outcome_col: str,
    alpha_grid: Sequence[float],
    l1_ratio_grid: Sequence[float],
    n_inner_splits: int,
    random_state: int,
    max_iter: int,
    selection_rule: str = "best",
    collinearity_cfg: Optional[dict] = None,
) -> tuple[float, float, float]:
    # Inner CV chooses the elastic-net strength used for the current binary endpoint.
    y = df_train[outcome_col].astype(int).to_numpy()
    resolved_splits = resolve_stratified_splits(y, n_inner_splits)
    if resolved_splits < 2:
        raise ValueError(f"Not enough class variation to tune logistic model for {outcome_col}.")

    splitter = StratifiedKFold(n_splits=resolved_splits, shuffle=True, random_state=random_state)
    candidate_rows: list[dict[str, float]] = []

    for alpha in alpha_grid:
        c_value = 1.0 / max(float(alpha), 1e-8)
        for l1_ratio in l1_ratio_grid:
            fold_scores: list[float] = []
            for inner_train_idx, inner_val_idx in splitter.split(df_train, y):
                inner_train = df_train.iloc[inner_train_idx].copy()
                inner_val = df_train.iloc[inner_val_idx].copy()
                x_train, x_val, _, _ = fit_transform_fold(
                    inner_train,
                    inner_val,
                    feature_cols,
                    clinical_cols,
                    collinearity_cfg=collinearity_cfg,
                )
                model = LogisticRegression(
                    penalty="elasticnet",
                    solver="saga",
                    C=c_value,
                    l1_ratio=float(l1_ratio),
                    max_iter=max_iter,
                    random_state=random_state,
                )
                try:
                    model.fit(x_train, inner_train[outcome_col].astype(int).to_numpy())
                    pred = model.predict_proba(x_val)[:, 1]
                    fold_scores.append(safe_auc(inner_val[outcome_col].astype(int).to_numpy(), pred))
                except Exception:
                    fold_scores.append(np.nan)

            candidate_rows.append(
                {
                    "alpha": float(alpha),
                    "l1_ratio": float(l1_ratio),
                    **summarize_resampling_scores(fold_scores),
                }
            )

    selected = select_candidate_by_rule(candidate_rows, strength_key="alpha", selection_rule=selection_rule)
    if selected is None:
        return float(alpha_grid[0]), float(l1_ratio_grid[0]), np.nan
    return float(selected["alpha"]), float(selected["l1_ratio"]), float(selected["mean_score"])


def fit_logistic_model(
    train_df: pd.DataFrame,
    test_df: pd.DataFrame,
    feature_cols: list[str],
    clinical_cols: list[str],
    outcome_col: str,
    alpha: float,
    l1_ratio: float,
    max_iter: int,
    random_state: int,
    collinearity_cfg: Optional[dict] = None,
) -> tuple[np.ndarray, pd.DataFrame]:
    x_train, x_test, feature_names, raw_names = fit_transform_fold(
        train_df,
        test_df,
        feature_cols,
        clinical_cols,
        collinearity_cfg=collinearity_cfg,
    )
    model = LogisticRegression(
        penalty="elasticnet",
        solver="saga",
        C=1.0 / max(float(alpha), 1e-8),
        l1_ratio=float(l1_ratio),
        max_iter=max_iter,
        random_state=random_state,
    )
    model.fit(x_train, train_df[outcome_col].astype(int).to_numpy())
    pred = model.predict_proba(x_test)[:, 1]
    coef_df = build_coef_df(feature_names, raw_names, model.coef_.ravel(), clinical_cols)
    return np.asarray(pred, dtype=float), coef_df


def evaluate_logistic_nested_cv(
    df: pd.DataFrame,
    feature_cols: list[str],
    clinical_cols: list[str],
    outcome_col: str,
    alpha_grid: Sequence[float],
    l1_ratio_grid: Sequence[float],
    n_outer_splits: int,
    n_outer_repeats: int,
    n_inner_splits: int,
    random_state: int,
    max_iter: int,
    selection_rule: str = "best",
    collinearity_cfg: Optional[dict] = None,
) -> tuple[pd.DataFrame, dict[str, float]]:
    # Repeated outer CV estimates binary-endpoint performance after inner tuning on each split.
    y = df[outcome_col].astype(int).to_numpy()
    resolved_splits = resolve_stratified_splits(y, n_outer_splits)
    if resolved_splits < 2:
        raise ValueError(f"Not enough class variation to run repeated CV for {outcome_col}.")

    splitter = RepeatedStratifiedKFold(
        n_splits=resolved_splits,
        n_repeats=n_outer_repeats,
        random_state=random_state,
    )
    rows: list[dict[str, float]] = []
    for fold_idx, (train_idx, test_idx) in enumerate(splitter.split(df, y), start=1):
        train_df = df.iloc[train_idx].copy()
        test_df = df.iloc[test_idx].copy()
        alpha, l1_ratio, inner_auc = tune_logistic_elastic_net(
            df_train=train_df,
            feature_cols=feature_cols,
            clinical_cols=clinical_cols,
            outcome_col=outcome_col,
            alpha_grid=alpha_grid,
            l1_ratio_grid=l1_ratio_grid,
            n_inner_splits=n_inner_splits,
            random_state=random_state + fold_idx,
            max_iter=max_iter,
            selection_rule=selection_rule,
            collinearity_cfg=collinearity_cfg,
        )
        pred, _ = fit_logistic_model(
            train_df=train_df,
            test_df=test_df,
            feature_cols=feature_cols,
            clinical_cols=clinical_cols,
            outcome_col=outcome_col,
            alpha=alpha,
            l1_ratio=l1_ratio,
            max_iter=max_iter,
            random_state=random_state + fold_idx,
            collinearity_cfg=collinearity_cfg,
        )
        y_test = test_df[outcome_col].astype(int).to_numpy()
        rows.append(
            {
                "fold": fold_idx,
                "inner_auc": inner_auc,
                "alpha": alpha,
                "l1_ratio": l1_ratio,
                "roc_auc": safe_auc(y_test, pred),
                "auprc": average_precision_score(y_test, pred) if len(np.unique(y_test)) > 1 else np.nan,
                "brier": brier_score_loss(y_test, pred),
                "n_train": len(train_df),
                "n_test": len(test_df),
                "events_test": int(y_test.sum()),
            }
        )

    fold_df = pd.DataFrame(rows)
    summary = {
        "roc_auc_mean": float(fold_df["roc_auc"].mean()),
        "roc_auc_sd": float(fold_df["roc_auc"].std(ddof=1)),
        "auprc_mean": float(fold_df["auprc"].mean()),
        "auprc_sd": float(fold_df["auprc"].std(ddof=1)),
        "brier_mean": float(fold_df["brier"].mean()),
        "brier_sd": float(fold_df["brier"].std(ddof=1)),
    }
    return fold_df, summary


def tune_on_full_data_logistic(
    df: pd.DataFrame,
    feature_cols: list[str],
    clinical_cols: list[str],
    outcome_col: str,
    alpha_grid: Sequence[float],
    l1_ratio_grid: Sequence[float],
    n_inner_splits: int,
    random_state: int,
    max_iter: int,
    selection_rule: str = "best",
    collinearity_cfg: Optional[dict] = None,
) -> tuple[float, float, float, pd.DataFrame]:
    alpha, l1_ratio, score = tune_logistic_elastic_net(
        df_train=df,
        feature_cols=feature_cols,
        clinical_cols=clinical_cols,
        outcome_col=outcome_col,
        alpha_grid=alpha_grid,
        l1_ratio_grid=l1_ratio_grid,
        n_inner_splits=n_inner_splits,
        random_state=random_state,
        max_iter=max_iter,
        selection_rule=selection_rule,
        collinearity_cfg=collinearity_cfg,
    )
    _, coef_df = fit_logistic_model(
        train_df=df,
        test_df=df,
        feature_cols=feature_cols,
        clinical_cols=clinical_cols,
        outcome_col=outcome_col,
        alpha=alpha,
        l1_ratio=l1_ratio,
        max_iter=max_iter,
        random_state=random_state,
        collinearity_cfg=collinearity_cfg,
    )
    return alpha, l1_ratio, score, coef_df


def bootstrap_stability_logistic(
    df: pd.DataFrame,
    feature_cols: list[str],
    clinical_cols: list[str],
    outcome_col: str,
    alpha: float,
    l1_ratio: float,
    n_bootstrap: int,
    selection_threshold: float,
    random_state: int,
    max_iter: int,
    selection_scope: str = "immune_only",
    collinearity_cfg: Optional[dict] = None,
) -> pd.DataFrame:
    # Stability selection counts how often immune features survive refitting across bootstrap samples.
    rng = np.random.default_rng(random_state)
    selection_counts: dict[tuple[str, str, str], int] = {}
    for iteration in range(n_bootstrap):
        sample_index = rng.choice(df.index.to_numpy(), size=len(df), replace=True)
        bootstrap_df = df.loc[sample_index].reset_index(drop=True)
        try:
            _, coef_df = fit_logistic_model(
                train_df=bootstrap_df,
                test_df=bootstrap_df,
                feature_cols=feature_cols,
                clinical_cols=clinical_cols,
                outcome_col=outcome_col,
                alpha=alpha,
                l1_ratio=l1_ratio,
                max_iter=max_iter,
                random_state=random_state + iteration,
                collinearity_cfg=collinearity_cfg,
            )
        except Exception:
            continue

        selected = coef_df[coef_df["abs_coef"] > selection_threshold].copy()
        if selection_scope == "immune_only":
            selected = selected[selected["feature_scope"] == "immune"]

        for row in selected.itertuples(index=False):
            key = (row.feature, row.raw_feature, row.feature_scope)
            selection_counts[key] = selection_counts.get(key, 0) + 1

    if not selection_counts:
        return pd.DataFrame(columns=["feature", "raw_feature", "feature_scope", "selection_count", "selection_frequency"])

    rows = []
    for (feature, raw_feature, feature_scope), count in selection_counts.items():
        rows.append(
            {
                "feature": feature,
                "raw_feature": raw_feature,
                "feature_scope": feature_scope,
                "selection_count": count,
                "selection_frequency": count / max(n_bootstrap, 1),
            }
        )
    return pd.DataFrame(rows).sort_values(["selection_frequency", "selection_count"], ascending=[False, False]).reset_index(drop=True)


def tune_cox_elastic_net(
    df_train: pd.DataFrame,
    feature_cols: list[str],
    clinical_cols: list[str],
    time_col: str,
    event_col: str,
    penalizer_grid: Sequence[float],
    l1_ratio_grid: Sequence[float],
    n_inner_splits: int,
    random_state: int,
    clinical_penalty_factor: float,
    selection_rule: str = "best",
    collinearity_cfg: Optional[dict] = None,
) -> tuple[float, float, float]:
    # Inner CV chooses the elastic-net setting for the PFS Cox model.
    y_event = df_train[event_col].astype(int).to_numpy()
    resolved_splits = resolve_stratified_splits(y_event, n_inner_splits)
    if resolved_splits < 2:
        raise ValueError(f"Not enough event variation to tune Cox model for {time_col}/{event_col}.")

    splitter = StratifiedKFold(n_splits=resolved_splits, shuffle=True, random_state=random_state)
    candidate_rows: list[dict[str, float]] = []

    for penalizer in penalizer_grid:
        for l1_ratio in l1_ratio_grid:
            fold_scores: list[float] = []
            for inner_train_idx, inner_val_idx in splitter.split(df_train, y_event):
                inner_train = df_train.iloc[inner_train_idx].copy()
                inner_val = df_train.iloc[inner_val_idx].copy()
                x_train, x_val, _, raw_names = fit_transform_fold(
                    inner_train,
                    inner_val,
                    feature_cols,
                    clinical_cols,
                    collinearity_cfg=collinearity_cfg,
                )
                weights = penalty_weights(raw_names, clinical_cols, penalizer, clinical_penalty_factor)
                fit_df = pd.concat(
                    [inner_train[[time_col, event_col]].reset_index(drop=True), x_train.reset_index(drop=True)],
                    axis=1,
                )
                try:
                    model = CoxPHFitter(penalizer=weights, l1_ratio=l1_ratio)
                    model.fit(fit_df, duration_col=time_col, event_col=event_col, show_progress=False)
                    risk = model.predict_partial_hazard(x_val).to_numpy().ravel()
                    fold_scores.append(
                        concordance_index(inner_val[time_col].to_numpy(), -risk, inner_val[event_col].to_numpy())
                    )
                except Exception:
                    fold_scores.append(np.nan)

            candidate_rows.append(
                {
                    "penalizer": float(penalizer),
                    "l1_ratio": float(l1_ratio),
                    **summarize_resampling_scores(fold_scores),
                }
            )

    selected = select_candidate_by_rule(candidate_rows, strength_key="penalizer", selection_rule=selection_rule)
    if selected is None:
        return float(penalizer_grid[0]), float(l1_ratio_grid[0]), np.nan
    return float(selected["penalizer"]), float(selected["l1_ratio"]), float(selected["mean_score"])


def fit_cox_model(
    train_df: pd.DataFrame,
    test_df: pd.DataFrame,
    feature_cols: list[str],
    clinical_cols: list[str],
    time_col: str,
    event_col: str,
    penalizer: float,
    l1_ratio: float,
    clinical_penalty_factor: float,
    collinearity_cfg: Optional[dict] = None,
) -> tuple[np.ndarray, pd.DataFrame]:
    # Clinical covariates stay in the combined Cox design but receive a weaker penalty than immune features.
    x_train, x_test, feature_names, raw_names = fit_transform_fold(
        train_df,
        test_df,
        feature_cols,
        clinical_cols,
        collinearity_cfg=collinearity_cfg,
    )
    weights = penalty_weights(raw_names, clinical_cols, penalizer, clinical_penalty_factor)
    fit_df = pd.concat(
        [train_df[[time_col, event_col]].reset_index(drop=True), x_train.reset_index(drop=True)],
        axis=1,
    )
    try:
        model = CoxPHFitter(penalizer=weights, l1_ratio=l1_ratio)
        model.fit(fit_df, duration_col=time_col, event_col=event_col, show_progress=False)
    except Exception:
        fallback_weights = np.maximum(weights, max(float(penalizer) * 0.1, 1e-4))
        model = CoxPHFitter(penalizer=fallback_weights, l1_ratio=0.0)
        model.fit(fit_df, duration_col=time_col, event_col=event_col, show_progress=False)

    risk = model.predict_partial_hazard(x_test).to_numpy().ravel()
    coef_df = build_coef_df(feature_names, raw_names, model.params_.to_numpy(), clinical_cols)
    return risk, coef_df


def evaluate_cox_nested_cv(
    df: pd.DataFrame,
    feature_cols: list[str],
    clinical_cols: list[str],
    time_col: str,
    event_col: str,
    penalizer_grid: Sequence[float],
    l1_ratio_grid: Sequence[float],
    n_outer_splits: int,
    n_outer_repeats: int,
    n_inner_splits: int,
    random_state: int,
    clinical_penalty_factor: float,
    selection_rule: str = "best",
    collinearity_cfg: Optional[dict] = None,
) -> tuple[pd.DataFrame, dict[str, float]]:
    # Repeated outer CV estimates PFS concordance after inner tuning on each split.
    y_event = df[event_col].astype(int).to_numpy()
    resolved_splits = resolve_stratified_splits(y_event, n_outer_splits)
    if resolved_splits < 2:
        raise ValueError(f"Not enough event variation to run repeated CV for {time_col}/{event_col}.")

    rows: list[dict[str, float]] = []
    for repeat_idx in range(n_outer_repeats):
        splitter = StratifiedKFold(
            n_splits=resolved_splits,
            shuffle=True,
            random_state=random_state + repeat_idx,
        )
        for fold_idx, (train_idx, test_idx) in enumerate(splitter.split(df, y_event), start=1):
            fold_id = repeat_idx * resolved_splits + fold_idx
            train_df = df.iloc[train_idx].copy()
            test_df = df.iloc[test_idx].copy()
            penalizer, l1_ratio, inner_cindex = tune_cox_elastic_net(
                df_train=train_df,
                feature_cols=feature_cols,
                clinical_cols=clinical_cols,
                time_col=time_col,
                event_col=event_col,
                penalizer_grid=penalizer_grid,
                l1_ratio_grid=l1_ratio_grid,
                n_inner_splits=n_inner_splits,
                random_state=random_state + fold_id,
                clinical_penalty_factor=clinical_penalty_factor,
                selection_rule=selection_rule,
                collinearity_cfg=collinearity_cfg,
            )
            risk, _ = fit_cox_model(
                train_df=train_df,
                test_df=test_df,
                feature_cols=feature_cols,
                clinical_cols=clinical_cols,
                time_col=time_col,
                event_col=event_col,
                penalizer=penalizer,
                l1_ratio=l1_ratio,
                clinical_penalty_factor=clinical_penalty_factor,
                collinearity_cfg=collinearity_cfg,
            )
            rows.append(
                {
                    "fold": fold_id,
                    "inner_cindex": inner_cindex,
                    "penalizer": penalizer,
                    "l1_ratio": l1_ratio,
                    "cindex": concordance_index(
                        test_df[time_col].to_numpy(),
                        -risk,
                        test_df[event_col].to_numpy(),
                    ),
                    "n_train": len(train_df),
                    "n_test": len(test_df),
                    "events_test": int(test_df[event_col].sum()),
                }
            )

    fold_df = pd.DataFrame(rows)
    summary = {
        "cindex_mean": float(fold_df["cindex"].mean()),
        "cindex_sd": float(fold_df["cindex"].std(ddof=1)),
    }
    return fold_df, summary


def tune_on_full_data_cox(
    df: pd.DataFrame,
    feature_cols: list[str],
    clinical_cols: list[str],
    time_col: str,
    event_col: str,
    penalizer_grid: Sequence[float],
    l1_ratio_grid: Sequence[float],
    n_inner_splits: int,
    random_state: int,
    clinical_penalty_factor: float,
    selection_rule: str = "best",
    collinearity_cfg: Optional[dict] = None,
) -> tuple[float, float, float, pd.DataFrame]:
    penalizer, l1_ratio, score = tune_cox_elastic_net(
        df_train=df,
        feature_cols=feature_cols,
        clinical_cols=clinical_cols,
        time_col=time_col,
        event_col=event_col,
        penalizer_grid=penalizer_grid,
        l1_ratio_grid=l1_ratio_grid,
        n_inner_splits=n_inner_splits,
        random_state=random_state,
        clinical_penalty_factor=clinical_penalty_factor,
        selection_rule=selection_rule,
        collinearity_cfg=collinearity_cfg,
    )
    _, coef_df = fit_cox_model(
        train_df=df,
        test_df=df,
        feature_cols=feature_cols,
        clinical_cols=clinical_cols,
        time_col=time_col,
        event_col=event_col,
        penalizer=penalizer,
        l1_ratio=l1_ratio,
        clinical_penalty_factor=clinical_penalty_factor,
        collinearity_cfg=collinearity_cfg,
    )
    return penalizer, l1_ratio, score, coef_df


def bootstrap_stability_cox(
    df: pd.DataFrame,
    feature_cols: list[str],
    clinical_cols: list[str],
    time_col: str,
    event_col: str,
    penalizer: float,
    l1_ratio: float,
    n_bootstrap: int,
    selection_threshold: float,
    random_state: int,
    clinical_penalty_factor: float,
    selection_scope: str = "immune_only",
    collinearity_cfg: Optional[dict] = None,
) -> pd.DataFrame:
    # Bootstrap stability for Cox mirrors the logistic path but uses penalized hazard model refits.
    rng = np.random.default_rng(random_state)
    selection_counts: dict[tuple[str, str, str], int] = {}
    for iteration in range(n_bootstrap):
        sample_index = rng.choice(df.index.to_numpy(), size=len(df), replace=True)
        bootstrap_df = df.loc[sample_index].reset_index(drop=True)
        try:
            _, coef_df = fit_cox_model(
                train_df=bootstrap_df,
                test_df=bootstrap_df,
                feature_cols=feature_cols,
                clinical_cols=clinical_cols,
                time_col=time_col,
                event_col=event_col,
                penalizer=penalizer,
                l1_ratio=l1_ratio,
                clinical_penalty_factor=clinical_penalty_factor,
                collinearity_cfg=collinearity_cfg,
            )
        except Exception:
            continue

        selected = coef_df[coef_df["abs_coef"] > selection_threshold].copy()
        if selection_scope == "immune_only":
            selected = selected[selected["feature_scope"] == "immune"]

        for row in selected.itertuples(index=False):
            key = (row.feature, row.raw_feature, row.feature_scope)
            selection_counts[key] = selection_counts.get(key, 0) + 1

    if not selection_counts:
        return pd.DataFrame(columns=["feature", "raw_feature", "feature_scope", "selection_count", "selection_frequency"])

    rows = []
    for (feature, raw_feature, feature_scope), count in selection_counts.items():
        rows.append(
            {
                "feature": feature,
                "raw_feature": raw_feature,
                "feature_scope": feature_scope,
                "selection_count": count,
                "selection_frequency": count / max(n_bootstrap, 1),
            }
        )
    return pd.DataFrame(rows).sort_values(["selection_frequency", "selection_count"], ascending=[False, False]).reset_index(drop=True)
