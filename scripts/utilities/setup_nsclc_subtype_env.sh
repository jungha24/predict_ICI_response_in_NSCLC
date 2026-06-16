#!/usr/bin/env bash

set -euo pipefail

ENV_NAME="${1:-nsclc-subtype}"
ENV_YAML="${2:-envs/nsclc-subtype.yml}"

if [[ ! -f "${ENV_YAML}" ]]; then
  echo "Environment file not found: ${ENV_YAML}" >&2
  exit 1
fi

if command -v mamba >/dev/null 2>&1; then
  SOLVER="mamba"
elif command -v conda >/dev/null 2>&1; then
  SOLVER="conda"
else
  echo "Neither mamba nor conda is available in this shell." >&2
  exit 1
fi

echo "[1/3] Removing an existing ${ENV_NAME} env if it already exists"
if conda env list | awk '{print $1}' | grep -Fxq "${ENV_NAME}"; then
  "${SOLVER}" env remove -y -n "${ENV_NAME}"
fi

echo "[2/3] Creating ${ENV_NAME} from ${ENV_YAML}"
"${SOLVER}" env create -n "${ENV_NAME}" -f "${ENV_YAML}"

echo "[3/3] Verifying key dependencies"
conda run -n "${ENV_NAME}" python -c "import scanpy, anndata, sklearn; print('scanpy', scanpy.__version__)"
conda run -n "${ENV_NAME}" Rscript -e "library(sva); cat('sva', as.character(packageVersion('sva')), '\n')"

cat <<EOF

Environment is ready.

Activate it with:
  conda activate ${ENV_NAME}

Or run the pipeline directly with:
  conda run -n ${ENV_NAME} python src/version1/run_subtype_pipeline.py --config configs/version1_subtype_base.yaml

EOF
