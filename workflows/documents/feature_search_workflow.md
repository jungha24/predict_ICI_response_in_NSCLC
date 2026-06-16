## Recommended structure

- `configs/`: versioned YAMLs only for this search pipeline
- `src/version2/feature_search/`: modular Python package for Stage 2-5
- `src/version2/run_feature_search.py`: thin CLI entrypoint


## What this v2 pipeline does

- **Stage 2**: single-feature add-on scan
  - baseline = `clinical_only`
  - test = `clinical + feature_j`
  - nested CV delta metric vs baseline per endpoint
- **Stage 3**: redundancy pruning
  - top-ranked features per endpoint
  - biological family cap
  - pairwise correlation pruning
  - VIF pruning
- **Stage 4**: subset search
  - exhaustive search for 2-3 feature subsets
  - heuristic beam search for 4+ feature subsets
- **Stage 5**: best-subset refits
  - full-data coefficients
  - bootstrap stability
  - coefficient stability across outer folds
- **Outer validation**: escape selection bias


## Multiple per-celltype feature CSVs (stage 1)

This pipeline can now read **multiple patient-level feature CSVs**, one per cell type or feature family.
Use `feature_library.tables` (preferred) or `feature_library.paths`. Each table must contain one sample/patient ID column and numeric feature columns.

When multiple tables are provided, feature names are automatically prefixed with the table alias/stem to avoid collisions, for example:
- `B__Centroid_mean_score`
- `CD14_Mono__curated_gene_hallmark_interferon_gamma_response__singscore`

Only samples listed in `feature_library.selected_patients_file` are kept before modeling.

## Run

```bash
PYTHONPATH=src conda run -n nsclc-subtype python src/version2/run_feature_search.py --config configs/version2_modeling_base.yaml
```
## script
1. entrypoint: run_feature_search.py 
2. 입출력 유틸: feature_search/io_utils.py
3. clinical/feature 결합: feature_search/data.py
4. 설계 정의: design.py
5. 모델 엔진: models.py
  - training-fold 기준 imputation/standardization/one-hot
  - zero-variance, optional corr/VIF pruning << pruning
  - nested CV elastic-net logistic/Cox
  - bootstrap stability
6. 실제 탐색 흐름: search.py
  - stage 2: baseline clinical model을 endpoint별로 돌리고, candidate immune feature을 하나씩 더한 걸 전부 평가해 baseline 대비 delta metric저장
  - stage 3: single-feature ranking 상위권에서 family cap을 걸고, pairwise corr와 VIF로 중복 feature 줄여 최종 후보 남기기 << pruning
  - stage 4: subset size 2,3은 exhaustive, size 4는 beam-search heuristic으로 탐색
  - stage 5: endpoint별 best subset 상위 3개를 full-data로 다시 튜닝/적압하고 coefficient, bootstrap sbtility, fold-wise coefficient stability저장


## version update
2026.04.03:
  - use pd-l1 tps as baselinie and change family cap
  - src/version2/run_feature_search_v2.py
  - configs/version2_modeling_base_v2.yaml
  - 기존의 stage2결과를 가져와서 stage 3부터 돌리기 가능하게 함. outer validation은 off, baseline: clinical_only

OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH=src conda run -n nsclc-subtype python src/version2/run_feature_search_v2.py --config configs/version2_modeling_base_v2.yaml


## history
2026.04.01 smoke test; OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH=src conda run -n nsclc-subtype python src/version2/run_feature_search.py --config configs/version2_modeling_smoketest.yaml
2026.04.01 full test; nohup env OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH=src conda run -n nsclc-subtype python src/version2/run_feature_search.py --config configs/version2_modeling_base.yaml > data/20260309_pilot/results/version2/feature_search_base_logs/stdout.log 2>&1 &

2026.04.03 20260401의 stage2 결과에서 family cap 정의만 바꾼 결과 확인; nohup env OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH=src conda run -n nsclc-subtype python src/version2/run_feature_search_v2.py --config configs/version2_modeling_base_v2.yaml --reuse-stage2-dir data/20260309_pilot/results/version2/feature_search_base > data/20260309_pilot/results/version2/feature_search_base_v2_logs/stdout.log 2>&1 &

2026.04.03 
mkdir -p data/20260309_pilot/results/version2/feature_search_base_v2_logs
nohup env OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH=src conda run -n nsclc-subtype python src/version2/run_feature_search_v2.py --config configs/version2_modeling_base_v2.yaml > data/20260309_pilot/results/version2/feature_search_base_v2_logs/stdout.log 2>&1 &

2026.04.06
nohup env OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH=src conda run -n nsclc-subtype python src/version2/run_feature_search_v2.py --config configs/version2_modeling_base_v2_single_feature_outer.yaml > data/20260309_pilot/results/version2/feature_search_base_v2_single_feature_outer_logs/stdout.log 2>&1 &

2026.04.08
conda run -n nsclc-subtype python src/version2/report_repeated_outer_top_features.py \
  --run-dir data/20260309_pilot/results/version2/feature_search_base_v2_single_feature_outer \
  --top-n 20

2026.04.29 local v2.0.1 archive integration:
  - source folder:
    analysis/n73_manual_trial_20260305/version2/v.2.0.1
  - 3-fold outer validation에서 각 fold의 Stage 2 single-feature scan 상위 30개를 요약
  - retained tables:
    workflows/version2/results/feature_search/base_v2_single_feature_outer/outer_search_validation/repeated_top_features__top30.csv
    workflows/version2/results/feature_search/base_v2_single_feature_outer/outer_search_validation/repeated_top_features__top30_details.csv
  - exploratory clustering branch:
    1. repeated_top_features__top30.csv를 mean/group rank 순서로 정렬
    2. representative top 30 feature를 patient-level matrix에서 선택
    3. missing value는 feature별 median으로 imputation
    4. clinical covariates와 immune features를 합쳐 FAMD 수행
    5. FAMD Dim2-Dim5를 embedding으로 사용
    6. k-nearest-neighbor graph 생성, k = 5
    7. igraph Louvain clustering 수행, resolution = 1
    8. UMAP 수행, n_neighbors = min(10, n - 1), min_dist = 0.3, metric = euclidean
  - retained exploratory figures:
    workflows/version2/figures/exploratory_rplots/20260427_version2_trial1_fold_group_rank.pdf
    workflows/version2/figures/exploratory_rplots/20260429_outervalidation_fold_rocauc.pdf
    workflows/version2/figures/exploratory_rplots/20260429_version2_trial1_FAMD_FD12.pdf
    workflows/version2/figures/exploratory_rplots/20260429_version2_trial1_FAMD_FD23.pdf
    workflows/version2/figures/exploratory_rplots/20260429_version2_trial1_FAMD_FD1_explain.pdf
    workflows/version2/figures/exploratory_rplots/20260429_version2_trial1_FAMD_FD2_explain.pdf
    workflows/version2/figures/exploratory_rplots/20260429_version2_trial1_FAMD_umap.pdf
    workflows/version2/figures/exploratory_rplots/20260429_version2_trial1_FAMD_umap_binarized_response.pdf
    workflows/version2/figures/exploratory_rplots/20260429_version2_trial1_FAMD_umap_PFS_6mo.pdf
    workflows/version2/figures/exploratory_rplots/20260429_version2_trial1_FAMD_umap_RECIST.pdf
    workflows/version2/figures/exploratory_rplots/20260429_version2_trial1_FAMD_umap_PDL1_TPS.pdf
  - raw patient-level feature matrix는 Git 포함 제외
