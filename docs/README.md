# file structure
1. notebook은 얇게
    - config load/override 여기에서
2. src는 라이브러리 함수로 순수하게 분석 로직만 갖게 하기
3. config 파일을 만들어 src내부의 parameter을 외부설명+실행 스냅샷으로 관리
    - RunConfig 기본값은 유지 (fallback)
    - 실제값은 configs/*.json에서 읽어 override
    - 실행할 때마다 최종설정을 logs/config_history.jsonl에 append
    - 저장하는 adata.uns['run_config']에도 같이 넣기

# parameter manipulation
1. cfg(기본 운영값, config 파일)
- 정의: RunConfig
- 파일 입력: configs/nsclc_v1.json
- 포함 항목:
    - min_counts, max_counts, min_genes, max_genes
    - max_pct_mt, max_pct_ribo, max_pct_hb
    - expected_doublet_rate
- 각 필터 on/off용 use_*_filter들

2. params(노트북 즉석 override)
- 정의: QCParams (in src/qc.py)
- 포함 항목:
    - 위와 동일한 threshold 항목들 전부
    - use_*_filter override도 가능

# things to consider for PBMC data
1. basic QC
- n_counts
- n_genes
- precentage_mt
- doublet score
2. comtamination
- RBC: hemoglobin gene; percentage_hb
- platelet: ambient RNA 형태로 섞일 수 있어 (platelet marker활용)
- ambient RNA
3. pre-analytical variable
- PBMC 채혈 후 isolation까지의 시간
4. 해석
- cell type proportion 신호와 state신호를 분리해서 볼 것

# version 2 protocol note
Source:
`analysis/n73_manual_trial_20260305/version2/protocol_v2.0.1.pages`

읽을 수 있는 본문 문자열과 v2 companion README 기준으로 정리:

1. Stage 1 feature library
    - patient-level feature library를 먼저 만든 뒤 Python search pipeline에 넣는 구조
    - feature type:
        - cell-level composition
        - potency/dynamics 또는 CellRank-like transition/priming surrogate
        - curated gene/pathway scores
        - pseudobulk program/module
        - interaction surrogate
        - selected latent axes, e.g. PC/FAMD summaries
    - 주요 cell-type block:
        - B lineage
        - monocyte
        - NK
        - CD4 T
        - CD8 T
        - nonconventional T
    - QC/filtering:
        - selected patient만 유지
        - min patients per feature, non-zero/detection 기준, unstable column 제거
        - 여러 CSV 병합 시 table alias/stem prefix로 feature collision 방지

2. Stage 2-5 feature search
    - Stage 2: clinical baseline + single immune feature add-on scan
    - Stage 3: family cap, pairwise correlation, VIF pruning
    - Stage 4: small subset exhaustive search and larger subset beam search
    - Stage 5: best subset refit, coefficient, bootstrap stability, fold-wise coefficient stability
    - outer validation: feature selection을 outer fold 안에서 반복해서 selection bias를 줄이는 목적

3. version 2.0.1 update
    - PD-L1 TPS를 clinical baseline에 포함
    - biological family cap 정의 변경
    - 기존 Stage 2 결과를 재사용해서 Stage 3부터 다시 실행 가능
    - 주요 config:
        - `configs/version2_modeling_base_v2.yaml`
        - `configs/version2_modeling_base_v2_single_feature_outer.yaml`
