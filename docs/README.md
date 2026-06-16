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
3. 해석
- cell type proportion 신호와 state신호를 분리해서 볼 것

