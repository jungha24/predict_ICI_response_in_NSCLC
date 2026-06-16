# 20260314
# in docker
# pip install pandas numpy scipy scikit-learn umap-learn matplotlib pyyaml pyarrow jupyter
# python src/version1/run_subtype_pipeline.py --config configs/version1_subtype_base.yaml


# 20260316
# bash scripts/setup_nsclc_subtype_env.sh
# conda activate nsclc-subtype
# python src/version1/run_subtype_pipeline.py --config configs/version1_subtype_base.yaml

# check feature importance
# python src/version1/permutation_feature_importance.py \
#  --run-dir data/20260309_pilot/results/version1/subtype_discovery/base_run_v7 \
#  --n-permutations 50
# conda env update -n nsclc-subtype -f envs/nsclc-subtype.yml 
# conda run -n nsclc-subtype python src/version1/permutation_feature_importance.py \
#  --run-dir data/20260309_pilot/results/version1/subtype_discovery/base_run_v7 \
#  --n-permutations 50

# check association with PD-1
# conda run -n nsclc-subtype python src/version1/cluster_clinical_association.py \
#  --run-dir data/20260309_pilot/results/version1/subtype_discovery/base_run_v7


# results/run_base
# - collection_batch, age, sex,
# - nPC: 1-5

# results/run_base_v2
# - collection_batch
# - nPC: 1-2

# results/run_base_v3
# - collection_batch
# - nPC: 1-2
# - use centered log-ration (CLR) for proportion
# - use highly variable gene set

# results/run_base_v4
# - pc_selection_file: run_base_v4/pca_pc_selection.json
# - umap_min_dist: 0.1
# - umap_n_neighbors: 10

# results/run_base_v5
# - use only proportion 
# - umap_min_dist: 0.1
# - umap_n_neighbors: 10

# results/run_base_v6
# - use only PCs
# - pc_selection_file: run_base_v4/pca_pc_selection.json
# - umap_min_dist: 0.1
# - umap_n_neighbors: 10

# results/run_base_v7
# - use proportion and PCs from naive/TCM CD4+T
# - pc_selection_file: run_base_v7/pca_pc_selection.json
# - umap_min_dist: 0.1
# - umap_n_neighbors: 10


# results/pc_only_celltype_sweep
# python scripts/run_celltype_pc_importance.py   --base-config configs/version1_subtype_base.yaml   --run-root data/20260309_pilot/results/version1/subtype_discovery/pc_only_celltype_sweep   --celltypes "effector CD8+T" "CD14 Mono" "NK" "naive/TCM CD4+T"   --npcs 10   --n-permutations 50   --seed 123
# use different nPC per cell type
# conda run -n nsclc-subtype python scripts/run_celltype_pc_importance.py \
#  --base-config configs/version1_subtype_base.yaml \
#  --run-root data/20260309_pilot/results/version1/subtype_discovery/pc_only_celltype_sweep \
#  --celltypes "effector CD8+T" "CD14 Mono" "NK" "naive/TCM CD4+T" \
#  --pc-selection-file data/20260309_pilot/results/version1/subtype_discovery/base_run_v4/pca_pc_selection.json \
#  --n-permutations 50 \
#  --seed 123
# write in commend
# conda run -n nsclc-subtype python scripts/run_celltype_pc_importance.py \
#  --base-config configs/version1_subtype_base.yaml \
#  --run-root data/20260309_pilot/results/version1/subtype_discovery/pc_only_celltype_sweep \
#  --celltypes "effector CD8+T" "CD14 Mono" "NK" "naive/TCM CD4+T" \
#  --celltype-npcs "effector CD8+T=10" "CD14 Mono=6" "NK=4" "naive/TCM CD4+T=8" \
#  --n-permutations 50 \
#  --seed 123
