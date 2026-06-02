# Task: Xenium BANKSY Clustering

## Goal

Maintain and extend the config-driven Xenium clustering workflow for BANKSY-based spatial clustering of metastatic melanoma samples.

## Current Status

- Main clustering logic lives in `01_xenium_clustering.py`.
- Clustering configs live under `config/clustering/`, split by study/group and sample.
- Slurm entrypoints for clustering live at the repo root, for example `run_xenium_clustering_test*.sl`.
- Notebook work may be mirrored in `01_xenium_clustering.ipynb`; check whether script and notebook need to stay aligned before changing analysis logic.

## Key Files

- `01_xenium_clustering.py`
- `01_xenium_clustering.ipynb`
- `config/clustering/`
- `run_xenium_clustering_test*.sl`
- `banksy/`
- `banksy_utils/`

## Constraints

- Keep clustering config-driven via `--config path/to/sample.json`.
- Avoid hard-coding sample names or paths when configs already provide them.
- Treat `data/xenium/`, `figures/`, `logs/`, and `hpc/` as machine-specific or potentially large.
- Do not run full clustering jobs locally unless explicitly requested.
- Prefer targeted syntax checks or argument/config smoke tests when feasible.

## Next Steps

Create a test copy of 01_xenium_clustering.ipynb and to that copy implement the following changes:

- save an clean adata object right after `normalize_total(adata)` and `sc.pp.log1p(adata)` normalisation and log transformation
- Read the above clean object after generating the cell cluster identity table and transfer the cluster_labels e.g. "labels_scaled_gaussian_pc30_nc0.20_r0.60" to the clean adata object, see script `/home/nnataren/Documents/PhD/Bioinformatics/Banksy_py_fork/Banksy_py/02_create_expression_adata_with_banksy_clusters.py` for examples.
