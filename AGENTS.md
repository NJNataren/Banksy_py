# Project Notes

## Project background
The aim of this project is to perform bioinformatic analysis of Xenium Spatial transcriptomic data of metastatic melanoma samples.The biological aim of the work is to identify the biological and mechanistic determinants of immune checkpoint inhibitor (ICI) response and to understand why some patients respond or do not respond to this treatment.

The lead human bioinformatician for this project is Nathalie.


## Repository Shape
- This is a research analysis fork of `Banksy_py` for spatial transcriptomics workflows, with core BANKSY code in `banksy/` and local analysis helpers in `banksy_utils/`.
- Xenium analysis notebooks and exported scripts live at the repo root, including `01_xenium_clustering.ipynb` and `01_xenium_clustering.py`.
- Runtime configuration is stored as JSON under `config/QC/` and `config/clustering/`, split by study/group and sample size.
- Slurm entrypoints live at the repo root, for example `run_xenium_clustering_test*.sl` and `run_xenium_QC_array_*.sl`.
- Data and generated outputs are expected under `data/xenium/`, `figures/`, `logs/`, and `hpc/`; treat these as potentially large or machine-specific.

## Environment
- Preferred environment files are `environment.yml` and `requirements.txt`.
- The conda environment name in `environment.yml` is `banksy-test`; the README refers to a generic `banksy` environment.
- Main scientific stack includes `scanpy`, `anndata`, `numpy`, `scipy`, `pandas`, `matplotlib`, `seaborn`, `igraph`, `leidenalg`, and `umap-learn`.

## Workflow Conventions
- Notebook work is often mirrored into `.py` exports. When changing analysis logic, check whether both notebook and script need to stay in sync.
- Clustering scripts are config-driven via `--config path/to/sample.json`.
- Avoid hard-coding sample names or paths when the existing script already reads from config.
- Existing scripts write directories with `os.makedirs`; preserve this style unless making a broader cleanup.
- Use `rg`/`rg --files` first for repo search.

## Git And File Safety
- The worktree may contain active notebook, script, and output changes. Do not revert or clean untracked files unless explicitly asked.
- Current active analysis files may include local edits to `01_xenium_clustering.ipynb`, QC notebooks/scripts, `banksy_utils/annotation_utils.py`, and Slurm scripts.
- Do not execute `git commit`, `git push` or `git pull` without permission.

## Validation
- There is no obvious test suite in this workspace. For code changes, prefer targeted syntax checks or small import checks when dependencies are available.
- For analysis script changes, a useful smoke test is usually argument parsing with one JSON config, but full execution may require local Xenium data and substantial compute.
