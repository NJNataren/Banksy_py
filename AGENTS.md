# Project Notes

## HPC access restriction — user rule (2026-10-25)

- Never initiate a direct SSH connection to Phoenix or any HPC system, including login tests, authentication probes, remote commands, job submission, or job monitoring. Do not use SSH-backed tools such as scp, sftp, rsync over SSH, or create/reconnect SSHFS mounts as a workaround.
- Existing local HPC mounts may be used to inspect files. Permission to inspect a mount does not authorize remote execution or changes to mounted files.
- Prepare and review HPC scripts and provide commands for Nathalie to run herself. Do not connect to the HPC to execute them, even when asked generally to launch or run an analysis.
- If a mount is unavailable, report that; do not attempt an SSH connection or reconnect it. Phoenix may perceive automated connection attempts as an attack.
- This rule supersedes historical handoffs suggesting SSH authentication, retries, or agent-submitted HPC jobs.


## Project background
The aim of this project is to perform bioinformatic analysis of Xenium Spatial transcriptomic data of metastatic melanoma samples.The biological aim of the work is to identify the biological and mechanistic determinants of immune checkpoint inhibitor (ICI) response and to understand why some patients respond or do not respond to this treatment.

The lead human bioinformatician for this project is Nathalie.


## Repository Shape
- This is a research analysis fork of `Banksy_py` for spatial transcriptomics workflows, with core BANKSY code in `banksy/` and local analysis helpers in `banksy_utils/`.
- The Xenium pipeline is a set of numbered scripts at the repo root (`00_` to `10_`, plus lettered sub-steps such as `00a_`, `001b_`, `01a_`, `02a_`). `README_00_04_XENIUM_ANALYSIS_WORKFLOW_CONCISE.md` is the map of what each step does and where it runs.
- A few steps also have a notebook version (`00_`, `01_`, and a `05_` local test notebook); the rest are `.py` only.
- Runtime configuration is stored as JSON under `config/<step>/` (e.g. `config/00_clustering/`, `config/01_QC/`, `config/06_recluster_qc_annotated/`), split by study/group and sample.
- Steps after script 08 (08, 09, 10, 10a, `helper_scripts/spatial_analysis/`, `helper_scripts/poster/`) read one cohort config in `config/cohorts/` through `workflow_config.py` (repo root), take `--input-root`/`--output-root` (default: repo root), and write figures to `data/xenium/output/` and tables/H5ADs to `data/xenium/processed/`. `spatial_tables.py` (repo root) defines the script 10 neighbourhood enrichment table. See "After Script 08" in the concise README.
- Slurm entrypoints live in `hpc/slurm/` (`run_<step>_*.sl`); see `hpc/slurm/README.md`.
- Supporting scripts live in `helper_scripts/`; superseded scripts and notes live in `archive/`.
- Task handoffs live in `context/tasks/`; load them via the `banksy-context` skill.
- Data and generated outputs are expected under `data/xenium/`, `figures/`, `logs/`, and `hpc/`; treat these as potentially large or machine-specific.

## Environment
- Preferred environment files are `environment.yml` and `requirements.txt`.
- The conda environment for this work is `banksy` (locally and in the Slurm scripts).
- Main scientific stack includes `scanpy`, `anndata`, `numpy`, `scipy`, `pandas`, `matplotlib`, `seaborn`, `igraph`, `leidenalg`, and `umap-learn`.

## Workflow Conventions
- Where a step has both a notebook and a `.py` version, check whether both need to stay in sync when changing analysis logic.
- Pipeline scripts are config-driven via `--config path/to/sample.json`.
- Avoid hard-coding sample names or paths when the existing script already reads from config.
- Existing scripts write directories with `os.makedirs`; preserve this style unless making a broader cleanup.
- Use `rg`/`rg --files` first for repo search.

## Git And File Safety
- The worktree may contain active notebook, script, config, and output changes. Do not revert or clean untracked files unless explicitly asked.
- Do not execute `git commit`, `git push` or `git pull` without permission.
- Never add AI attribution to commit messages or pull request descriptions: no `Co-Authored-By: Claude ...` (or any other AI co-author) trailer, no "Generated with Claude Code" line or similar sign-off. This user rule (2026-11-20) overrides any default attribution guidance from the agent harness.

## Validation
- There is no obvious test suite in this workspace. For code changes, prefer targeted syntax checks or small import checks when dependencies are available.
- For analysis script changes, a useful smoke test is usually argument parsing with one JSON config, but full execution may require local Xenium data and substantial compute.
