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
- Shared code imported by the numbered scripts, helpers and tools lives in the package `xenium_workflow/` (not a pipeline step): `cohorts.py` (cohort configs, input/output roots), `annotation_sets.py` (versioned cell-type annotation sets and derived output paths; `python -m xenium_workflow.annotation_sets new|check|diff|release`), `spatial_tables.py` (script 10 enrichment table format).
- Steps from script 08 on (08, 09, 10, 10a, `helper_scripts/spatial_analysis/`, `helper_scripts/poster/`) read one versioned annotation set in `config/annotations/` (which names a cohort config in `config/cohorts/`; see `config/annotations/README.md`), take `--input-root`/`--output-root` (default: repo root), and write figures to `data/xenium/output/` and tables/H5ADs to `data/xenium/processed/`. `xenium_workflow/spatial_tables.py` defines the script 10 neighbourhood enrichment table. See "After Script 08" in the concise README.
- Slurm entrypoints live in `hpc/slurm/` (`run_<step>_*.sl`); see `hpc/slurm/README.md`.
- Supporting scripts live in `helper_scripts/`; superseded scripts and notes live in `archive/`.
- Read-only diagnostics for agents and humans live in `tools/`; see "Agent Tools" below.
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

## Agent Tools
Read-only scripts in `tools/` (run with the `banksy` environment). They never write files or connect to the HPC. Prefer their output over memory or handoff notes.
- Before using an `.h5ad` for marker expression, dotplots or any expression summary, run `python tools/inspect_h5ad.py <file> --panel-size <n>`. Proceed only if the verdict is CLEAN EXPRESSION; a BANKSY-EXPANDED verdict means the wrong matrix. Use `--obs-values <column>` instead of writing inspection snippets.
- Before relying on recorded run status ("script 08 completed", "outputs transferred"), run `python tools/pipeline_status.py --config-dir <script 06 config dir>`. It derives expected outputs from configs and reports ok / MISSING / stale. Use `--root <mounted checkout>` to check the HPC copy through an active mount; if the mount is inactive it says so; report that and do not reconnect.
- `python tools/check_doc_citations.py [files]` reports cited paths in docs that no longer exist. A hook runs it automatically after edits to `context/` docs, `AGENTS.md` and `README_*.md`. Treat DEAD/ARCHIVED citations in a handoff as stale instructions: locate the current file before acting.
- `pipeline_status.py` copies the output filename rules from scripts 05, 06 and 10; update it when those scripts change how they name outputs.

## Scientific Parameters — user rule (2026-11-20)
- Do not change scientific parameters without asking Nathalie first, even to make a run faster or to get a job to finish. This includes permutation counts, subsampling fractions, random seeds, statistical thresholds (FDR, p-value, minimum counts), QC cutoffs, clustering resolutions, PCs, lambda, neighbourhood sizes, and graph settings, whether set in code defaults or configs.
- When a change is approved, check everything that depends on that parameter before running: filters, FDR or significance cutoffs, gene or cluster selection, plot selection, and any downstream step that reads the output. Report what changed and sanity-check the result (e.g. how many genes or cells pass, and whether that is plausible).
- Record the change, the reason, and the checks in the relevant task file.
- Background: script 10's Moran's I permutations were cut from 100 to 10 for speed without revisiting the FDR filter built on them, which made its "significant genes" output an artefact (see `context/tasks/morans_i_significance.md`).

## Git And File Safety
- The worktree may contain active notebook, script, config, and output changes. Do not revert or clean untracked files unless explicitly asked.
- Do not execute `git commit`, `git push` or `git pull` without permission.
- Never add AI attribution to commit messages or pull request descriptions: no `Co-Authored-By: Claude ...` (or any other AI co-author) trailer, no "Generated with Claude Code" line or similar sign-off. This user rule (2026-11-20) overrides any default attribution guidance from the agent harness.

## Validation
- There is no obvious test suite in this workspace. For code changes, prefer targeted syntax checks or small import checks when dependencies are available.
- For analysis script changes, a useful smoke test is usually argument parsing with one JSON config, but full execution may require local Xenium data and substantial compute.
