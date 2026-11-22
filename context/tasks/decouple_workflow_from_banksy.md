# Future Work: Decouple Xenium Workflow From Banksy_py Fork

## Goal

Separate the project-specific Xenium analysis workflow from the `Banksy_py` source tree, so BANKSY can be used as an installed dependency and the analysis pipeline can be maintained as its own workflow package/repository.

This should be done as a dedicated future branch/project, not mixed into active QC filtering, reclustering, or config cleanup work.

Suggested branch names:

```text
refactor/decouple-banksy-workflow
workflow-package-split
xenium-workflow-package
```

## Motivation

The current repository mixes several layers:

```text
upstream BANKSY code
local BANKSY fork changes
Xenium workflow scripts
sample/run configs
HPC Slurm wrappers
notebooks
analysis outputs/context
```

This made sense during early workflow development, but it now makes ownership blurry:

- Which files are upstream BANKSY internals?
- Which files are local workflow logic?
- Which utilities are reusable project code?
- Which changes should be proposed upstream versus kept in the workflow?
- How should BANKSY itself be version-pinned for reproducibility?

Long term, the Xenium workflow should be independent of the BANKSY source tree.

## Target Model

Use BANKSY as a dependency:

```bash
pip install banksy
```

or, if a fork/commit is still required:

```bash
pip install git+https://github.com/<owner>/Banksy_py.git@<commit-or-branch>
```

Maintain the Xenium workflow separately:

```text
xenium-banksy-workflow/
  pyproject.toml
  src/xenium_banksy_workflow/
    qc/
    filtering/
    reclustering/
    plotting/
    config/
    io/
  scripts/ or console entry points
  config/
    samples/
    runs/
  hpc/
  context/
  notebooks/
```

The intended boundary:

```text
BANKSY package      = spatial clustering dependency
workflow package    = Xenium QC/filter/recluster/plot orchestration
configs             = project/sample/run metadata
HPC wrappers        = execution layer
outputs             = generated artifacts and provenance
```

## Candidate Package Structure

Possible Python package layout:

```text
src/xenium_banksy_workflow/
  __init__.py
  config/
    resolver.py
  io/
    paths.py
    anndata.py
  qc/
    script01_qc.py
    masks.py
  filtering/
    apply_filters.py
    area_percentile.py
    provenance.py
  reclustering/
    banksy_recluster.py
  plotting/
    qc_plots.py
    dotplots.py
  cli/
    qc.py
    apply_filters.py
    recluster.py
    export_dotplot.py
```

Current numbered scripts could become thin CLI entry points or wrappers around package functions.

Example commands:

```bash
xenium-qc --config resolved_script01_config.json
xenium-apply-qc-filters --config resolved_script05_config.json
xenium-recluster --config resolved_script06_config.json
xenium-export-dotplot --config resolved_script03_config.json
```

## Migration Plan

Do this incrementally.

1. Inventory imports from `banksy` and `banksy_utils` in the workflow scripts.
2. Identify which imports are true BANKSY API calls and which are local workflow conveniences.
3. Move project-specific helper logic out of the BANKSY tree into a new workflow package namespace.
4. Add a `pyproject.toml` for the workflow package.
5. Convert one low-risk script first, likely script 05, into a package-backed CLI while preserving the existing script entry point.
6. Pin BANKSY as an environment dependency using either a release version or a specific git commit.
7. Update Slurm wrappers to call installed workflow commands or `python -m xenium_banksy_workflow...`.
8. Validate one known sample/branch against the current repo-local script behavior.
9. Gradually migrate scripts 01, 06, 03/04 dotplot tooling, 08, 09, 10, and the post-08 helpers (`helper_scripts/spatial_analysis/`, `helper_scripts/poster/`, `helper_scripts/conference_figures/`). (Script 07 was renamed to 10 on 2026-10-25.)
10. Once stable, keep the BANKSY fork only for actual BANKSY source modifications.

## Relationship To Config Refactor

This is related to the sample/run config refactor, but should remain a separate project.

Config refactor goal:

```text
reduce duplicate workflow configs
```

Decoupling goal:

```text
separate workflow code from BANKSY source code
```

The config resolver could become part of the new workflow package, but the package split should not be required before improving the config model.

Update (2026-11-19): the current lean is to treat the two as one project on a
separate branch, building the layered config (proposed `config/resolver.py`) inside the
new package rather than refactoring the fork and then moving it. Do not add new
workflow code to `banksy_utils/`; it would deepen the coupling this task removes.

Notes for migration:

- `helper_scripts/spatial_analysis/spatial_common.py` has no BANKSY imports. It
  maps onto `config/` (sample loading), `io/anndata.py` (h5py label and
  coordinate reads), and `plotting/` (response colours, `save_figure`).
- `resolve_path` is duplicated in script 10 and `spatial_common.py`;
  `save_figure` closely matches script 10's `save_spatial_figure`. Consolidate
  them in the package.
- The cross-sample spatial helpers should become numbered sibling steps in the
  package once the spatial analysis extensions are finished; neither depends on
  the other (see `spatial_analysis.md`).
- Post-08 cleanup (2026-11-20, branch `cleanup/post08-flow`) added two
  BANKSY-free shared modules at the repo root that map straight onto the
  package: `workflow_config.py` (cohort loading, sample selection, input/output
  roots) -> `config/` and `io/paths.py`; `spatial_tables.py` (script 10
  enrichment table format) -> `io/`. `spatial_common.py` now imports both; its
  `resolve_path` duplicate is gone. Script 10 still has its own
  `save_spatial_figure`.

## Coupling Inventory and Effort Estimate (2026-11-21)

A quick inventory (discussion only; nothing implemented) showed that the coupling is
narrower than first assumed. This covers migration plan steps 1–2.

**Which scripts import BANKSY.** Only 3 of the ~18 workflow scripts (~14k
lines) import from `banksy/` or `banksy_utils/`:

- `00_xenium_clustering_clean_adata.py`
- `01_QC_xenium_spatial_clean_clustered.py`
- `06_recluster_qc_annotated_with_banksy.py`

Scripts 02–10, `workflow_config.py`, `spatial_tables.py` and `helper_scripts/`
have no BANKSY imports. (`banksy_tutorials/` also imports BANKSY but is not
workflow code.)

**What they import.**

- True BANKSY API (`banksy/`): `initialize_banksy`, `generate_banksy_matrix`,
  `run_Leiden_partition_parallel`, `concatenate_all`,
  `median_dist_to_nearest_neighbour`, `plot_results`.
- `banksy_utils/` helpers: `pca_umap`, `create_spatial_nonspatial_adata`,
  `pad_clusters`, `load_adata`, `filter_utils`, `object_downcasting_utils`,
  `plot_utils`. Most are workflow conveniences that belong in the new package.
  Decide case by case which stay with BANKSY.

**Catch: the fork changes BANKSY's own code.** `git log` shows local commits
to `banksy/` and `banksy_utils/` that upstream does not have:

- parallel Leiden in `banksy/cluster_methods.py`, used by scripts 00 and 06
  through `run_Leiden_partition_parallel`
- the `initialize_banksy` fix (`plt_edge_hist=False`, non-`scaled_gaussian`
  weights)
- a `weights_scheme` parameter on `create_spatial_nonspatial_adata`
- cluster colour fixes and configurable UMAP transparency
- a zero-count QC helper and the annotation helpers in `banksy_utils/`

The package therefore cannot depend on the published `banksy` release without
risking different results. Depend on this fork at a pinned commit (step 6), or
port the workflow-only changes into the package. Keep the fork only for the
genuine BANKSY changes (step 10). This partly answers the first open question
below.

**Revised effort estimate.** The decoupling itself is mostly mechanical:

- package layout and `pyproject.toml`
- move the scripts, keeping the numbered entry points as wrappers
- move the `banksy_utils` workflow helpers
- pin the fork commit
- switch the Slurm wrappers to `python -m`

That is about 3–6 focused days, or roughly one extra week on top of the config
refactor. Combined estimate: about 3–4 weeks of calendar time, down from an
earlier 4–8 week guess. Most of that time goes on HPC validation runs and review,
not coding.

**Validation shortcut.** Moving code without changing it should give
identical outputs. One CK_skin_res run through scripts 00 and 06, compared with
the current outputs, is the main check. First confirm that BANKSY clustering is
deterministic with a fixed seed, because the parallel Leiden path may not be.

**Branching.** Doing this together with the config refactor in one project
remains reasonable. Merge phase by phase (as with `cleanup/post08-flow`) to
avoid long-lived conflicts with `config/` changes on main.

## Guardrails

- Do not break existing numbered scripts during migration.
- Keep old script paths available as compatibility wrappers until HPC workflows are updated.
- Pin the BANKSY dependency by version or commit for reproducibility.
- Keep generated outputs and large data outside the package source tree.
- Keep project-specific sample/run configs outside the installed Python package when possible.
- Avoid modifying BANKSY internals unless the change genuinely belongs upstream or in a maintained fork.
- Validate each migrated script against a known sample before moving the next one.

## Open Questions

- Is BANKSY available as a stable pip-installable package for the exact functionality needed, or do we need to pin a git commit/fork? (Lean, 2026-11-21: pin the fork, because it modifies `banksy/`; see "Coupling Inventory".)
- Which current `banksy_utils` functions are upstream utilities versus local workflow helpers? (The imported set is listed under "Coupling Inventory".)
- Should the workflow live in a new repository or remain in this repo while treating BANKSY as an external dependency?
- Should numbered script names be preserved as wrappers for continuity?
- Should HPC environments install the workflow package editable from source or from a tagged release?
- How should notebooks reference the package during transition?

## Acceptance Criteria

- Xenium workflow scripts live outside the BANKSY source tree or are clearly namespaced as workflow code.
- BANKSY is installed as a dependency rather than used implicitly because scripts are inside the fork.
- Slurm wrappers can run the workflow through installed commands or module entry points.
- Existing QC filtering/reclustering behavior is reproduced for at least one validated sample branch.
- BANKSY version or commit is pinned in the environment for reproducibility.
- The old repo-local workflow remains available until the new package path is validated.
