# Future Work: Sample/Run Config Refactor

## Goal

Reduce the volume and duplication of hand-maintained JSON configs across the Xenium workflow by separating stable sample metadata from run-specific workflow parameters.

This should be done on a dedicated branch, not mixed into active QC filtering or reclustering changes.

Suggested branch names:

```text
refactor/config-layering
config-sample-run-split
```

## Motivation

The current config tree has grown as workflow branches have multiplied. Many configs duplicate the same sample facts while changing only a few run-policy fields.

The duplication pattern is roughly:

```text
sample x script x area mode x percentile x transcript threshold x reclustering settings
```

This makes it harder to tell which configs are canonical, which are exploratory, and which are stale. It also increases the chance of copying an outdated `pc_label`, `cluster_col`, `coord_keys`, `group_map`, or run label into a new branch.

## Proposed Model

Split human-authored configuration into two layers:

```text
sample config = stable facts about a sample
run config    = parameters for an analysis branch or workflow policy
```

Resolved per-script configs can still exist, but they should be generated as provenance/audit artifacts rather than hand-maintained source files.

## Scope Update (2026-11-19)

The model below was written before scripts 08–10 and the post-08 helpers
existed. It should span the whole pipeline, 00 through 10, including the
spatial, poster, and conference helpers. Post-08 steps currently have the
deepest config chains (infiltration -> comparison -> script 10 -> script 08 ->
palette/mapping) and three different root overrides (`input_root`,
`--annotation-root`, `result_root`).

Layering decisions agreed so far:

- **Sample facts** include `response` and `tissue`. These currently live only
  in `config/10_squidpy_recluster_analysis/vbct/cell_type_annotation/nhood_enrichment_comparison.json`.
- **Selected analysis resolution is a branch/run setting, not a sample fact.**
  The post-08 resolutions (e.g. CK_skin_res 1.00, MF_skin_non_res 0.50) were
  chosen for the `min15_area95` recluster branch and differ by branch.
- **Reviewed artefacts** (cluster-to-cell-type mapping CSV, shared palette,
  area group map) are referenced by path from the run config, not copied.
- **Data roots:** one local and one HPC root, overridable the same way in every
  step, replacing the per-script root keys.
- The resolver belongs in the standalone workflow package
  (`config/resolver.py` in `decouple_workflow_from_banksy.md`), not in
  `banksy_utils/`, because the goal is to move workflow code out of this fork.
  Prefer building the layered config inside the new package over refactoring
  the fork first and moving it afterwards. Do this on a separate branch.

The interim post-08 cleanup is done on branch `cleanup/post08-flow` (see
`context/tasks/post08_flow_cleanup/post08_flow_cleanup.md`): `config/cohorts/vbct/min15_area95_selected_resolutions.json`
holds sample facts plus the selected per-branch resolution, read by scripts 08,
09, 10, 10a and the post-08 helpers through `workflow_config.py`, which also
provides the shared `--input-root`/`--output-root` options (replacing
`input_root`, `--annotation-root`, `result_root`, and
`project_root_parent_depth` for those steps). In this refactor, split the
cohort file into the sample and run layers, and fold `workflow_config.py` into
the package's `config/resolver.py`.

## Sample Configs

One canonical sample config per sample, shared across scripts 00–10, dotplot export, the post-08 helpers, and related workflow stages.

Possible location:

```text
config/samples/vbct/<sample>.json
```

Example contents:

```json
{
  "project": "vbct",
  "raw_subdir": "vbct",
  "dataset_name": "CK_skin_res",
  "pc_label": "30",
  "lambda_label": "0.20",
  "res_label": ["0.50", "0.60", "0.70", "0.80", "0.90", "1.00", "1.10"],
  "plot_res_label": "0.70",
  "nbr_weight_decay": "scaled_gaussian",
  "coord_keys": ["x", "y", "xy"],
  "cluster_col": "labels_scaled_gaussian_pc30_nc0.20_r0.70",
  "new_labels": {},
  "area_group_map": {}
}
```

Stable sample facts include:

- `project`
- `raw_subdir`
- `dataset_name`
- `pc_label`
- `lambda_label`
- default clustering resolutions from script 00/01
- plotting resolution
- coordinate keys
- selected cluster column
- reviewed annotations
- reviewed area group map, when available

## Run Configs

One run config per analysis branch or workflow policy.

Possible location:

```text
config/runs/vbct/<run_id>.json
```

Proposed examples (not yet created):

```text
grouped_area95_min_transcript_gte15.json
grouped_area99_default_min_transcript.json
unmerged_area95_default_min_transcript.json
unmerged_area99_default_min_transcript.json
```

Run configs should contain only branch-specific settings, split by script where needed.

Example:

```json
{
  "run_id": "grouped_area95_min_transcript_gte15",
  "script01": {
    "qc_output_label": "min_transcript_gte15",
    "min_transcript_threshold": 15,
    "min_transcript_pass_operator": "gte"
  },
  "script05": {
    "input_qc_label": "min_transcript_gte15",
    "output_label": "filtered_qc_v1_grouped_area95_min_transcript_gte15",
    "area_percentile_filter": {
      "mode": "cell_type_group",
      "percentile": 0.95,
      "groupby_key": "cluster_col",
      "group_map_key": "area_group_map"
    }
  },
  "script06": {
    "input_label": "filtered_qc_v1_grouped_area95_min_transcript_gte15",
    "run_label": "filtered_qc_v1_grouped_area95_min_transcript_gte15_qc_pass_only",
    "recluster_inclusion": "qc_pass_only",
    "res_label": ["0.50", "0.60", "0.70", "0.80", "0.90", "1.00", "1.10", "1.20", "1.30", "1.40", "1.50"]
  }
}
```

## Resolved Configs

Each script should still receive or write a fully resolved config containing all values needed to reproduce that exact run.

However, these resolved configs should be generated outputs rather than manually duplicated inputs.

Possible resolved config path pattern:

```text
data/xenium/output/<project>/<stage>/<sample>/<sample>_<script>_<run_id>_resolved_config.json
```

The resolved config should include:

- merged sample config
- selected run config section
- derived labels and paths
- exact area grouping used
- exact threshold settings used
- provenance copied from upstream objects where applicable

## Config Resolver

Add a small shared utility that can merge:

```text
sample config + run config + script-specific block
```

The resolver should support:

- loading a sample config by sample name
- loading a run config by run ID/path
- selecting a script block such as `script01`, `script05`, or `script06`
- resolving references like `groupby_key = "cluster_col"`
- resolving `group_map_key = "area_group_map"`
- writing the resolved config for audit/provenance

Utility location (updated 2026-11-19): the workflow package's
`config/resolver.py` (proposed; see `decouple_workflow_from_banksy.md`). The earlier
suggestions, `banksy_utils/config_resolver.py` and
`helper_scripts/config/resolve_workflow_config.py`, are superseded because
they keep workflow code inside the BANKSY fork.

## Migration Plan

Do this incrementally on a dedicated branch.

1. Keep all existing script-specific configs working.
2. Add the config resolver utility without changing existing scripts.
3. Add canonical sample configs for VBCT.
4. Add one run config for a known branch, such as `grouped_area95_min_transcript_gte15`.
5. Generate resolved configs for one sample and compare them to the current hand-written configs.
6. Update script 05 first, because it currently carries the most config branching complexity.
7. Validate one small sample end to end against the current config behavior.
8. Extend the pattern to scripts 01 and 06.
9. Update Slurm wrappers to accept `SAMPLE_CONFIG_DIR` plus `RUN_CONFIG` or a batch manifest.
10. Only after validation, archive or regenerate old duplicated config directories.
11. Update `tools/pipeline_status.py` to the new config layout (see "Follow-up: Pipeline Status Tool" below).

## Follow-up: Pipeline Status Tool (added 2026-11-20)

Task: once the refactor is complete, update `tools/pipeline_status.py`.
Decision (2026-11-20): do the refactor first, then extend the tool, rather
than adding per-step adapters to the current configs.

The tool reports, per sample in a branch, which pipeline outputs exist
(ok / MISSING / stale). It currently:

- starts from a script 06 config directory, links to script 05 configs through
  `output_label`, and to scripts 08/10 through `config/cohorts/` and the
  script 10 config that names the cohort;
- covers scripts 01, 05, 06, 08 and 10 only (not 00, 02, 03, 04, 09, 10a or
  the post-08 helpers);
- duplicates the output filename rules of scripts 05, 06 and 10, so it must
  be edited whenever those scripts rename outputs.

After the refactor:

1. Read the new sample and run configs (or resolved configs) instead of
   per-script config directories; take a run config or manifest as the
   branch argument.
2. Import output paths from the package's shared path module (planned
   `io/paths.py` in `decouple_workflow_from_banksy.md`) so the scripts and the
   tool use one set of filename rules; remove the tool's own copy.
3. Extend coverage to the remaining steps (00, 03, 04, 09, 10a, and the
   helpers that write outputs).
4. Keep its behaviour: read-only, `--root` for checking a mounted HPC
   checkout, refuses to run against an inactive mount, `--json` output.
5. Update the "Agent Tools" section of `AGENTS.md` if its usage changes.

## Slurm/Batches

Instead of using one directory of per-sample script configs, a future wrapper could use:

```bash
sbatch --array=0-7%2   --export=ALL,SAMPLE_CONFIG_DIR=config/samples/vbct,RUN_CONFIG=config/runs/vbct/grouped_area95_min_transcript_gte15.json,SCRIPT_BLOCK=script05   hpc/slurm/run_05_apply_qc_filters_for_reclustering.sl
```

Or use a proposed manifest:

```json
{
  "run_config": "config/runs/vbct/grouped_area95_min_transcript_gte15.json",
  "samples": [
    "BE_brain_non_res",
    "CK_bowel_res",
    "CK_skin_res",
    "EL_skin_res",
    "GR_lung_non_res_roi",
    "HW_brain_res",
    "MF_skin_non_res",
    "MG_gastric_non_res"
  ]
}
```

## Guardrails

- Do not remove or rewrite existing configs until replacement runs have been validated.
- Preserve current config loading as the fallback path during migration.
- Resolved configs should be written for every run so provenance remains easy to inspect without opening `.h5ad` objects.
- Keep biological/sample facts out of run configs where possible.
- Keep branch/run policy out of sample configs where possible.
- Do not mix this refactor with active changes to filtering semantics.

## Open Questions

- Should scripts accept `--config` only, with resolved configs generated before submission, or should scripts accept `--sample-config` plus `--run-config` directly?
- Should config inheritance be supported for run variants?
- Where should generated resolved configs live: under `data/xenium/output`, `logs`, or a dedicated `resolved_config` directory?
- Should existing script-specific config directories be archived after migration or retained as examples?
- Should dotplot/export configs join this system immediately or remain separate until scripts 01/05/06 are stable?
- Should the config refactor and the package split (`decouple_workflow_from_banksy.md`) be one branch/project? Current lean (2026-11-19): yes, build the layered config inside the new package.
