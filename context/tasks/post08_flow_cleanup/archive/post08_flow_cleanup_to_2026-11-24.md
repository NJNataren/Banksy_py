> **Archived 2026-11-24.** Verbatim copy of `context/tasks/post08_flow_cleanup.md`
> (now `context/tasks/post08_flow_cleanup/post08_flow_cleanup.md`) before it was
> converted to the current-state layout. Historical only: paths, run status and
> next steps here may be stale (shared modules moved to `xenium_workflow/` and
> outputs gained annotation-set folders on 2026-11-22).

# Task: Untangle the Pipeline Flow After Script 08

## Goal

Simplify control flow, configs, data roots, and output locations for every step
downstream of script 08 (08 itself, 09, 10, the cross-sample spatial helpers,
and the poster/conference helpers) inside this fork, without changing any
scientific result.

This is an interim cleanup. The full sample/run config refactor and the package
split stay separate future work (`config_refactor.md`,
`decouple_workflow_from_banksy.md`); this cleanup must stay compatible with
them, so its cohort config later becomes their sample + run layers.

Branch: `cleanup/post08-flow` (not `main`).

Related context: `spatial_analysis.md`, `cell_type_annotation.md`,
`hca_poster.md`, `conference_figures.md`.

## Current Status (2026-11-20)

All phases (0–5) are done on branch `cleanup/post08-flow`. Next: push, then
Nathalie's one-sample HPC test, then merge to `main`. Do not merge before the
HPC test passes (Nathalie, 2026-11-20).

Commits on the branch (no AI attribution trailers; see AGENTS.md rule):
`adced39` phase 1, `9cb6b47` phase 2, `7c79e34` phase 3, `c269407` Moran seed,
`74a6b4d` phase 4, `0d9e8c2` AGENTS.md attribution rule, `2b34089` phase 5
docs. The branch history was rewritten once (2026-11-20) to strip
`Co-Authored-By` trailers and force-pushed with `--force-with-lease`.
`a6f3c4e` on `main` still has a trailer; left deliberately (rewriting `main`
would desync the HPC checkout).

Phase 5 as built: concise README "After Script 08" section (plain-text
diagram, run table, commands, history) and rewritten "Where Things Run"; full
README sections for 08, 09, 10, 10a; `config/spatial_analysis/README.md` and
`config/poster/README.md`; AGENTS.md repository shape. Task files:
`spatial_analysis.md` current sections rewritten in place; dated "Path
Update" notes in `cell_type_annotation.md`, `hca_poster.md`,
`conference_figures.md`, `filtering.md`; `config_refactor.md` and
`decouple_workflow_from_banksy.md` updated.

Old figures in `figures/` (cell_type_annotation, cell_type_proportions,
spatial_analysis, conference poster) stay where they are (decision
2026-11-20). After the reruns refill the new locations, move the non-poster
ones to a dated archive such as `figures/_pre_post08_cleanup_2026-11-20/`; leave
the poster folder until it is confirmed no design file links into it.

Follow-ups outside this cleanup: conference-figure helpers onto the cohort
config and output rule; promote spatial helpers in the workflow package.

HPC test plan (Nathalie runs): on Phoenix, `git fetch && git checkout
cleanup/post08-flow`, then `sbatch --array=0 hpc/slurm/run_10_cell_type_spatial_analysis.sl`
(CK_skin_res). Check: job succeeds; figures under
`data/xenium/output/vbct/spatial_analysis/...`; tables, result H5AD, and
`_nhood_enrichment.csv` under `data/xenium/processed/vbct/spatial_analysis/...`.
Locally: compare its CSVs with a local run (should be byte-identical given
both seeds), then `10a --overwrite` and rerun the comparison helper. Note the
HPC run overwrites the existing CK_skin_res result H5AD (unseeded) on Phoenix.

Moran seed (2026-11-20, Nathalie's request): script 10 passes `seed=1234` to
`sq.gr.spatial_autocorr`. Two local CK_skin_res runs are byte-identical.
Stored HPC results were unseeded, so simulated p-values from a rerun will
differ from them.

Phase 4 as built:

- Renamed `10_squidpy_recluster_spatial_analysis.py` -> `10_cell_type_spatial_analysis.py`;
  renamed `config/10_squidpy_recluster_analysis/vbct/cell_type_annotation/` ->
  `config/10_cell_type_spatial_analysis/vbct/`. Output file names unchanged
  (still `..._pcannotated_nccell_type_r<res>`).
- Legacy cluster mode, its 8 + 2 smoke configs, and
  `run_10_squidpy_recluster_spatial_analysis.sl` removed from the branch;
  local snapshot with README in `archive/2026-11-20_script10_legacy_cluster_mode/`
  (archive/ is gitignored). Script 10 without a `cohort` key raises a pointer
  to the archive.
- Docs updated: concise and full workflow READMEs (full README script 10
  section rewritten), `hpc/slurm/README.md`, script 10 config README.
- Validation: validate-only 8/8; wrapper tasks 0/7 -> CK_skin_res/CK_bowel_res;
  CK_skin_res local run byte-identical to pre-rename run; 10a validates 8/8
  stored results on Phoenix.

Phase 3 as built:

- Figures: `data/xenium/output/vbct/<step>/min15_area95_selected_resolutions/`
  with steps `cell_type_annotation` (08, plus `clean_scanpy/`),
  `cell_type_proportions` (09), `spatial_analysis` (10 and helpers), `poster`.
- Tables/H5ADs: same pattern under `data/xenium/processed/vbct/`. Script 09
  has a new `table_dir`. Script 10 writes centrality/Moran CSVs to
  `<result_dir>/<sample>/<branch>_r<res>/`.
- Configs moved: poster -> `config/poster/vbct/`; spatial helpers ->
  `config/spatial_analysis/vbct/`.
- HPC staging list now includes `workflow_config.py` (renamed to `xenium_workflow/cohorts.py` 2026-11-22), `spatial_tables.py`, and
  the cohort config.
- Existing figures in `figures/` were not moved in phase 3 (archived 2026-11-21; see Session History).
- Out of scope: conference-figure helpers (separate tool, also used by pre-08
  configs, mixes summary CSVs with figures; own sample list and palette copy).
- Validation: all 14 comparable tables byte-identical to baseline; output
  split clean. Local script 10 run on CK_skin_res (scratch) reproduced stored
  centrality CSV, FDR-significant Moran genes, and enrichment edge counts.

Mount note (2026-11-20): `mountpoint` reports the Phoenix sshfs mount as "not a
mountpoint" from inside the agent sandbox even when it is connected. Use
`findmnt -T ~/mnt/phoenix_scratch` or list a known file. Earlier notes in this
file saying the mount was disconnected on 2026-11-19 were wrong.

Phase 2 as built (commit hashes in older notes predate the 2026-11-20 rewrite; see the list above):

- `spatial_tables.py` (repo root): table format, `enrichment_from_labels`,
  `write_enrichment_table`, `read_enrichment_table` (round-trip float parsing).
- Script 10 writes `<result_dir>/<sample>/<stem>_nhood_enrichment.csv` after
  each result H5AD.
- `10a_export_nhood_enrichment_tables.py`: backfill from existing result H5ADs
  (`--input-root ~/mnt/phoenix_scratch/Banksy_py`), local writes only. Run on
  2026-11-20: all eight tables written locally under
  `data/xenium/processed/vbct/spatial_analysis/min15_area95_selected_resolutions/<sample>/`.
- Comparison helper reads the tables under `--input-root`; `result_root`,
  `--recompute-missing`, `--force-recompute`, and the recompute cache path are
  removed (the old `recomputed/*.npz` files remain on disk, unused). `source`
  column still reads `stored_script10`.
- Validation: all eight `cross_sample_nhood_enrichment` CSVs byte-identical to
  the phase 0 baseline (z-scores included). Script 08 annotation validated for
  CK_skin_res against Phoenix script 06 inputs (closes the phase 0/1 gap for
  one sample).
- After the next HPC script 10 run, run 10a with `--overwrite` to refresh the
  local tables (or copy the new CSVs).

Phase 1 as built:

- `config/cohorts/vbct/min15_area95_selected_resolutions.json` + README;
  loader `workflow_config.py` (repo root; renamed to `xenium_workflow/cohorts.py` 2026-11-22): `load_cohort`, `select_samples`,
  `add_root_arguments`, `resolve_path`, `REPO_ROOT`.
- Cohort sample keys: `sample`, `response`, `tissue`, `resolution`,
  `recluster_h5ad` (06), `annotated_h5ad` (08), `groupby`, `umap_key`;
  cohort-level `project`, `branch`, `branch_label`, `palette_json`,
  `mapping_csv`, `coord_keys`.
- Step configs: `"cohort"` plus optional `samples` (names or objects adding
  plot-only fields; redefining a cohort fact raises). Script 08 replot config
  uses `"input_key": "annotated_h5ad"`. Script 09 facets renamed `tissue` ->
  `facet`. Comparison config: `samples` is a name list; `script10_result_dir`
  replaces the chain to the script 10 config; `result_root` stays until phase 2.
- `--input-root`/`--output-root` on 08, 09, 10, spatial and poster helpers.
  Local script 08 now needs `--input-root ~/mnt/phoenix_scratch/Banksy_py`.
- Script 10 provenance key `annotation_config` -> `cohort_config`.
- Conference-figure helpers unchanged: they never chained to the 08 config;
  their config keeps its own sample list (absolute Phoenix paths) and palette
  copy. Candidate for phase 3 or later.
- Validation: identical sample records old vs new for all consumers; 09,
  infiltration, and poster dotplot reruns byte-identical to baseline; Slurm
  index -> sample mapping unchanged. Not run: script 08 annotate (Phoenix
  inputs), full script 10, full comparison helper (mount down).
- Before the next HPC run: pull this branch (or merged main); `run_08` and
  `run_10` wrappers work unchanged otherwise.

Phase 0 baseline: `data/xenium/baselines/post08_cleanup_2026-11-19/` (README,
`SHA256SUMS` for 36 CSVs, `--validate-only` logs). Script 08 could not be
validated locally because its script 06 inputs are only on Phoenix; validate it
when the mount is connected or on the HPC. `context/` is gitignored, so task
files are local only.

## Problems Being Fixed

Current flow after 08:

```text
06 clean recluster objects (Phoenix)
 └─ 08 annotate   [HPC batch or local; input_root defaults to the Phoenix mount, HPC overrides with --input-root]
     ├─ annotated H5ADs ──────────────┬─ 10 Squidpy   [HPC array, per sample; --annotation-root override]
     │                                │    ├─ _squidpy.h5ad → processed/  (Phoenix only)
     │                                │    └─ figures/CSVs → output/
     │                                ├─ infiltration helper  [local; reads 08 objects only]
     │                                ├─ poster helpers       [local; configs live in config/08_.../]
     │                                └─ conference figures   [local; config/conference_figures/]
     ├─ reviewed counts CSV ─── 09 proportions  [local]
     └─ 10 results + 08 objects ── nhood comparison helper  [local; needs mount via result_root; recompute fallback]
```

1. Config chains up to four hops (infiltration -> comparison -> script 10 ->
   script 08 -> palette/mapping), each resolving its root via its own
   `project_root_parent_depth` (3 or 4). Sample definitions live in 08 but
   response/tissue live only in the comparison config.
2. Three root mechanisms: `input_root`/`--input-root` (08),
   `--annotation-root` (10), `result_root` (comparison helper). The comparison
   helper reads local 08 objects and mounted 10 results at once.
3. Configs filed under the step they read from (poster configs in
   `config/08_cell_type_annotation/`, spatial helper configs in
   `config/10_squidpy_recluster_analysis/`; both since renamed in phase 3/4).
4. Figures spread across `figures/` (08, 09, poster) and `data/xenium/output/`
   (10 and spatial helpers since 2026-11-19).
5. Leftovers: "recluster" in script 10's name and config folder; a legacy
   cluster mode in script 10; a recompute fallback in the comparison helper
   that can produce the same result with different z-scores.

## Decisions (2026-11-19)

- **Figures for 08 and 09 move to `data/xenium/output/`**, matching 10 and the
  spatial helpers. Tables and H5ADs stay under `data/xenium/processed/`.
- **Archive script 10's legacy cluster mode.** It groups by numbered BANKSY
  recluster clusters from script 06 objects; its 8 per-sample configs point at
  the superseded `filtered_qc_v1_qc_pass_only` branch with a fixed resolution
  of 1.50. Archive `run_legacy_config`, `config/10_squidpy_recluster_analysis/vbct/*.json`,
  `config/10_squidpy_recluster_analysis/vbct/smoke/` (archived), and
  `hpc/slurm/run_10_squidpy_recluster_spatial_analysis.sl` (archived; done in
  phase 4, see `archive/2026-11-20_script10_legacy_cluster_mode/`). Existing outputs
  are not deleted. If pre-annotation cluster-level spatial analysis is ever
  needed (e.g. PTMT), make the grouping column configurable in the cell-type
  mode instead of restoring a second code path.
- **Rename script 10** to `10_cell_type_spatial_analysis.py` and its config
  folder to `config/10_cell_type_spatial_analysis/`. Nathalie pulls on the HPC
  before the next run.
- The cross-sample spatial helpers stay in `helper_scripts/spatial_analysis/`
  for now; `spatial_common.py` stays there too (not `banksy_utils/`).

## Plan

### Phase 0: Baseline

Record checksums and copies of key current outputs, which every later phase
must reproduce exactly (only paths may change):

- `nhood_enrichment_long.csv` and the other comparison tables
- `infiltration_distance_summary.csv`, `infiltration_distance_histogram.csv`
- script 08 `reviewed_cell_type_counts.csv`; script 09 plotted values
- `--validate-only` output for 08, 10, and both helpers

### Phase 1: One cohort config for post-08 steps

- Add `config/cohorts/vbct/min15_area95_selected_resolutions.json`: samples
  with response, tissue, branch, selected resolution, and annotated H5AD path;
  mapping CSV and palette paths; `data_root` with local and HPC values.
- 08, 09, 10, spatial helpers, and poster/conference helpers read it directly.
  Each step config keeps only its own parameters and output dirs, plus
  `"cohort": <path>`. One hop instead of four.
- Replace `input_root`, `--annotation-root`, and `result_root` with one
  `--data-root` override that behaves the same everywhere. Update
  `run_08_cell_type_annotation.sl` and `run_10_cell_type_spatial_analysis.sl`.

### Phase 2: Remove the Phoenix-mount dependency

- Script 10 also writes neighbourhood enrichment as small tables (counts,
  z-scores, cell-type order) beside its other CSVs.
- The comparison helper reads those tables instead of `_squidpy.h5ad` on the
  mount. Remove the recompute fallback.
- One-off backfill (like `00a`) extracts the tables from the existing stored
  H5ADs while the mount is connected, so script 10 need not be rerun.

### Phase 3: One output rule

- Figures: `data/xenium/output/<project>/<step>/<branch>/` for 08, 09, 10,
  helpers.
- Tables and H5ADs: `data/xenium/processed/...`.
- Configs live under the step that uses them (poster configs move out of
  `config/08_cell_type_annotation/`).

### Phase 4: Remove leftovers

- Rename script 10 and its config folder (see Decisions).
- Archive the legacy cluster mode, its configs, and the old wrapper.

### Phase 5: Documentation

- Post-08 diagram in `README_00_04_XENIUM_ANALYSIS_WORKFLOW_CONCISE.md`
  pipeline map; rewrite "Where Things Run".
- Update `spatial_analysis.md`, `cell_type_annotation.md`, `hca_poster.md`,
  `conference_figures.md`, and the config READMEs.

## Validation

- `--validate-only` for every step after each phase.
- Phase 0 baseline outputs identical after each phase.
- Phases 2 and 4: one small HPC script 10 run (`CK_skin_res`), run by
  Nathalie; its new enrichment tables must match the stored H5AD values.

## Guardrails

- Never SSH to or submit on the HPC; prepare commands for Nathalie.
- No scientific setting changes (graphs, seeds, permutations, thresholds).
- Do not modify or move files on the Phoenix mount.
- Commit per phase on the branch; merge to `main` only after validation.

## Session History

- **2026-11-19:** Plan drafted and decisions recorded. Script 10 and spatial
  helper figure output already moved to `data/xenium/output/` on `main`
  (commit `a6f3c4e`).
- **2026-11-20:** Phases 2–5 done; Moran's I seeded; legacy script 10 mode
  archived; branch history rewritten to drop AI co-author trailers and
  force-pushed; AGENTS.md now forbids AI attribution in commits/PRs. Mount
  found to be connected throughout (`mountpoint` misreport). Old `figures/`
  outputs left in place until reruns. Next: HPC one-sample test, then merge.
- **2026-11-20, HPC test (job 16757390, array task 0, CK_skin_res):** passed.
  No errors; figures only under `data/xenium/output/...` (21 files, same names
  as local); result H5AD, enrichment table, centrality and Moran CSVs under
  `processed/`. Enrichment table, centrality CSV, and FDR Moran CSV
  byte-identical to a local run (same package versions: scanpy 1.9.5, anndata
  0.10.2, numpy 1.24.4, squidpy 1.6.5). `moran_scores_all` differs only in
  simulated p-values: I and pval_norm agree to ~1e-14, but floating-point
  near-ties flip permutation counts (seeded RNG is identical; results also
  differ locally between 8 and 16 numba threads). With n_perms=10 the minimum
  pval_sim is 1/11; for CK_skin_res no gene passes FDR 0.1 (0 in both runs), but
  other samples flip to nearly all genes passing. Follow-up task:
  `morans_i_significance.md`. Local CK_skin_res
  enrichment table refreshed with 10a; identical to the HPC-written table.
  The stored CK_skin_res result on Phoenix is now the seeded one.
- **2026-11-21:** HPC rerun job 16758519 (tasks 1–7) finished without errors,
  so all eight samples have seeded script 10 results on Phoenix. Local
  enrichment tables refreshed with `10a --overwrite`; comparison and
  infiltration helpers rerun to their normal locations. Against the phase 0
  baseline, every table is byte-identical except the `zscore` column of
  `nhood_enrichment_long.csv` and `nhood_enrichment_focus_pairs.csv` (all 342
  primary-graph rows; max change 37.7, the unseeded-to-seeded shift). Log2 O/E,
  edge counts, differences, combined skin, tiles, and infiltration unchanged.
- **2026-11-21:** `AGENTS.md` Scientific Parameters rule committed (`00fb996`).
  Script 09 and script 08's clean-UMAP replot regenerated in the new output
  locations (09 CSVs byte-identical to baseline). Local `figures/spatial_analysis/`,
  `figures/cell_type_annotation/`, and `figures/cell_type_proportions/` moved to
  `figures/_pre_post08_cleanup_2026-11-20/` (161 MB, README there); docs updated
  (`9b00d44`). Script 08's titled UMAPs and combined panels exist only in the
  archive until script 08 is rerun. Poster folder and Phoenix copies untouched.
  Commits `00fb996` and `9b00d44` are local, not yet pushed.
