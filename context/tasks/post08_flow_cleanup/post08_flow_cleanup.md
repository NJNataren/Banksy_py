# Task: Untangle the Pipeline Flow After Script 08

<!-- Task file rules (see Banksy_py/AGENTS.md, "Task files"): edit Current
State, Next Step and Open Questions in place; add one dated line per decision;
no "Session Update" sections; move superseded text to archive/ beside this file;
keep this file under ~200 lines. -->

## Goal

Keep the steps after script 08 (08, 09, 10, 10a, the spatial helpers, the
poster helpers) simple to follow: one sample definition, one rule for data
roots, one rule for output locations, and a step order that matches how data
flows. Change no scientific result.

The first cleanup (2026-11-19 to 2026-11-21) is done and merged. The remaining
work is the restructure of the spatial steps described under Next Step. It is
not a small job and is not scheduled yet.

Related: `context/tasks/spatial_analysis.md` (the analysis itself),
`context/tasks/config_refactor.md` and `context/tasks/decouple_workflow_from_banksy.md`
(later config and package work; this task must stay compatible with them).

## Current State (updated 2026-11-24)

**Done and on `main`** (branch `cleanup/post08-flow`, merged after the HPC
test; details in the archive):
- One cohort config, `config/cohorts/vbct/min15_area95_selected_resolutions.json`.
  Since 2026-11-22 steps name a versioned annotation set in
  `config/annotations/` that points at it (`config/annotations/README.md`).
- Shared code in `xenium_workflow/`: `cohorts.py` (cohort and
  `--input-root`/`--output-root`), `annotation_sets.py`, `spatial_tables.py`.
- Script 10 writes a small enrichment CSV per sample;
  `10a_export_nhood_enrichment_tables.py` backfills it from result H5ADs. The
  comparison helper reads only those CSVs for enrichment.
- Output rule: figures under `data/xenium/output/`, tables and H5ADs under
  `data/xenium/processed/`, by project / step / cohort / annotation version /
  label level (`xenium_workflow/annotation_sets.py`, `step_dir`).
- Script 10 renamed to `10_cell_type_spatial_analysis.py`; legacy cluster mode
  archived (`archive/2026-11-20_script10_legacy_cluster_mode/`). Moran's I
  seeded (1234); its significance output is still unusable
  (`context/tasks/morans_i_significance.md`).
- Validation: phase 0 baseline in `data/xenium/baselines/post08_cleanup_2026-11-19/`.
  After the seeded HPC rerun (jobs 16757390, 16758519) every table matched it
  except the Squidpy `zscore` column (seed change); log2 O/E unchanged.

**How the spatial steps are laid out now** (the reason for Next Step):

```text
08 annotate ─┬─ 09 composition                      [local]
             ├─ 10 per-sample Squidpy (HPC) ─ 10a backfill tables
             │     └─ helper: compare_nhood_enrichment_from_config.py   [local]
             │          reads 10's tables, but also rebuilds per-sample graphs
             │          from coordinates (knn15, radius 25/50/100 µm, 500 µm tiles)
             └─ helper: plot_immune_infiltration_from_config.py         [local]
                  computes per-sample immune-to-tumour distances and plots
                  across samples in one script
```

- Both helpers live in `helper_scripts/spatial_analysis/` with
  `spatial_common.py`, although they produce the main spatial results.
- Each helper mixes per-sample computation with the cross-sample comparison,
  so local runs still need the annotated H5ADs (coordinates and labels).
- The two helpers no longer share a config: each names the annotation set.
  The infiltration helper's docstring (lines 15–17) still says its samples come
  from the comparison config; that text is stale.
- The planned tumour core vs border mask and TLS check
  (`context/tasks/spatial_analysis.md`) have no place in the current layout.

**Loose ends from the first cleanup:**
- Conference-figure helpers still keep their own sample list and palette copy
  and write to `figures/conference_min15_area95_selected_resolutions_poster/`
  (`context/tasks/conference_figures.md`).
- That poster folder is not archived: confirm no design file links into it.
- Script 08's titled UMAPs and combined panels exist only in
  `figures/_pre_post08_cleanup_2026-11-20/` until script 08 is rerun.

## Decisions

- 2026-11-19 (Nathalie): one cohort config for every step after 08; figures for
  08 and 09 move to `data/xenium/output/`; archive script 10's legacy cluster
  mode; rename script 10.
- 2026-11-19 (Nathalie): keep the cross-sample spatial helpers as helpers until
  the analysis extensions settle; `spatial_common.py` not into `banksy_utils/`.
- 2026-11-20 (Nathalie): merge to `main` only after a one-sample HPC test; old
  `figures/` outputs stay until reruns refill the new locations (archived
  2026-11-21).
- 2026-11-24 (Nathalie): the spatial restructure below is future work, not
  started; numbering the two helpers by dependency alone (e.g. 08a / 10b) was
  discussed but not chosen.

## Next Step

Restructure the spatial steps so each number is one stage and each analysis
is split into per-sample computation and cross-sample comparison. Not
scheduled; do it on its own branch with a baseline, as for the first cleanup.

Proposed layout (2026-11-24, for discussion):

| Step | Does | Reads | Runs |
|---|---|---|---|
| 08 | annotate cell types | recluster H5AD | HPC |
| 09 | cell-type composition | 08 | local |
| 10 | per-sample neighbourhood: enrichment for all graph settings and tiles, centrality, Moran's I | 08 | HPC |
| 10a | backfill enrichment tables (exists) | 10 | local |
| 11 | per-sample regions: tumour mask, core vs border (later TLS) | 08 | HPC or local |
| 12 | per-sample infiltration distances | 08 (+ 11 once the mask exists) | HPC or local |
| 13 | cross-sample comparison and plots (13a neighbourhood, 13b infiltration) | tables from 10 and 12 | local |

Rules behind it:
- Per-sample steps write small CSVs; step 13 reads only those, so local work
  never needs the H5ADs or the Phoenix mount (the pattern already used for
  enrichment tables).
- Step order is run order; a new per-sample analysis takes the next number
  before the cross-sample step, a new cross-sample view becomes 13c.
- Letters keep their current meaning: small sub-steps (10a) or parallel views
  of one stage (13a/b).
- Each step has its own config naming the annotation set; no step reads its
  sample list from another step's config.
- `spatial_common.py` moves into `xenium_workflow/`; poster helpers stay
  helpers.

Work involved: split each helper into compute and compare, define table
formats in `xenium_workflow/` (like `spatial_tables.py`), move the sensitivity
graphs and tiles into script 10 (longer HPC job, same settings), new configs
and Slurm wrapper changes, docs (`README_00_04_XENIUM_ANALYSIS_WORKFLOW_CONCISE.md`,
config READMEs, `AGENTS.md`), and a byte-for-byte check of every table against a
baseline taken first. No scientific parameter changes; if any setting has to
move, confirm with Nathalie first.

## Open Questions

- Will the tumour mask be its own step used by several analyses, or live inside
  the infiltration step? Decides whether slot 11 exists; settle before
  restructuring.
- Full compute/compare split, or rename the two helpers first and split later?
- Move Moran's I out of script 10 into its own step, to keep it apart from the
  enrichment results until `context/tasks/morans_i_significance.md` is fixed?
- Accept the longer script 10 HPC job for sensitivity graphs and tiles, or
  keep them in 13a and let 13a read coordinates?

## History

- To 2026-11-24: problems found, decisions, phases 0–5 as built, HPC test and
  seeded reruns, figure archiving, commit list. All in
  `context/tasks/post08_flow_cleanup/archive/post08_flow_cleanup_to_2026-11-24.md`.
