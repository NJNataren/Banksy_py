# Task: VBCT Cell-Type Annotation (Annotation Sets)

<!-- Task file rules (see Banksy_py/AGENTS.md, "Task files"): edit Current
State, Next Step and Open Questions in place; add one dated line per decision;
no "Session Update" sections; move superseded text to archive/ beside this file;
keep this file under ~200 lines. -->

## Goal

Give VBCT clusters shared, reviewed cell-type labels with one fixed palette, and
carry them through script 08 (annotation, UMAPs), 09 (proportions), 10 (spatial
analysis) and the spatial/poster helpers. Labels are revised by manual cluster
merging in numbered annotation sets, so every iteration is recorded and
reproducible.

**How-to:** `config/annotations/README.md` (layout, tiers, `new`/`check`/`diff`/
`release`, outputs); `config/08_cell_type_annotation/README.md` (script 08, HPC);
`config/10_cell_type_spatial_analysis/vbct/README.md` (script 10). Code:
`xenium_workflow/annotation_sets.py`. Marker dot plots for reviewing clusters:
`context/tasks/dotplot/dotplot.md`.

## Current State (updated 2026-11-24)

**Code.** Branch `feature/annotation-sets`, pushed to `origin`, not merged into
`main`. The Phoenix checkout is on this branch. Shared modules live in
`xenium_workflow/` (`cohorts.py`, `annotation_sets.py`, `spatial_tables.py`).

**Cohort.** `config/cohorts/vbct/min15_area95_selected_resolutions.json`:
min15_area95 QC branch, 127 clusters, 2,083,396 cells. Selected resolutions:
CK_skin_res 1.00, EL_skin_res 0.80, GR_lung_non_res_roi 1.00, BE_brain_non_res
0.90, HW_brain_res 0.80, MF_skin_non_res 0.50, MG_gastric_non_res 0.70,
CK_bowel_res 0.80.

**Annotation sets** (`config/annotations/vbct/min15_area95_selected_resolutions/`):

| Set | Status | Content | Script 08 (Phoenix) | Script 10 (Phoenix) |
|---|---|---|---|---|
| v1 | released 2026-11-22 (`f60c669e…`) | first reviewed labelling, tier 1 only, 8 types | done, job 16939096; counts identical to the pre-set labelling | done 8/8, job 16941700, tier 1; logs show no errors |
| v2 | released 2026-11-24 (`9d6d99cc…`) | parent v1; BE c5 (6,652 cells) and HW c15 (431 cells) moved from Epithelial cell to Melanoma cell | done, job 17027103 (run while v2 was a draft, same hash); `changes_vs_v1.csv` lists exactly those 2 clusters | done 8/8, tier 1 (verified through the mount with `tools/pipeline_status.py` 2026-11-25) |

Script 08 status was verified through the mount on 2026-11-24 with
`tools/pipeline_status.py --annotation <set> --root <mount>`.

**Script 09 (local, reads counts from the mount):** v1 done 2026-11-24, both
configs (standard and by tissue). Proportions are identical to the
pre-annotation-set version. Outputs:
`data/xenium/{output,processed}/vbct/cell_type_proportions/min15_area95_selected_resolutions/v1/tier_1/`.
v2 done 2026-11-24 (same folder under `v2/`). Only BE and HW change from v1:
BE melanoma 24.4% → 33.8% and epithelial 9.4% → 0%; HW melanoma 71.5% → 71.6%
and epithelial 0.1% → 0%. Epithelial cells remain only in samples whose
epithelial clusters v2 did not change.

**Palette (tier 1, unchanged since v1):** Melanoma cell `#AA4499`, Macrophage
`#EE7733`, T cell `#4477AA`, B cell `#66CCEE`, Fibroblast `#228833`,
Endothelial cell `#CC3311`, Epithelial cell `#CCBB44`, Pericyte `#777777`.

**Outputs** (paths derived from set and tier; see `config/annotations/README.md`):
- Phoenix: script 08 objects, tables and UMAPs under
  `data/xenium/{processed,output}/vbct/cell_type_annotation/min15_area95_selected_resolutions/<version>/`.
  Script 10 results under `.../spatial_analysis/min15_area95_selected_resolutions/v1/tier_1/<sample>/`.
- Local: cluster marker dot plots for v1 (8 samples) and v2 (BE, HW) under
  `data/xenium/output/vbct/cell_type_annotation/min15_area95_selected_resolutions/<version>/tier_1/cluster_marker_dotplots/`.
  Cross-sample spatial helpers: run 2026-11-25 for v1 (neighbourhood comparison)
  and v2 (both helpers; the helper configs now name v2). Results and v1 → v2
  changes: `context/tasks/spatial_analysis/spatial_analysis.md`.
- Archived: outputs from before annotation sets (old `cell_type` column, no
  stamp) are in `data/xenium/_pre_annotation_sets_2026-11-22/`, locally and on
  Phoenix (`helper_scripts/archive_pre_annotation_sets_outputs.sh`, README
  inside). The conference-figure config reads the archived objects there.

**Caveat for this week's meeting:** script 10 ran with the current Moran's I
settings. Do not present its significant-gene CSVs or top-Moran plots
(`context/tasks/morans_i_significance.md`). Neighbourhood enrichment,
centrality, co-occurrence and the cross-sample comparisons are unaffected.

## Decisions

- 2026-10-24: Never transfer labels between samples, branches or resolutions by
  matching numeric cluster IDs; label each clustering from its own evidence.
- 2026-10-24 (Nathalie): Annotate all eight VBCT samples at the resolutions
  listed above.
- 2026-10-25 (Nathalie): Use the min15_area95 branch and the eight-type tier-1
  vocabulary and colours above; label spellings normalised (e.g. `Melanoma` →
  `Melanoma cell`).
- 2026-10-25 (Nathalie): Scanpy UMAPs; PNG, SVG and PDF at 300 DPI with
  editable text and rasterised points.
- 2026-11-22 (Nathalie): Labels live in versioned annotation sets (plain `vN`,
  each with a `CHANGES.md`); the cohort keeps clustering facts only. Script 08
  stays the annotation step.
- 2026-11-22 (Nathalie): Tiers are `cell_type_tier_1` (broad) and optional
  `cell_type_tier_2` (finer, nested under tier 1). The old `cell_type` column is
  retired. A blank tier 2 inherits its tier-1 label and is flagged
  (`cell_type_tier_2_inherited`).
- 2026-11-22 (Nathalie): Released sets are locked by hash; any change means a
  new version. Scripts 09 and 10 and the helpers need a released set
  (`--allow-draft` for local tests only).
- 2026-11-22 (Nathalie): Regenerate v1 under the new layout rather than move old
  outputs; archive the unversioned outputs. Results so far are preliminary.
- 2026-11-22 (Nathalie): Shared modules moved into the `xenium_workflow/`
  package, not `banksy_utils/` (keeps upstream BANKSY separate).
- 2026-11-22 (Nathalie): Run script 10 for v1 and v2 before the Moran's I fix,
  with unchanged settings, for this week's meeting.
- 2026-11-23 (Nathalie): v2 moves BE c5 and HW c15 from Epithelial cell to
  Melanoma cell (melanocytic markers with KRT18, reviewed on the v1 cluster dot
  plots).

## Next Step

Decide when to merge `feature/annotation-sets` into `main` (v2 script 10 and
the spatial helpers are done as of 2026-11-25; see Open Questions).

## Open Questions

- **Merging the branch:** when should `feature/annotation-sets` be merged into
  `main`?
- **Moran's I:** script 10 results need a rerun after the fix
  (`context/tasks/morans_i_significance.md`).
- **Tier 2:** not used yet. A tier-2 script 10 run would repeat Moran's I,
  which does not depend on labels; decide whether to skip it then.
- **Small clusters seen in the v1 dot plots:** CK_bowel_res c16, c17 and c18,
  and GR_lung_non_res_roi c13. Candidates for a later version; details in
  `context/tasks/dotplot/dotplot.md`.
- **Label consistency across samples:** a cross-sample label dot plot is an open
  dot plot question (`context/tasks/dotplot/dotplot.md`).
- **v2 dot plots** for the other six samples, if a full v2 set is wanted.

## History

- To 2026-11-24: the initial annotation plan, the selected-resolution review,
  script 08/09/10 build-out, the annotation-set design and its implementation
  and verification, v1 regeneration, archiving, and the v2 draft and dot plots.
  All in `context/tasks/cell_type_annotation/archive/cell_type_annotation_to_2026-11-24.md`.
