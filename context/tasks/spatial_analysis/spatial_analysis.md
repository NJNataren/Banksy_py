# Task: Cross-Sample Spatial Analysis of Response

<!-- Task file rules (see Banksy_py/AGENTS.md, "Task files"): edit Current
State, Next Step and Open Questions in place; add one dated line per decision;
no "Session Update" sections; move superseded text to archive/ beside this file;
keep this file under ~200 lines. -->

## Goal

Identify spatial features of the tumour microenvironment that differ between
ICI responders (`res`) and non-responders (`non_res`), from the script 10
cell-type spatial outputs for the VBCT `min15_area95` selected-resolution
objects. Current focus: whether tumour cells exclude or admit immune cells, and
how tumour, immune, and endothelial cells are arranged relative to each other.

**Tools** (local, `banksy` env, from the `Banksy_py` root):
`helper_scripts/spatial_analysis/compare_nhood_enrichment_from_config.py`
(neighbourhood log2 O/E, sensitivity graphs, tiles) and
`helper_scripts/spatial_analysis/plot_immune_infiltration_from_config.py`
(distance to nearest tumour cell), sharing `spatial_common.py`. Usage, config
keys and output paths: `config/spatial_analysis/README.md`; pipeline position
and commands: `README_00_04_XENIUM_ANALYSIS_WORKFLOW_CONCISE.md`, "After
Script 08"; script 10 run history: `config/10_cell_type_spatial_analysis/vbct/README.md`.

Related: `context/tasks/cell_type_annotation/cell_type_annotation.md`
(annotation sets, palette), `context/tasks/post08_flow_cleanup/post08_flow_cleanup.md`
(the 2026-11-19/20 restructure), `context/tasks/hca_poster.md` (these results
were kept off the poster), `context/tasks/morans_i_significance.md`.

## Current State (updated 2026-11-25)

**Samples** (response, tissue and resolution come from
`config/cohorts/vbct/min15_area95_selected_resolutions.json`; never parse them
from names, as `non_res` contains `res`):

| Sample | Response | Tissue | Resolution |
|---|---|---|---|
| `CK_skin_res` | res | skin | 1.00 |
| `EL_skin_res` | res | skin | 0.80 |
| `CK_bowel_res` | res | bowel | 0.80 |
| `HW_brain_res` | res | brain | 0.80 |
| `MF_skin_non_res` | non_res | skin | 0.50 |
| `MG_gastric_non_res` | non_res | gastric | 0.70 |
| `BE_brain_non_res` | non_res | brain | 0.90 |
| `GR_lung_non_res_roi` | non_res | lung | 1.00 |

**Inputs.** Both helper configs name annotation set **v2** (`tier_1`) since
2026-11-25. v2 script 08 and script 10 outputs are complete for all eight
samples on Phoenix (`tools/pipeline_status.py --annotation .../v2 --root <mount>`,
2026-11-25). The enrichment tables and annotated H5ADs exist only on the mount,
so run the helpers with `--input-root ~/mnt/phoenix_scratch/Banksy_py`
(read-only); outputs are written locally.

**Outputs (2026-11-25), both versions local:**
`data/xenium/{output,processed}/vbct/spatial_analysis/min15_area95_selected_resolutions/<v1|v2>/tier_1/`
with `cross_sample_nhood_enrichment/` and (v2 only) `immune_infiltration/`.
v1 figures are kept for reference; v2 is current.

**Kanban #65 (triangle heatmaps), done 2026-11-25.**
`compare_nhood_enrichment_from_config.py` has a config key `heatmap_layout`
(`"triangle"` or `"full"`, default `"full"`); the VBCT config sets
`"triangle"`. Every heatmap (all-sample grid, per-comparison, difference,
combined skin) shows the lower triangle with the diagonal, labelled on both
axes. The all-sample grid places responders (left, blue heading) and
non-responders (right, orange heading) side by side as 2×2 blocks with one
shared colour scale. A cell type is kept in a figure if any of its samples has
it (blank row/column where a sample lacks it); it is omitted only when no
sample in that figure has it (in v2, Epithelial cell is absent from both brain
samples, so the brain comparison shows six types).

Checks on the v1 run: difference, combined-skin and tile tables byte-identical
to `data/xenium/baselines/post08_cleanup_2026-11-19/`; long and focus-pair
tables differ only in reference z-scores (max 37.7, seeded reruns).

**v1 → v2 changes** (v2 relabels BE c5 and HW c15 from Epithelial to Melanoma
cell). Only BE_brain_non_res and HW_brain_res change in any table; skin and GI
results are identical.
- HW: epithelial row gone; all other values unchanged (431 cells).
- BE (knn6 log2 O/E, v1 → v2): melanoma–melanoma 1.49 → 1.26, melanoma–T
  −0.95 → −1.41, melanoma–macrophage −1.26 → −1.72, melanoma–B −0.17 → 0.46,
  melanoma–fibroblast −3.09 → −3.53.
- Brain difference (res − non_res): melanoma–T 0.43 → 0.89, melanoma–B
  −0.38 → −1.02. All-sample difference: melanoma–T 1.69 → 1.81.
- BE infiltration: contact ratio T 1.23 → 1.33, B 1.53 → 2.53, macrophage
  1.06 → 1.14; median B-cell distance 19.2 → 9.7 µm. HW unchanged.

**Findings so far** (v2; descriptive only: 4 vs 4 samples, five tissue sites,
different resolutions):
- Melanoma–fibroblast segregation and strong same-type clustering appear in
  every sample and do not separate the response groups.
- Melanoma–T exclusion is strong in `MG_gastric_non_res` (log2 O/E −6.9),
  `MF_skin_non_res` (−5.4) and `EL_skin_res` (−5.1); the others range from
  −0.5 to −1.4. It persists at a 100 µm radius (−3 to −5) and in 500 µm tile
  medians (MG −5.7, EL −4.0, MF −3.8), so exclusion does not separate `res`
  from `non_res`. The earlier reading of raw z-score heatmaps was a scale
  artefact.
- The GI pair is the clearest contrast: T-cell contact ratio 1.36 in
  `CK_bowel_res` vs 0.10 in `MG_gastric_non_res` (median distance to tumour
  12 vs 142 µm); B-cell contact ratios 0.48 vs 0.04.
- T-cell contact ratios: CK_bowel 1.36, BE 1.33, HW 1.10, GR 0.96, CK_skin
  0.74, MF 0.63, EL 0.18, MG 0.10.
- `EL_skin_res` is excluded but has a dense T-cell rim and a large nearby
  T-cell aggregate; `MG_gastric_non_res` has a thinner rim with most T cells
  distant. In excluded samples, macrophages are the immune cells in contact
  with tumour (contact ratio about 2.4–2.6 in EL and MG).
- Skin, responder mean (CK_skin, EL_skin) vs MF_skin_non_res: almost the whole
  melanoma row is red in the difference (melanoma–T +2.2, –macrophage +2.0,
  –fibroblast +2.3), so MF's tumour is separated from every cell type.
  Macrophage–T is −2.6. The responder mean hides a split (melanoma–T −1.3 in
  CK_skin, −5.1 in EL_skin); show the side-by-side figure wherever that pair
  matters.
- Endothelial focus pairs (2026-10-28) are exploratory; nothing recorded yet.

**How to read the figures.** Red (>0) = neighbours more often than chance,
blue = less; each step of 1 is a doubling. In difference panels (res − non_res)
red means more together, or less apart, in responders; both groups can be
below 0. Values describe arrangement relative to each sample's own cell mix,
not cell-type amounts (script 09 has those).

## Decisions

- 2026-10-26: compare samples with log2(observed/expected) neighbour edges,
  not Squidpy z-scores (which scale with cell number); directed kNN counts
  symmetrised, expected from edge shares, pseudocount 1.
- 2026-10-26: primary graph is script 10's `knn6`; other graphs use KD-tree
  counts (match Squidpy edge-for-edge) and have no z-scores.
- 2026-10-26: group means weight each sample equally; pooled edge counts
  rejected (EL_skin ~236,000 cells vs CK_skin ~6,000). Differences are
  descriptive only.
- 2026-10-26: tiles are 500 µm squares; a pair is reported only where both
  types have ≥20 cells.
- 2026-10-27: infiltration bands ≤20 µm (contact), 20–100 µm, >100 µm;
  contact ratio is relative to all non-tumour cells.
- 2026-10-28: focus pairs changed to the endothelial set (old outputs in
  `previous_focus_pairs_2026-10-26/`).
- 2026-11-15: combined figures show all types; values from fewer samples than
  the group are marked †.
- 2026-11-19 (Nathalie, with Kat): responder and non-responder O/E heatmaps go
  side by side as triangles with labels on both axes; do not flip and merge
  the triangles into one plot (kanban #65).
- 2026-11-25 (Nathalie): spatial helpers use annotation set v2; triangle
  layout for all heatmap figures; keep blank rows/columns for cell types a
  sample lacks, so samples stay comparable.
- Figure colours: res `#2a78d6`, non_res `#eb6834`; tissue by marker shape;
  heatmaps blue–grey–red centred at 0.

## Next Step

Review the endothelial focus-pair outputs (v2:
`nhood_enrichment_focus_pairs`, `nhood_enrichment_graph_sensitivity`,
`nhood_enrichment_tiles_knn6_500um`) and record findings above.

## Open Questions

1. Should a cell type absent from every sample in a figure (Epithelial cell in
   the v2 brain comparison) still get a blank row, so all figures share the
   eight-type layout?
2. Distance to nearest tumour cell cannot tell intratumoural from border
   cells; add an inside-tumour vs border mask (e.g. local tumour density).
3. Does the `EL_skin_res` T-cell aggregate have tertiary lymphoid structure
   features?
4. Will PTMT use the same helpers? Needs only a PTMT cohort config and helper
   configs.
5. Formal res vs non_res testing is not meaningful at 4 vs 4; per-sample
   robustness (graphs, tiles) is the available evidence.
6. Moran's I significance in script 10 is unusable as is; do not use the Moran
   FDR CSVs or plots (`context/tasks/morans_i_significance.md`).
7. Promoting the helpers to numbered steps is part of a proposed spatial-step
   restructure (`context/tasks/post08_flow_cleanup/post08_flow_cleanup.md`).

## Pitfalls

- Never compare Squidpy z-score heatmaps across samples; use log2 O/E.
- Annotated objects list all palette types as categories, including absent
  ones; drop absent types in palette order, as script 10 does.
- After any HPC script 10 run, refresh tables with script 10a (`--overwrite`).
- Project-wide gotchas (mount check, `conda run` stdin, fonts) are under
  "Environment" in `AGENTS.md`.

## History

- 2026-10-26 to 2026-11-21: helpers built, post-08 cleanup, seeded HPC reruns,
  script-by-script description, configs and session history:
  `archive/spatial_analysis_to_2026-11-25.md`.
- 2026-11-25: file converted to the task-file layout and moved into
  `context/tasks/spatial_analysis/`; #65 triangle layout implemented; helpers
  switched to v2 and rerun (v1 outputs kept).
