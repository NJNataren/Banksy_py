# Task: Cross-Sample Spatial Analysis of Response

## Goal

Identify spatial features of the tumour microenvironment that differ between
ICI responders (`res`) and non-responders (`non_res`), starting from the script
10 cell-type spatial outputs for the VBCT `min15_area95` selected-resolution
objects. Current focus: whether tumour cells exclude or admit immune cells, and
how tumour, immune, and endothelial cells are arranged relative to each other.

Related context:

- `README_00_04_XENIUM_ANALYSIS_WORKFLOW_CONCISE.md`, section "After Script 08":
  post-08 flow diagram, run table, and commands.
- `config/spatial_analysis/README.md`: helper usage; `config/10_cell_type_spatial_analysis/vbct/README.md`:
  script 10 HPC run history.
- `context/tasks/post08_flow_cleanup.md`: the 2026-11-19/20 restructuring
  behind the paths below (branch `cleanup/post08-flow`).
- `context/tasks/cell_type_annotation/cell_type_annotation.md`: script 08 annotations and shared palette.
- `context/tasks/hca_poster.md`: these results were kept off the HCA poster.

## Current Status (2026-11-21)

All analyses below are built and have been run on all eight samples. This file
describes `main`, which includes the post-08 cleanup (merged 2026-11-20 after
the HPC test passed). Before that merge, script 10 was
`10_squidpy_recluster_spatial_analysis.py` (now archived) and the helpers read
the result H5ADs on the mount.

- All eight samples were rerun on the HPC with seeds on 2026-11-20/21 (jobs
  16757390 and 16758519) from merged `main`; these replace the unseeded
  October results on Phoenix. Local enrichment tables were refreshed with 10a
  and both helpers rerun on 2026-11-21: only the reference z-scores changed
  (max 37.7); log2 O/E, all difference/combined/tile tables, infiltration, and
  the findings below are unchanged.
- Script 10 seeds neighbourhood enrichment and Moran's I (1234) and writes a
  small enrichment table per sample; the comparison helper reads these local
  tables and no longer needs the mount.
- Before the reruns, a comparison from tables backfilled from the old unseeded
  results reproduced every baseline table byte for byte (z-scores included).
- Current outputs are in the new locations
  (`data/xenium/output/...` and `data/xenium/processed/...`); pre-cleanup
  figures were archived locally in `figures/_pre_post08_cleanup_2026-11-20/spatial_analysis/` (2026-11-21).

## Samples

| Sample | Response | Tissue | Resolution | Script 10 result |
|---|---|---|---|---|
| `CK_skin_res` | res | skin | 1.00 | stored |
| `EL_skin_res` | res | skin | 0.80 | stored |
| `CK_bowel_res` | res | bowel | 0.80 | stored (HPC rerun) |
| `HW_brain_res` | res | brain | 0.80 | stored |
| `MF_skin_non_res` | non_res | skin | 0.50 | stored (HPC rerun) |
| `MG_gastric_non_res` | non_res | gastric | 0.70 | stored |
| `BE_brain_non_res` | non_res | brain | 0.90 | stored |
| `GR_lung_non_res_roi` | non_res | lung | 1.00 | stored |

Response and tissue labels are set explicitly in the cohort config
(`config/cohorts/vbct/min15_area95_selected_resolutions.json`); do not parse
them from sample names (`non_res` contains `res`). Sample order in figures is
set by each helper config's `samples` list (responders first).

## Scripts

Run from the `Banksy_py` root with the `banksy` conda environment.

- `10_cell_type_spatial_analysis.py` (repo root, HPC): per-sample Squidpy
  analyses. Writes the result H5AD, `<stem>_nhood_enrichment.csv`, and
  centrality/Moran CSVs to `processed/`, figures to `output/`.
- `10a_export_nhood_enrichment_tables.py` (repo root, local): rebuilds the
  enrichment tables from stored result H5ADs (reads only the enrichment and
  `cell_type`; verifies the round trip). Run after every HPC script 10 run.
- `helper_scripts/spatial_analysis/` (local). Both analysis scripts accept
  `--config`, `--sample` (repeatable), `--validate-only`, `--input-root`, and
  `--output-root` (both roots default to the repository root).

### `spatial_common.py` (shared module)

- `load_sample_records`: reads the cohort via `workflow_config.py` (renamed to `xenium_workflow/cohorts.py` 2026-11-22; now reads an annotation set) and returns
  one record per sample (response, tissue, branch, resolution, annotated H5AD,
  enrichment table path), in the helper config's order.
- `read_labels_and_coords`: reads only `obs/cell_type`, `obs/x`, and `obs/y`
  from an H5AD and drops absent cell types in palette order, as script 10 does.
- Shared response colours, tissue markers, panel titles, label spacing
  (`spread_label_positions`), and figure export (`save_figure`: png/svg/pdf
  with editable text).

### `compare_nhood_enrichment_from_config.py` (neighbourhood enrichment)

- Loads each sample's script 10 enrichment from its table (`load_enrichment`,
  via `spatial_tables.read_enrichment_table`). A missing table is an error that
  points to script 10a. The local recompute fallback was removed on 2026-11-20.
- Converts edge counts to log2 observed/expected (`log2_obs_exp`,
  `enrichment_to_long`).
- Neighbourhood-size sensitivity: repeats the edge counts for each configured
  graph with KD-tree counting (`neighbour_type_counts`).
- Tile analysis: focus-pair log2 O/E inside 500 µm squares
  (`tile_focus_pairs`).
- Figures:
  - `plot_sample_grid`: all eight samples on one shared colour scale,
    responders on the top row.
  - `plot_comparison`: per-sample panels for one comparison plus a
    res − non_res difference panel.
  - `plot_combined_comparison`: for comparisons with `combine_by_response`, one
    mean panel per response group plus their difference (`group_mean_matrix`).
    Each sample counts equally. A cell type absent from some samples in a group
    is averaged over the samples that have it and marked with †, explained in a
    footnote.
  - `plot_focus_pairs`, `plot_graph_sensitivity`, `plot_tile_distributions`:
    focus-pair strip plot, values across neighbourhood sizes, and tile
    distributions.

### `plot_immune_infiltration_from_config.py` (distance to tumour)

- For each sample, measures the distance from every T cell, B cell, and
  macrophage to its nearest melanoma cell (`tumour_distances`), plus all
  non-tumour cells as a geometry reference.
- `summarise`: median and quartile distances, the share of cells in each
  distance band, and the contact ratio.
- Figures: `plot_bin_fractions` (stacked bars of distance bands),
  `plot_distance_ecdf` (cumulative distance curves; samples told apart by
  colour plus line style), and `plot_spatial_map` (T cells coloured by distance
  band over tumour). `--skip-maps` omits the maps.

## Configs

- Cohort: `config/cohorts/vbct/min15_area95_selected_resolutions.json`
  (samples, response, tissue, resolution, object paths, palette, mapping,
  coordinate keys).
- Script 10: `config/10_cell_type_spatial_analysis/vbct/min15_area95_selected_resolutions.json`
  (`figure_dir`, `result_dir`, Moran settings).
- `config/spatial_analysis/vbct/nhood_enrichment_comparison.json`:
  - `samples` (order), `script10_result_dir` (where enrichment tables live,
    under `--input-root`)
  - pseudocount, `primary_graph` (`knn6`) and `graph_settings` (`knn6`,
    `knn15`, `radius25um`, `radius50um`, `radius100um`)
  - `tiling`: `knn6`, 500 µm tiles, at least 20 cells of each type
  - `focus_pairs` (currently Endothelial–Endothelial, Endothelial–Melanoma,
    Endothelial–T, Endothelial–Macrophage, and Melanoma–T; see Session History)
  - `comparisons`: `all_samples`, `gastrointestinal` (CK_bowel_res vs
    MG_gastric_non_res), `skin` (CK_skin_res and EL_skin_res vs
    MF_skin_non_res, with `"combine_by_response": true`), and `brain`
    (HW_brain_res vs BE_brain_non_res)
- `config/spatial_analysis/vbct/immune_infiltration.json`: `samples` (order),
  tumour cell type, query cell types (T, B, macrophage), distance bands (20 and
  100 µm), plot range, and map settings.

Each helper config points straight at the cohort; there is no config chain.

## Running

```bash
# Only after a new HPC script 10 run (reads the mount, writes locally):
python 10a_export_nhood_enrichment_tables.py \
  --config config/10_cell_type_spatial_analysis/vbct/min15_area95_selected_resolutions.json \
  --input-root ~/mnt/phoenix_scratch/Banksy_py --overwrite

S=config/spatial_analysis/vbct
python helper_scripts/spatial_analysis/compare_nhood_enrichment_from_config.py --config $S/nhood_enrichment_comparison.json
python helper_scripts/spatial_analysis/plot_immune_infiltration_from_config.py --config $S/immune_infiltration.json
```

Each helper takes about 30 seconds and needs only local files: script 08
annotated objects and the enrichment tables.

## Inputs and Outputs

Inputs (local):

- Script 08 annotated objects:
  `data/xenium/processed/vbct/<sample>/cell_type_annotation/adata_expression_clean_<sample>_min15_area95_cell_type_annotated_r<res>.h5ad`
- Enrichment tables:
  `data/xenium/processed/vbct/spatial_analysis/min15_area95_selected_resolutions/<sample>/<stem>_nhood_enrichment.csv`
  (backfilled 2026-11-20 from the Phoenix result H5ADs)
- Shared palette: `config/annotations/vbct/min15_area95_selected_resolutions/v1/palette.json` (moved there from the script 08 config folder
  on 2026-11-22; named by the annotation set)

Outputs (figures as png/svg/pdf):

```text
data/xenium/output/vbct/spatial_analysis/min15_area95_selected_resolutions/cross_sample_nhood_enrichment/
  nhood_enrichment_log2_obs_exp_all_samples
  nhood_enrichment_comparison_{all_samples,gastrointestinal,skin,brain}
  nhood_enrichment_comparison_skin_combined
  nhood_enrichment_focus_pairs, nhood_enrichment_graph_sensitivity, nhood_enrichment_tiles_knn6_500um
data/xenium/output/vbct/spatial_analysis/min15_area95_selected_resolutions/immune_infiltration/
  infiltration_distance_bin_fractions, infiltration_distance_ecdf
  maps/infiltration_map_<sample>_T_cell
data/xenium/processed/vbct/spatial_analysis/min15_area95_selected_resolutions/cross_sample_nhood_enrichment/
  nhood_enrichment_long.csv               all samples, all pairs, all graphs
  nhood_enrichment_focus_pairs.csv
  nhood_enrichment_difference_<comparison>.csv
  nhood_enrichment_combined_skin.csv      group means and samples per value
  nhood_enrichment_tiles_knn6_500um.csv   focus pairs only
  previous_focus_pairs_2026-10-26/        tables from the original focus pairs
  recomputed/                             old recompute cache, no longer read
data/xenium/processed/vbct/spatial_analysis/min15_area95_selected_resolutions/immune_infiltration/
  infiltration_distance_summary.csv, infiltration_distance_histogram.csv
```

Figures made before 2026-11-19 (including `previous_focus_pairs_2026-10-26/`)
were archived locally on 2026-11-21 in
`figures/_pre_post08_cleanup_2026-11-20/spatial_analysis/vbct_min15_area95_selected_resolutions/`. Phase 0 reference
tables: `data/xenium/baselines/post08_cleanup_2026-11-19/`.

## Decisions Made

- **Metric:** Squidpy z-scores grow with cell number (CK_skin −29 to 70 versus
  MF_skin −750 to 1179), so samples are compared with log2(observed / expected)
  neighbour edges. Directed kNN counts are symmetrised, and expected edges come
  from each type's share of edges. Pseudocount 1. Z-scores remain in
  `nhood_enrichment_long.csv` for reference only.
- **Primary graph:** script 10's 6-nearest-neighbour graph (`knn6`). Heatmaps
  and comparison figures use it. Coordinates are in µm (nearest cell about
  7–8 µm; the 6th neighbour about 15 µm).
- **Edge counting:** non-primary graphs use KD-tree counts, which match
  Squidpy's `knn6` counts edge-for-edge for all eight samples. They have no
  z-scores.
- **Group means:** within a comparison, each sample counts equally (mean of
  per-sample log2 O/E). Pooling edge counts was rejected because EL_skin_res
  (about 236,000 cells) would outweigh CK_skin_res (about 6,000 cells) and it
  merges different patients' tissue. Differences are labelled descriptive only.
- **Missing cell types in combined figures:** all types are shown; values from
  fewer samples than the group are marked †. In the side-by-side difference
  panel, a pair is shown only when every sample has both types.
- **Tiles:** 500 µm squares; a pair is reported only where both types have at
  least 20 cells. Edges come from the whole-section graph and are assigned to
  the tile of their source cell.
- **Infiltration bands:** ≤20 µm (contact), 20–100 µm (margin), >100 µm
  (distant). The contact ratio is a population's ≤20 µm share divided by the
  share for all non-tumour cells (1 = what tumour shape alone gives).
- **Figures:** blue `#2a78d6` = res and orange `#eb6834` = non_res (both
  validated for colour-vision deficiency); tissue is marked by marker shape.
  Heatmaps use a blue–grey–red diverging scale centred at 0. Distance bands use
  a three-step blue ramp; a fourth, lighter step fails contrast, so at most two
  bin edges are allowed.

## How to Read the Figures

- **Per-sample or group heatmaps:** red (above 0) means two cell types are
  neighbours more often than chance; blue (below 0) means less often. Each step
  of 1 is a doubling (−3 = one eighth of chance).
- **Difference panels (res − non_res):** red means the pair is more together,
  or less apart, in responders than in non-responders; blue means the reverse.
  Red does not mean the cells cluster in responders: both groups can be below
  0 (e.g. melanoma–T in skin, −3.2 vs −5.4).
- Values describe arrangement relative to each sample's own cell mix, not the
  amounts of each cell type (see script 09 proportions for those).

## Findings So Far

Descriptive only: 4 versus 4 samples, five tissue sites, different resolutions.

- Melanoma–fibroblast segregation and strong same-type clustering appear in
  every sample and do not separate the response groups.
- Melanoma–T cell exclusion is strong in `MG_gastric_non_res` (log2 O/E −6.9),
  `MF_skin_non_res` (−5.4), and `EL_skin_res` (−5.1); the other samples range
  from −0.5 to −1.4. It persists at a 100 µm radius (−3 to −5) and in 500 µm
  tile medians (MG −5.7, EL −4.0, MF −3.8). Exclusion therefore does not
  separate `res` from `non_res`. The earlier visual reading of raw z-score
  heatmaps suggested otherwise; that was an artefact of scale.
- The GI pair is the clearest contrast. The T-cell contact ratio is 1.36 in
  `CK_bowel_res` and 0.10 in `MG_gastric_non_res`, with median T-cell distance
  to tumour of 12 µm and 142 µm. B-cell contact ratios are 0.48 and 0.04.
- T-cell contact ratios for all samples: CK_bowel 1.36, BE 1.23, HW 1.10,
  GR 0.96, CK_skin 0.74, MF 0.63, EL 0.18, MG 0.10.
- `EL_skin_res` is excluded but has a dense T-cell rim along the tumour border
  and a large nearby T-cell aggregate. `MG_gastric_non_res` has a thinner rim,
  and most T cells are distant in the surrounding tissue.
- In the excluded samples, macrophages are the immune cells in contact with
  tumour (contact ratio about 2.4–2.6 in EL and MG).
- Skin, combined responders (mean of CK_skin_res and EL_skin_res) versus
  MF_skin_non_res: almost the whole melanoma row is red in the difference
  panel (melanoma–T +2.2, melanoma–macrophage +2.0, melanoma–fibroblast +2.3),
  so MF's tumour is separated from every cell type, not specifically from T
  cells. Macrophage–T is −2.6: in responders they sit apart, in MF they mix
  about randomly. The responder mean hides a split: melanoma–T is −1.3 in
  CK_skin_res but −5.1 in EL_skin_res, so the +2.2 comes almost entirely from
  CK. Show the side-by-side figure wherever that pair matters.
- Endothelial focus pairs (2026-10-28) are exploratory; no findings recorded yet.

## Open Questions and Next Steps

1. HPC test of branch `cleanup/post08-flow` (Nathalie runs it): one-sample
   script 10 run, then script 10a and a local comparison rerun. Steps in
   `post08_flow_cleanup.md`. A full seeded script 10 rerun would replace the
   unseeded stored results; log2 O/E findings will not change.
2. Distance to the nearest tumour cell cannot distinguish intratumoural cells
   from border cells. Add an inside-tumour versus border mask, for example from
   local tumour-cell density.
3. Check whether the `EL_skin_res` aggregate has tertiary lymphoid structure
   features (B-cell and T-cell composition, organisation).
4. Review the endothelial focus-pair outputs.
5. Decide whether the PTMT cohort will use the same helpers: add a PTMT cohort
   config and helper configs pointing at it; no code changes needed. Add
   `"combine_by_response": true` to any comparison that needs a combined figure.
6. Formal res versus non_res testing is not meaningful with the current numbers
   of samples; per-sample robustness (graph sizes, tiles) is the available
   evidence.
7. Moran's I significance in script 10 is unusable as is (permutation FDR
   threshold artefact, all-NaN analytic FDR when a gene is constant, arbitrary
   top-Moran plots). Do not use the Moran FDR CSVs or plots; fix planned in
   `morans_i_significance.md`.
8. Promote the two helpers to numbered pipeline steps in the planned workflow
   package once the analysis extensions are finished (siblings; neither depends
   on the other). See `decouple_workflow_from_banksy.md`.

## Session History

- **2026-10-26/27:** seed added to script 10; comparison helper, sensitivity
  graphs, tiles, infiltration script, and shared module built and run.
  CK_bowel_res and MF_skin_non_res recomputed locally while their script 10
  results were missing.
- **2026-10-28:** `focus_pairs` changed from the original Melanoma–T,
  Melanoma–Macrophage, Melanoma–Fibroblast, T–Fibroblast, and T–Endothelial set
  to the current endothelial set. Previous outputs are backed up in
  `previous_focus_pairs_2026-10-26/` in the figure and table directories. The
  tile table and focus-pair figures now cover only the current pairs;
  `nhood_enrichment_long.csv` still holds every pair for every graph.
- **2026-11-15:** comparison figure layout fixed (title height, panel spacing);
  combined skin comparison added (`combine_by_response`). The full rerun read
  stored script 10 results for all eight samples.
- **2026-11-19:** `figure_dir` in the script 10 config and both helper configs
  moved from `figures/spatial_analysis/...` to
  `data/xenium/output/vbct/spatial_analysis/min15_area95_selected_resolutions/...`.
  Config-only change; existing figures were not moved. The HPC checkout needs
  the updated script 10 config before the next run.
- **2026-11-19/20 (post-08 cleanup, `cleanup/post08-flow`):** helpers read the
  cohort config directly; comparison helper reads per-sample enrichment tables
  (backfilled with script 10a) instead of result H5ADs, and its recompute
  fallback, `result_root`, `--recompute-missing`, and `--force-recompute` were
  removed; helper configs moved to `config/spatial_analysis/vbct/`; script 10
  renamed to `10_cell_type_spatial_analysis.py` and seeds Moran's I. All
  comparison and infiltration tables reproduced the baseline byte for byte.


## Pitfalls

- Do not compare Squidpy z-score heatmaps across samples; use log2 O/E.
- Script 08 annotated objects list all eight palette cell types as categories,
  including absent ones. Drop absent types in palette order, as script 10 does,
  before computing anything.
- After any HPC script 10 run, refresh the local tables with script 10a
  (`--overwrite`); otherwise the comparison uses the old ones.
- The Phoenix mount is only needed for script 10a (and local script 08). Check
  it with `findmnt -T ~/mnt/phoenix_scratch`; `mountpoint` reports "not a
  mountpoint" from the agent sandbox even when it is connected. Never write to
  the mount; if it is truly unavailable, report it and never reconnect or use
  SSH (see `AGENTS.md`).
- `conda run` does not forward stdin, so heredoc Python scripts silently do
  nothing. Write a script file or use `python -c`.
- Arial is not installed locally; the helpers fall back to DejaVu Sans without
  warnings. SVG text remains editable.
- Squidpy 1.6 is installed locally (`banksy` conda env); `anndata.io.read_elem`
  falls back to `anndata.experimental.read_elem` on anndata 0.10.
