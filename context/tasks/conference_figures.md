# Task: Conference Spatial Figures

## Goal

Create presentation-ready spatial figures from unfiltered Xenium/BANKSY objects for a conference poster, with enough config-driven structure that the same helper can be reused across selected exemplar samples locally or on the HPC.

## Path Update (2026-11-20, post-08 cleanup)

Paths below are historical. On branch `cleanup/post08-flow` (see
`post08_flow_cleanup.md`; merge to `main` pending an HPC test):

- The conference-figure helpers and their configs are unchanged: their config
  keeps its own sample list (absolute Phoenix paths) and palette copy, and
  writes to `figures/conference_min15_area95_selected_resolutions_poster/`.
  Moving them onto the cohort config and the output rule is a noted follow-up.
- The cohort UMAP grid and marker dot plot (`helper_scripts/poster/`) now use
  `config/poster/vbct/` and write to `data/xenium/output/vbct/poster/...`.

## Background

The conference focus is not a full cell-typing analysis. The immediate aim is to make clear spatial views that show tissue structure and highlight endothelial-rich regions. Melanoma clusters can be collapsed into a single presentation label because they are not the central focus of the poster.

The current exemplar samples are:

```text
CK_bowel_res      resolution r1p10
MG_gastric_non_res resolution r0p70
```

For MG_gastric_non_res, resolution `0.70` is the supervisor-reviewed focus resolution for now.

## Implemented In This Session

- Created and iterated on a config-driven helper:

```text
helper_scripts/conference_figures/plot_spatial_conference_figures_from_config.py
```

- Added a conference config with inline cluster-to-label maps:

```text
config/conference_figures/vbct_exemplar_ck_bowel_mg_gastric.json
```

- Added a helper README:

```text
helper_scripts/conference_figures/README.md
```

- Generated per-sample figure families into organized output subdirectories:

```text
summary/
spatial_clusters/
cell_type_labels/
umap/
marker_spatial/
dotplots/
crops/
cluster_abundance/
```

- Added endothelial marker plots and dotplots using:

```text
AQP1, CALCRL, CDH5, ECSCR, PLVAP, SELP, VWF
```

`TFPI` was deliberately left out of the current marker set.

- Added manual crop windows for MG_gastric_non_res:

```json
{
  "name": "view 1",
  "xlim": [3000, 4200],
  "ylim": [1500, 2500]
}
```

```json
{
  "name": "view 2",
  "xlim": [4550, 5500],
  "ylim": [1000, 1450]
}
```

- Added cluster crop plots with endothelial cells redrawn on top via `cluster_highlight`.
- Current endothelial highlight point size is `15`.
- Current MG view 1 cluster base point size is `15`; current MG view 2 cluster base point size is `10`.
- Added explicit configured colors through `cell_type_palette`.
- Rebuilt cropped cluster legends so they show cell-type labels, not cluster-number prefixes, with duplicate labels collapsed and sorted alphabetically.
- DAPI/OME-TIFF overlay code was explored but commented out because the DAPI-only crop was not visually convincing for this conference use case. Xenium Explorer screenshots may be better for side-by-side histology-like context.

## Current Run Command

Run MG only from inside `Banksy_py`:

```bash
conda run -n banksy python helper_scripts/conference_figures/plot_spatial_conference_figures_from_config.py \
  --config config/conference_figures/vbct_exemplar_ck_bowel_mg_gastric.json \
  --only-sample MG_gastric_non_res
```

Run all samples in the config:

```bash
conda run -n banksy python helper_scripts/conference_figures/plot_spatial_conference_figures_from_config.py \
  --config config/conference_figures/vbct_exemplar_ck_bowel_mg_gastric.json
```

## Key MG Outputs

The most actively reviewed outputs are:

```text
figures/conference/MG_gastric_non_res_unfiltered_r0p70/crops/MG_gastric_non_res_spatial_clusters_r0p70_view_1.png
figures/conference/MG_gastric_non_res_unfiltered_r0p70/crops/MG_gastric_non_res_spatial_clusters_r0p70_view_2.png
figures/conference/MG_gastric_non_res_unfiltered_r0p70/crops/MG_gastric_non_res_spatial_cell_type_labels_r0p70_view_1.png
figures/conference/MG_gastric_non_res_unfiltered_r0p70/crops/MG_gastric_non_res_spatial_cell_type_labels_r0p70_view_2.png
```

Last successful MG run reported:

```text
MG_gastric_non_res: 143,512 cells, 304 genes
Crop view 1: 11,124 cells, 324 endothelial-highlighted cells
Crop view 2: 4,011 cells, 191 endothelial-highlighted cells
```

## Next Step

Identify a good CK_bowel_res crop window that contains visually useful endothelial cells.

Suggested workflow:

1. Run or inspect the CK_bowel_res full spatial cell-type and cluster axis-scout plots.
2. Use the numeric axis-scout plots to choose candidate `xlim` and `ylim` windows.
3. Add one or more CK_bowel_res `crop_windows` entries to the conference config.
4. Re-run only CK_bowel_res with `--only-sample CK_bowel_res`.
5. Review the CK cropped cluster and cell-type plots, then tune crop coordinates and point sizes.

Useful expected CK scout outputs:

```text
figures/conference/CK_bowel_res_unfiltered_r1p10/spatial_clusters/CK_bowel_res_spatial_clusters_r1p10_axis_scout.png
figures/conference/CK_bowel_res_unfiltered_r1p10/cell_type_labels/CK_bowel_res_spatial_cell_type_labels_r1p10_axis_scout.png
```

The immediate thing to do next is pick a CK_bowel_res view with a good local endothelial signal, then add it as `view 1` for CK in the conference config.

## Plot Tuning Still Needed

Marker spatial plots such as:

```text
figures/conference/MG_gastric_non_res_unfiltered_r0p70/crops/MG_gastric_non_res_spatial_endothelial_markers_view_1.png
```

need colour-scale tuning so endothelial marker-gene expression is easier to see. The current default marker-expression scale is useful for a first pass, but for conference figures it may need per-marker or percentile-based `vmin`/`vmax` settings, a different cmap, or expression clipping to make local positive cells stand out more clearly.

## Handoff Notes: 2026-08-29 Poster Figure Polish

Current exemplar figure generation is centered on:

```text
config/conference_figures/vbct_exemplar_ck_bowel_mg_gastric.json
helper_scripts/conference_figures/plot_spatial_conference_figures_from_config.py
```

Run from inside `Banksy_py`:

```bash
conda run --no-capture-output -n banksy python helper_scripts/conference_figures/plot_spatial_conference_figures_from_config.py \
  --config config/conference_figures/vbct_exemplar_ck_bowel_mg_gastric.json
```

Current state:

- Both `MG_gastric_non_res` and `CK_bowel_res` regenerate successfully.
- CK cluster `27` / `Unassigned_27` is presentation-masked with `exclude_clusters: ["27"]`; last run excluded 60 cells.
- CK crop windows now include:

```json
{
  "name": "view 1",
  "xlim": [6700, 7900],
  "ylim": [7100, 8100],
  "point_size": 10,
  "cluster_point_size": 10
}
```

```json
{
  "name": "view 2",
  "xlim": [1950, 3150],
  "ylim": [850, 1850],
  "point_size": 10,
  "cluster_point_size": 10
}
```

- MG crop windows now include views 1, 2, and 3. View 3 is:

```json
{
  "name": "view 3",
  "xlim": [3000, 4000],
  "ylim": [2250, 2750],
  "point_size": 10,
  "cluster_point_size": 10
}
```

- Cell-type palette matching now normalizes case, spaces, underscores, and punctuation so shared biological labels use consistent colours across CK and MG.
- Marker expression plots use `Reds` with `vmin: 0` and `vmax: "p99"` instead of viridis, to avoid dark low-expression backgrounds.
- SVG output is generated for editable/poster use, but dense scatter point layers are rasterized inside SVGs so PowerPoint/Inkscape remain usable.
- Axis-scout plots intentionally do not save SVGs; they remain PNG/PDF investigation outputs only.
- Individual marker plots are written under `marker_spatial/individual_markers/` and crop-specific individual marker plots under `crops/individual_markers/`.
- Poster-scale sizing is currently configured as:

```json
{
  "dpi": 300,
  "full_spatial_figsize": [12, 14],
  "full_cell_type_figsize": [12, 14],
  "crop_figsize": [9, 7],
  "dotplot_figsize": [10, 6],
  "cluster_abundance_figsize": [12, 5],
  "marker_figsize": [20, 10],
  "individual_marker_figsize": [7, 7]
}
```

Note: `marker_figsize` and `individual_marker_figsize` are stored inside the endothelial marker-set config as `figsize` and `individual_figsize`.

Recommended poster workflow:

- Use PDFs in Inkscape for cleanup and composition of figure elements.
- Use SVGs when PowerPoint editability is useful and import behaves well.
- Use PNGs as reliable preview/fallback exports.

Next thing to fix:

- Dotplots became too spread out after increasing `dotplot_figsize` to `[10, 6]`. Next time, make dotplots more compact so the dots sit closer together, while still saving PDF outputs for Inkscape/poster use. The helper currently saves dotplots as PNG and SVG only, so add PDF saving to `plot_marker_dotplot()` and either reduce `dotplot_figsize` or add dotplot-specific spacing controls. A likely first try is `dotplot_figsize: [7, 4.5]` or `[8, 4.5]`, then regenerate only the dotplots or rerun the full config.

## Conference and Poster Output Map

The directory inventory below was verified locally on 2026-11-18. Paths are
relative to the `Banksy_py/` project root. Individual files may be regenerated
or added later.

### Quick Navigation

The main output roots are:

```text
figures/conference/
figures/conference_filtered_reclustered_qc_v1/
figures/conference_filtered_reclustered_qc_v1_broad_cell_types/
figures/conference_min15_area95_selected_resolutions_poster/
figures/dotplots/vbct/recluster_filtered_qc_pass_only/
```

Use
`figures/conference_filtered_reclustered_qc_v1_broad_cell_types/Final_poster_figure/`
for the hand-selected and edited assets assembled for the final poster. The
other directories contain reproducible figure families generated from configs.

### Unfiltered CK/MG Exemplar Figures

Generated by
`config/conference_figures/vbct_exemplar_ck_bowel_mg_gastric.json`:

```text
figures/conference/CK_bowel_res_unfiltered_r1p10/
figures/conference/MG_gastric_non_res_unfiltered_r0p70/
```

An earlier single-sample CK skin figure set is retained at:

```text
figures/conference/CK_skin_res_unfiltered_r1p10/
```

The CK bowel and MG gastric directories contain the standard output families
described below, plus `umap_rework/` and `umap_outlined/` variants.

### Filtered and Reclustered CK/MG Figures

The reviewed post-QC recluster outputs are:

```text
figures/conference_filtered_reclustered_qc_v1/
├── CK_bowel_res_filtered_reclustered_qc_v1_r1p10/
└── MG_gastric_non_res_filtered_reclustered_qc_v1_r0p70/
```

These are generated from
`config/conference_figures/vbct_exemplar_ck_bowel_mg_gastric_filtered_reclustered_qc_v1.json`.
In addition to the standard figure families, these sample directories contain
`individual_clusters/`, `umap_rework/`, and `umap_outlined/` where those helper
runs have been performed.

The broad presentation-label version is separate:

```text
figures/conference_filtered_reclustered_qc_v1_broad_cell_types/
├── CK_bowel_res_filtered_reclustered_qc_v1_r1p10/
├── MG_gastric_non_res_filtered_reclustered_qc_v1_r0p70/
└── Final_poster_figure/
```

It is generated from
`config/conference_figures/vbct_exemplar_ck_bowel_mg_gastric_filtered_reclustered_qc_v1_broad_cell_types.json`.

#### Final Poster Staging Area

`figures/conference_filtered_reclustered_qc_v1_broad_cell_types/Final_poster_figure/`
contains the selected, copied, and manually edited poster assets. Its organized
subdirectories are:

```text
cell_types/
crops/
dotplot/
individual_cluster/
marker_spatial/
multimodal/
umap/
workflow/
```

The same directory also contains the assembled ICI cohort figure and poster
logos/assets. Treat this as the presentation assembly area rather than a fully
reproducible script output directory.

### Eight-Sample min15/area95 Poster Figures

The current eight-sample selected-resolution output root is:

```text
figures/conference_min15_area95_selected_resolutions_poster/
```

It is configured by
`config/conference_figures/vbct_min15_area95_selected_resolutions.json` and
contains:

```text
BE_brain_non_res_min15_area95_r0p90/
CK_bowel_res_min15_area95_r0p80/
CK_skin_res_min15_area95_r1p00/
EL_skin_res_min15_area95_r0p80/
GR_lung_non_res_roi_min15_area95_r1p00/
HW_brain_res_min15_area95_r0p80/
MF_skin_non_res_min15_area95_r0p50/
MG_gastric_non_res_min15_area95_r0p70/
cohort_marker_dotplot/
cohort_umap_grid/
```

`cohort_marker_dotplot/` contains the combined canonical-marker dotplots.
`cohort_umap_grid/` contains the combined 2-by-4 and 4-by-2 UMAP panels, with
older variants retained under `cohort_umap_grid/umaps_old/`.

There is also an older/partial output root at
`figures/conference_min15_area95_selected_resolutions/`. The poster-suffixed
directory above is the current configured destination.

### Dotplots Stored Outside the Conference Roots

Post-QC, all-resolution dotplots and source plots used by the conference configs
live under:

```text
figures/dotplots/vbct/recluster_filtered_qc_pass_only/
```

The focused annotated CK `r1.10` and MG `r0.70` comparison is:

```text
figures/dotplots/vbct/recluster_filtered_qc_pass_only/CK_bowel_res_MG_gastric_non_res_selected_res_all_genes_conference_annotated_zscore_dotplot.png
```

The filtered CK/MG conference configs also refer to the corresponding
sample-level `*_all_genes_recluster_qc_pass_only_all_res_zscore_dotplot.png`
files in this directory.

### Standard Per-Sample Directory Layout

The main config-driven helper organizes each sample output as follows:

| Directory | Contents |
|---|---|
| `summary/` | Spatial-coordinate summary CSVs. |
| `spatial_clusters/` | Full-tissue BANKSY cluster maps, on-data labels, and axis-scout plots. |
| `cell_type_labels/` | Full-tissue reviewed/presentation cell-type maps and axis-scout plots. |
| `umap/` | UMAPs generated during the full figure-pack run. |
| `umap_rework/` | Lighter UMAP variants generated by the dedicated UMAP helper. |
| `umap_outlined/` | Cartoon-style outlined UMAP variants. |
| `marker_spatial/` | Combined endothelial marker-expression panels. |
| `marker_spatial/individual_markers/` | One full-tissue spatial plot per marker. |
| `dotplots/` | Endothelial marker dotplots for that sample and selected resolution. |
| `crops/` | Configured spatial close-ups for clusters, cell types, and marker panels. |
| `crops/individual_markers/` | Crop-specific plots for individual marker genes. |
| `cluster_abundance/` | Cluster abundance bar plots. |
| `individual_clusters/` | One highlighted spatial map per selected BANKSY cluster. |

Most generated figures are available as PNG plus PDF and/or SVG. Dense cell
scatter layers are rasterized inside PDF/SVG files so labels, axes, legends,
and colour bars remain editable without making the vector files impractically
large. Axis-scout plots are investigation outputs and are intentionally not
saved as SVG.

### Code and Configuration Locations

The output roots above are controlled by:

```text
config/conference_figures/
helper_scripts/conference_figures/
```

The main generator is
`helper_scripts/conference_figures/plot_spatial_conference_figures_from_config.py`.
The helper README contains run commands and config-field documentation. These
are source/control files, not figure outputs.
