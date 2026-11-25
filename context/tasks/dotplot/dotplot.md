# Task: VBCT Marker Dot Plots

<!-- Task file rules (see Banksy_py/AGENTS.md, "Task files"): edit Current
State, Next Step and Open Questions in place; add one dated line per decision;
no "Session Update" sections; move superseded text to archive/ beside this file;
keep this file under ~200 lines. -->

## Goal

Marker-expression dot plots for the VBCT cohort, used to review and justify
cluster annotations. PTMT dot plots are in `ptmt_dotplot.md`.

**Core rule:** summarise clean biological expression only (`.raw`, or the
clean `.X`), never BANKSY-expanded `adata_spatial.X`. Run
`tools/inspect_h5ad.py` on any object before summarising it.

**Tools:** script 03 exports per-cluster summaries (HPC, or locally from the
Phoenix mount); script 04 renders them locally. The cohort poster dot plot is
`helper_scripts/poster/plot_cohort_marker_dotplot_from_config.py`. Options and
config layout: `config/04_plot_dotplot/README.md`; full description:
`README_00_04_XENIUM_ANALYSIS_WORKFLOW.md`. Use script 02 only for legacy or
archived-label objects.

## Current State (updated 2026-11-24)

**Gene sources** (`data/xenium/raw_data/gene_markers/`):
- `01_xenium_gene_list_annotation_master_v1_May_2026_Tier1_2026-06-15.csv`:
  304 panel genes, grouped by `Tier_1_annotation`.
- `02_xenium_melanoma_canonical_genes_2026-05-22.csv`: canonical genes,
  highlighted red.
- `03_xenium_gene_review_for_dotplot_manual_annotation_2026-06-11.csv`:
  review filter (`manual_review_True_if_keep`), keeps 89 genes. It drops
  several subtype/state markers (e.g. TOP2A, MZB1, CLEC9A, CCR7, LYVE1, PROX1,
  SELE, S100B).

**1. Per-sample cluster plots labelled by annotation set** (made 2026-11-23,
to review clusters before the v2 merges)
- Configs: `config/04_plot_dotplot/active/vbct_annotation_sets/min15_area95_selected_resolutions/`
  (one per sample; README there). Uses script 04's `annotation` block.
- Rows: one sample's clusters at its selected resolution (min15_area95
  branch), labelled `cluster <id> · <label> (n=<cells>)`, grouped by label.
- Input: the `min15_area95` rows of the filter-sensitivity exports (below); MF
  uses `MF_skin_non_res_min15_area95_all_res_dotplot_summary.csv`. Checked
  against v1's `cluster_annotation_counts.csv` on Phoenix: all 127 clusters and
  cell counts identical, 304 genes from `.raw`.
- Outputs: `data/xenium/output/vbct/cell_type_annotation/min15_area95_selected_resolutions/<version>/tier_1/cluster_marker_dotplots/`.
  v1 and v2: all eight samples. v2 (released 2026-11-24) changes only BE c5
  and HW c15 (Epithelial → Melanoma); all eight v2 plots were re-rendered on
  2026-11-24 after release (127 clusters; titles read "released").
- Under v2, epithelial cells remain only in MG_gastric_non_res (17.0%),
  CK_bowel_res (4.3%) and EL_skin_res (0.6%) (script 09 v2 proportions).
- CDH1 is not in the 304-gene panel. The panel's epithelial genes are
  skin-derived (keratinocyte/eccrine: KRT5, KRT15, KRT18, CLDN4, TP63, TACSTD2,
  AQP3, DSC1, LY6D, S100A14, PTGR1; LYZ is mis-tiered as epithelial); no
  EPCAM, KRT8 or KRT19.

```bash
C=config/04_plot_dotplot/active/vbct_annotation_sets/min15_area95_selected_resolutions
for f in $C/*.json; do python 04_plot_multi_sample_dotplot_from_config_local.py --config "$f"; done
# other version/tier: add --annotation config/annotations/vbct/min15_area95_selected_resolutions/v2 [--label-level tier_2]
```

**2. QC filter-sensitivity plots** (2026-10)
- Each sample at its selected resolution, five conditions stacked as
  pseudo-samples: `prefiltered` (script 00 object and clusters), `min10_area99`,
  `min10_area95`, `min15_area99`, `min15_area95`.
- Configs: `config/03_export_summary/active/vbct_filter_sensitivity/selected_resolutions/`
  and `config/04_plot_dotplot/active/vbct_filter_sensitivity/selected_resolutions/`
  (README there: resolutions, commands, mount-based local export). Extra
  comparisons (EL 1.50, GR/MG 1.00, MG 0.70) sit in the parent folders.
- Local exports: `data/xenium/processed/vbct/cross_sample_dotplot_exports/filter_sensitivity/`;
  figures: `figures/dotplots/filter_sensitivity/`. Seven samples have
  five-condition plots; MF has only a min15_area95 all-resolution plot.

**3. Paper cell-type marker dot plots** (issue #64, rendered 2026-11-24)
- Script `helper_scripts/paper_figures/plot_cell_type_marker_dotplot_from_config.py`;
  configs `config/paper_figures/vbct/cell_type_marker_dotplot.json` (by sample,
  z-score within sample) and `..._cohort.json` (cell-type rows, equal sample
  weight, z-score across cell types); a third, `..._by_sample_global_zscore.json`,
  z-scores across all 50 rows as a supplementary detection panel (MF near white;
  HW immune rows carry melanoma markers); README in `config/paper_figures/`. v2,
  tier 1; reads the same script 03 exports as (1).
- 127 clusters merged into 50 sample x cell type rows; counts match script 09 v2.
- Outputs: `data/xenium/{output,processed}/vbct/paper_figures/min15_area95_selected_resolutions/v2/tier_1/cell_type_marker_dotplot/`
  (figure; merged expression and marker specificity tables).
- Markers are provisional (poster set; epithelial KRT18, CLDN4, TACSTD2).
  Specificity summary flags: B-cell markers MZB1 and TNFRSF17 are plasma-cell
  genes and rank top with CD79A (MS4A1/CD19/CD79B absent); CLDN4 is weak in
  epithelial cells (12% expressing), S100A14 and MAL2 are more specific; Pericyte
  occurs only in GR (207 cells) and RGS5 is weak there while MYH11 is specific
  (smooth-muscle-like); C1QA is the weakest macrophage pick (FCER1G, AIF1 rank
  higher). EL macrophages show MITF/MLANA signal. Endothelial rows carry
  RGS5/NOTCH3/ACTA2 in most samples (pericytes likely inside endothelial
  clusters where no Pericyte label exists).

**4. Cohort poster dot plot** (one row per label, samples averaged with equal
weight, genes min-max scaled)
- Config: `config/poster/vbct/poster_cohort_marker_dotplot_rdbu_large_font.json`
  (v1, tier 1, three markers per label; epithelial uses KRT18 and KRT5 only).
  Poster decisions are in `hca_poster.md`.

**Older configs** (`vbct_clean_script00/`, `vbct_recluster_qc_pass_only/`,
`vbct_recluster_grouped_area95_trial/` under `config/03_export_summary/active/`
and `config/04_plot_dotplot/active/`) predate the QC branches; not in current use.

## Decisions

- 2026-10-23: Filter branches are reclustered independently. Cluster IDs only
  identify clusters within their own branch; never transfer labels between
  branches by matching IDs.
- 2026-10-23: Filter-sensitivity colour is z-scored per gene across all
  conditions within one sample, clipped at 2.5; it is not comparable across
  samples.
- 2026-10-24: Added the `prefiltered` baseline as the first condition.
- 2026-11-23 (Nathalie): Annotation-labelled plots use the 89 review-filtered
  genes, not all 304; a 304-gene version stays an option if clusters need
  splitting.
- 2026-11-23 (Nathalie): Annotation-labelled plots z-score per gene across the
  plotted sample's clusters, clipped at ±2.5; dot size is percent expressing.
- 2026-11-24 (Nathalie): The paper's main marker figure merges clusters per
  cell type within each sample, rows grouped by sample (skin, brain, GI, lung;
  res before non-res), three markers per cell type, no minimum-cell cutoff;
  per-cluster plots are supporting.
- 2026-11-24 (Nathalie): Make both a by-sample plot and a cohort plot (equal
  sample weight). Colour is a Seurat-style per-gene z-score (clip ±2.5): within
  each sample for the by-sample plot, across cell types for the cohort plot.
  Replaces per-gene min-max across all rows, which made low-detection samples
  (MF) and non-maximal samples look pale.
- 2026-11-24 (Nathalie): Paper figure scripts live in `helper_scripts/paper_figures/`
  (post-08, named by content), not as a numbered pipeline step.

## Next Step

Issue #64: Nathalie reviews `marker_specificity_summary.csv` and the two paper
dot plots (State 3), then finalises three markers per cell type (same list in
both configs) and decides how to handle the B-cell (plasma-cell) and Pericyte
labels.

Related: #63 (markers for myofibroblasts, neurons, dendritic cells).

## Open Questions

- **B-cell label:** the "B cell" rows are dominated by plasma-cell genes
  (MZB1, TNFRSF17); relabel (e.g. "B/plasma cell") in a new annotation set?
- **Pericyte:** only in GR (207 cells), ACTA2/MYH11-high; pericyte, smooth
  muscle or myofibroblast? Overlaps #63.
- **Which plot is main:** by-sample or cohort plot as the main figure, the
  other as supplementary.
- **BE melanoma clusters 9 and 11:** little MITF/MLANA/PMEL signal; check the
  label (annotation question, not plotting).
- **Small clusters seen in the v1 plots** (to judge, not conclusions):
  CK_bowel_res c16 (n=376) and c17 (n=238), labelled melanoma, each with one
  strongly expressed atypical gene; CK_bowel_res c18, a 28-cell T-cell cluster;
  GR_lung_non_res_roi c13 (n=403), labelled melanoma, B-cell-like signal.

## History

- To 2026-11-24, including script 02/03/04 notes, PTMT sections and the
  superseded next steps: `context/tasks/dotplot/archive/dotplot_to_2026-11-24.md`.
