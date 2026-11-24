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
  v1: all eight samples. v2 (draft): BE_brain_non_res and HW_brain_res only
  (the samples v2 changes: BE c5 and HW c15 Epithelial → Melanoma). Rerun the
  other six once v2 is final.

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

**3. Cohort poster dot plot** (one row per label, samples averaged with equal
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

## Next Step

Issue #64 (supervisor meeting 2026-11-19): per-sample cell-type marker dot
plots instead of the pooled cohort plot.
1. Per-sample rows rather than pooled.
2. Drop cell types that are implausible for the tissue (e.g. epithelial in
   brain). Nathalie decides the list.
3. Use the most recognisable marker per cell type.
4. Check whether CDH1 is in the panel and whether it is a good epithelial marker.

Related: #63 (markers for myofibroblasts, neurons, dendritic cells).

## Open Questions

- **Annotation version:** v1 (released) or v2 once final?
- **Build route for per-sample label rows:**
  - (b1) Script 04 from the existing exports plus the set's mapping, no H5ADs.
    Label values must be cell-count-weighted over their clusters (`n_cells`),
    not averaged equally. Needs the `annotation` block to accept several samples.
  - (b2) Add a per-sample row option to the poster helper, which computes label
    values from the annotated H5ADs.
- **z-score scope** for any cross-sample figure: within each sample (contrast
  between clusters) or across the cohort (comparable colours, but high-expressing
  samples dominate).
- **All 127 clusters in one figure** (stacked by sample) is the alternative
  to label rows; it is tall.
- **Small clusters seen in the v1 plots** (to judge, not conclusions):
  CK_bowel_res c16 (n=376) and c17 (n=238), labelled melanoma, each with one
  strongly expressed atypical gene; CK_bowel_res c18, a 28-cell T-cell cluster;
  GR_lung_non_res_roi c13 (n=403), labelled melanoma, B-cell-like signal.

## History

- To 2026-11-24, including script 02/03/04 notes, PTMT sections and the
  superseded next steps: `context/tasks/dotplot/archive/dotplot_to_2026-11-24.md`.
