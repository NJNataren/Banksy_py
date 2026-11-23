# Task: Xenium Dotplot Export And Plotting

## Goal

Maintain and extend the Xenium dotplot workflow for clean expression summaries and multi-sample plotting. Dotplots should use biological expression values from clean expression AnnData objects, not BANKSY-expanded `adata_spatial.X` matrices.

## Current Status

- Script 00 now creates the preferred clean clustered AnnData handoff:
  - `data/xenium/processed/<project>/<sample>/adata_expression_clean_<sample>_with_banksy_clusters_<resolutions>.h5ad`
- Current BANKSY cluster dotplots should usually skip script 02 and use the clean script 00 object directly with:
  - `03_export_dotplot_data_from_clean_script00_config.py`
- Script 02 remains useful for archived-label augmentation or legacy regeneration from separate expression and `adata_spatial_*` files.
- Script 03 exports long-format dotplot summary CSVs.
- Script 04 renders multi-sample dotplot PNGs from script 03 CSVs.
- PTMT panel marker/review files now exist under `data/xenium/raw_data/gene_markers/ptmt/`.

## Script 02: Clean Expression Object Creation Or Archived-Label Augmentation

Main script:

```text
02_create_expression_adata_with_banksy_clusters.py
```

Script 02 reads a clean/count-like expression object, normalizes/log-transforms it, and copies label metadata from one or more `adata_spatial_*` objects into `.obs`. It can also attach archived cell-level labels.

Current config directories:

```text
config/02_create_expression/vbct/
config/02_create_expression/vbct_archive/
config/02_create_expression/testing/
```

Use script 02 when:

- cluster labels only exist in separate `adata_spatial_*` objects;
- archived labels need to be joined onto clean expression objects;
- a legacy workflow expects `*_normalised_log1p_with_banksy_clusters.h5ad` outputs.

Do not use script 02 just to make current BANKSY cluster dotplots if script 00 already saved the clean clustered object with all needed `.obs` label columns.

## Script 03: Dotplot Summary Export

Legacy/general script:

```text
03_export_dotplot_data_from_config.py
```

Current clean script 00 export script:

```text
03_export_dotplot_data_from_clean_script00_config.py
```

The clean script 00 exporter reads one clean clustered AnnData object per sample and summarizes clean expression from `.raw` when available, otherwise `.X`. It groups cells using BANKSY cluster columns already present in `.obs`.

For cross-resolution per-sample review, script 03 can now read script 01 QC configs and attach human-readable `new_labels` to the exported `group_label` for the configured `cluster_col` resolution. This is enabled by default when `project` is present and `config/01_QC/<project>/<sample>.json` exists; set `use_qc_annotations: false`, `qc_config_dir`, or per-sample/per-grouping `qc_config` to override that behavior. Numeric `cluster_id`/`group_id` values are still preserved for ordering and traceability.

It validates that:

- the configured AnnData object has cells and genes;
- each requested `.obs` grouping column exists;
- each requested grouping column has at least one non-null label.

Current clean script 00 config directories:

```text
config/03_export_summary/active/vbct_clean_script00/all_genes/
config/03_export_summary/active/vbct_clean_script00/canonical_markers/
config/03_export_summary/active/ptmt_clean_script00/ptmt_panel_pc35/
```

PTMT script 03 configs use:

```text
data/xenium/raw_data/gene_markers/ptmt/ptmt_combined_xenium_panel_markers.csv
```

with:

```json
"gene_column": "Gene",
"marker_group_column": "Tier_1_annotation"
```

PTMT script 03 HPC submission:

```bash
sbatch --array=0-23%3 \
  --export=CONFIG_DIR=config/03_export_summary/active/ptmt_clean_script00/ptmt_panel_pc35 \
  run_03_xenium_dotplot_export_from_clean_script00_config.sl
```

Expected PTMT script 03 output directory:

```text
data/xenium/processed/ptmt/cross_sample_dotplot_exports/
```

## Script 04: Multi-Sample Dotplot Rendering

Local plotting script:

```text
04_plot_multi_sample_dotplot_from_config_local.py
```

Script 04 reads one or more script 03 summary CSVs and creates a combined dotplot PNG. It can:

For per-sample all-resolution plots, `sample_filters` now accepts `"resolution": "all"`, `"include_all_resolutions": true`, or a `"resolutions"` list. Resolution filtering and ordering compare numeric-equivalent labels such as `0.7` and `0.70`, so older and newer CSVs can be mixed more safely.

- select samples and resolutions with `sample_filters`;
- use a marker/gene file to choose and group genes;
- apply a gene review file to remove genes from plotting;
- z-score mean expression within each gene across plotted rows;
- scale dot size using `percent_expressing`;
- write a configured PNG output.

Current local configs include:

```text
config/04_plot_dotplot/active/vbct_clean_script00/per_sample_canonical_markers_all_res/CK_skin_res_canonical_markers_all_res_local.json
config/04_plot_dotplot/active/vbct_clean_script00/multi_sample/all_genes_multi_sample_local.json
config/04_plot_dotplot/active/vbct_clean_script00/multi_sample/canonical_markers_multi_sample_local.json
config/04_plot_dotplot/archive/vbct_archive_labels/all_genes_archive_cell_type_labels_multi_sample_local.json
config/04_plot_dotplot/archive/vbct_archive_labels/canonical_markers_archive_cell_type_labels_multi_sample_local.json
config/04_plot_dotplot/active/ptmt_clean_script00/multi_sample/ptmt_panel_multi_sample_local.json
```

The PTMT plotting config reads the 24 PTMT panel summary CSVs from:

```text
data/xenium/processed/ptmt/cross_sample_dotplot_exports/*_ptmt_panel_dotplot_summary.csv
```

It uses:

```text
data/xenium/raw_data/gene_markers/ptmt/ptmt_combined_xenium_panel_markers.csv
data/xenium/raw_data/gene_markers/ptmt/ptmt_combined_xenium_panel_gene_review_for_dotplot.csv
```

The initial PTMT plotting config filters every sample to resolution `0.7` and writes:

```text
figures/dotplots/local_ptmt_panel_multi_sample_r0p70_zscore_dotplot.png
```

Run locally after copying script 03 PTMT summary CSVs back from the HPC:

```bash
conda run -n banksy python 04_plot_multi_sample_dotplot_from_config_local.py \
  --config config/04_plot_dotplot/active/ptmt_clean_script00/multi_sample/ptmt_panel_multi_sample_local.json
```

## Key Output Columns

Script 03 summary CSVs should contain columns including:

```text
sample
resolution
cluster_id
group_id
group_label
sample_cluster
sample_group
groupby
groupby_label
gene
marker_group
mean_expression
percent_expressing
n_cells
expression_source
adata_path
```

Script 04 uses `mean_expression` or `mean_expression_zscore` for dot colour and `percent_expressing` for dot size.

## PTMT Marker Panel Files

Combined PTMT panel marker file:

```text
data/xenium/raw_data/gene_markers/ptmt/ptmt_combined_xenium_panel_markers.csv
```

Review scaffold:

```text
data/xenium/raw_data/gene_markers/ptmt/ptmt_combined_xenium_panel_gene_review_for_dotplot.csv
```

The combined marker file was built from:

```text
hImmuno_submitted_to_designer_tier1.csv
hProstate_100g_with_Tier_1_annotations.csv
```

It contains one row per unique gene. Genes present in both source panels have `panel_source` and `Tier_1_annotation` values joined by `;` for human review.

## Validation

Useful checks:

```bash
python -m py_compile 03_export_dotplot_data_from_clean_script00_config.py \
  04_plot_multi_sample_dotplot_from_config_local.py

python3 -m json.tool config/04_plot_dotplot/active/ptmt_clean_script00/multi_sample/ptmt_panel_multi_sample_local.json
```

Full script 03 and script 04 runs need the relevant `.h5ad` or summary CSV files to exist on the machine running the command.

## Next Steps

The next task is to adjust the dot spacing and related scale settings in `04_plot_multi_sample_dotplot_from_config_local.py` so dots sit closer together on the plot, possibly overlapping.

Likely areas to inspect:

- figure sizing parameters in the config, especially `width_per_gene`, `height_per_cluster`, `min_width`, and `min_height`;
- dot size scaling via `max_dot_size` and the code path that maps `percent_expressing` to marker area;
- axis/grid layout in the plotting code, including how categorical x/y positions are assigned;
- whether a config option should allow tighter row/column spacing independently from figure size;
- whether overlapping dots remain readable after z-scoring and gene group labelling.

Start with the PTMT config because it has many samples and many genes, so it exposes the spacing problem clearly.

## Session Update: Config Cleanup and PTMT PC55 Dotplots

Script 03 export configs have been reorganised under:

```text
config/03_export_summary/active/
config/03_export_summary/local_helpers/
config/03_export_summary/archive/
```

Current PTMT config locations:

```text
config/03_export_summary/active/ptmt_clean_script00/ptmt_panel_pc35/
config/03_export_summary/active/ptmt_pc55_clean_script00/ptmt_panel/
config/03_export_summary/active/ptmt_pc55_clean_script00/all_genes/
```

PTMT PC55 panel exports write to:

```text
data/xenium/processed/ptmt_pc55/cross_sample_dotplot_exports/*_ptmt_panel_dotplot_summary.csv
```

Script 04 plot configs have been reorganised under:

```text
config/04_plot_dotplot/active/
config/04_plot_dotplot/archive/
```

The config layout is documented in:

```text
config/04_plot_dotplot/README.md
```

Current PTMT PC55 plotting configs:

```text
config/04_plot_dotplot/active/ptmt_pc55_clean_script00/multi_sample/ptmt_pc55_panel_multi_sample_local.json
config/04_plot_dotplot/active/ptmt_pc55_clean_script00/multi_sample/ptmt_pc55_panel_run_1_local.json
config/04_plot_dotplot/active/ptmt_pc55_clean_script00/per_sample_panel_all_res/run_1/*.json
config/04_plot_dotplot/active/ptmt_pc55_clean_script00/per_sample_all_genes_all_res/*.json
```

The temporary `available` PC55 plot configs were removed after `10708_run_1_1644` script 03 results were copied back. The canonical run 1 PC55 plot config now includes all eight run 1 samples.

Run the PTMT PC55 run 1 panel plot with:

```bash
conda run -n banksy python 04_plot_multi_sample_dotplot_from_config_local.py \
  --config config/04_plot_dotplot/active/ptmt_pc55_clean_script00/multi_sample/ptmt_pc55_panel_run_1_local.json
```

Run one per-sample all-resolution PTMT PC55 panel dotplot per run 1 sample with:

```bash
for cfg in config/04_plot_dotplot/active/ptmt_pc55_clean_script00/per_sample_panel_all_res/run_1/*.json; do
  conda run -n banksy python 04_plot_multi_sample_dotplot_from_config_local.py \
    --config "$cfg"
done
```

Important distinction: existing PTMT panel script 03 exports already contain all available resolutions for the PTMT panel genes. New script 03 `all_genes` configs are only needed if the intended plot is literally every gene in the clean AnnData object, not just every gene in the PTMT panel/review file.

Red gene-label highlighting in PTMT plots is controlled by:

```json
"highlight_gene_file": "data/xenium/raw_data/gene_markers/ptmt/ptmt_combined_xenium_panel_gene_review_for_dotplot.csv",
"highlight_gene_column": "colour_canonical"
```

with the visual style set in the `figure` block by `highlight_label_color` and `highlight_label_weight`.

## Session Update: PTMT PC55 scRNA-seq Marker Dotplots

The PTMT PC55 run 1 dotplots were extended to use the harmonised scRNA-seq
marker file:

```text
data/xenium/raw_data/gene_markers/ptmt/2026-08-14_dotplot_marker_genes_harmonised_scRNAseq_analysis_Kat.csv
```

The plots use the eight `run_1` PC55 samples. Marker genes are grouped by
`minor_class` and ordered using `plot_order`. Gene clustering and the gene
dendrogram are disabled, so the plotted gene order follows the marker CSV
rather than a hierarchical clustering result.

Resolution-specific configs were created for BANKSY resolutions 0.7, 1.0,
and 1.5. The current resolution 1.0 config displays marker labels by joining
the `gene` and `cell_types` columns with ` - `, for example:

```text
ANPEP - LE-KLK3; Mono
```

Script 04 was refactored to introduce a `GeneMetadata` dataclass. Gene order,
groups, display labels, highlighted genes, review-filter state, and optional
clustering metadata now travel together through the plotting workflow. The
script also supports the config options:

```json
"gene_order_column": "plot_order",
"gene_label_columns": ["gene", "cell_types"],
"gene_label_separator": " - "
```

The active PTMT PC55 multi-sample configs are now:

```text
config/04_plot_dotplot/active/ptmt_pc55_clean_script00/multi_sample/ptmt_pc55_panel_multi_sample_local.json
config/04_plot_dotplot/active/ptmt_pc55_clean_script00/multi_sample/ptmt_pc55_minor_class_run_1_r_0p70_local.json
config/04_plot_dotplot/active/ptmt_pc55_clean_script00/multi_sample/ptmt_pc55_minor_class_run_1_r1p00_gene_cell_types_labels_local.json
config/04_plot_dotplot/active/ptmt_pc55_clean_script00/multi_sample/ptmt_pc55_minor_class_run_1_r1p50_local.json
```

The all-sample Xenium panel config remains active. Superseded major-class,
run-1 panel, plain resolution 1.0, and cell-type-sorted experimental configs
were moved to:

```text
config/04_plot_dotplot/archive/ptmt_pc55_config_trials/
```

The latest plot was generated with:

```bash
conda run -n banksy python 04_plot_multi_sample_dotplot_from_config_local.py \
  --config config/04_plot_dotplot/active/ptmt_pc55_clean_script00/multi_sample/ptmt_pc55_minor_class_run_1_r1p00_gene_cell_types_labels_local.json
```

It writes:

```text
figures/dotplots/local_ptmt_pc55_scRNA-seq_markers_minor_class_run_1_r1p00_plot_order_gene_cell_types_labels_zscore_dotplot.png
```

All four active configs render at 150 DPI.

## Current Next Step

The next task is to review the `cell_types` column in
`2026-08-14_dotplot_marker_genes_harmonised_scRNAseq_analysis_Kat.csv` and
choose which cell-type label should be displayed for each marker gene. Some
rows contain multiple candidate labels separated by semicolons. Once those
preferred labels have been selected, store them in a dedicated display-label
column and update `gene_label_columns` in the active resolution 1.0 config to
use that curated column.

This current next step supersedes the older dot-spacing task described above.

## Next Task: VBCT Filter-Sensitivity All-Gene Dotplots

Create all-gene dotplots per VBCT biological sample across the full available
BANKSY reclustering range (`0.50` through `1.50`), rather than starting at
`1.00`. Use per-resolution figures for cross-condition comparisons when needed
and all-resolution single-branch plots for resolution review.

### Design And Interpretation

Represent filter branches as pseudo-samples and stack them on the y-axis in
this fixed order:

1. `min10_area99`: minimum 10 transcripts and per-cluster 99th
   area-percentile removal.
2. `min10_area95`: minimum 10 transcripts and per-cluster 95th
   area-percentile removal.
3. `min15_area99`: minimum 15 transcripts and per-cluster 99th
   area-percentile removal.
4. `min15_area95`: minimum 15 transcripts and per-cluster 95th
   area-percentile removal.

Use the existing multi-sample workflow: script 03 exports these four branches
as pseudo-samples, and script 04 stacks them using `sample_order`. Use a
row-label template equivalent to `"{sample} | cluster {group_label}"`.

Each branch is reclustered independently. Numeric cluster IDs only identify
clusters within their own branch: matching IDs across filter conditions do not
establish biological equivalence. Do not propagate annotations or merge maps
between conditions merely by matching the numeric IDs.

### Inputs And Settings

Script 03 requires one clean script 06 reclustering H5AD per branch and sample.
That H5AD already contains all reclustering resolutions; do not create or copy
an H5AD per resolution. At 1.00, use the grouping column from the corresponding
script 06 config, conventionally:

```text
labels_scaled_gaussian_pc<pcs>_nc0.20_r1.00_recluster_<run_label>
```

Take `<pcs>` and `<run_label>` from that sample's script 06 config. In
particular, `MF_skin_non_res` uses 20 PCs; settings are not universal.

Export all genes from `.raw`, excluding `excluded_by_qc`. Retain the existing
gene ordering, review, and highlight sources:

```text
data/xenium/raw_data/gene_markers/01_xenium_gene_list_annotation_master_v1_May_2026_Tier1_2026-06-15.csv
data/xenium/raw_data/gene_markers/02_xenium_melanoma_canonical_genes_2026-05-22.csv
data/xenium/raw_data/gene_markers/03_xenium_gene_review_for_dotplot_manual_annotation_2026-06-11.csv
```

Keep shared expression z-scoring across all four branches, clipped at 2.5, so
colour is comparable within a biological sample.

### Templates And Run Status

GR and MG are the working templates:

```text
config/03_export_summary/active/vbct_filter_sensitivity/GR_lung_non_res_roi_r1.00_all_genes.json
config/03_export_summary/active/vbct_filter_sensitivity/MG_gastric_non_res_r1.00_all_genes.json
config/04_plot_dotplot/active/vbct_filter_sensitivity/GR_lung_non_res_roi_r1.00_filter_sensitivity_local.json
config/04_plot_dotplot/active/vbct_filter_sensitivity/MG_gastric_non_res_r1.00_filter_sensitivity_local.json
```

As of 2026-10-23, all four clean script 06 branches existed for
`BE_brain_non_res`, `CK_skin_res`, `EL_skin_res`, `GR_lung_non_res_roi`, and
`MG_gastric_non_res`. The following were being run on HPC and must be checked
before their dotplot configs are made: `CK_bowel_res` min10/area95,
`HW_brain_res` min10/area95, and `MF_skin_non_res` min10/area95 plus
min10/area99.

### Execution Sequence

1. Confirm all four clean script 06 H5AD inputs exist for a sample.
2. Copy a GR or MG script 03 config and update its four branch paths,
   resolution-specific grouping columns, and output stem.
3. Run script 03 on HPC, where the H5ADs reside:

   ```bash
   conda run -n banksy python 03_export_dotplot_data_from_clean_script00_config.py \
     --config config/03_export_summary/active/vbct_filter_sensitivity/<sample>_r1.00_all_genes.json
   ```

4. Transfer only the produced summary CSV locally, avoiding four large H5AD
   transfers per sample.
5. Copy the matching script 04 config, update its local CSV path and figure
   stem, then render locally:

   ```bash
   conda run -n banksy python 04_plot_multi_sample_dotplot_from_config_local.py \
     --config config/04_plot_dotplot/active/vbct_filter_sensitivity/<sample>_r1.00_filter_sensitivity_local.json
   ```

6. Validate the CSV has exactly the four pseudo-samples in the prescribed
   order and all y-axis labels include their condition.

Expected outputs:

```text
data/xenium/processed/vbct/cross_sample_dotplot_exports/filter_sensitivity/<sample>_r1.00_filter_sensitivity_dotplot_summary.csv
figures/dotplots/filter_sensitivity/<sample>_r1.00_filter_sensitivity_zscore_dotplot.png
```

After review, decide whether evidence supports condition-specific biological
annotations or explicit cross-condition cluster mappings. Such mappings must
be based on the expression and marker evidence, not reused numeric IDs.

This VBCT task supersedes the PTMT marker-label curation task above as the
active next dotplot task.

## Session Update: Prefiltered Baseline in VBCT Sensitivity Plots (2026-10-24)

The active GR and MG filter-sensitivity export and plotting configs now include
five conditions at resolution `1.00`, ordered as `prefiltered`, `min10_area99`,
`min10_area95`, `min15_area99`, and `min15_area95`. This supersedes the four-only
ordering and validation requirement above; retain the four QC branches unchanged.

`prefiltered` means the original script 00 clean expression object and original
BANKSY clusters before the additional script 05 QC filtering and script 06
reclustering. It still includes script 00's initial preprocessing. Use the
original script 00 config's H5AD path and grouping column, not a recluster object
with original labels restricted to QC-passing cells.

For GR and MG, the local five-condition summary CSVs reuse the resolution 1.00
rows from `cross_sample_dotplot_exports/all_genes_all_res/<sample>_all_genes_dotplot_summary.csv`.
Both baseline exports contain all 304 genes and report `expression_source=raw`
with the expected clean script 00 H5AD provenance. Only condition/label metadata
is changed; expression values, percentages, and cell counts are preserved.
GR has 18 prefiltered clusters; MG has 24. These are existing export summaries,
not fresh exports from H5ADs verified on HPC in this session.

Script 03 configs now export all five conditions on a future HPC rerun and
explicitly disable automatic QC annotations to retain branch-specific numeric
cluster labels. Script 04 computes shared gene-wise z-scores over all five
conditions (clipped at 2.5), so colours change relative to the four-condition plots.
Cluster IDs remain condition-specific, including the prefiltered baseline.

## Future Work: Combined Selected-Resolution Dot Plot Across Samples (noted 2026-11-23)

Not started; options to discuss with Nathalie before implementing.

### Already done

Per-sample cluster marker dot plots for annotation set v1 exist (one per VBCT
sample). Each shows only that sample's clusters at its selected resolution
(min15_area95 branch), labelled `cluster <id> · <tier-1 label> (n=<cells>)`
and grouped by label. They use the 89 genes kept by
`03_xenium_gene_review_for_dotplot_manual_annotation_2026-06-11.csv`, z-scored
per gene within the sample, clipped at ±2.5.

- Configs: `config/04_plot_dotplot/active/vbct_annotation_sets/min15_area95_selected_resolutions/`
  (README there).
- Outputs: `data/xenium/output/vbct/cell_type_annotation/min15_area95_selected_resolutions/v1/tier_1/cluster_marker_dotplots/`.
- Mechanism: the `annotation` block in script 04 (`apply_annotation_labels`).
- Inputs: the `min15_area95` rows of the local filter-sensitivity exports, verified
  identical to v1's clusters and cell counts. The same configs replot for another
  set with `--annotation <set dir>` (drafts allowed).
- Context: `cell_type_annotation.md`.

### Goal

One figure covering all eight samples, each at only the resolution selected for
script 08.

### Options

**(a) All 127 clusters, stacked by sample.**
- **Layout.** One block per sample, rows labelled
  `<sample> · cluster <id> · <label>`, grouped by label within each sample and
  separated by sample.
- **Pros and cons.** Shows every cluster in one place, but the figure is tall
  (127 rows).
- **Decision needed: z-score scope.**
  - *Within each sample:* matches the per-sample plots and keeps contrast
    between clusters, but colours are not comparable between samples.
  - *Across the whole cohort:* colours are comparable between samples, but
    samples with high overall expression dominate.

**(b) One row per label per sample (cross-sample consistency check).**
- **Layout.** Rows like `T cell · CK_skin_res`, `T cell · EL_skin_res`, ... (about 50
  rows), against the same marker genes.
- **Purpose.** Shows whether a label means the same thing in every sample. For
  example: does every `T cell` group express CD3E/CD2, or is one sample's
  `Fibroblast` pericyte-like (RGS5)? Do labels used in only one or two samples
  (`Epithelial cell`, `Pericyte`) have convincing markers? Inconsistent labels
  would distort script 10's responder vs non-responder comparisons.
- **Method.** Each label's values must be recomputed from its clusters weighted
  by cell count (mean expression and percent expressing weighted by `n_cells`),
  not averaged equally across clusters.
  - Can be built from the same exports and the set's mapping, so no H5ADs are
    needed.
  - Alternatively, the poster marker dot plot helper
    (`helper_scripts/poster/plot_cohort_marker_dotplot_from_config.py`) computes
    per-label values from the annotated H5ADs, but it averages samples into one
    cohort row per label. It would need a per-sample row option.
- **Comparing versions.** Rerun after v2 to see whether the labels became more
  consistent.

### Implementation notes (either option)

- Extend script 04's `annotation` block to accept several samples (one block per
  sample, or the whole cohort), rather than one `sample`.
- Keep the check that every exported cluster is mapped and every mapped cluster
  is exported, for each sample.
- Keep the review filter, gene grouping, and canonical-gene highlighting from the
  per-sample configs unless Nathalie chooses otherwise.
- Write outputs next to the per-sample plots under the set's `<version>/<tier>/`
  folder.

### Observations from the v1 per-sample plots (to judge, not conclusions)

Some very small clusters stand out and may be worth deciding on in v2:

- CK_bowel_res clusters 16 (n=376) and 17 (n=238), labelled melanoma, each with
  one strongly expressed atypical gene.
- CK_bowel_res cluster 18, a 28-cell T-cell cluster.
- GR_lung_non_res_roi cluster 13 (n=403), labelled melanoma, with a B-cell-like
  marker signal.
