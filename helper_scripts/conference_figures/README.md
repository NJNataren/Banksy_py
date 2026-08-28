# Conference Figure Helpers

These helpers make quick, reproducible figure sets for talks and posters from
existing Xenium BANKSY outputs. They are intentionally outside the formal
numbered workflow: the goal is to support curated presentation figures without
changing the main analysis pipeline or overwriting biological data objects.

## Background

This directory was added while preparing conference figures from VBCT Xenium
samples. The first use case was a `scanpy` spatial embedding plot for
`CK_skin_res`, similar to:

```python
sc.pl.embedding(adata, basis="spatial", color="clusters")
```

The scope then expanded to presentation-ready figure packs for exemplar samples,
especially:

- `CK_bowel_res` at BANKSY resolution `1.10`
- `MG_gastric_non_res` at BANKSY resolution `0.70`

For the current poster use case, melanoma substructure is not the focus, so the
config supports collapsing melanoma-like cluster labels into a simpler
`Melanoma` category while preserving the detailed labels in memory.

Important workflow rule: use clean expression AnnData objects for biological
marker-expression plots. Do not use `adata_spatial_*.h5ad` BANKSY spatial objects
for marker expression, because their `.X` matrices contain BANKSY-expanded and
scaled neighbour features. These helpers are designed to read clean script 00
objects such as:

```text
data/xenium/processed/<project>/<sample>/adata_expression_clean_<sample>_with_banksy_clusters_<resolutions>.h5ad
```

## Main Script

Use the config-driven helper for new work:

```text
plot_spatial_conference_figures_from_config.py
```

It can generate, per sample:

- spatial cluster map
- spatial cluster map with labels placed on the tissue
- UMAP colored by the selected BANKSY cluster column, if a copied UMAP is present
- spatial marker-expression panels
- marker dotplot by cell-type labels when provided, otherwise by cluster ID
- cluster abundance bar plot
- spatial coordinate summaries and axis-scout plots for choosing manual closeups
- optional cropped closeup plots from config-defined `crop_windows`

The older CK-only script is kept as a simple starting example:

```text
plot_ck_skin_res_unfiltered_r1p10.py
```

Prefer the config-driven script for anything involving more than one sample or
manual label edits.

## Example Configs

Configs live under:

```text
config/conference_figures/
```

Current examples:

```text
config/conference_figures/ck_skin_res_unfiltered_r1p10.json
config/conference_figures/vbct_exemplar_ck_bowel_mg_gastric.json
```

The exemplar config currently includes:

- `CK_bowel_res`, `labels_scaled_gaussian_pc30_nc0.20_r1.10`
- `MG_gastric_non_res`, `labels_scaled_gaussian_pc35_nc0.20_r0.70`
- endothelial marker genes: `AQP1`, `CALCRL`, `CDH5`, `ECSCR`, `PLVAP`,
  `SELP`, `VWF`, `TFPI`

## Running Locally

From the repository root:

```bash
conda run -n banksy python helper_scripts/conference_figures/plot_spatial_conference_figures_from_config.py \
  --config config/conference_figures/vbct_exemplar_ck_bowel_mg_gastric.json
```

Outputs are written under the configured `output_dir`, usually:

```text
figures/conference/<sample_output_subdir>/
```

For example:

```text
figures/conference/CK_bowel_res_unfiltered_r1p10/
figures/conference/MG_gastric_non_res_unfiltered_r0p70/
```

Local plotting has successfully run for:

- `CK_bowel_res`: 683,656 cells, input `.h5ad` about 979 MB
- `MG_gastric_non_res`: 143,512 cells, input `.h5ad` about 213 MB

This is much lighter than rerunning BANKSY clustering, but large spatial plots
can still take several minutes to render.

## Running On HPC

Use the Slurm wrapper when local memory is tight or when the input `.h5ad` files
only exist on Phoenix:

```bash
cd /scratchdata1/users/a1210419/Banksy_py

sbatch --export=ALL,CONFIG=config/conference_figures/vbct_exemplar_ck_bowel_mg_gastric.json \
  hpc/slurm/run_conference_spatial_figures.sl
```

The wrapper activates the `banksy` conda environment and runs the same
config-driven helper.

## Config Fields

Top-level fields:

- `output_dir`: root folder for figure outputs.
- `dpi`: raster output resolution.
- `point_size`: default point size for spatial and UMAP plots.
- `marker_sets`: named marker panels to plot.
- `label_collapse_rules`: optional presentation-only relabeling rules.
- `write_coordinate_summary`: write min/max and percentile x/y ranges; defaults to `true`.
- `write_axis_scout_plots`: write full-tissue spatial plots with axes; defaults to `true`.
- `crop_windows`: optional global manual closeup windows.
- `crop_point_size`: point size for cropped spatial plots; defaults to `8.0`.
- `crop_show_axes`: show axes on cropped cluster/cell-type plots; defaults to `true`.
- `samples`: list of samples to process.

Each sample needs:

- `sample`: sample name used in output filenames.
- `resolution`: label used in plot filenames, for example `r1p10` or `r0p70`.
- `adata_path`: clean expression AnnData with BANKSY labels copied into `.obs`.
- `cluster_col`: selected cluster column in `.obs`.
- `output_subdir`: sample-specific output folder under `output_dir`.
- `crop_windows`: optional sample-specific manual closeup windows.

Optional label sources:

- `label_map`: inline mapping from cluster ID to display label.
- `label_file`: CSV mapping cluster IDs to labels.
- `label_obs_col`: existing `.obs` column to use as labels.
- `label_config`: JSON config containing a label map, usually `new_labels`.

For conference figures, prefer `label_map` because it is explicit and easy to
manually edit without piggy-backing on another workflow config.

Example inline label map:

```json
"label_map": {
    "0": "T_Cell_0",
    "1": "Melanoma_Cell_or_Melanocyte_1",
    "2": "Macrophage_2"
}
```

## Manual Label Tweaks

Edit the `label_map` directly in the conference config when you want nicer
poster labels. For example:

```json
"3": "Endothelial cells"
```

The helper stores detailed labels in memory as:

```text
conference_cell_type_detail
```

It then plots:

```text
conference_cell_type
```

That plotted label can be simplified by `label_collapse_rules`.

## Collapsing Melanoma Labels

The VBCT exemplar config currently collapses labels containing `Melanoma` or
`melanocyte` into one display category:

```json
"label_collapse_rules": [
    {
        "label": "Melanoma",
        "contains": [
            "Melanoma",
            "melanocyte"
        ]
    }
]
```

This is presentation-only. It does not modify the AnnData on disk or the source
QC labels. If you manually rename a cluster to a label that no longer contains
`Melanoma` or `melanocyte`, it will no longer be collapsed by this rule.

## Manual Closeups

The helper supports less-opinionated closeups through config-defined
`crop_windows`. The script does not choose endothelial regions for you. Instead,
it writes two scouting aids per sample:

```text
<sample>_spatial_coordinate_summary.csv
<sample>_spatial_clusters_<resolution>_axis_scout.png
<sample>_spatial_cell_type_labels_<resolution>_axis_scout.png
```

Use those outputs to choose a region, then add a crop window to either the
top-level config or one sample block. A crop can be specified with explicit
limits:

```json
"crop_windows": [
    {
        "name": "endothelial_zoom_1",
        "xlim": [2500, 3500],
        "ylim": [1500, 2300],
        "point_size": 10
    }
]
```

or with a center and dimensions:

```json
"crop_windows": [
    {
        "name": "endothelial_zoom_1",
        "center": [3000, 1900],
        "width": 1000,
        "height": 800,
        "point_size": 10
    }
]
```

For each crop, the helper writes closeup versions of the cluster plot, cell-type
label plot when labels are available, and marker spatial panels. Crop windows are
applied to `obsm["X_spatial"]`, which is copied from `obsm["spatial"]` or
`obsm["xy"]` when needed.

## Validation

Before running after edits:

```bash
python3 -m py_compile helper_scripts/conference_figures/plot_spatial_conference_figures_from_config.py
python3 -m json.tool config/conference_figures/vbct_exemplar_ck_bowel_mg_gastric.json
```

The full plotting run requires the `banksy` environment because it depends on
Scanpy and Matplotlib.
