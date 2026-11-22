# Task: Xenium QC Script Development

## Goal

Develop the Xenium QC workflow around `00_QC_testing_xenium_spatial.py` (since archived to `archive/`; the active QC script is `01_QC_xenium_spatial_clean_clustered.py`), preserving the existing plotting/output style while adding sample-aware QC metrics and masks that can be used before reclustering.

## Current Direction

- Prefer storing QC masks/flags in `.obs` so cells can be colored as pass/fail and excluded downstream without immediately dropping them.
- It is acceptable to also save a filtered AnnData object, but the main downstream pattern should be to use masks to exclude failing cells from clustering.
- For now, color QC plots by failing cells before filtering.
- Add onto the existing QC plotting/output style already present in the QC scripts rather than replacing the workflow wholesale.
- `CK_skin_res` is the preferred small local smoke-test sample because it is tiny.
- The expected analysis rhythm is: dirty clustering first to get a feel for the data, then clustering, then QC to mask cells, then reclustering with filters.
- Never overwrite anything in `data/xenium/raw_data/`.

## Current QC Script Comparison

- `00_QC_testing_xenium_spatial.py` (since archived to `archive/`) is the more standard config-driven QC script candidate. It reads a JSON config, loads `<dataset_name>_raw.h5ad`, adds `xy`, removes zero-count cells, downcasts and saves a float32 object, calculates Scanpy QC metrics, plots count/gene QC, applies minimum transcript and top-2-percent maximum transcript masks, colors spatial plots by pass/fail, exports inspection CSVs, and then continues into cluster/marker/GOI spatial analysis.
- `00_QC_testing_xenium_spatial_with_cell_subset_check.py` (no longer in the repository) is explicitly labelled as a local subset-inspection/debugging version, "Not for use on HPC." It contains extra commented blocks for checking individual cell IDs, comparing raw matrix/count slots, and investigating count discrepancies. It otherwise overlaps heavily with the standard QC script and also continues into clustering/marker/GOI plotting.
- Nathalie will decide which QC script should be treated as active/canonical after reviewing this comparison.

## Recent Notebook Work

Recent QC development has been happening in `01_QC_xenium_spatial_clean_clustered.ipynb`, using `CK_skin_res` as the local test sample.

- The active notebook workflow reads the clean clustered AnnData object from script 00 rather than re-reading raw `.h5ad` inputs directly.
- Local CK testing has focused on `labels_scaled_gaussian_pc30_nc0.20_r0.70` and updated manual cluster labels for nine clusters.
- The notebook and script were cleaned up to use `adata` consistently instead of mixed `adata_clean` / `adata_clean_clean` object names.
- The three helper plotting functions under the notebook `FUNCTIONS` section now have docstrings for quicker reuse and review.
- Negative-control probe QC has been explored from two angles:
  - aggregate per-cell `.obs` columns such as `control_probe_counts`;
  - transcript-level `neg_tx` data for identifying which specific `NegControlProbe_*` features contribute to the negative-control signal.
- Useful negative-control summaries include:
  - percentage of cells with at least one negative-control probe count;
  - global percentage contribution of each negative-control probe to all negative-control counts;
  - within-cell proportion of negative-control counts attributable to each probe.
- For plotting negative-control probe dominance, a global bar plot is the clearest first diagnostic. Violin plots can be useful for within-cell proportions, but `cut=0` and `plt.ylim(0, 1)` are needed because kernel smoothing can visually extend beyond the valid proportion range.
- Spatial and UMAP QC plots should use categorical pass/fail or descriptive columns rather than raw boolean columns. Scanpy can error when plotting boolean columns directly.
- `cluster_ann_col` is defined early as a Python variable name, but the actual `.obs` column only exists after the later cluster-annotation assignment step. This was important for debugging KeyErrors during plotting.
- Negative-control-positive cells are currently best treated as an inspection signal rather than automatic failures until a filtering threshold is chosen.
- The first area filter under development is a top-percentile maximum `cell_area` filter, stored as an `.obs` mask such as `max_area_threshold_99` plus a categorical plotting column.
- `plot_umap_qc_metric()` was revised so it builds output filenames internally from `dataset_name` and `qc_metric` rather than requiring a manual filename argument.
- The same function now accepts lightweight tags that can be appended to the saved filename and plot title, which made it easier to support looping over multiple clustering resolutions without making the call sites messy.
- The QC UMAP plotting workflow is now resolution-aware:
  - it loops across `res_label_list`;
  - it constructs the matching cluster column for each resolution;
  - it colors the cluster UMAP using that resolution-specific cluster column;
  - it saves files with the resolution in the filename;
  - it includes the resolution in the plot title.
- A valid small-sample QC config was prepared for `config/01_QC/vbct/CK_skin_res.json` so the clean clustered QC script can be smoke tested locally or on HPC.
- The `run_01_xenium_QC_array_vbct.sl` wrapper was updated so it points at the small QC config directory, changes into the repo before running, and calls `01_QC_xenium_spatial_clean_clustered.py`.

## Negative-Control Probe Plot Options

- Global contribution bar plot: percent of all negative-control counts contributed by each probe.
- Detection-rate bar plot: percent of cells with each probe detected at least once.
- Raw-count violin plot: per-cell distribution of counts for each probe.
- Cluster-by-probe heatmap: cluster-level mean counts or detection rates for each probe.
- Spatial burden plot: tissue plot of total negative-control probe counts or a pass/fail category.
- QC scatter plots: total negative-control burden versus other QC metrics such as total counts or cell area.
- Within-cell proportion violin plot: per-cell fraction of total negative-control signal attributable to each probe.

## Next Work
- In the next work window, develop the clustree workflow so the clustree script can be executed across multiple samples rather than one sample at a time.

## Session Update: PTMT PC55 QC and Clustree

PTMT PC55 QC configs live under:

```text
config/01_QC/ptmt_pc55/
```

They must mirror the full script 00 PC55 resolution list through `1.5`; otherwise script 01 tries to open the old 7-resolution clean clustered object and fails with `FileNotFoundError`.

Run all eight run 1 PTMT PC55 QC jobs on the HPC with:

```bash
cd /scratchdata1/users/a1210419/Banksy_py

sbatch --array=0-7%2 \
  --export=ALL,CONFIG_DIR=config/01_QC/ptmt_pc55 \
  run_01_xenium_QC_array_ptmt.sl
```

Array order is:

```text
0  10329_run_1_DMSO
1  10331_run_1_1644
2  10707_run_1_DMSO
3  10708_run_1_1644
4  10722_run_1_DMSO
5  10723_run_1_1644
6  10733_run_1_DMSO
7  10734_run_1_1644
```

Script 01 QC outputs are project-scoped:

```text
data/xenium/output/ptmt_pc55/<sample>/
```

A local PTMT PC55 clustree runner was added:

```text
run_01a_clustree_ptmt_pc55_local.sh
```

Before running clustree locally, copy script 00 cluster assignment CSVs from the HPC:

```bash
rsync -av --include='*/' --include='*_cell_cluster_id_across_clustering_res_*.csv' --exclude='*' \
  a1210419@phoenix:/scratchdata1/users/a1210419/Banksy_py/data/xenium/processed/ptmt_pc55/ \
  data/xenium/processed/ptmt_pc55/
```

Then run:

```bash
./run_01a_clustree_ptmt_pc55_local.sh
```

The runner uses prefix `labels_scaled_gaussian_pc55_nc0.20_r` and writes clustree outputs to:

```text
data/xenium/output/ptmt_pc55/<sample>/clustree_qc/
```

Annotated clustree outputs are only produced when the matching script 01 QC config contains both `cluster_col` and `new_labels`. Standard topology and SC3-stability clustree plots are produced without annotation labels.

