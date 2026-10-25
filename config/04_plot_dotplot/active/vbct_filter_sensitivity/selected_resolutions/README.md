# Selected-resolution VBCT annotation review

One plot per biological sample, following the GR filter-sensitivity template.
Each plot includes `prefiltered`, `min10_area99`, `min10_area95`,
`min15_area99`, and `min15_area95`, in that order.

| Sample | Resolution |
|---|---|
| CK_skin_res | 1.00 |
| EL_skin_res | 0.80 |
| GR_lung_non_res_roi | 1.00 |
| BE_brain_non_res | 0.90 |
| HW_brain_res | 0.80 |
| MF_skin_non_res | 0.50 |
| MG_gastric_non_res | 0.70 |
| CK_bowel_res | 0.80 |

The paired export configs are in
`config/03_export_summary/active/vbct_filter_sensitivity/selected_resolutions/`.
This directory contains exactly eight JSON configs, avoiding the additional
MG resolution 1.00 comparison in the parent directory.

## Export on HPC

From the HPC repository root, after syncing the configs:

```bash
mkdir -p logs
sbatch --array=0-7%3 \
  --export=ALL,CONFIG_DIR=config/03_export_summary/active/vbct_filter_sensitivity/selected_resolutions \
  hpc/slurm/run_03_xenium_dotplot_export_from_clean_script00_config.sl
```

Array order is alphabetical: BE, CK bowel, CK skin, EL, GR, HW, MF, MG.
Required source objects must exist on HPC; their current HPC availability has
not been verified. The baseline for HW, MF, and CK bowel follows the established
clean normalized-expression plus external cluster-label CSV workflow. Other
baselines use the existing clean clustered script 00 objects. All four filtered
conditions use the corresponding script 06 `cells_used_for_banksy` objects.
MF's baseline uses 30 PCs and its filtered branches use 20 PCs, as specified
by the source configs. Automatic QC annotation labels are disabled.

## Plot locally

Copy the exported summaries from
`data/xenium/processed/vbct/cross_sample_dotplot_exports/filter_sensitivity/`
to the same local relative directory, then run from the local repository root:

```bash
for cfg in config/04_plot_dotplot/active/vbct_filter_sensitivity/selected_resolutions/*.json; do
  conda run -n banksy python 04_plot_multi_sample_dotplot_from_config_local.py \
    --config "$cfg" || break
done
```

Outputs are under `figures/dotplots/filter_sensitivity/`, named by sample and
resolution. GR 1.00 and MG 0.70 reuse their existing output names; the plotting
script archives earlier images automatically.

All genes are exported from `.raw`; plotting retains the template's manual
gene review, canonical-gene highlighting, Tier 1 ordering, and figure settings.
Thus the final plot is the reviewed gene set, not literally all exported genes.
Expression is z-scored per gene across all five conditions within each sample,
clipped at 2.5. Colours are not a shared z-score across biological samples.
Match each condition to its own spatial plot using the exact grouping column
in the export config; matching numeric cluster IDs across conditions does not
imply matching cell types.

As of 2026-10-24, GR 1.00 and MG 0.70 summaries were available locally and their
plots were rendered. The other six summaries require HPC export and transfer.


## Local export through mounted HPC storage

Local copies of the five pending sample export configs are in
`config/03_export_summary/local_helpers/vbct_filter_sensitivity/selected_resolutions/`.
Their H5AD and baseline label inputs use absolute paths under
`/home/nnataren/mnt/phoenix_scratch/Banksy_py`; summary outputs remain local.
The mount must be connected. MF is excluded until its reclustering is ready;
GR and MG already have summaries. Run sequentially from the local repository:

```bash
for cfg in config/03_export_summary/local_helpers/vbct_filter_sensitivity/selected_resolutions/*.json; do
  conda run -n banksy python 03_export_dotplot_data_from_clean_script00_config.py \
    --config "$cfg" || break
done
```

Use the corresponding script 04 config after each successful export. Local
execution logs for the initial run are under `logs/local_filter_sensitivity/`.


## Local Export Completed (2026-10-24)

The five remaining available samples were exported and plotted locally using
HPC inputs under `/home/nnataren/mnt/phoenix_scratch/Banksy_py`:
CK_skin_res 1.00, BE_brain_non_res 0.90, EL_skin_res 0.80,
HW_brain_res 0.80, and CK_bowel_res 0.80. The sequential batch completed
successfully. Local input configs are under
`config/03_export_summary/local_helpers/vbct_filter_sensitivity/selected_resolutions/`.
All five summary CSVs passed condition, resolution, grouping, provenance,
304-gene coverage, duplicate-row, expression-percentage, and cell-count checks.
All five PNG files passed image-integrity checks; visual layout inspection was
blocked by the image tool's sandbox startup error.

Together with GR 1.00 and MG 0.70, seven samples now have local summaries and
plots in the standard `filter_sensitivity` output directories. MF 0.50 remains
pending because its new clustering results are not yet ready. This supersedes
the earlier statement that the five samples require HPC export and transfer.
No biological annotations were applied. Execution logs are under
`logs/local_filter_sensitivity/`.
