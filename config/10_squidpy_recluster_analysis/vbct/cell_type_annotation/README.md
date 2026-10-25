# Spatial analysis using reviewed cell types

This configuration references script 08's configuration directly. Sample names,
branches, resolutions, original grouping columns, reviewed mapping, annotated
output paths, and shared palette are taken from that source of truth.
Script 10 reads the `output_h5ad` objects written by script 08, not its unannotated
HPC inputs. It groups spatial analyses by `cell_type` and processes one sample
at a time. Original annotation objects remain unchanged.

From the Banksy_py project root:

```bash
# Read-only metadata validation for all eight annotated objects.
python 10_squidpy_recluster_spatial_analysis.py --config config/10_squidpy_recluster_analysis/vbct/cell_type_annotation/min15_area95_selected_resolutions.json --validate-only

# Run one sample; repeat --sample to select additional samples.
python 10_squidpy_recluster_spatial_analysis.py --config config/10_squidpy_recluster_analysis/vbct/cell_type_annotation/min15_area95_selected_resolutions.json --sample CK_skin_res
```

Without `--sample`, a full run processes all eight sequentially. Use the HPC
workflow for full spatial calculations. Annotated objects must be available
there at script 08's relative output paths. `--annotation-root` overrides their
root directory (including when script 08 previously used --output-root).
Run script 08 first with `hpc/slurm/run_08_cell_type_annotation.sl`; it creates
the annotated H5AD inputs in their configured HPC-relative paths. The adjacent
`hpc_staging_files.txt` lists the code and configuration files needed for the
08 → 10 workflow; generated annotated H5ADs are no longer listed as staging
inputs.

The dedicated `hpc/slurm/run_10_cell_type_spatial_analysis.sl` wrapper uses one
task per sample (`0-7%2`) and selects samples from script 08's config. Submit it
yourself after script 08 has completed on HPC (or annotated files are transferred):

```bash
mkdir -p logs
sbatch hpc/slurm/run_10_cell_type_spatial_analysis.sl
```

Do not use the legacy config-per-task wrapper for this batch config. If using
that legacy wrapper with this directory, only one task is appropriate because
the config itself contains all eight samples. See the workflow README for the
script 08 HPC input-root override and full 08 → 10 sequence.

Figures and CSVs are written under
`figures/spatial_analysis/vbct_min15_area95_selected_resolutions/<sample>/<branch>_r<resolution>/squidpy/`.
Spatial H5AD results are written separately under
`data/xenium/processed/vbct/spatial_analysis/min15_area95_selected_resolutions/<sample>/`.
Re-running a full analysis replaces its generated results using existing script
10 behavior. The annotated source objects are not overwritten.

Validation checks unique cell IDs, cell-type vocabulary, finite coordinates,
and agreement with the reviewed mapping for the exact original cluster column,
branch, and resolution. Absent cell types are excluded from calculations while
present types retain canonical palette order and colours.

The scientific analyses are unchanged and all four still run. Moran settings
retain the earlier lightweight values: 25% subsampling, 10 permutations, FDR 0.10,
and five gene plots. These settings have not been upgraded for final inference.
All present cell types receive co-occurrence plots. Legacy single-sample cluster
configs still work and now also accept --validate-only.

Validation completed on all eight objects (2,083,396 cells). A CK dispatch check
verified grouping, palette, output separation, and provenance with scientific
calculations and H5AD writes mocked. No full spatial analysis has been launched.
