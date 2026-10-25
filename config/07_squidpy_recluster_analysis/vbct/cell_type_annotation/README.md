# Spatial analysis using reviewed cell types

This configuration references script 08's configuration directly. Sample names,
branches, resolutions, original grouping columns, reviewed mapping, annotated
output paths, and shared palette are taken from that source of truth.
Script 07 reads the `output_h5ad` objects written by script 08, not its unannotated
HPC inputs. It groups spatial analyses by `cell_type` and processes one sample
at a time. Original annotation objects remain unchanged.

From the Banksy_py project root:

```bash
# Read-only metadata validation for all eight annotated objects.
python 07_squidpy_recluster_spatial_analysis.py --config config/07_squidpy_recluster_analysis/vbct/cell_type_annotation/min15_area95_selected_resolutions.json --validate-only

# Run one sample; repeat --sample to select additional samples.
python 07_squidpy_recluster_spatial_analysis.py --config config/07_squidpy_recluster_analysis/vbct/cell_type_annotation/min15_area95_selected_resolutions.json --sample CK_skin_res
```

Without `--sample`, a full run processes all eight sequentially. Use the HPC
workflow for full spatial calculations. Annotated objects must be available
there at script 08's relative output paths. `--annotation-root` overrides their
root directory (including when script 08 previously used --output-root).
The supplied config uses one shared batch definition, so an HPC array pointed at
this directory must have just one task, not eight; per-sample job submissions can
instead call Python with --sample.

Figures and CSVs are written under
`figures/spatial_analysis/vbct_min15_area95_selected_resolutions/<sample>/<branch>_r<resolution>/squidpy/`.
Spatial H5AD results are written separately under
`data/xenium/processed/vbct/spatial_analysis/min15_area95_selected_resolutions/<sample>/`.
Re-running a full analysis replaces its generated results using existing script
07 behavior. The annotated source objects are not overwritten.

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
