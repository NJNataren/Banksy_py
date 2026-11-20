# Spatial analysis using reviewed cell types

This configuration names the cohort config
(`config/cohorts/vbct/min15_area95_selected_resolutions.json`). Sample names,
branches, resolutions, original grouping columns, reviewed mapping, annotated
object paths, coordinate keys, and shared palette are taken from it.
Script 10 reads the `annotated_h5ad` objects written by script 08, not its unannotated
HPC inputs. It groups spatial analyses by `cell_type` and processes one sample
at a time. Original annotation objects remain unchanged.

From the Banksy_py project root:

```bash
# Read-only metadata validation for all eight annotated objects.
python 10_cell_type_spatial_analysis.py --config config/10_cell_type_spatial_analysis/vbct/min15_area95_selected_resolutions.json --validate-only

# Run one sample; repeat --sample to select additional samples.
python 10_cell_type_spatial_analysis.py --config config/10_cell_type_spatial_analysis/vbct/min15_area95_selected_resolutions.json --sample CK_skin_res
```

Without `--sample`, a full run processes all eight sequentially. Use the HPC
workflow for full spatial calculations. Annotated objects must be available
there at the cohort's relative `annotated_h5ad` paths. `--input-root` overrides
their root and `--output-root` the root for figures and result H5ADs; both
default to the repository root (`--annotation-root` was replaced on 2026-11-19).
Run script 08 first with `hpc/slurm/run_08_cell_type_annotation.sl`; it creates
the annotated H5AD inputs in their configured HPC-relative paths. The adjacent
`hpc_staging_files.txt` lists the code and configuration files needed for the
08 → 10 workflow; generated annotated H5ADs are no longer listed as staging
inputs.

The dedicated `hpc/slurm/run_10_cell_type_spatial_analysis.sl` wrapper uses one
task per sample (`0-7%2`) and selects samples in cohort order. Each task
requests 8 CPUs, 180 GB RAM, and 24 hours. Submit it
yourself after script 08 has completed on HPC (or annotated files are transferred):

```bash
mkdir -p logs
sbatch hpc/slurm/run_10_cell_type_spatial_analysis.sl
```

The legacy config-per-task wrapper and cluster-mode configs were archived on
2026-11-20 (`archive/2026-11-20_script10_legacy_cluster_mode/`); this folder was
renamed from `config/10_squidpy_recluster_analysis/vbct/cell_type_annotation/`.
See the workflow README for the full 08 → 10 sequence.

Figures are written under
`data/xenium/output/vbct/spatial_analysis/min15_area95_selected_resolutions/<sample>/<branch>_r<resolution>/squidpy/`
(moved from `figures/spatial_analysis/vbct_min15_area95_selected_resolutions/`
on 2026-11-19; earlier runs remain at the old path). Under
`data/xenium/processed/vbct/spatial_analysis/min15_area95_selected_resolutions/<sample>/`
go the result H5AD, `<stem>_nhood_enrichment.csv`, and (since 2026-11-20) the
centrality and Moran CSVs in `<branch>_r<resolution>/`.
Re-running a full analysis replaces its generated results using existing script
10 behavior. The annotated source objects are not overwritten.

Validation checks unique cell IDs, cell-type vocabulary, finite coordinates,
and agreement with the reviewed mapping for the exact original cluster column,
branch, and resolution. Absent cell types are excluded from calculations while
present types retain canonical palette order and colours.

The scientific analyses are unchanged and all four still run. Moran settings
retain the earlier lightweight values: 25% subsampling, 10 permutations, FDR 0.10,
and five gene plots. These settings have not been upgraded for final inference.
Since 2026-11-20 the Moran permutations are seeded (1234), so simulated p-values
and the FDR gene list are reproducible; runs before then were unseeded.
All present cell types receive co-occurrence plots.

Validation completed on all eight objects (2,083,396 cells). A CK dispatch check
verified grouping, palette, output separation, and provenance with scientific
calculations and H5AD writes mocked.

The first full HPC array run on 2026-10-25 completed tasks 0, 2, and 3 (CK skin,
GR lung, and BE brain). Tasks 1, 4, and 6 reached top-Moran spatial plotting but
failed because Squidpy 1.6 forwarded the unsupported `show` keyword to
Matplotlib. Tasks 5 and 7 (MF skin and CK bowel) reached co-occurrence but hit
the original eight-hour limit. Script 10 now omits that keyword and the active
wrapper allows 24 hours. After updating the HPC checkout, rerun only failed task
indices with `sbatch --array=1,4-7%2 hpc/slurm/run_10_cell_type_spatial_analysis.sl`.
Existing annotated script 08 inputs remain unchanged.

## Enrichment tables and cross-sample helpers

Script 10 seeds neighbourhood enrichment (1234), so z-scores are reproducible;
earlier HPC results were unseeded and differ by permutation noise only (edge
counts are deterministic). Since 2026-11-20 it also writes
`<result_dir>/<sample>/<stem>_nhood_enrichment.csv` (edge counts, z-scores, cell
counts; format in `spatial_tables.py`) beside each result H5AD. After an HPC run,
refresh the local copies with `10a_export_nhood_enrichment_tables.py`.

The cross-sample helpers that read these tables (neighbourhood enrichment
comparison, immune infiltration) are documented in
`config/spatial_analysis/README.md`.
