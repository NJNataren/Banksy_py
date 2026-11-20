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
python 10_squidpy_recluster_spatial_analysis.py --config config/10_squidpy_recluster_analysis/vbct/cell_type_annotation/min15_area95_selected_resolutions.json --validate-only

# Run one sample; repeat --sample to select additional samples.
python 10_squidpy_recluster_spatial_analysis.py --config config/10_squidpy_recluster_analysis/vbct/cell_type_annotation/min15_area95_selected_resolutions.json --sample CK_skin_res
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

Do not use the legacy config-per-task wrapper for this batch config. If using
that legacy wrapper with this directory, only one task is appropriate because
the config itself contains all eight samples. See the workflow README for the
script 08 HPC input-root override and full 08 → 10 sequence.

Figures and CSVs are written under
`data/xenium/output/vbct/spatial_analysis/min15_area95_selected_resolutions/<sample>/<branch>_r<resolution>/squidpy/`
(moved from `figures/spatial_analysis/vbct_min15_area95_selected_resolutions/`
on 2026-11-19; earlier runs remain at the old path).
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
calculations and H5AD writes mocked.

The first full HPC array run on 2026-10-25 completed tasks 0, 2, and 3 (CK skin,
GR lung, and BE brain). Tasks 1, 4, and 6 reached top-Moran spatial plotting but
failed because Squidpy 1.6 forwarded the unsupported `show` keyword to
Matplotlib. Tasks 5 and 7 (MF skin and CK bowel) reached co-occurrence but hit
the original eight-hour limit. Script 10 now omits that keyword and the active
wrapper allows 24 hours. After updating the HPC checkout, rerun only failed task
indices with `sbatch --array=1,4-7%2 hpc/slurm/run_10_cell_type_spatial_analysis.sl`.
Existing annotated script 08 inputs remain unchanged.

## Cross-sample neighbourhood enrichment comparison

Script 10 now fixes the `nhood_enrichment` permutation seed (1234), so z-scores
are reproducible; earlier HPC results were unseeded and differ by permutation
noise only (edge counts are deterministic).

Since 2026-11-20, script 10 also writes a small per-sample table,
`<result_dir>/<sample>/<stem>_nhood_enrichment.csv` (edge counts, z-scores,
cell counts; format in `spatial_tables.py`), beside each result H5AD.
`helper_scripts/spatial_analysis/compare_nhood_enrichment_from_config.py` reads
those tables and compares samples using log2(observed / expected) neighbour
edges, which, unlike Squidpy z-scores, does not scale with cell number. It no
longer reads the result H5ADs or needs the Phoenix mount, and the local
recompute fallback (`--recompute-missing`, `result_root`) was removed.

For results written before that change, or after a new HPC script 10 run,
bring the tables to this repository with the backfill script. It reads only the
enrichment and `cell_type` from each result H5AD and verifies the table round
trip exactly:

```bash
python 10a_export_nhood_enrichment_tables.py --config config/10_squidpy_recluster_analysis/vbct/cell_type_annotation/min15_area95_selected_resolutions.json --input-root ~/mnt/phoenix_scratch/Banksy_py
python helper_scripts/spatial_analysis/compare_nhood_enrichment_from_config.py --config config/10_squidpy_recluster_analysis/vbct/cell_type_annotation/nhood_enrichment_comparison.json --validate-only
python helper_scripts/spatial_analysis/compare_nhood_enrichment_from_config.py --config config/10_squidpy_recluster_analysis/vbct/cell_type_annotation/nhood_enrichment_comparison.json
```

Existing tables are kept unless `--overwrite` is passed to script 10a.

Response and tissue labels come from the cohort config; sample order, focus
pairs, and matched-tissue comparisons are set in this JSON. `script10_result_dir`
under `--input-root` locates the enrichment tables. Difference panels are descriptive (1–4 samples per group).

Sensitivity options in the same config: `graph_settings` repeats edge counts
for wider neighbourhoods (`n_neighs` or `radius_um`; coordinates are µm), and
`tiling` computes focus-pair log2 O/E inside square tiles. Non-primary graphs
use KD-tree edge counts, verified edge-for-edge against Squidpy's 6-NN counts
for all eight samples; they have no z-scores. Comparison heatmaps use the
primary (script 10) graph; difference CSVs cover every graph.

## Immune infiltration distances

`helper_scripts/spatial_analysis/plot_immune_infiltration_from_config.py`
(config `immune_infiltration.json`) measures each T/B/macrophage cell's
distance to the nearest melanoma cell and reports the share within ≤20 µm,
20–100 µm, and >100 µm, against all non-tumour cells as a geometry reference.
It reads samples and response labels from the cohort config directly (sample
order is listed in its own JSON).
Both helpers import `helper_scripts/spatial_analysis/spatial_common.py`.

```bash
python helper_scripts/spatial_analysis/plot_immune_infiltration_from_config.py --config config/10_squidpy_recluster_analysis/vbct/cell_type_annotation/immune_infiltration.json
```

Use `--skip-maps` to omit the per-sample spatial maps.
