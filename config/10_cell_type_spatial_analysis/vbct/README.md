# Spatial analysis using reviewed cell types

Since 2026-11-22 this configuration names an annotation set (`"annotation"`)
and a tier (`"label_level"`, default `tier_1`); see `config/annotations/README.md`.
Sample names, branches, resolutions, grouping columns, the mapping, annotated
object paths, coordinate keys, and palette come from the set and its cohort.
Script 10 reads the set's script 08 objects, checks their annotation stamp,
groups spatial analyses by the tier's column (e.g. `cell_type_tier_1`), and
processes one sample at a time. It needs a released set (`--allow-draft` for
local tests); `--annotation` and `--label-level` override the config.
Original annotation objects remain unchanged.

From the Banksy_py project root:

```bash
# Read-only metadata validation for all eight annotated objects.
python 10_cell_type_spatial_analysis.py --config config/10_cell_type_spatial_analysis/vbct/min15_area95_selected_resolutions.json --validate-only

# Run one sample; repeat --sample to select additional samples.
python 10_cell_type_spatial_analysis.py --config config/10_cell_type_spatial_analysis/vbct/min15_area95_selected_resolutions.json --sample CK_skin_res
```

Without `--sample`, a full run processes all eight sequentially. Use the HPC
workflow for full spatial calculations. Annotated objects must be available
there at the set's derived paths (written by script 08 on the HPC). `--input-root` overrides
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
# another version or tier with the same config:
sbatch --export=ALL,ANNOTATION=config/annotations/vbct/min15_area95_selected_resolutions/<version>,LABEL_LEVEL=tier_2 \
  hpc/slurm/run_10_cell_type_spatial_analysis.sl
```

A tier-2 run repeats Moran's I, which does not depend on labels (open point in
`context/tasks/cell_type_annotation/cell_type_annotation.md`).

The legacy config-per-task wrapper and cluster-mode configs were archived on
2026-11-20 (`archive/2026-11-20_script10_legacy_cluster_mode/`); this folder was
renamed from `config/10_squidpy_recluster_analysis/vbct/cell_type_annotation/`.
See the workflow README for the full 08 → 10 sequence.

Figures are written under
`data/xenium/output/<project>/spatial_analysis/<cohort>/<version>/<tier>/<sample>/<branch>_r<resolution>/squidpy/`.
Under `data/xenium/processed/<project>/spatial_analysis/<cohort>/<version>/<tier>/<sample>/`
go the result H5AD, `<stem>_nhood_enrichment.csv`, and the centrality and Moran
CSVs in `<branch>_r<resolution>/`. Unversioned folders
(`.../spatial_analysis/min15_area95_selected_resolutions/<sample>/`) hold results
from before 2026-11-22 (v1 labels, `cell_type` column); earlier figures are in
`figures/_pre_post08_cleanup_2026-11-20/spatial_analysis/`.
Re-running a full analysis replaces its generated results using existing script
10 behavior. The annotated source objects are not overwritten.

Validation checks the annotation stamp, unique cell IDs, the tier's vocabulary,
finite coordinates, and agreement with the set's mapping for the exact original cluster column,
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
`<result folder>/<sample>/<stem>_nhood_enrichment.csv` (edge counts, z-scores, cell
counts; format in `xenium_workflow/spatial_tables.py`) beside each result H5AD. After an HPC run,
refresh the local copies with `10a_export_nhood_enrichment_tables.py`.

The cross-sample helpers that read these tables (neighbourhood enrichment
comparison, immune infiltration) are documented in
`config/spatial_analysis/README.md`.
