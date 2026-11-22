# Cross-sample spatial helpers

Local helpers in `helper_scripts/spatial_analysis/` compare script 10's
cell-type spatial results across the cohort (`res` vs `non_res`). Configs in
`vbct/` moved here from `config/10_squidpy_recluster_analysis/vbct/cell_type_annotation/`
on 2026-11-20 (post-08 cleanup, phase 3). Since 2026-11-22 each config names an
annotation set (`"annotation"`), a tier (`"label_level"`, default `tier_1`) and an
`output_subdir`; samples, response, tissue, palette, and coordinate keys come
from the set and its cohort (`config/annotations/README.md`). They list their
own sample order and take `--input-root`/`--output-root` (default: repository
root), `--annotation`/`--label-level` overrides, and `--allow-draft`. Cell-type
names in these configs (focus pairs, tumour and query types) must be labels of
the chosen tier. Figures go to
`data/xenium/output/<project>/spatial_analysis/<cohort>/<version>/<tier>/<output_subdir>/`,
tables to the same path under `data/xenium/processed/`. Analysis status and
findings: `context/tasks/spatial_analysis.md`.

## Neighbourhood enrichment comparison

Script 10 seeds neighbourhood enrichment (1234); earlier HPC results were
unseeded and differ by permutation noise only. Since 2026-11-20, script 10 also writes a small per-sample table,
`<script 10 result folder>/<sample>/<stem>_nhood_enrichment.csv` (edge counts, z-scores,
cell counts; format in `xenium_workflow/spatial_tables.py`), beside each result H5AD.
`helper_scripts/spatial_analysis/compare_nhood_enrichment_from_config.py` reads
those tables and compares samples using log2(observed / expected) neighbour
edges, which, unlike Squidpy z-scores, does not scale with cell number. It no
longer reads the result H5ADs or needs the Phoenix mount, and the local
recompute fallback (`--recompute-missing`, `result_root`) was removed.

For results written before that change, or after a new HPC script 10 run,
bring the tables to this repository with the backfill script. It reads only the
enrichment and the tier's label column from each result H5AD and verifies the
table round trip exactly:

```bash
python 10a_export_nhood_enrichment_tables.py --config config/10_cell_type_spatial_analysis/vbct/min15_area95_selected_resolutions.json --input-root ~/mnt/phoenix_scratch/Banksy_py
python helper_scripts/spatial_analysis/compare_nhood_enrichment_from_config.py --config config/spatial_analysis/vbct/nhood_enrichment_comparison.json --validate-only
python helper_scripts/spatial_analysis/compare_nhood_enrichment_from_config.py --config config/spatial_analysis/vbct/nhood_enrichment_comparison.json
```

Existing tables are kept unless `--overwrite` is passed to script 10a.

Response and tissue labels come from the cohort config; sample order, focus
pairs, and matched-tissue comparisons are set in this JSON. The enrichment
tables are found in script 10's derived result folder for the same set and tier
under `--input-root`. Difference panels are descriptive (1–4 samples per group).

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
python helper_scripts/spatial_analysis/plot_immune_infiltration_from_config.py --config config/spatial_analysis/vbct/immune_infiltration.json
```

Use `--skip-maps` to omit the per-sample spatial maps.
