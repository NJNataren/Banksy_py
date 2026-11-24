# Task: Shared Cell-Type Annotation and Colours Across Samples

## Goal

Identify clusters representing the same biological cell type, group them under
shared cell-type labels, and apply one consistent colour palette across samples.
Use these annotations to begin comparing cell populations, expression, and
spatial distributions across VBCT samples.

## Path Update (2026-11-20, post-08 cleanup)

Superseded again on 2026-11-22: mapping, palette, and annotated-object paths
now live in versioned annotation sets; see "Design: Versioned Annotation Sets"
and "Session Update: Annotation Sets Implemented" at the end of this file.

Paths below are historical. On branch `cleanup/post08-flow` (see
`post08_flow_cleanup.md`; merge to `main` pending an HPC test):

- Samples, mapping CSV, and palette are defined once in
  `config/cohorts/vbct/min15_area95_selected_resolutions.json`; the script 08
  and 09 configs point to it. Script 08 takes `--input-root`/`--output-root`
  (default: repo root); locally it needs `--input-root ~/mnt/phoenix_scratch/Banksy_py`.
- Script 08 figures: `figures/cell_type_annotation/vbct_min15_area95_selected_resolutions/`
  -> `data/xenium/output/vbct/cell_type_annotation/min15_area95_selected_resolutions/`.
- Script 09 figures: `figures/cell_type_proportions/...` ->
  `data/xenium/output/vbct/cell_type_proportions/min15_area95_selected_resolutions/`;
  proportion CSV -> `data/xenium/processed/vbct/cell_type_proportions/...` (`table_dir`).
  The tissue-facet field is now `facet`.
- Script 10: `10_squidpy_recluster_spatial_analysis.py` -> `10_cell_type_spatial_analysis.py`;
  config folder `config/10_squidpy_recluster_analysis/vbct/cell_type_annotation/`
  -> `config/10_cell_type_spatial_analysis/vbct/`; `--annotation-root` ->
  `--input-root`; figures under `data/xenium/output/vbct/spatial_analysis/...`.
- Local copies of the old `figures/spatial_analysis/`, `figures/cell_type_annotation/`, and `figures/cell_type_proportions/` folders were moved to `figures/_pre_post08_cleanup_2026-11-20/` on 2026-11-21 (README there); copies on Phoenix are untouched. Script 08's titled UMAPs, combined panels, `per_sample_comparison/`,
  and `cell_type_composition/` exist only in the archive until script 08 is
  rerun; `clean_scanpy/` and script 09 were regenerated in the new locations.

## Current Status

Planning stage, 2026-10-24. No annotation mappings, merged cell-type objects, or
shared palette have been created for this task. The workflow below is the
proposed approach; biological labels require review before application.

Related context:

- `README_00_04_XENIUM_ANALYSIS_WORKFLOW_CONCISE.md`: pipeline and object conventions.
- `context/tasks/dotplot/dotplot.md`: marker-expression plots and filter-sensitivity comparisons.

GR and MG have filter-sensitivity plots at resolution `1.00`, including a
prefiltered baseline and four QC branches. MG also has a five-condition plot at
resolution `0.70`. These are candidate review inputs, not final selections of
branch or resolution for annotation. Check current files before relying on them.

## Decisions to Agree Before Annotation

1. Start with MG and GR, or include all VBCT samples?
2. Annotate one selected QC branch per sample, or all four branches plus the
   prefiltered baseline?
3. Which resolution should be used for each sample and branch? A shared numeric
   resolution does not guarantee equivalent biological granularity.
4. What broad cell types and finer subtypes or states should the shared
   vocabulary include?
5. Which plots and comparisons should be produced first after annotation?

Do not assume these decisions have been made from the existing dotplot settings.

## Biological Annotation and Review

Review clean biological marker expression together with spatial context and
available prior annotations. Record the evidence supporting each proposed label,
including ambiguous or conflicting markers. Existing annotations are supporting
context and must be checked against the actual clustering being annotated.

Clusters are identified within a particular sample, analysis branch, and
resolution. Numeric cluster IDs do not establish biological equivalence across
samples, branches, or resolutions. Never transfer annotations by matching cluster
numbers alone.

Maintain a shared vocabulary with separate broad cell-type and finer subtype or
state labels. This allows broad comparisons while preserving potentially
meaningful differences within a cell type. Do not force uncertain populations
into a known type; retain explicit unknown, ambiguous, or mixed labels for review.
Distinguish an unreviewed cluster from a reviewed but unresolved population.

## Group Clusters Through Cell-Level Annotation

The proposed meaning of "merge" is to assign the same cell-type label to cells
from multiple biologically matching clusters within an object. Preserve the
original cluster columns and cell identifiers. This step does not require
rerunning BANKSY or physically combining all sample objects.

Attach reviewed annotation columns to clean expression AnnData objects while
retaining sample, branch, resolution, and original-cluster provenance. Keep
existing expression values unchanged and write annotated outputs separately from
the source objects.

Possible annotation columns, subject to agreement:

- `cell_type`: shared broad biological identity.
- `cell_subtype`: finer identity or state, where supported.
- `annotation_status`: proposed, reviewed, or unresolved.

A broad cell type can include multiple original clusters while those clusters
remain available for subtype review and sensitivity comparisons.

## Mapping File Proposal

Use a reviewable table with one row per original cluster in each selected
sample/branch/grouping. Suggested fields:

| Field | Purpose |
|---|---|
| `sample` | Biological sample identifier |
| `branch` | Prefiltered baseline or specific QC branch |
| `resolution` | BANKSY resolution associated with the source labels |
| `adata_path` | Source clean expression object |
| `groupby` | Exact original `.obs` cluster column |
| `cluster_id` | Original cluster ID within this grouping |
| `cell_type` | Proposed or reviewed shared broad label |
| `cell_subtype` | Optional supported subtype or state |
| `annotation_status` | Review state |
| `evidence` | Markers, spatial observations, and relevant references |
| `notes` | Uncertainty and decisions requiring follow-up |

Choose the final table location and schema when implementation begins. Check for
duplicate or conflicting mappings and report unmapped clusters explicitly.
Apply only mappings matching the actual source object and grouping column.

## Shared Colour Palette

Maintain one project-level mapping from canonical cell-type labels to fixed
colours. Reuse it across every sample and condition, even when some cell types
are absent. Do not assign colours independently from category order in each plot.

Apply the palette to categorical cell-type displays, including spatial plots,
UMAPs, and annotation legends. Reserve distinguishable colours for unresolved
categories. Agree on subtype colours separately if subtype plots are needed.

For expression dotplots, retain the continuous colour scale for expression;
show cell-type colours through row annotations or labels rather than replacing
expression values with categorical colours.

## Summaries and Comparisons After Grouping

Recalculate cell-type expression summaries from the constituent cells in the
clean expression object. Do not average cluster-level means or percentages
equally: differently sized clusters would contribute incorrectly. Use clean
expression `.raw` or another explicitly verified expression source, never
BANKSY-expanded `adata_spatial.X`.

Keep summaries separate by biological sample and analysis condition. Specify the
cell inclusion rules and denominator for cell-type proportions. Treat QC branches
as alternative analyses of the same biological sample, not independent samples.

Use shared gene selection and plot settings for comparisons. Consistent labels
and colours aid comparison but do not alone resolve differences in sampling,
normalisation, tissue context, or annotation confidence.

## Proposed Outputs

- A reviewed cluster-to-cell-type mapping table with supporting evidence.
- A shared cell-type vocabulary and colour palette.
- Clean expression objects with added cell-type annotations and preserved
  original cluster labels.
- Spatial plots and UMAPs coloured consistently by cell type.
- Cell-type expression dotplots and per-sample cell counts/proportions, according
  to the agreed comparison scope.

Resolve output paths and naming before generating these artifacts. Preserve
source objects and distinguish outputs by sample, branch, and grouping.

## Validation

- Every selected source cluster has a mapping or an explicit unresolved status.
- Each source cluster maps to at most one broad cell type within its grouping.
- Adding annotations preserves cell IDs, cell counts, original labels, and
  expression values.
- Counts summed over cell types reconcile with the selected cell population;
  exclusions and denominators are recorded.
- The same cell-type label has the same colour in every plot and legend.
- No cross-sample or cross-branch mapping relies only on numeric cluster IDs.
- Expression summaries use the intended clean expression source and cells.

## Next Step

Agree on the starting samples, QC branch scope, and resolutions. Then review the
relevant marker-expression and spatial plots and prepare proposed cluster labels
for biological review before merging annotations or generating comparisons.

To load this task in a new chat:

> Use $banksy-context for cell_type_annotation.md.


## Session Update: Selected-Resolution Filter Review (2026-10-24)

The user selected all eight VBCT samples for initial marker/spatial review:
CK_skin_res 1.00, EL_skin_res 0.80, GR_lung_non_res_roi 1.00,
BE_brain_non_res 0.90, HW_brain_res 0.80, MF_skin_non_res 0.50,
MG_gastric_non_res 0.70, and CK_bowel_res 0.80. Compare the prefiltered baseline
and all four QC branches for each sample before selecting annotation mappings.
This resolves initial review scope; it does not approve a final filtering branch
or biological labels.

Eight paired configs reside under `selected_resolutions/` within both
`config/03_export_summary/active/vbct_filter_sensitivity/` and
`config/04_plot_dotplot/active/vbct_filter_sensitivity/`. They follow the GR
1.00 template, including gene review/highlighting and five-condition z-scoring.
The plotting directory README documents HPC and local commands.

All 16 configs were validated. GR 1.00 and MG 0.70 existing summaries passed
resolution/grouping/provenance checks and their plots were rendered. The other
six lack complete local branch inputs and summaries; export on HPC and transfer
their CSVs before rendering. Current HPC availability remains unverified.
HW, MF, and CK bowel baselines reuse clean normalized expression plus original
cluster-label CSVs from the existing helper workflow. MF baseline and recluster
labels use 30 and 20 PCs respectively. No cell-type annotations or shared palette
have been applied.


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


## Additional EL Resolution 1.50 Review

Generated `figures/dotplots/filter_sensitivity/EL_skin_res_r1.50_filter_sensitivity_zscore_dotplot.png`
from four filtered branches at resolution 1.50 using the existing exporter and
plotting settings. The prefiltered object has no 1.50 labels (its original
clustering stops at 1.10), so the baseline is omitted rather than mixing
resolutions. Branch cluster counts are min10_area99: 32, min10_area95: 33,
min15_area99: 32, min15_area95: 29. The summary includes all 304 genes;
the plot applies the existing manual gene review and four-condition z-scoring.
Summary provenance/resolution/value checks and PNG integrity checks passed.
Configs are in the parent `vbct_filter_sensitivity/` directories for active
script 03, local-helper script 03, and active script 04, named with
`EL_skin_res_r1.50`. The selected-resolution eight-sample batch is unchanged.


## Selected-Resolution Cell-Type Annotation and Scanpy UMAPs (2026-10-25)

Implemented `08_annotate_cell_types_and_plot_umaps.py` with a reviewed 127-row
cluster mapping, a shared eight-cell-type palette, and the min15/area95 selected
resolution config under `config/08_cell_type_annotation/vbct/`. All eight source
objects were annotated successfully (2,083,396 cells total), with one categorical
`cell_type` column and consistent `cell_type_colors`; per-cluster and per-cell-type
count tables reconcile with the source objects.

The active figure workflow now uses Scanpy only and writes PNG, SVG, and PDF at 300 DPI.
Dense UMAP points are rasterized in SVG/PDF while titles and on-data labels remain
editable text for Inkscape. The completed output contains eight per-sample figures
under `figures/cell_type_annotation/vbct_min15_area95_selected_resolutions/per_sample_scanpy/`
and one combined panel in the parent directory, for 27 files total. All PNGs,
SVGs, and PDFs passed format/integrity checks; every SVG contains editable text
and an embedded raster point layer. Previously generated Matplotlib comparison
PNGs remain as reference outputs.


## Handoff: Cell-Type UMAPs Complete; Proportion Plot Next (2026-10-25)

### Work completed in this session

The user supplied reviewed cluster-to-cell-type labels for all eight VBCT samples.
Label spelling was normalised to the canonical broad cell types in
`config/08_cell_type_annotation/vbct/shared_cell_type_palette.json`; for example,
`Melanoma`, `Melanoma cells`, and `Melanoma cell` are represented as
`Melanoma cell`. The canonical cell-type order and colours are:

1. Melanoma cell — `#AA4499`
2. Macrophage — `#EE7733`
3. T cell — `#4477AA`
4. B cell — `#66CCEE`
5. Fibroblast — `#228833`
6. Endothelial cell — `#CC3311`
7. Epithelial cell — `#CCBB44`
8. Pericyte — `#777777`

The chosen source branch is the clean, post-filtered min15/area95 population:
minimum 15 transcripts and grouped cell-area 95th-percentile filtering. The
selected BANKSY resolutions, in the intended sample order, are:

- CK_skin_res: 1.00
- EL_skin_res: 0.80
- GR_lung_non_res_roi: 1.00
- BE_brain_non_res: 0.90
- HW_brain_res: 0.80
- MF_skin_non_res: 0.50
- MG_gastric_non_res: 0.70
- CK_bowel_res: 0.80

`08_annotate_cell_types_and_plot_umaps.py` applies the reviewed mapping, validates
complete and unique cluster coverage, writes categorical `cell_type` metadata and
the shared `cell_type_colors`, produces reconciled count summaries, and renders
the saved BANKSY UMAP embeddings through Scanpy. Its active run configuration is
`config/08_cell_type_annotation/vbct/min15_area95_selected_resolutions.json` and
the reviewed mapping is
`config/08_cell_type_annotation/vbct/min15_area95_selected_resolutions_cell_type_mapping.csv`.
The mapping contains 127 source clusters across 2,083,396 cells. All eight
annotated H5AD objects were written successfully to the sample-specific
`data/xenium/processed/vbct/<sample>/cell_type_annotation/` paths recorded in the
run config. The latest plotting rerun used `--no-write-annotated`, so it replaced
figures without unnecessarily rewriting those validated H5AD objects.

The user preferred the Scanpy rendering over the direct Matplotlib comparison.
The active Scanpy outputs are under
`figures/cell_type_annotation/vbct_min15_area95_selected_resolutions/`:

- eight per-sample UMAPs in `per_sample_scanpy/`;
- one combined eight-panel figure named `combined_scanpy_cell_type_umaps`;
- PNG, SVG, and PDF for every figure, giving 27 active files.

The config now defaults to 300 DPI for poster use. Individual PNGs are
3287 × 2148 pixels and the combined PNG is 4160 × 5802 pixels. SVG/PDF exports
use hybrid rasterisation: the dense cell point layers are embedded at 300 DPI,
while titles and on-data cell-type labels remain editable vector text in
Inkscape. All 27 files passed format and integrity checks, PNG DPI checks, and
SVG checks for both editable text and embedded point layers. Scanpy emitted only
known deprecation/colour-mapping warnings. Older Matplotlib comparison PNGs were
left in place as reference outputs and are no longer the active figure style.

The reconciled summary tables are:

- `data/xenium/processed/vbct/cell_type_annotation/min15_area95_selected_resolutions/reviewed_cluster_to_cell_type_counts.csv`
- `data/xenium/processed/vbct/cell_type_annotation/min15_area95_selected_resolutions/reviewed_cell_type_counts.csv`

### Next step: stacked cell-proportion bar plot

Create a 100% stacked bar plot showing broad cell-type proportions for all eight
samples. Use `reviewed_cell_type_counts.csv` directly; it already contains
`sample`, `branch`, `resolution`, `cell_type`, and `n_cells`, so the annotated
H5AD files do not need to be loaded again. For each sample, calculate
`proportion = n_cells / sum(n_cells)` and verify that proportions sum to 1. Fill
cell types absent from a sample with zero before plotting.

Preserve the sample order and selected resolutions listed above. Preserve the
cell-type stacking order and exact colours from the shared palette (superseded:
`shared_cell_type_palette.json` moved to `config/annotations/vbct/min15_area95_selected_resolutions/v1/palette.json` on 2026-11-22) so the bar plot matches the UMAPs. A sensible x-axis label is the sample name,
with its selected resolution either on a second line or in the caption. Label
the y-axis `Cell proportion` and format it as 0–100%. Use one shared legend and
avoid labels inside very small segments unless specifically requested.

Implement this as a reusable, config-driven script rather than an ad hoc notebook,
following the existing numbered workflow (the next script can be
`09_plot_cell_type_proportions.py`). Write PNG, SVG, and PDF at 300 DPI to a new
cell-type-composition subdirectory without overwriting the UMAPs. Keep SVG/PDF
text editable for Inkscape. Validate the input total of 2,083,396 cells, eight
unique samples, per-sample proportion sums, the expected cell-type vocabulary,
output file integrity, and 300-DPI PNG metadata.


## Cell-Proportion Plot Implemented (2026-10-25)

Implemented the separate, config-driven `09_plot_cell_type_proportions.py` with
`config/09_cell_type_proportions/vbct/min15_area95_selected_resolutions.json`.
This completes the stacked-bar next step above. The script reads the reviewed
count table, preserves sample/resolution order and the exact shared palette,
fills absent cell types with zero, and uses each sample's annotated post-filtered
cell total as its denominator. No H5AD objects are loaded or changed.

Outputs are `cell_type_proportions.{png,svg,pdf,csv}` under
`figures/cell_type_annotation/vbct_min15_area95_selected_resolutions/cell_type_composition/`.
Checks passed for eight samples, 2,083,396 cells, 64 zero-filled sample/type rows,
per-sample proportions summing to one, syntax validity, PNG integrity and 300-DPI
metadata (4500 x 2400 pixels), editable SVG text, and PDF header/end markers.
Visual layout inspection remains pending because the image viewer failed to
initialize its filesystem sandbox; the PDF was not visually rendered.
Usage is documented in `config/09_cell_type_proportions/README.md`.


## Dedicated Proportion Figure Directory and Rerun (2026-10-25)

At the user’s request, script 09 was rerun with its own top-level figure directory:
`figures/cell_type_proportions/vbct_min15_area95_selected_resolutions/`.
This supersedes the previous nested UMAP-directory output path. PNG, SVG, PDF,
and CSV were generated successfully. Checks passed for PNG integrity, 4500 x
2400 dimensions and 300 DPI, editable SVG text, PDF signatures, 64 sample/type
rows, eight samples, 2,083,396 cells, and per-sample proportions summing to one.
Earlier outputs remain in their original location.


## Tissue-Faceted Proportions (2026-10-25)

Implemented user-approved Skin, Brain, Lung, and Gastrointestinal facets in
script 09. The new config is
`config/09_cell_type_proportions/vbct/min15_area95_selected_resolutions_by_tissue.json`.
Facets have 3, 2, 1, and 2 bars with proportional panel widths and equal bar widths.
Gastric and bowel remain distinct samples within Gastrointestinal. Labels show
patient and response; gastrointestinal labels retain the anatomical site.
Shared scale, palette, legend, and per-sample denominators are preserved.
Outputs use `cell_type_proportions_by_tissue.{png,svg,pdf,csv}` in the dedicated
`figures/cell_type_proportions/vbct_min15_area95_selected_resolutions/` directory.
Original outputs are retained. Counts/proportions exactly match the original CSV.
Checks passed for panel membership, equal bar widths, axis-label canvas bounds,
PNG integrity (5100 x 2400 at 300 DPI), SVG text, and PDF signatures.


## Script 10 Aligned With Script 08 (2026-10-25)

Script 10 now accepts a configuration referencing script 08's selected-resolution
configuration. The entrypoint is now
`config/10_cell_type_spatial_analysis/vbct/min15_area95_selected_resolutions.json`
(moved there from the original script 10 config root on 2026-11-20).
It reads script 08's annotated output H5ADs, uses `cell_type`, and reuses the
reviewed mapping and shared palette. Sample/branch/resolution definitions are not
duplicated. Legacy cluster configs remain supported. CLI options include
`--sample`, `--annotation-root`, and read-only `--validate-only`.

All eight annotated inputs passed backed metadata validation (2,083,396 cells),
including mapping agreement and finite coordinates. A CK dispatch smoke test
passed with scientific analyses and output H5AD writes mocked. Full spatial
analyses were not run. Moran settings retain 25% subsampling and 10 permutations;
review these before final inference. All present types receive co-occurrence
plots. New figures/tables go under `figures/spatial_analysis/`; spatial H5ADs go
under `data/xenium/processed/vbct/spatial_analysis/`, preserving annotation inputs.
Usage and HPC batch considerations are in the new config directory README.


## Spatial Launch Requested; SSH Blocked (2026-10-25)

User authorized launch. Phoenix SSH authentication failed (Permission denied);
no jobs were submitted. Prepared `hpc/slurm/run_10_cell_type_spatial_analysis.sl`
with one task per sample in script 08 order (0-7%2; 180G; 8 hours).
The script 10 cell_type_annotation config directory now contains
`hpc_staging_files.txt`, listing code/configs and eight annotated H5ADs
(about 2.0 GiB). Checked mounted CK annotation directory was absent.
Superseded by the user’s HPC access rule below: do not retry SSH or submit jobs
from the agent. Nathalie handles any remote connection and job submission.


## HPC Access Rule (2026-10-25)

The user prohibits direct SSH connections to Phoenix/HPC, including connection
tests and remote job execution. Existing local HPC mounts may be inspected;
do not establish or reconnect mounts or use SSH-backed transfer tools. Prepare
scripts and commands for Nathalie to run herself. This supersedes prior launch
and SSH-retry instructions. The persistent rule is recorded in both workspace
and analysis-project AGENTS.md files.


## Spatial Step Renamed to 10 (2026-10-25)

At the user’s request, renamed the spatial script, config root, and both HPC
wrappers from 07 to 10; references and relevant READMEs were updated.
The active dependency is 06 → 08 → 10, with 09 optional after annotation.
Use `10_squidpy_recluster_spatial_analysis.py`,
`config/10_squidpy_recluster_analysis/`, and
`hpc/slurm/run_10_cell_type_spatial_analysis.sl`. No old-name alias is retained.
Superseded 2026-11-20: script 10 is now `10_cell_type_spatial_analysis.py` with
its config in `config/10_cell_type_spatial_analysis/` (see `post08_flow_cleanup.md`).
Existing annotated data and scientific settings are unchanged. HPC commands
remain user-run only; no HPC connection or job submission was attempted.


## Script 08 HPC Wrapper Added (2026-10-25)

Added `hpc/slurm/run_08_cell_type_annotation.sl` so annotation is performed in a
Slurm allocation rather than on an HPC login node. It runs all eight configured
samples sequentially in one job (8 CPUs, 180 GB, 8 hours), with the HPC project
directory passed as both input and output root. A Slurm array is intentionally
not used because script 08 writes shared count tables and a combined UMAP after
processing every sample. The wrapper refuses to replace existing outputs by
default; `OVERWRITE=1` is an explicit recovery/rerun option. Relevant workflow,
config, and Slurm READMEs were updated. No HPC connection or submission was made.


## Script 10 First HPC Run and Corrections (2026-10-26)

Read-only inspection of mounted logs for array job 16021129 confirmed that tasks
0, 2, and 3 completed (CK skin, GR lung, and BE brain). Tasks 1, 4, and 6 (EL
skin, HW brain, and MG gastric) completed the analyses through Moran score export
but failed when `sq.pl.spatial_scatter()` forwarded unsupported `show=False` to
Matplotlib under Squidpy 1.6.5. Tasks 5 and 7 (MF skin and CK bowel) timed out
during co-occurrence at the original eight-hour limit.

Removed the incompatible `show` argument, updated the pandas categorical dtype
check, removed a duplicated docstring line, and increased the active script 10
wrapper to 24 hours per task. The targeted rerun is array indices `1,4-7%2` after
the updated commit reaches HPC. Completed indices 0, 2, and 3 do not need to be
rerun. No HPC connection, file modification, monitoring, or submission occurred.

## Clean UMAPs and Larger Figure Typography (2026-10-26)

Added
`config/08_cell_type_annotation/vbct/min15_area95_selected_resolutions_clean_umaps.json`
for presentation-ready per-sample Scanpy UMAPs. This config reads the eight
existing annotated H5ADs and disables H5AD and summary-table writes. The figures
have no sample title, figure title, axis labels, or external legend; broad cell
types are labelled directly on the embedding with 12-point text. The original
titled figures remain unchanged. The 24 generated PNG, SVG, and PDF files are in
`figures/cell_type_annotation/vbct_min15_area95_selected_resolutions/clean_scanpy/per_sample_scanpy/`.

Script 08 now exposes title, legend, axis-label, combined-figure, and
summary-table controls while preserving the previous defaults. Script 09 now
reads its typography sizes from configuration. Both proportion configs were
updated with larger text and canvases, and both figure variants were regenerated
under `figures/cell_type_proportions/vbct_min15_area95_selected_resolutions/`.
The standard plot is 5100 x 2700 pixels and the tissue-faceted plot is
5700 x 2700 pixels, both at 300 DPI. Their overall titles and captions are
removed; the shared legends and tissue facet headings are retained. Input totals
remain eight samples and 2,083,396 cells, and the regenerated CSV values match
the previous results.

## Consolidated Session Handoff (2026-10-26)

### Current workflow and HPC boundary

The active order is 06 reclustering -> 08 reviewed annotation -> 10 Squidpy
spatial analysis. Script 09 is an optional summary-figure step after script 08.
The persistent user rule remains absolute: never connect to Phoenix or another
HPC with SSH, submit or monitor jobs remotely, use SSH-backed transfer tools, or
reconnect a mount. Existing local HPC mounts may only be inspected read-only.
Nathalie runs all HPC commands.

### Script 10 fixes and rerun state

Read-only log inspection for Slurm array job 16021129 established:

- completed: index 0 CK_skin_res, index 2 GR_lung_non_res_roi, and index 3
  BE_brain_non_res;
- reached Moran export but failed at top-Moran plotting: index 1 EL_skin_res,
  index 4 HW_brain_res, and index 6 MG_gastric_non_res;
- timed out during co-occurrence: index 5 MF_skin_non_res and index 7
  CK_bowel_res.

The plotting failure came from passing `show=False` to
`sq.pl.spatial_scatter()`; Squidpy 1.6.5 forwarded that keyword to Matplotlib,
which rejected it. The argument was removed from
`10_squidpy_recluster_spatial_analysis.py` (since renamed `10_cell_type_spatial_analysis.py`). The script also now uses
`isinstance(dtype, pd.CategoricalDtype)` and no longer has a duplicated
docstring line. The Slurm wrapper
`hpc/slurm/run_10_cell_type_spatial_analysis.sl` now requests 24 hours instead
of 8 hours. Local syntax and a small Squidpy 1.6.5 plotting smoke test passed.

After the updated commit is available on the HPC, the outstanding tasks can be
rerun from the HPC repository root with:

```bash
sbatch --array=1,4-7%2 hpc/slurm/run_10_cell_type_spatial_analysis.sl
```

Script 10 plots and result CSVs are written to:

```text
figures/spatial_analysis/vbct_min15_area95_selected_resolutions/
  <sample>/min15_area95_r<resolution>/squidpy/
```

This includes neighbourhood-enrichment, centrality, per-cell-type
co-occurrence, and top-Moran spatial-expression PNGs, plus centrality and Moran
CSVs. Spatially augmented H5AD files are written separately to:

```text
data/xenium/processed/vbct/spatial_analysis/
  min15_area95_selected_resolutions/<sample>/
```

### Clean Scanpy UMAPs

`08_annotate_cell_types_and_plot_umaps.py` now has config-driven controls for
sample titles, figure titles, external legends, axis labels, combined panels,
and summary-table writes. All new controls default to the previous behaviour, so
the reviewed annotation workflow remains compatible.

The new config
`config/08_cell_type_annotation/vbct/min15_area95_selected_resolutions_clean_umaps.json`
reads the existing script 08 annotated H5ADs, does not rewrite H5ADs or summary
tables, and produces eight clean per-sample Scanpy UMAPs. They have no titles,
axis labels, or external legends. Cell-type names are placed directly on the
embedding at 12 points, increased from 8 points. PNG, SVG, and PDF were generated
for every sample here:

```text
figures/cell_type_annotation/vbct_min15_area95_selected_resolutions/
  clean_scanpy/per_sample_scanpy/
```

The original titled per-sample and combined UMAP figures were retained.

### Clean cell-type composition figures

`09_plot_cell_type_proportions.py` now reads base, title, facet, legend,
caption, axis-label, and tick-label sizes from configuration. The supplied
configs use 14-point legend text, a 15-point legend title and y-axis label,
13-point tick labels, and 16-point tissue headings. Overall titles and captions
are disabled for these clean exports, while the shared cell-type legend remains.
The tissue plot retains Skin, Brain, Lung, and Gastrointestinal headings.

Both standard and tissue-faceted variants were regenerated as PNG, SVG, PDF,
and CSV under:

```text
figures/cell_type_proportions/vbct_min15_area95_selected_resolutions/
```

The standard PNG is 5100 x 2700 and the tissue-faceted PNG is 5700 x 2700, both
at 300 DPI. All eight samples and 2,083,396 cells were retained; per-sample
proportions sum to one and the regenerated values match the earlier results.

### Validation and repository state

Python compilation, JSON parsing, Git whitespace checks, output counts, PNG
integrity and DPI, PDF signatures, and SVG text-presence/absence checks passed.
The local image-viewer sandbox was unavailable, so the new figures have not had
a final visual inspection inside Codex.

Files changed for the clean UMAP/composition work are:

- `08_annotate_cell_types_and_plot_umaps.py`
- `09_plot_cell_type_proportions.py`
- `config/08_cell_type_annotation/README.md`
- `config/08_cell_type_annotation/vbct/min15_area95_selected_resolutions_clean_umaps.json`
- `config/09_cell_type_proportions/README.md`
- both active JSON configs under `config/09_cell_type_proportions/vbct/`

The worktree also contains conference-figure helper/config changes that were
already present and were deliberately left untouched. Review and stage changes
selectively before committing.

## Next Session Priority: Script 10 Vector Exports (2026-10-26)

The next session will update script 10 spatial-analysis plotting so its plots are
exported as SVG and PDF as well as PNG, with raster components rendered at
300 DPI. Preserve editable vector text and layout elements where the plotting
libraries allow it, while rasterizing dense point layers when needed to keep
files practical. Apply this consistently to the script 10 neighbourhood
enrichment, centrality, co-occurrence, and top-Moran spatial-expression plots
under:

```text
figures/spatial_analysis/vbct_min15_area95_selected_resolutions/
  <sample>/min15_area95_r<resolution>/squidpy/
```

Before rerunning the complete HPC analysis, verify locally that all three formats
are produced with stable filenames, PDFs are valid, SVG text remains editable,
and embedded raster layers use 300 DPI. Update the script 10 config/README and
HPC handoff instructions as part of that work. The Phoenix no-SSH rule continues
to apply; Nathalie will run any required HPC jobs.


## Design: Versioned Annotation Sets (2026-11-22)

Agreed with Nathalie in discussion on 2026-11-22. Implementation is on branch
`feature/annotation-sets`. This section supersedes the single-mapping layout
recorded above. The poster run is not special: it becomes annotation `v1`.

### Concept

Clustering facts and labels change at different rates, so they are split:

- **Cohort** (`config/cohorts/<project>/<cohort>.json`) holds the clustering facts: samples,
  response, tissue, branch, resolution, script 06 object, `groupby`, `umap_key`, and
  coordinate keys. It no longer holds the mapping, the palette or the annotated paths.
- An **annotation set** (`config/annotations/<project>/<cohort>/<version>/`) holds the labels.
  One cohort can have several sets. Versions are plain `v1`, `v2`, ..., and each
  folder's `CHANGES.md` records what the version represents and why.

```text
config/annotations/vbct/min15_area95_selected_resolutions/
  v1/
    annotation.json   # version, cohort, parent, status, date, release hash
    mapping.csv       # one row per source cluster
    palette.json      # per-tier order, colours, tier hierarchy
    CHANGES.md        # what this version is and why (human-written)
  v2/ ...             # parent: v1
```

Script 08 remains the annotation step; no new numbered script.

### Tiers

- `cell_type_tier_1` is the broad type (the current 8 types). `cell_type_tier_2` is an
  optional finer level. Further tiers follow the same pattern.
- Mapping columns are `cell_type_tier_1` and optionally `cell_type_tier_2`. Several
  clusters sharing a label is how clusters are merged, at any tier.
- **Hierarchy.** `palette.json` declares each tier-2 label's tier-1 parent. Script 08
  rejects any mapping row whose tier-2 label has a different parent, so tier-1 results
  stay comparable across versions.
- **Blank tier 2.** The cell inherits its tier-1 label and colour, and
  `obs["cell_type_tier_2_inherited"]` is True, so it is not mistaken for a real split.
  (Implemented as this flag rather than an `annotation_status` value;
  `annotation_status` in the mapping is `reviewed` or `proposed`.)
- **H5AD.** Script 08 writes one obs column per tier (`cell_type_tier_1`,
  `cell_type_tier_2`), each with `uns["<column>_colors"]`.
- **Downstream.** Steps choose `label_level: tier_1 | tier_2` (default `tier_1`). The
  level is part of the output path, so the two levels never overwrite each other.
- The old `cell_type` column is retired. There is no alias, to avoid duplicate label
  columns.

### Draft and released status

- **Draft.** A draft set can be edited and rerun freely (`--overwrite`).
- **Releasing.** `python -m xenium_workflow.annotation_sets release <set dir>` validates the set, then
  records a SHA-256 of `mapping.csv` and `palette.json` in `annotation.json`.
- **Edits after release.** Script 08 refuses to run a released set whose files no longer
  match the hash. Changes therefore always mean a new version
  (`python -m xenium_workflow.annotation_sets new <parent dir>` copies the parent as the next draft).
- **Script 10 needs a released set.** It is the expensive HPC step. Pass `--allow-draft`
  only for local smoke tests.

### Provenance and change records

- **In the H5AD.** `uns["cell_type_annotation"]` records the annotation version, cohort,
  parent, status, hash, tiers, mapping path, `groupby`, resolution and UMAP key. Script
  10 and the helpers check the stamp against the set they were asked to use.
- **Change table.** When a set has a parent on the same clustering, script 08 writes
  `changes_vs_<parent>.csv`. It has one row per source cluster whose label changed at
  any tier, with old and new labels and `n_cells`. `xenium_workflow/annotation_sets.py diff` prints the
  same comparison without H5ADs (no cell counts).
- **Run manifest.** Script 08 also writes a run manifest (JSON): inputs, outputs,
  set hash, and cell totals.

### Derived paths

Step configs no longer contain output folders. They name an annotation set and a
level, and both can be overridden on the command line (`--annotation`,
`--label-level`), so v2 reuses the same step configs.

```text
annotated H5AD   data/xenium/processed/<project>/cell_type_annotation/<cohort>/<version>/<sample>/
                   adata_expression_clean_<sample>_<branch>_cell_type_annotated_<version>_r<res>.h5ad
08 tables        data/xenium/processed/<project>/cell_type_annotation/<cohort>/<version>/
08 figures       data/xenium/output/<project>/cell_type_annotation/<cohort>/<version>/<level>/[<subdir>/]
other steps      data/xenium/{output,processed}/<project>/<step>/<cohort>/<version>/<level>/[<subdir>/]
```

Steps: `cell_type_proportions` (09), `spatial_analysis` (10 and its helpers), and
`poster`. An optional `output_subdir` in a step config separates helper outputs, for
example `immune_infiltration`.

### Validation in script 08

- Every source cluster is mapped exactly once, and every label is in its tier palette.
- The tier hierarchy is consistent.
- Cell IDs, shape and expression are unchanged.
- Per-tier counts sum to `n_obs`.
- A released set matches its hash.
- Tier-1 colours shared with the parent are unchanged; differences are reported.

### Migration

1. Build `v1` from the current mapping (`cell_type` -> `cell_type_tier_1`) and palette.
   Its `CHANGES.md` notes it was the labelling used for the HCA poster. Release it.
2. Rerun script 08 for v1 on the HPC to write stamped objects under the new paths.
3. Script 10 is not rerun just for migration. Its next run (pending the Moran's I
   review, `morans_i_significance.md`) reads v1 under the new layout.
4. Archive the unversioned annotation, proportion, spatial and poster outputs in one
   dated folder with a README, locally and on the HPC (Nathalie runs the HPC moves).

### Open points

- Moran's I does not depend on labels, so a tier-2 script 10 run repeats it. Revisit
  this when tier 2 is first used.
- Conference-figure helpers keep their own sample list and are out of scope.


## Session Update: Annotation Sets Implemented (2026-11-22)

Implemented the design above on branch `feature/annotation-sets` (not
committed; created from `main` at `bf83ff1`). Usage is documented in
`config/annotations/README.md`.

### What changed

- New `xenium_workflow/annotation_sets.py` (written at the repository root; moved into the package the same day). It loads and validates sets
  (palette tiers and hierarchy, mapping coverage against the cohort, release hash)
  and derives every output path. Its command line provides `new`, `check`, `diff`
  and `release`.
- `xenium_workflow/cohorts.py`: cohorts hold clustering facts only. It rejects the retired
  `palette_json`, `mapping_csv` and `annotated_h5ad` keys.
- v1 was created at `config/annotations/vbct/min15_area95_selected_resolutions/v1/`
  from the old mapping and palette (`git mv`; column `cell_type` became
  `cell_type_tier_1`; the 127 rows are otherwise byte-identical). It was released
  on 2026-11-22 with SHA-256 `f60c669e...540ef`.
- Script 08 writes the tier columns, stamp, per-tier UMAPs, `cluster_annotation_counts.csv`,
  `cell_type_counts_<tier>.csv`, `changes_vs_<parent>.csv` and the run manifest.
- Scripts 09, 10 and 10a, both spatial helpers and both poster helpers now take the
  set and level from `annotation` and `label_level` in the step config, with
  `--annotation` and `--label-level` overrides. They check the stamp and derive
  their paths. Scripts 09 and 10 and the helpers need a released set
  (`--allow-draft` overrides this).
- Step configs dropped `cohort` and all output folders. The helpers use
  `output_subdir`.
- Slurm wrappers accept `ANNOTATION=`, and script 10 also accepts `LABEL_LEVEL=` and `CONFIG=`.
- `tools/pipeline_status.py` checks the highest released set by default
  (`--annotation`, `--label-level`).
- Docs updated: concise and full workflow READMEs; the cohorts, 08, 09, 10,
  spatial_analysis and poster config READMEs; `hpc/slurm/README.md`.
- No scientific parameters changed. Script 10's Moran settings are untouched.
  The only change is in Squidpy figure filenames: `nccell_type` became
  `nctier_1`, for example.

### Verification (local, outputs written to the session scratchpad, not the repo)

- **GR v1 script 08.** Run on the local GR script 06 object (clean expression,
  24,040 cells). Labels, `.X`, cell order and colours are identical to the old
  annotated object, and the stamp is present.
- **GR v1 script 10.** Completed in about 1 minute. Its enrichment table is
  bit-identical to the existing GR table. Script 10a reproduces it exactly.
- **Scripts 09 and helpers.** Script 09 (both configs, all 8 samples, built from
  the existing count table) matches the existing proportion CSVs. The
  infiltration and enrichment-comparison helpers' GR rows match existing outputs.
  Both poster helpers run.
- **Throwaway v2.** Built with `new`, it contained two GR melanoma tier-2 splits
  and one Pericyte-to-Fibroblast merge.
  - `diff` and `changes_vs_v1.csv` listed exactly those 3 clusters (5,391 cells).
  - Tier counts reconciled (12,753 melanoma; 3,551 fibroblast; 18,856 inherited).
  - The tier-2 UMAP rendered correctly, and script 10 ran at tier 2.
  - v2 was then deleted.
- **Refusals confirmed.**
  - An object stamped by a different set.
  - Editing a released set (v1 was restored afterwards and checks OK).
  - A tier-2 label under the wrong parent.
  - Releasing with TODOs left in `CHANGES.md`.
  - Retired step-config keys.
  - Script 10 on a draft set without `--allow-draft`.
- Script 08 `--validate-only` passed on GR and MG, the only script 06 objects
  available locally. The Phoenix mount was not active, so the other six were not
  checked locally.

### Not done / next steps (for Nathalie)

1. Review the branch diff and commit if happy. Nothing has been committed.
2. Push the branch to the HPC checkout. Then, from the HPC repository root, regenerate v1:
   `sbatch hpc/slurm/run_08_cell_type_annotation.sl`. This writes the stamped v1
   objects under `data/xenium/processed/vbct/cell_type_annotation/min15_area95_selected_resolutions/v1/`.
3. Script 10 is not rerun just for migration. Its next run (after the Moran's I
   review in `morans_i_significance.md`) reads v1. Until then, the old unversioned
   script 10 results cannot be read by the updated helpers, which expect v1 paths
   and stamps.
4. After step 2, archive the old unversioned outputs in one dated folder with a
   README, locally and on the HPC (Nathalie runs the HPC moves). This was not done
   yet because the old outputs are the only complete script 08/10 results until v1
   is regenerated.
5. For the second labelling: `python -m xenium_workflow.annotation_sets new <v1 dir>`, edit
   `v2/mapping.csv`, `palette.json` and `CHANGES.md`, run script 08 with
   `ANNOTATION=<v2 dir>`, review, then `release`.


## Session Update: Shared Modules Moved Into a Package (2026-11-22)

At Nathalie's request, the three unnumbered shared modules at the repository
root moved into the new package `xenium_workflow/`, so the root holds only
pipeline steps:

- `workflow_config.py` -> `xenium_workflow/cohorts.py` (renamed; `REPO_ROOT` now
  resolves one level up)
- `annotation_sets.py` -> `xenium_workflow/annotation_sets.py`; its command line is
  now `python -m xenium_workflow.annotation_sets new|check|diff|release`
- `spatial_tables.py` -> `xenium_workflow/spatial_tables.py`

Imports were updated in scripts 08, 09, 10 and 10a, both spatial helpers and both poster helpers,
`tools/pipeline_status.py`, and the script 10 Slurm wrapper. Current docs
(`AGENTS.md`, both workflow READMEs, config READMEs) were updated. Older task
handoffs keep the old names as history. `banksy_utils/` was deliberately not
used, to keep upstream BANKSY code separate from pipeline plumbing.

## Session Update: v1 Regenerated on Phoenix (2026-11-22)

Nathalie ran script 08 for v1 on Phoenix from branch `feature/annotation-sets`
(commit `36529a2`, pushed to `origin`), Slurm job 16939096. It finished at 13:43.
It was checked read-only through the mount:

- **Samples.** All 8 validated and written under
  `data/xenium/processed/vbct/cell_type_annotation/min15_area95_selected_resolutions/v1/`.
  `tools/pipeline_status.py` reports `08-annot` 8/8.
- **Manifest and totals.** The run manifest shows v1 released (hash `f60c669e...`),
  2,083,396 cells and 8 H5ADs. The 127 clusters account for all 2,083,396 cells.
- **Counts.** `cell_type_counts_tier_1.csv` (52 rows) is identical to the old
  `reviewed_cell_type_counts.csv`.
- **Objects.** Every object's stamp matches released v1. Each has
  `cell_type_tier_1` and no `cell_type`. CK_skin_res `inspect_h5ad`: CLEAN EXPRESSION.
- **Warnings.** The `.err` contains only Scanpy deprecation and colormapping warnings.

Remaining from the migration plan: archive the unversioned outputs, and run
script 10 on v1 after the Moran's I review.

## Session Update: Old Outputs Archived Locally; v2 Draft Created (2026-11-22)

- **Archive script.** New `helper_scripts/archive_pre_annotation_sets_outputs.sh` is a
  one-off. It runs dry by default and moves only with `--apply`. It moves the
  unversioned script 08/09/10, helper and poster outputs into
  `data/xenium/_pre_annotation_sets_2026-11-22/` (same relative paths, README
  included) and never touches `vN/` folders.
- **Local run.** Applied locally: 32 paths, 2.1 GB. Not yet run on Phoenix
  (Nathalie runs it there after pulling).
- **Conference figures.** The config `config/conference_figures/vbct_min15_area95_selected_resolutions.json`
  read the old annotated objects on the Phoenix mount by absolute path. Its 8
  `adata_path` values now point at the archive location, which is valid once the
  script has run on Phoenix.
- **v2 draft.** `python -m xenium_workflow.annotation_sets new` created
  `config/annotations/vbct/min15_area95_selected_resolutions/v2/` (draft, parent
  v1, identical copy). Nathalie is editing it for the second round of cluster
  merges.
- **Script 10 for this week's meeting.** Nathalie wants script 10 for v1 and v2
  before the Moran's I fix, using the current settings with no parameter changes.
  Do not present the Moran FDR gene lists or the top-Moran plots from these runs
  (`morans_i_significance.md`). Neighbourhood enrichment, centrality and
  co-occurrence are unaffected.

## Session Update: v1 Cluster Marker Dot Plots (2026-11-23)

Per-sample dot plots of v1's clusters at each sample's selected resolution, to
review markers before the v2 merges. No new script 03 export was needed:

- **Input check.** The `min15_area95` rows of the local filter-sensitivity exports
  match v1 exactly for all eight samples: same clusters, same cell counts as v1's
  `cluster_annotation_counts.csv` on Phoenix (127 clusters), 304 genes from `.raw`,
  and the cohort's script 06 object.
- **Script 04 change.** It gained an optional `annotation` block, read by
  `apply_annotation_labels`. The block keeps one sample's selected-resolution
  clusters and labels them from a set. It requires every exported cluster to be
  mapped and every mapped cluster to be exported.
- **Script 04 rendering.** It also orders rows by label, can derive the output
  path, and accepts `--annotation`/`--label-level` overrides, so the same configs
  replot for v2. A new `rows_top_to_bottom` option is on by default only for
  annotation plots; existing configs are unchanged.
- **Configs and plot settings.** The configs are in
  `config/04_plot_dotplot/active/vbct_annotation_sets/min15_area95_selected_resolutions/`
  (README there). Nathalie chose the reviewed 89-gene filter and z-scoring within
  each sample, clipped at ±2.5.
- **Outputs.** All eight rendered locally (127 clusters labelled) to
  `data/xenium/output/vbct/cell_type_annotation/min15_area95_selected_resolutions/v1/tier_1/cluster_marker_dotplots/`.
- **Still to discuss.** Nathalie wants a combined dot plot that shows only each
  sample's selected resolution.


## Session Update: v2 Draft Committed; v2 Dot Plots (2026-11-23)

- **Edits.** Nathalie edited v2: BE_brain_non_res cluster 5 and HW_brain_res
  cluster 15 moved from Epithelial cell to Melanoma cell (melanocytic markers
  with KRT18). `CHANGES.md` was written; both rows are `reviewed`.
- **Loader fix.** Saving `mapping.csv` from a spreadsheet re-sorted the rows and
  rewrote resolutions (`1.00` became `1`), which broke `diff`.
  `load_mapping` now writes each row's resolution in the cohort's spelling. v2's
  file was re-sorted into v1 order, so it differs from v1 only in those two rows.
  `diff` reports 2 clusters changed.
- **Committed.** `7363034` (v2 still a draft) and pushed.
- **Dot plots.** BE and HW cluster dot plots were replotted with v2 labels
  (`--annotation .../v2`) to
  `data/xenium/output/vbct/cell_type_annotation/min15_area95_selected_resolutions/v2/tier_1/cluster_marker_dotplots/`.
  Method: `context/tasks/dotplot/dotplot.md` (Current State, item 1).
- **Next.** Script 08 for v2 on Phoenix
  (`ANNOTATION=config/annotations/vbct/min15_area95_selected_resolutions/v2`).
  Then check its output (expect `changes_vs_v1.csv` with the two clusters),
  release v2, commit and push, and run script 10 for v2.
