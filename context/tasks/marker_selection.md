# Task: PTMT PC55 Run 1 Cluster Annotation And Shared Marker Selection

## Objective

Annotate every BANKSY resolution 1.5 cluster in the eight PTMT PC55 run 1
samples and identify a reproducible shared marker set for each resulting cell
type, aiming for five markers per cell type where the data support that many.

This is a review-first analysis. Do not write labels into the script 01 QC
configs or overwrite source marker files until the proposed assignments have
been reviewed and approved.

## Agreed Annotation Rules

- Work at resolution 1.5 (`r1p5`).
- Use the detailed `minor_class` vocabulary rather than only broad major
  classes.
- Assign every cluster a proposed label, but use `Unresolved` or a mixed label
  when the evidence is insufficient or conflicting.
- Include a confidence rating for every cluster.
- Record the positive markers supporting each assignment.
- Record a secondary candidate for ambiguous clusters.
- Select shared markers using recurrence across independent samples, not just
  repeated clusters within one sample.
- Aim for five shared markers per cell type, but do not force five when fewer
  genes are specific and consistently supported.

## Samples And Cluster Counts

The eight run 1 samples contain 230 r1.5 clusters in total:

| Sample | Treatment | r1.5 clusters |
| --- | --- | ---: |
| `10329_run_1_DMSO` | DMSO | 19 |
| `10331_run_1_1644` | 1644 | 18 |
| `10707_run_1_DMSO` | DMSO | 29 |
| `10708_run_1_1644` | 1644 | 52 |
| `10722_run_1_DMSO` | DMSO | 27 |
| `10723_run_1_1644` | 1644 | 31 |
| `10733_run_1_DMSO` | DMSO | 26 |
| `10734_run_1_1644` | 1644 | 28 |

Cluster IDs run from zero to the cluster count minus one in each sample.

## Primary Inputs

### Combined r1.5 dot plot

Use this exact plot:

```text
figures/dotplots/local_ptmt_pc55_scRNA-seq_markers_minor_class_run_1_r1p50_plot_order_gene_cell_types_labels_zscore_dotplot.png
```

It was generated with:

```text
config/04_plot_dotplot/active/ptmt_pc55_clean_script00/multi_sample/ptmt_pc55_minor_class_run_1_r1p50_gene_cell_types_labels_local.json
```

The plot contains all eight run 1 samples at resolution 1.5. Genes are grouped
by `minor_class`, ordered by `plot_order`, and displayed as
`gene - cell_types`. Expression colour is a per-gene z-score and dot size is
percent expressing.

### Harmonised scRNA-seq marker reference

```text
data/xenium/raw_data/gene_markers/ptmt/2026-08-14_dotplot_marker_genes_harmonised_scRNAseq_analysis_Kat.csv
```

Relevant columns include:

```text
plot_order
gene
major_class
minor_class
cell_types
```

The file contains 195 unique marker genes, five major classes, and 38 minor
classes. Treat it as the primary annotation vocabulary and marker reference.
Some genes have multiple candidate cell types in `cell_types`; do not assume
the first semicolon-separated label is automatically correct.

### Per-sample r1.5 top-marker plots

Each sample has 11 plots, one per resolution. For this task, use only the
files ending in `r1p5_wilcoxon.png` under:

```text
/home/nnataren/mnt/phoenix_scratch/Banksy_py/data/xenium/output/ptmt_pc55/<sample>/top_marker_plot/
```

Example:

```text
/home/nnataren/mnt/phoenix_scratch/Banksy_py/data/xenium/output/ptmt_pc55/10733_run_1_DMSO/top_marker_plot/rank_genes_groups_10733_run_1_DMSO_pc55_nc0.20_r1p5_wilcoxon.png
```

### Per-sample r1.5 top-marker tables

Use the corresponding CSVs ending in `r1p5_wilcoxon.csv` under:

```text
/home/nnataren/mnt/phoenix_scratch/Banksy_py/data/xenium/output/ptmt_pc55/<sample>/top_marker_tables/
```

Each table has:

```text
group
names
scores
logfoldchanges
pvals
pvals_adj
```

These tables are the exact ranked-marker data behind the top-marker plots and
should be used for calculations instead of trying to transcribe gene names
from the PNGs.

### Script 03 dot-plot summaries

The plot-ready expression summaries are under:

```text
data/xenium/processed/ptmt_pc55/cross_sample_dotplot_exports/
```

There is one `*_ptmt_panel_dotplot_summary.csv` per sample. Each contains all
11 resolutions and 474 panel genes. Filter to resolution 1.5 for this task.
Useful columns include mean expression, percent expressing, cluster ID, and
sample/cluster identifiers. These data provide complementary evidence when a
gene ranks highly but is expressed in very few cells, or when two candidate
cell types share markers.

## Available Annotation Vocabulary

The harmonised marker file includes these broad families and detailed types:

- Epithelial: `LE-KLK3`, `LE-KLK4`, `Luminal`, `Basal`, `Club`, `Hillock`,
  `LPCs`, `NE`
- T/NK: `T cell`, `CD4T`, `CD4 Trm`, `CD4 naive-cm`, `Treg`, `CD8T`,
  `CD8 Trm`, `CD8 cytotoxic`, `T_MKI67`, `NK`, `NK cell`, `NK CD16pos`,
  `NK CD16neg`
- Myeloid: `MNP`, `Mono`, `Mph`, `Mac1`, `Mac2`, `Mac-MT1`, `Mac-cycling`,
  `Mye_cycling`, `DCs`
- B/mast: `B cell`, `Plasma_B`, `Mast`
- Stromal: `Endo`, `Fibroblast`, `CAFs`, `Mural`, `Pericytes`

Use these names where supported. A broader parent label is acceptable when
the detailed subtype is not distinguishable.

## Required Analysis

### 1. Build the cluster evidence table

For every sample and r1.5 cluster, collect:

- ranked positive markers;
- Wilcoxon score, log fold-change, and adjusted p-value;
- matching `major_class`, `minor_class`, and `cell_types` annotations;
- mean expression and percent expressing from the script 03 summary;
- the strongest and second-strongest candidate cell types;
- markers that contradict the proposed assignment.

Do not let marker classes with more genes win merely because they have a
larger reference set. Normalise enrichment scores for marker-set size and
inspect the actual genes driving each score.

### 2. Assign cluster labels

Use convergent evidence from:

1. specific high-ranking positive markers with positive log fold-change;
2. several markers from the same cell type rather than a single shared gene;
3. expression level and percent expressing in the combined dot plot;
4. recurrence of the same signature in other run 1 samples;
5. broad-lineage consistency before selecting a detailed subtype.

Distinguish biological state from lineage. For example, interferon response,
stress, cycling, or immediate-early genes should not replace the underlying
cell lineage unless the intended vocabulary explicitly defines that state.

Suggested confidence scheme:

- `High`: several specific markers agree, with good expression/coverage and
  little contradictory evidence.
- `Medium`: lineage is clear but subtype evidence is incomplete, or a smaller
  number of specific markers support the label.
- `Low`: assignment relies on one marker, shared markers, sparse expression,
  or conflicting signatures.
- `Unresolved`: no defensible lineage or a genuinely mixed signature.

### 3. Select shared markers by cell type

After cluster labels are assigned, rank candidate markers for each cell type
using:

- number of independent samples in which the marker is positively ranked;
- number and fraction of matching clusters containing the marker;
- median or robust average rank;
- positive log fold-change and adjusted significance;
- percent expressing and mean expression in matching clusters;
- specificity against clusters assigned to other cell types;
- agreement with the harmonised scRNA-seq marker reference.

Prioritise genes observed across both DMSO and 1644 samples where possible.
Avoid selecting generic epithelial, immune, housekeeping, stress, or cycling
genes when more specific markers are available. Record fewer than five
markers when five defensible shared markers do not exist.

## Deliverables

Create reviewable outputs without modifying the source marker CSV:

1. `ptmt_pc55_run1_r1p5_cluster_label_review.csv`
   - one row per sample/cluster;
   - proposed detailed label and broad class;
   - confidence;
   - supporting markers;
   - secondary candidate;
   - contradictory markers;
   - concise rationale and review notes.
2. `ptmt_pc55_run1_r1p5_shared_markers_by_cell_type.csv`
   - one row per selected cell-type/marker pair;
   - marker rank within cell type;
   - sample recurrence;
   - cluster recurrence;
   - summary rank/log-fold-change/expression statistics;
   - whether it is present in the harmonised marker reference;
   - caveats.
3. A concise Markdown summary of strong assignments, ambiguous groups, and
   cell types for which five shared markers could not be supported.

Proposed output directory:

```text
data/xenium/processed/ptmt_pc55/cluster_annotation_review/run_1_r1p5/
```

Keep all outputs provisional until reviewed by the user. Do not update
`config/01_QC/ptmt_pc55/*.json` as part of the first pass.

## Validation

- Confirm exactly 230 cluster rows are present in the cluster review table.
- Confirm all eight samples are represented with the expected cluster counts.
- Confirm cluster IDs are unique within each sample.
- Confirm every proposed marker exists in the available panel/export data.
- Confirm selected shared markers are supported in more than one sample when
  claiming that they are shared; otherwise label them sample-specific.
- Confirm no cell type is given five markers merely by duplicating aliases or
  selecting broadly nonspecific genes.
- Spot-check high-, medium-, low-confidence, and unresolved assignments
  against the dot plot and per-sample top-marker plots.

## Current Analysis Status (2026-09-23)

Both the detailed and broad analyses have been implemented and provisional
results have been generated. The source marker reference and script 01 QC
configs remain unchanged.

### Reproducible scripts

Detailed annotation and shared-marker analysis:

```text
analysis/ptmt_pc55_run1_r1p5_marker_triage/ptmt_pc55_run1_r1p5_marker_selection.py
```

Broad annotation, conditional recurrence, and sample-presence analysis:

```text
analysis/ptmt_pc55_run1_r1p5_marker_triage/ptmt_pc55_run1_r1p5_broad_cell_types.py
```

Broad scoring sensitivity analysis:

```text
analysis/ptmt_pc55_run1_r1p5_marker_triage/ptmt_pc55_run1_r1p5_broad_sensitivity.py
```

The detailed, broad, and sensitivity scripts passed `py_compile` on 2026-09-23.
The broad and sensitivity scripts import the validated input-loading and
per-gene contribution functions from the detailed script.

Run the broad workflow from the `Banksy_py` directory with:

```bash
/home/nnataren/miniforge3/bin/python3.12 \
  analysis/ptmt_pc55_run1_r1p5_marker_triage/ptmt_pc55_run1_r1p5_broad_cell_types.py
```

Use `--scratch-root /path/to/Banksy_py` if the Phoenix project is mounted at a
different location.

### Consolidated methodology

The complete methodology, formulas, origins of the numerical constants,
presence definitions, recurrence tiers, validation, and limitations are in:

```text
data/xenium/processed/ptmt_pc55/cluster_annotation_review/methods.md
```

The scoring coefficients, caps, floors, confidence cutoffs, and tier cutoffs
are exploratory heuristics introduced for this workflow. They were not fitted
to PTMT labels, supplied by BANKSY or Scanpy, taken from a published annotation
model, or calibrated as probabilities. This limitation must remain explicit.

### Detailed results

Directory:

```text
data/xenium/processed/ptmt_pc55/cluster_annotation_review/run_1_r1p5/
```

Files:

```text
ptmt_pc55_run1_r1p5_cluster_label_review.csv
ptmt_pc55_run1_r1p5_shared_markers_by_cell_type.csv
ptmt_pc55_run1_r1p5_annotation_summary.md
```

The detailed review has exactly 230 cluster rows. Its provisional confidence
distribution is:

- High: 58
- Medium: 54
- Low: 81
- Unresolved: 37

### Broad and subset-aware results

Directory:

```text
data/xenium/processed/ptmt_pc55/cluster_annotation_review/run_1_r1p5_broad_cell_types/
```

Files:

```text
ptmt_pc55_run1_r1p5_broad_cluster_label_review.csv
ptmt_pc55_run1_r1p5_broad_shared_markers.csv
ptmt_pc55_run1_r1p5_broad_cell_type_presence.csv
ptmt_pc55_run1_r1p5_broad_analysis_summary.md
```

The broad analysis re-scored the original marker evidence; it did not merely
rename detailed labels. Its provisional confidence distribution is:

- High: 80
- Medium: 54
- Low: 58
- Unresolved/Mixed: 38

Broad vocabulary:

- Luminal epithelial
- Basal epithelial
- Club/Hillock epithelial
- Cycling epithelial
- Neuroendocrine epithelial
- T cell
- NK
- B cell
- Plasma cell
- Mast
- Monocyte
- Macrophage
- Dendritic cell
- Cycling myeloid
- Endothelial
- Fibroblast/CAF
- Mural/Pericyte
- Unresolved/Mixed

The current broad shared-marker table contains 63 selected markers classified
as:

- Pan-cohort: 30
- Core conditional: 12
- Recurrent conditional: 0
- Partial support: 3
- Sample-specific: 18

The sample-by-cell-type presence table contains 136 combinations:

- Present: 66
- Possibly present: 16
- Not detected: 54

Presence definitions:

- `Present`: at least one High/Medium cluster in that sample.
- `Possibly present`: no High/Medium cluster but at least one Low cluster.
- `Not detected`: no assigned cluster; this is not proof of biological absence.

Marker-support tiers:

- `Pan-cohort`: marker positive in all eight samples.
- `Core conditional`: cell type strictly present in at least three samples and
  marker positive in every strict-positive sample.
- `Recurrent conditional`: cell type strictly present in at least three samples
  and marker positive in at least 75% of strict-positive samples.
- `Partial support`: marker positive in at least two inclusive-positive samples
  and at least 50% of inclusive-positive samples.
- `Sample-specific`: does not meet a shared-marker tier.

The shared-marker CSV reports strict and inclusive recurrence denominators,
DMSO and 1644 recurrence separately, supporting sample names, cluster
recurrence, specificity, expression summaries, and treatment-associated
caveats.

### Pan-cohort broad markers

Six broad cell types have selected markers supported in all eight run 1
samples:

- Basal epithelial: `KRT5`, `KRT14`, `SDC1`, `ID1`, `TMEM173`
- Fibroblast/CAF: `C1S`, `MMP2`, `DCN`, `LUM`, `PDGFRA`
- Luminal epithelial: `TMPRSS2`, `FOLH1`, `HOXB13`, `SLC45A3`, `KLK2`
- Macrophage: `FCGR2A`, `CD163`, `AIF1`, `CXCL8`, `CSF1R`
- Mast: `CPA3`, `HPGDS`, `PTGS1`, `KIT`, `MS4A2`
- T cell: `IL7R`, `CD2`, `CD3E`, `CD3D`, `PRF1`

`IGKC` is present in all seven samples with an assigned B-cell cluster, but B
cell does not meet the literal eight-sample criterion because no B-cell cluster
was assigned in the remaining sample.

### Sensitivity analysis

Directory:

```text
data/xenium/processed/ptmt_pc55/cluster_annotation_review/run_1_r1p5_broad_sensitivity/
```

Files:

```text
ptmt_pc55_run1_r1p5_broad_sensitivity_assignments.csv
ptmt_pc55_run1_r1p5_broad_sensitivity_cluster_stability.csv
ptmt_pc55_run1_r1p5_broad_sensitivity_summary.csv
ptmt_pc55_run1_r1p5_broad_sensitivity_summary.md
```

The sensitivity workflow compared the current composite score with rank-only,
unweighted-hit, rank-plus-significance, and no-coverage variants. Stable labels
across all five variants occurred in 115 of 230 clusters. Among clusters
resolved by the current composite score, 96 of 192 kept the same resolved label
across all variants. Rank-only, rank-plus-significance, and no-coverage variants
each matched the current labels for 207 of 230 clusters; unweighted-hit scoring
matched 120 of 230 clusters.

### Validation completed

- Exactly 230 rows are present in both cluster review tables.
- All eight samples have the expected cluster counts.
- Sample/cluster identifiers are unique.
- Selected markers occur in the panel exports.
- Recurrence fractions are bounded from zero to one.
- All broad markers have a recognized support tier.
- No source marker CSV or script 01 QC config was modified.

### Comparable minor-class workflow

New analysis directories:

```text
data/xenium/processed/ptmt_pc55/cluster_annotation_review/run_1_r1p5_minor_class/
data/xenium/processed/ptmt_pc55/cluster_annotation_review/run_1_r1p5_minor_class_sensitivity/
```

The current composite assignments reproduce the earlier detailed results:
High 58, Medium 54, Low 81, and Unresolved 37. The comparable marker table has
104 rows across 28 represented resolved labels: 18 Pan-cohort, 13 Core
conditional, 9 Recurrent conditional, 32 Partial support, and 32
Sample-specific.

Minor-class labels are stable across all five scoring variants for 103 of 230
clusters. Among the 193 current-composite resolved clusters, 83 remain resolved
with the same label across every variant. Rank-only, rank-plus-significance,
no-coverage, and unweighted-hit scoring match 203, 201, 206, and 108 current
labels respectively.

The all-cluster and priority-review dotplots are:

```text
figures/dotplots/automated_marker_selection/local_ptmt_pc55_run1_r1p5_minor_class_marker_review_zscore_dotplot.png
figures/dotplots/automated_marker_selection/local_ptmt_pc55_run1_r1p5_minor_class_marker_priority_review_zscore_dotplot.png
```

The priority plot includes Low, Unresolved, or sensitivity-unstable clusters
and is ordered by proposed minor class. It was visually inspected through a
reduced preview and showed nonblank, ordered lineage blocks without title or
axis overlap.

## Next Goal: Wet-Lab Expert Marker Augmentation

### Objective

Evaluate whether provisional gene-to-cell-type labels supplied by the wet-lab
team can strengthen existing assignments, resolve ambiguous clusters, or
recover additional cell types such as dendritic cells. This is a reference
augmentation task, not a new differential-expression analysis.

All candidate genes are already among the 474 Xenium targets. Their per-cluster
Wilcoxon rank, score, log fold-change, adjusted p-value, mean expression, and
percent expressing are already available. The missing information is the
gene-to-cell-type mapping that determines which candidate labels those genes
may support.

### Required input audit

First locate the wet-lab annotation file and establish its schema without
modifying it. Determine which columns contain gene symbol, provisional cell
type, comments, confidence, positive/negative marker role, and source. Check
gene-symbol normalization and intersect the proposed genes with the 474-gene
panel.

Classify every wet-lab gene-label pair as:

- concordant: the scRNA reference assigns the same gene to the same label;
- conflicting: the scRNA reference assigns the gene to another label;
- secondary-only: the label appears only in the existing multi-label
  `cell_types` field, which current automated scoring ignores;
- expert-only: the gene has no automated reference label; or
- absent/invalid: the gene cannot be matched to the 474-gene panel.

Produce a normalized, read-only-derived marker table with at least `gene`,
`proposed_cell_type`, `source`, `evidence_status`, `marker_role`, and `notes`.
Preserve the original wet-lab wording in a separate column.

### Three-reference comparison

Run the same eligibility and cluster-scoring workflow in three parallel modes:

1. Baseline scRNA reference only, reproducing the existing current-composite
   assignments.
2. Wet-lab expert reference only, showing what the provisional labels support
   independently.
3. Combined scRNA plus wet-lab reference, using each unique gene-label pair
   once.

Do not double-count a gene because two sources agree. Source agreement should
be reported as stronger provenance, not as two expression contributions. If a
gene is assigned to multiple candidate labels, retain the conflict and apply
ambiguity normalization across the unique labels rather than silently choosing
one.

The initial combined analysis should not give an automatic numerical bonus to
either source. Any later source weighting must be introduced only as a named
sensitivity variant.

### Assignment comparison

For every sample/cluster, compare baseline, expert-only, and combined labels
and classify the outcome as:

- unchanged;
- strengthened, where the label is unchanged but gains independent expert
  marker support;
- newly resolved;
- relabelled within the same broad lineage;
- relabelled across broad lineages;
- new cell type introduced; or
- conflicting/less resolved.

Report the supporting genes separately by source, together with their rank,
Wilcoxon score, log fold-change, adjusted p-value, mean expression, and percent
expressing. Retain the existing close-runner-up and confidence logic for the
first comparison so changes can be attributed to the augmented reference.

### New and rare cell-type safeguards

Treat a label supported by one expert-only gene as a weak candidate, even when
its numerical score is high. A new resolved label should normally require at
least two independent compatible markers; one-marker cases should be reported
as `<cell type>-like candidate` or Low confidence for visual review.

For dendritic cells specifically, first identify all expert-labelled
dendritic genes in the panel and determine whether they are concordant,
conflicting, secondary-only, or expert-only. Compare their evidence against
Monocyte, Macrophage, Myeloid, and Unresolved clusters. Require co-occurrence
of compatible markers and inspect contradictory myeloid evidence.

A rare population may be absorbed into a larger BANKSY cluster. If no distinct
cluster-level dendritic call appears, inspect individual cells within candidate
myeloid and unresolved clusters, provided the cell-level expression object is
available. Test whether multiple dendritic markers occur in the same cells and
form a coherent spatial subgroup before considering selective reclustering or
a higher-resolution BANKSY analysis.

### Validation and interpretation

- Keep the original scRNA-only assignments as an immutable baseline.
- Preserve source provenance in every augmented output.
- Report DMSO and 1644 support separately and retain sample-specific evidence.
- Recalculate broad and minor-class sensitivity after augmentation.
- Generate all-cluster and priority-review dotplots using the same display
  rules as the current workflow.
- Treat dotplots as review interfaces rather than independent validation.
- Do not use a gene to define a new label and then describe the same expression
  pattern as independent validation.
- Where available, compare candidate labels with morphology, spatial context,
  pathology review, protein staining, or an external prostate reference.
- Do not modify the source marker CSV, wet-lab file, or script 01 QC configs.

### Proposed outputs

Write new results under
`data/xenium/processed/ptmt_pc55/cluster_annotation_review/run_1_r1p5_expert_marker_augmentation/`
without overwriting current broad or minor-class outputs. Include:

- normalized expert marker reference and reference-overlap audit;
- baseline, expert-only, and combined per-cluster assignments;
- assignment-change classification table;
- source-separated supporting-marker evidence;
- updated cell-type presence and tiered marker tables;
- sensitivity and stability summaries;
- dendritic-candidate review table, including cell-level evidence if available;
- reproducible methods/validation summary; and
- overview and priority-review dotplots under
  `figures/dotplots/automated_marker_selection/expert_marker_augmentation/`.

### Completion criteria

The goal is complete when every wet-lab gene-label pair has an audit status,
all three reference modes have reproducible assignments for exactly 230
clusters, label changes are explicitly classified, new cell types are
supported by transparent multi-gene evidence or clearly marked as weak
candidates, and no existing analysis or QC configuration has been overwritten.

## Environment Notes

- The local image viewer failed with
  `bwrap: loopback: Failed RTM_NEWADDR: Operation not permitted` for these
  files, and no browser surface was available through computer use.
- If that failure persists, use the exact marker tables and script 03 summary
  CSVs as the primary evidence, and use local cropping/OCR only as a secondary
  visual aid. Do not claim that a PNG was visually inspected when it was not.
- Commands reading the Phoenix scratch paths may require escalated sandbox
  permission.
- On 2026-09-23, the Phoenix path
  `/home/nnataren/mnt/phoenix_scratch/Banksy_py` was available again. The broad
  script was rerun from the original Wilcoxon tables and passed validation. The
  cluster-label and presence outputs remained structurally unchanged, while the
  shared-marker table increased from the previously reconstructed 46 rows to 63
  rows because source-table recurrence restored 17 additional sample-specific
  markers. No markers were removed.
- The working tree contained an existing modification to
  `04_plot_multi_sample_dotplot_from_config_local.py`. Preserve it and do not
  revert unrelated changes.
- The new r1.5 labelled config was uncommitted at the end of the previous
  session.

## Recommended Starting Point For The Next Window

1. Read this document completely.
2. Locate and inspect the wet-lab gene-label file; do not assume its schema or
   modify the source.
3. Build the normalized expert marker table and reference-overlap audit.
4. Run baseline scRNA-only, expert-only, and combined comparisons without
   double-counting concordant gene-label pairs.
5. Review newly resolved, relabelled, new-cell-type, and conflicting clusters;
   prioritize dendritic evidence within Myeloid and Unresolved populations.
6. If cluster-level evidence is insufficient for a rare population, inspect
   cell-level co-expression and spatial localization before reclustering.
7. Regenerate sensitivity summaries and marker-review dotplots for augmented
   labels.
8. Do not transfer labels into script 01 QC configs until the user reviews and
   approves the augmented assignments.
