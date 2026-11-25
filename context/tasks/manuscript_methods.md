# Task: Manuscript Methods

<!-- Task file rules (see Banksy_py/AGENTS.md, "Task files"): edit Current
State, Next Step and Open Questions in place; add one dated line per decision;
no "Session Update" sections; keep this file under ~200 lines. -->

## Goal

Write the paper's Methods section in the repo, next to the code, one file
per subsection, updated as each analysis is finalised so that every number
traces to a config or script. The cell-type marker dot plots (issue #64) are
the first worked example.

## Current State (2026-11-24)

Nothing has been created yet. The structure below was proposed on 2026-11-24,
and Nathalie agreed with the idea; location and rules are still to confirm.

**Proposed layout** (`Banksy_py/manuscript/`):

```
manuscript/
  README.md                       conventions, section order, export (pandoc -> docx)
  methods/
    01_samples_and_xenium.md      cohort, panel, runs
    02_preprocessing_qc.md        scripts 01/05: filters, thresholds
    03_clustering_banksy.md       scripts 00/06: lambda, PCs, resolutions
    04_cell_type_annotation.md    script 08, annotation sets v1 -> v2
    05_marker_visualisation.md    first example (issue #64)
    06_spatial_analysis.md        script 10: neighbourhood enrichment, Moran's I
    07_software.md                package versions from the banksy env
    references.md
  figure_legends/
```

**Each methods file:**
1. Paper-ready text (past tense, journal level of detail).
2. A provenance block (not for the paper): scripts, configs, annotation
   version, git commit, and every parameter value with the config it came from.
3. `[verify]` markers on anything not yet checked against code.

**Proposed rules:**
- Write or update the section in the same session an analysis is finalised.
- Numbers come from configs or code, never from memory or handoff notes.
- A scientific parameter change updates its methods section in the same
  change (extends the Scientific Parameters rule in `AGENTS.md`).
- Run `tools/check_doc_citations.py` on `manuscript/` so stale script/config
  references are flagged.
- Methods = what and why; READMEs = how to run; task files = status. Link,
  don't duplicate.
- The repo is the source of truth; export to docx or a shared doc for
  supervisor review.

**First example: marker visualisation (facts established 2026-11-24)**
Sources: `config/paper_figures/README.md`,
`helper_scripts/paper_figures/plot_cell_type_marker_dotplot_from_config.py`,
`context/tasks/dotplot/dotplot.md` (Decisions).
- Per cluster (script 03): mean log-normalised expression and percent of cells
  with expression > 0, over all cells of the cluster (zeros included), from
  clean expression (`.raw`), never BANKSY-expanded features.
- Clusters labelled with annotation v2, tier 1; clusters of one cell type
  within a sample merged by cell-count-weighted averaging (exact, because
  script 03 averages over all cells). Merged counts match script 09.
- Dot size: percent expressing (comparable across all rows and plots).
- Colour: per-gene z-score of mean expression, sample SD, clipped at +/-2.5
  (Seurat DotPlot scaling convention):
  - by-sample plot: within each sample;
  - cohort plot: per-sample values averaged with equal sample weight, then
    z-scored across cell types;
  - supplementary panel: across all 50 sample x cell type rows.
- Wording: describe colour as "mean log-normalised expression, z-scored per
  gene (Seurat scaling convention)", not "Seurat DotPlot". Seurat averages
  expm1 values before log; we average log values (Scanpy convention). This
  affects absolute values, not row order, so plots are near-identical.
- Specificity table `diff_vs_rest` is a difference of mean log values
  (approximately, not exactly, a log fold change).
- Marker list is not final yet (dotplot task Next Step); write the section
  with the list as a placeholder until Nathalie confirms it.
- Upstream steps (QC branch min15_area95, BANKSY reclustering, annotation)
  have no sections yet: cross-reference them and mark claims `[verify]`.

## Decisions

- 2026-11-24 (Nathalie): Start a project-wide Methods draft in the repo, with
  the marker dot plot work as the first example.

## Next Step

Confirm location and rules with Nathalie, then create `manuscript/README.md`,
`methods/05_marker_visualisation.md` (text + provenance block) and empty
stubs for the other sections.

## Open Questions

- **Location/name:** `Banksy_py/manuscript/` (proposed) or elsewhere?
- **AGENTS.md rule:** add "parameter changes update the methods section" now?
- **Figure legends:** start `figure_legends/` with the marker dot plots, or
  wait until figures are final?
- **Supervisor review format:** docx export, or a shared doc?

## History

- 2026-11-24: Proposed during issue #64 (marker dot plots); see
  `context/tasks/dotplot/dotplot.md`.
