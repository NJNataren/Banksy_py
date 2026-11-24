# Task: HCA Scientific Conference Poster

Started: 2026-10-25. Latest update: 2026-10-28.

## Path Update (2026-11-20, post-08 cleanup)

Paths below are historical. On branch `cleanup/post08-flow` (see
`post08_flow_cleanup.md`; merge to `main` pending an HPC test):

- Poster configs moved from `config/08_cell_type_annotation/vbct/poster_*.json`
  to `config/poster/vbct/` (README there). They name the cohort config instead
  of the script 08 config.
- New poster helper outputs: figures to
  `data/xenium/output/vbct/poster/min15_area95_selected_resolutions/`, dot plot
  tables to `data/xenium/processed/vbct/poster/min15_area95_selected_resolutions/cohort_marker_dotplot/`.
  The poster's existing figures stay in
  `figures/conference_min15_area95_selected_resolutions_poster/`; do not move
  them while a poster design file may link to them.
- Script 08/09/10 paths changed as listed in `cell_type_annotation.md`
  ("Path Update"); script 10 is now `10_cell_type_spatial_analysis.py`.
- Local copies of the old `figures/spatial_analysis/`, `figures/cell_type_annotation/`, and `figures/cell_type_proportions/` folders were moved to `figures/_pre_post08_cleanup_2026-11-20/` on 2026-11-21 (README there); copies on Phoenix are untouched. The poster folder itself was not moved.

## Purpose and current status

Nathalie requested a new poster-editing project based on the Xenium/BANKSY work,
then requested a scientific-conference critique incorporating the newer
cell-type annotation and spatial analyses. A simple A0 portrait PDF mock-up
was created with draft text and seven labelled figure placeholders. Nathalie
responded positively and requested this task handoff.

The mock-up is a proposed structure, not a completed scientific poster or
approval of every scientific statement. No original poster, analysis input,
figure asset, or mounted HPC file was modified. No analyses were rerun and no
HPC connection was initiated. The editable Canva source has not been supplied.

This handoff supersedes the initial poster-project orientation notes in
`../poster/README.md` (relative to the analysis root), which still say the
poster source and dimensions are unknown. The supplied source is now known,
and its text and visual layout have been inspected.

## Path conventions and access rules

Analysis root:

`/home/nnataren/Documents/PhD/Bioinformatics/Banksy_py_fork/Banksy_py`

Workspace root:

`/home/nnataren/Documents/PhD/Bioinformatics/Banksy_py_fork`

Unless marked otherwise, paths below are relative to the analysis root.

Existing local HPC mount inspected read-only:

`/home/nnataren/mnt/phoenix_scratch/Banksy_py`

Never initiate SSH, scp, sftp, SSH-backed rsync, remote jobs, connection probes,
or mount reconnection. Existing mounted files may be inspected, not modified.
Nathalie handles all remote connections and execution. A directory listing is
not proof of job completion; the output inventory below is a dated snapshot.

## Context loaded before making recommendations

Read the following fully:

1. Workspace `AGENTS.md` and analysis `AGENTS.md`.
2. `README_00_04_XENIUM_ANALYSIS_WORKFLOW_CONCISE.md`.
3. `context/tasks/conference_figures.md`, including the August 29 poster-polish
   update superseding earlier crop-selection notes.
4. `context/tasks/cell_type_annotation/cell_type_annotation.md` (moved and converted
   2026-11-24; the October 24-25 updates are in its `archive/`).

Skills applied:

- `/home/nnataren/.codex/skills/banksy-context/SKILL.md`
- `/home/nnataren/.codex/plugins/cache/openai-primary-runtime/pdf/26.1022.11436/skills/pdf/SKILL.md` (Codex plugin cache; no longer present)
- `/home/nnataren/.codex/skills/python-script-documenter/SKILL.md`

The conference handoff describes older unfiltered CK bowel (resolution 1.10)
and MG gastric (0.70) exemplars. It is historical background, not evidence that
the supplied PDF uses those exact objects. It records these source locations:

- `config/conference_figures/vbct_exemplar_ck_bowel_mg_gastric.json`
- `helper_scripts/conference_figures/plot_spatial_conference_figures_from_config.py`
- `helper_scripts/conference_figures/README.md`
- `figures/conference/`

Those files were referenced through the handoff; their complete contents were
not independently read for the mock-up. Do not substitute this old branch for
the newer selected-resolution min15/area95 branch.

## Original poster: inspected source and findings

Source PDF:

`figures/conference_filtered_reclustered_qc_v1_broad_cell_types/HCA_research_poster_2026_Nathalie_Nataren_v1.pdf`

Verified using `pdfinfo`, `pdftotext -layout`, and a rendered image:

- One page, A0 portrait (2384.25 x 3370.5 PDF points).
- Canva is the PDF creator; author is Nathalie Nataren.
- Title: **Spatial profiling of melanoma vasculature in immune checkpoint inhibitor response**.
- Authors: Nathalie Nataren, Lisa Ebert, Katherine Pillman, Claudine S Bonder.
- Three columns with introduction, cohort/panel, aims and methods on the left;
  GI exemplar UMAPs, spatial maps and endothelial marker dotplots in the middle;
  CDH5 expression, endothelial distribution, multimodal morphology, discussion
  and future directions on the right.
- Full-width take-home statement, references, sponsors and contact information.

Original result panels:

1. Figure 1: GI responder/non-responder UMAPs and spatial cell-type maps.
2. Figure 2: endothelial-marker expression dotplots.
3. Figure 3: CDH5 expression and endothelial cluster location.
4. Figure 4: spatial crop and multimodal vessel morphology comparison.

Visual assessment:

- Clear vascular focus and useful tissue images; retain the complementary
  marker/morphology evidence as the centrepiece.
- Workflow claims repeat across introduction, aim, methods and conclusion.
- Two GI exemplars dominate; quantitative evidence occupies relatively little
  space, and dotplots/labels need more room.
- New results justify a cohort-wide story rather than only workflow development.
- Use a single shared cell-type palette across retained and new panels.

The original PDF was initially rendered but could not be displayed because the
image tool's filesystem sandbox failed. A later smaller JPEG preview was
successfully displayed and inspected. Visual review of the original poster is
therefore complete at overview scale, not still blocked.

Temporary original-poster previews (not durable deliverables):

- `/tmp/hca_poster_preview.png`
- `/tmp/hca_poster_review.jpg`

### Nearby editable figure assets found

The supplied PDF's parent directory was searched for SVG/PPTX/ODP files.
Individual SVG assets were found; no editable whole-poster PPTX/ODP was returned.
These were inventoried, not individually visually reviewed:

`figures/conference_filtered_reclustered_qc_v1_broad_cell_types/Final_poster_figure/`

Examples under that directory:

- `ici_cohort_figure.svg`
- `dotplot/dotplot_edit.svg`
- `umap/CK_bowel_res_umap.svg`
- `umap/MG_gastric_non_res_umap_cell_type_labels_r0p70_filtered_reclustered_on_data_edit.svg`
- `cell_types/CK_bowel_res/CK_bowel_res_spatial_cell_type_labels_r1p10_filtered_reclustered_edits.svg`
- `cell_types/MG_gastric_non_res/MG_gastric_non_res_spatial_cell_type_labels_r0p70_filtered_reclustered_edits.svg`
- `marker_spatial/MG_gastric_non_res_spatial_endothelial_CDH5_edit.svg`
- `individual_cluster/MG_gastric_non_res_spatial_individual_cluster_r0p70_filtered_reclustered_cluster_7_Endothelial.svg`
- `multimodal/crop_cluster_CK_bowel_res.svg`
- `multimodal/crop_multimodal_CK_bowel_res.svg`
- `multimodal/crop_multimodal_CK_bowel_res_no_arrows.svg`
- `crops/crop_multimodal_CK_bowel_res_no_arrow.svg`
- `ccb_logo (1).svg` and `au_logo (1).svg`

## New annotation and composition: evidence examined

The complete annotation handoff records:

- All eight VBCT samples annotated using reviewed mappings.
- 127 source clusters, 2,083,396 post-QC cells.
- Branch: min15/area95 (at least 15 transcripts; grouped cell-area p95 filtering).
- Shared broad labels: Melanoma cell, Macrophage, T cell, B cell, Fibroblast,
  Endothelial cell, Epithelial cell, Pericyte.
- Cell IDs, clean expression and original clustering provenance preserved.
- Annotation and composition completed; these are no longer future tasks.

Selected resolutions from the handoff:

| Sample | Resolution |
|---|---:|
| CK_skin_res | 1.00 |
| EL_skin_res | 0.80 |
| GR_lung_non_res_roi | 1.00 |
| BE_brain_non_res | 0.90 |
| HW_brain_res | 0.80 |
| MF_skin_non_res | 0.50 |
| MG_gastric_non_res | 0.70 |
| CK_bowel_res | 0.80 |

Mounted annotation directory actually listed:

`/home/nnataren/mnt/phoenix_scratch/Banksy_py/figures/cell_type_annotation/vbct_min15_area95_selected_resolutions/`

It contains `combined_scanpy_cell_type_umaps.{png,svg,pdf}` and
`per_sample_scanpy/<sample>_cell_type_umap_scanpy.{png,svg,pdf}` for all eight
samples. Existence was checked; the new eight-panel UMAP montage was not
visually inspected in this session and annotated H5ADs were not loaded.

Local composition files located:

- `figures/cell_type_proportions/vbct_min15_area95_selected_resolutions/cell_type_proportions.csv`
- `figures/cell_type_proportions/vbct_min15_area95_selected_resolutions/cell_type_proportions_by_tissue.csv`

The complete `cell_type_proportions_by_tissue.csv` was read. It contains cell
type, count, sample, branch, resolution, sample total and proportion. Directly
observed endothelial fractions (rounded):

| Sample | Endothelial cells | Fraction of retained annotated cells |
|---|---:|---:|
| CK_skin_res | 542 | 8.88% |
| EL_skin_res | 18,518 | 7.86% |
| GR_lung_non_res_roi | 1,329 | 5.53% |
| BE_brain_non_res | 4,666 | 6.61% |
| HW_brain_res | 13,495 | 3.71% |
| MF_skin_non_res | 35,206 | 5.36% |
| MG_gastric_non_res | 7,307 | 6.00% |
| CK_bowel_res | 57,485 | 9.51% |

Thus the mock-up uses the descriptive range **3.7-9.5%**. These are fractions
of retained annotated cells, not vascular density or whole-tissue fractions.
The table also contains 207 pericyte-labelled cells in GR lung and zero
pericyte-labelled cells in the other samples. This is an annotation result,
not proof of biological absence elsewhere or independent validation of identity.

The handoff reports composition figures with the same stems in PNG/SVG/PDF,
and a tissue-faceted version (Skin 3, Brain 2, Lung 1, GI 2). These figures were
recommended for the poster but were not visually reviewed in this session.
The corresponding mounted `figures/cell_type_proportions/` directory was absent
when checked; the CSV used here was local.

Related paths recorded in the handoff, not fully reread for this mock-up:

- `08_annotate_cell_types_and_plot_umaps.py`
- `09_plot_cell_type_proportions.py`
- `config/08_cell_type_annotation/vbct/min15_area95_selected_resolutions.json`
- `config/annotations/vbct/min15_area95_selected_resolutions/v1/mapping.csv` (moved 2026-11-22 from `config/08_cell_type_annotation/vbct/min15_area95_selected_resolutions_cell_type_mapping.csv`; column `cell_type` renamed `cell_type_tier_1`)
- `config/annotations/vbct/min15_area95_selected_resolutions/v1/palette.json` (moved there from the script 08 config folder on 2026-11-22)
- `config/09_cell_type_proportions/vbct/min15_area95_selected_resolutions_by_tissue.json`
- `data/xenium/processed/vbct/cell_type_annotation/min15_area95_selected_resolutions/reviewed_cell_type_counts.csv`
- `data/xenium/processed/vbct/cell_type_annotation/min15_area95_selected_resolutions/reviewed_cluster_to_cell_type_counts.csv`

## Spatial analysis: actual inventory versus historical handoff

Mounted directory inventoried recursively:

`/home/nnataren/mnt/phoenix_scratch/Banksy_py/figures/spatial_analysis/vbct_min15_area95_selected_resolutions/`

At inspection, five sample directories existed:

| Sample/subdirectory | Outputs observed |
|---|---|
| `CK_skin_res/min15_area95_r1.00/squidpy/` | Neighbourhood heatmap; centrality PNG/CSV; co-occurrence PNGs for five present types; all/FDR Moran CSVs |
| `EL_skin_res/min15_area95_r0.80/squidpy/` | Neighbourhood heatmap; centrality PNG/CSV |
| `GR_lung_non_res_roi/min15_area95_r1.00/squidpy/` | Neighbourhood heatmap; centrality PNG/CSV; co-occurrence PNGs for seven present types including pericyte; all/FDR Moran CSVs |
| `BE_brain_non_res/min15_area95_r0.90/squidpy/` | Neighbourhood heatmap; centrality PNG/CSV; co-occurrence PNGs for seven present types; all/FDR Moran CSVs |
| `HW_brain_res/min15_area95_r0.80/squidpy/` | Neighbourhood heatmap; centrality PNG/CSV |

No MF skin, MG gastric or CK bowel spatial directory was returned. This shows
progress beyond the older handoff's unrun spatial stage, but does not prove
that any job finished successfully or that missing outputs failed. Recheck
before updating the final poster.

Read the full local spatial config:

`config/10_squidpy_recluster_analysis/vbct/cell_type_annotation/min15_area95_selected_resolutions.json`

It references script 08's config and uses `cell_type` annotations, with:

- Coordinates `x`, `y`.
- Moran subsampling fraction 0.25.
- Moran permutations 10, jobs 1, FDR threshold 0.1, top genes 5.
- Figures under `figures/spatial_analysis/vbct_min15_area95_selected_resolutions`.
- Results under `data/xenium/processed/vbct/spatial_analysis/min15_area95_selected_resolutions`.

Inspected relevant sections of `10_squidpy_recluster_spatial_analysis.py`:
spatial graph creation, neighbourhood enrichment, centrality CSV/PNG saving,
co-occurrence calls, and searched Moran settings/output logic. The initial
graph uses `sq.gr.spatial_neighbors(..., coord_type="generic")`;
neighbourhood enrichment and co-occurrence use library defaults for settings
not explicitly passed. The Moran subsample graph uses `delaunay=True`.
No installed Squidpy version/defaults or mounted execution logs were checked,
so the local config is not proof of the exact completed remote run settings.

### Spatial figures visually inspected

Under the mounted spatial root, displayed these two actual PNGs:

- `CK_skin_res/min15_area95_r1.00/squidpy/co_occurrence_cluster_Endothelial cell_CK_skin_res_pcannotated_nccell_type_r1p00.png`
- `CK_skin_res/min15_area95_r1.00/squidpy/neighbourhood_enrichment_CK_skin_res_pcannotated_nccell_type_r1p00.png`

The co-occurrence plot has five curves, generic distance/value labels and an
oversized mathematical title. The heatmap has strong diagonal structure and
a wide continuous scale. Neither was used to claim a response-associated
mechanism. Suggested poster improvements:

- Choose an explicit sample-specific biological observation before selecting a panel.
- Add verified distance units and a co-occurrence reference line at ratio 1.
- Label the enrichment score explicitly and use a clear colour scale.
- Use selected cell-type relationships and pair the quantitative panel with
  tissue context where useful.
- Distinguish spatial association from direct cell-cell interaction.

Temporary inspection files: `/tmp/endothelial_cooccurrence_review.jpg` and
`/tmp/neighbour_review.jpg`. The first JPEG conversion gave transparency a dark
background; that is a preview-conversion issue, not evidence of a source defect.

Other spatial PNGs, centrality CSV values and Moran CSV values were inventoried
but not analysed or visually reviewed. Moran's I was not recommended as a
headline result: it is peripheral to the vascular story and the current local
config's 10 permutations warrants review before final inference.

## Scientific and editorial recommendations

Proposed story:

> We annotated major cell populations across eight melanoma metastases,
> identified endothelial populations using molecular and morphological evidence,
> and began characterising their spatial relationships with the tumour
> microenvironment.

Three connected results:

1. Cohort-wide annotation and cell composition.
2. Endothelial identity supported by marker expression and matched morphology.
3. Exploratory vascular neighbourhood analysis.

Editorial choices used in the mock-up:

- Proposed title: **Spatial characterisation of vascular populations in metastatic melanoma**.
- Treat ICI response as the biological motivation; no response-associated
  vascular signature has been established by this review.
- Shorten repeated workflow prose; retain a compact workflow schematic.
- Move completed annotation and composition out of Future Directions.
- Retain and enlarge endothelial evidence, especially matched tissue morphology.
- Prioritise tissue-faceted composition over a full eight-UMAP montage.
  A small UMAP is optional if it answers a clear question; the full montage
  can be supplementary. No UMAP placeholder is required in the current mock-up.
- Use common labels/colours; old and new figures may use different palettes.
- Prefer "supports endothelial identity" to blanket "validation" claims.
  Marker evidence used for labelling is not independent validation.
- Describe BANKSY embeddings as spatially informed, rather than attributing
  all UMAP separation to transcription alone.

Outstanding scientific checks:

- Verify unique patient count. The original cohort schematic says pictograms
  represent patients, while CK appears in skin and bowel sample names. Do not
  equate eight samples with eight independent patients without metadata.
- Treat tissue site, sample size, ROI sampling and patient dependence explicitly
  in future response comparisons; cells are not independent patient replicates.
- Reconcile each retained panel's branch and resolution. Earlier CK bowel
  assets carry r1p10, while the new selected annotation is r0.80. MG remains
  r0.70 but branch provenance still requires checking.
- Confirm clean expression sources for all marker panels; never interpret
  BANKSY-expanded/scaled matrices as biological expression.
- Review the 207-cell pericyte-labelled population before elevating it to a result.
- Verify literature claims, references, contact spelling and final captions
  during full poster editing. No external literature search or reference
  validation was performed for this mock-up.

## Mock-up deliverables and layout

All paths here are absolute and local to the workspace:

- PDF: `/home/nnataren/Documents/PhD/Bioinformatics/Banksy_py_fork/poster/output/pdf/melanoma_poster_mockup.pdf`
- JPEG preview: `/home/nnataren/Documents/PhD/Bioinformatics/Banksy_py_fork/poster/output/melanoma_poster_mockup.jpg`
- Editable generation script: `/home/nnataren/Documents/PhD/Bioinformatics/Banksy_py_fork/poster/build_mockup.py`
- Intermediate PNG: `/home/nnataren/Documents/PhD/Bioinformatics/Banksy_py_fork/poster/tmp/pdfs/melanoma_poster_mockup.png`
- Initial project brief: `/home/nnataren/Documents/PhD/Bioinformatics/Banksy_py_fork/poster/README.md`

The PDF is A0 portrait (841 x 1189 mm), generated using ReportLab, with navy
headings, teal accents, three columns and labelled light-grey figure boxes.
It contains draft prose and real descriptive numbers from the reviewed context,
but **no actual data figures or simulated scientific plots**. All seven image
positions are explicit placeholders:

| Label | Placeholder | Purpose |
|---|---|---|
| A | Eight pretreatment metastases | Tissue/response cohort schematic; patient count pending |
| B | Cell-type composition across samples | Tissue-faceted 100% stacked bars, shared palette |
| C | Canonical marker-expression dotplot | Endothelial annotation with other cell types as context |
| D | Endothelial spatial map | Tissue location, matched to morphology |
| E | Matched multimodal tissue crop | Vessel morphology, scale bar, channels; optional CDH5 inset |
| F | Endothelial neighbourhood enrichment | Selected sample-specific spatial result |
| G | Endothelial co-occurrence by distance | Selected relationships, distance units, ratio-1 baseline |

Left column: rationale; cohort/question; five-step workflow.
Centre: result 1 (cohort annotation/composition); result 2 (endothelial identity).
Right: result 3 (vascular neighbourhoods); interpretation/limitations; next step.
Header: title/authors/affiliation and eight samples / 304 genes / 2.08 million
post-QC cells. Footer: take-home message, references and acknowledgement spaces,
and an explicit layout-mock-up/date label.

Draft take-home wording in the mock-up:

> Cohort-wide annotation identifies endothelial populations across eight
> melanoma metastases, supported by marker expression and vessel morphology,
> enabling exploratory study of vascular organisation.

This is draft framing; do not imply morphology was independently checked for
every sample. Full panel-level evidence review remains outstanding.

## Generation and validation

The environment's sandbox helper repeatedly failed with
`bwrap: loopback: Failed RTM_NEWADDR: Operation not permitted`. Local shell
commands were run through approved escalation; no SSH workaround was used.
The default `python` command was absent. `python3` was available but lacked
ReportLab, so dependencies were installed into `/tmp/poster_pdf_dependencies`.
That temporary dependency directory may not persist.

From the workspace root, the successful render commands were:

```bash
PYTHONPATH=/tmp/poster_pdf_dependencies python3 poster/build_mockup.py
pdftoppm -scale-to 1800 -singlefile -png poster/output/pdf/melanoma_poster_mockup.pdf poster/tmp/pdfs/melanoma_poster_mockup
convert poster/tmp/pdfs/melanoma_poster_mockup.png -background white -alpha remove -quality 75 poster/output/melanoma_poster_mockup.jpg
python3 -m py_compile poster/build_mockup.py
```

The PDF skill's artifact-start marker was run successfully once before authoring.
The script syntax check passed. The PDF was rendered and the complete JPEG was
successfully displayed and visually inspected: no obvious clipping, overlaps or
text overflow were observed at overview scale. The final PDF and JPEG were
linked to Nathalie. The source poster remains unchanged.

## Next session

Load this file with the banksy-context skill. Review the mock-up with Nathalie
before replacing placeholders or changing the scientific emphasis. Then:

1. Confirm the preferred layout and draft title/story.
2. Recheck current spatial output availability on the existing mount only.
3. Select exact figure panels and verify branch, resolution, expression source,
   sample identity, colour palette and scientific interpretation.
4. Confirm patient metadata and response labels.
5. Replace placeholders with selected, legible figures and concise captions.
6. Verify references and wording, render again and inspect at poster scale.

These are handoff suggestions, not authorization to rerun analyses, connect to
HPC, modify mounted files, or edit the original Canva design.

## Conference-figure regeneration from min15/area95 objects (2026-10-26)

The reusable conference helper can consume the reviewed script 08 annotated
objects directly. A new eight-sample config is available at
`config/conference_figures/vbct_min15_area95_selected_resolutions.json`. It reads
the existing mounted H5ADs under
`/home/nnataren/mnt/phoenix_scratch/Banksy_py/data/xenium/processed/vbct/` and
writes figures locally to
`figures/conference_min15_area95_selected_resolutions_poster/`; mounted inputs remain
read-only.

Each sample uses its selected min15/area95 BANKSY cluster column, saved UMAP key,
and existing `cell_type` annotation. The helper now honours the shared biological
cell-type order as well as the exact script 08 palette, so legends and colours
remain consistent with the cell-type UMAPs and proportion figures. The configured
endothelial marker panel is AQP1, CALCRL, CDH5, ECSCR, PLVAP, SELP and VWF.

The poster-focused config generates cell-type tissue maps, endothelial dotplots, and cell-type crops; UMAPs are excluded because script 08 already produced them. The exact established CK bowel and MG gastric crop windows were copied from the earlier broad-cell-type conference config. The full eight-sample local batch completed successfully using read-only backed access to the mounted H5ADs. Full-tissue maps are available as PNG/PDF, dotplots as PNG/PDF/SVG, and all three MG plus both CK bowel crops as PNG/PDF/SVG under `figures/conference_min15_area95_selected_resolutions_poster/`. No UMAP files were written. All expected files were present and non-empty; 29 PNGs passed Pillow verification, 29 PDFs passed header/EOF checks, and 29 SVGs parsed as XML. Exact palette RGB values from scripts 08/09 were confirmed in the tissue-map PNGs.

### Cell-type legend correction (2026-10-26)

The first regenerated `cell_type_labels` and crop figures used correct point
colours but an incorrectly paired custom legend: labels were alphabetically
sorted after the palette had been assembled in the shared biological order.
`force_categorical_legend()` now preserves the categorical order, keeping each
label paired with its script 08/09 colour. All eight full-tissue maps and all five
CK bowel/MG gastric crops were regenerated in PNG, PDF and SVG. Automated SVG
inspection checked 13 legends against the exact shared label/order/hex mapping
and found zero errors.

### Endothelial marker spatial maps (2026-10-26)

The same eight-sample config now enables `plot_marker_spatial`. A marker-only local run generated a seven-gene overview and seven individual whole-tissue expression maps for every sample from clean biological expression in the reviewed min15/area95 annotated H5ADs. The genes are AQP1, CALCRL, CDH5, ECSCR, PLVAP, SELP and VWF. Outputs are stored per sample under `figures/conference_min15_area95_selected_resolutions_poster/<sample>/marker_spatial/`, with single-gene maps in `individual_markers/`. Each figure is available as PNG, PDF and editable SVG. The completed batch contains 8 combined panels and 56 individual maps in each format (64 files per format; 192 total). Pillow, PDF header/EOF and XML parsing checks found zero integrity errors and all 56 expected sample-gene SVG combinations were present. The CK bowel combined panel was visually inspected and rendered cleanly. Expression uses a red sequential scale from zero to each gene's 99th percentile within each sample, so colour-bar values must be read per gene and sample rather than compared as a shared absolute scale.

### Endothelial BANKSY single-cluster spatial plots (2026-10-26)

`helper_scripts/conference_figures/plot_individual_cluster_spatial_from_config.py` now supports repeatable `--only-cell-type` filtering and honours the config's read-only backed H5AD mode. A focused local run generated one grey-context spatial plot for each reviewed BANKSY cluster annotated as `Endothelial cell`, highlighting that cluster with the exact script 08/09 colour `#CC3311`. Ten clusters were exported: CK skin 5; EL skin 8, 12 and 16; GR lung 9; BE brain 8; HW brain 10; MF skin 8; MG gastric 6; and CK bowel 2. Outputs are stored under each sample's `individual_clusters/` directory in PNG, PDF and editable SVG formats. All 30 files passed format-integrity checks, all expected sample-cluster combinations were present, and direct PNG inspection confirmed exact `(204, 51, 17)` pixels in every plot. The CK bowel plot was visually inspected after export. These plots preserve BANKSY cluster identity; the shared red indicates the reviewed broad endothelial annotation and should not be presented as evidence that different endothelial clusters are equivalent subtypes.

## Poster development and workflow handoff (2026-10-28)

This section supersedes the older `## Next session` list above. The poster,
figure selection and workflow graphics have progressed substantially since the
first mock-up.

### Current scientific framing

Present a cohort-wide spatial-transcriptomic survey of eight metastatic melanoma
lesions, with endothelial populations as a biologically motivated case study.
Do not claim a demonstrated vascular mechanism of ICI response. The defensible
sequence is:

1. Spatially informed BANKSY clustering and cohort-wide annotation identify
   major tumour, immune and stromal populations.
2. Cohort-consistent marker expression and vessel-like spatial morphology
   support the endothelial annotations.
3. One representative lesion illustrates preliminary cell-type neighbourhood
   structure with explicit exploratory caveats.

The current working poster was visually reviewed at:

`/home/nnataren/Documents/PhD/Conferences/HCA_2026/Poster/HCA_research_poster_2026_Nathalie_Nataren_v1 (1).pdf`

Its overall flow is coherent: biological question -> cohort/workflow -> cellular
mapping -> endothelial evidence -> preliminary spatial analysis -> future
directions. The endothelial marker/morphology panel is the strongest result.
Enlarge the neighbourhood grid and endothelial dotplots enough to read at poster
distance. Use “supports endothelial identity” rather than an unqualified
“validates”. Resolve the remaining figure-number sequence and proofreading in
the final poster.

Recommended future directions:

- Validate endothelial annotations and confirm spatial patterns across different
  neighbourhood distances.
- Resolve candidate endothelial states and test whether vascular regions are
  associated with immune access to melanoma compartments.
- Extend to a larger tissue-matched cohort to separate tissue-site effects from
  possible associations with ICI response.

### Workflow content and QC wording

The workflow now explains why discovery clustering preceded final QC:
provisional biological groups allowed transcript and cell-area distributions to
be reviewed in context, reducing disproportionate removal of rare or
morphologically distinct populations.

Concise QC wording:

- `>=15 transcripts`
- `Transcript count <99th percentile`
- `Cell area <95th percentile within each cluster`
- `<2 distinct negative-control probes`

The negative-control rule counts distinct probe identities, not total
negative-control transcripts. Cells detecting at least two distinct
negative-control probes fail this filter.

Clustree-guided resolution selection and cohort-wide marker annotation belong
inside the iterative QC/filter-selection loop. Clustree assessed stability
across resolutions and identified the final biologically meaningful split.
Annotation reviewed the full 304-gene panel across samples and prioritised
markers with consistent cluster-enriched expression across the cohort.

### Editable workflow and icon assets

Current editable workflow:

`/home/nnataren/Documents/PhD/Conferences/HCA_2026/Poster/Image_editing_for_poster/workflow/xenium_banksy_workflow_vertical_v3.svg`

The notes below were written for v2, which has since been archived to `workflow/archive/`.

It was rendered in Chrome and checked for clipping and overlaps. Its stages are:

1. H5AD counts/coordinates plus JSON configuration/metadata.
2. Pre-QC discovery BANKSY clustering.
3. QC review, Clustree resolution selection and cohort-wide annotation inside
   an iterative review loop.
4. Cell-aware QC filtering.
5. Final BANKSY reclustering at selected resolutions.
6. A standalone `Squidpy spatial analysis — Script 10` node.
7. Cohort landscape, vascular evidence and spatial-organisation outputs.

The script 10 node lists neighbourhood enrichment, centrality, co-occurrence and
Moran's I. Treat the SVG as the source of truth; confirm any nearby PNG reflects
the latest SVG before using it.

Editable vector icons in the same folder:

- `downstream_interpretation_icons.svg` — combined sheet.
- `cohort_landscape_icon.svg` — tissue map and composition bar.
- `vascular_evidence_icon.svg` — marker dotplot and vessel morphology.
- `spatial_organisation_icon.svg` — neighbourhood-enrichment grid.

Each icon retains grouped and individually editable Inkscape shapes. Final icons
represent scientific outputs rather than file formats.

### Spatial analysis selected for the poster

Script 10 produces per-sample Squidpy neighbourhood-enrichment z-score heatmaps,
centrality figures/CSVs, co-occurrence plots for present cell types, Moran's I
tables and top-gene spatial maps, plus a separate spatial-results H5AD. Moran's I
runs on a subsampled copy and is exported as tables/maps rather than retained in
the main written result object.

The cross-sample `log2(observed/expected)` grid is a later helper output, not the
direct script 10 z-score heatmap. It compares observed spatial-neighbour graph
edges with the number expected from each type's share of graph-edge ends. The
grid is symmetric; neither axis is a special index or neighbour direction.

Selected representative MG gastric panel:

- Cell-type map:
  `figures/conference_min15_area95_selected_resolutions_poster/MG_gastric_non_res_min15_area95_r0p70/cell_type_labels/MG_gastric_non_res_spatial_cell_type_labels_r0p70.svg`
- Neighbourhood montage containing the MG panel:
  `figures/spatial_analysis/vbct_min15_area95_selected_resolutions/cross_sample_nhood_enrichment/nhood_enrichment_log2_obs_exp_all_samples.svg`

Preferred heading:

> Preliminary neighbourhood-enrichment analysis of a representative melanoma
> metastasis

Preferred caption:

> **Figure X. Preliminary spatial neighbourhood analysis of a representative
> melanoma metastasis.** (A) Spatial distribution of broad cell types in the MG
> gastric lesion. (B) Cell-type neighbourhood enrichment, expressed as
> log2(observed/expected edges); red indicates enrichment and blue indicates
> depletion. Endothelial cells occupied spatially coherent regions and shared
> fewer neighbour edges with melanoma cells than expected. This exploratory
> analysis does not imply direct cellular interaction or an association with
> treatment response.

BANKSY already incorporates spatial information, so same-type enrichment is not
an independent discovery. Keep this panel explicitly preliminary.

### Repository snapshot

On 2026-10-28, `git status --short` under `Banksy_py/` showed modifications to
`10_squidpy_recluster_spatial_analysis.py`, the concise workflow README, and the
script 10 cell-type spatial-analysis README/config. It showed untracked
`immune_infiltration.json`, `nhood_enrichment_comparison.json`, and
`helper_scripts/spatial_analysis/`. Conference-figure helper/config files were
not reported as modified or staged. Recheck status before committing; do not
rely on earlier conversational staging state.

## Immediate next task: enlarge MG gastric cell-type points

Increase the plotted cell size in the MG gastric full-tissue cell-type map. It
currently looks washed out because `point_size: 2.0` makes the tissue appear
visually sparse.

Target directory:

`figures/conference_min15_area95_selected_resolutions_poster/MG_gastric_non_res_min15_area95_r0p70/cell_type_labels/`

The principal poster source is
`MG_gastric_non_res_spatial_cell_type_labels_r0p70.svg`; PNG/PDF and `_on_data`
variants are also present.

Relevant code:

- `config/conference_figures/vbct_min15_area95_selected_resolutions.json`
- `helper_scripts/conference_figures/plot_spatial_conference_figures_from_config.py`
- Shared full-plot `point_size` is currently `2.0`.
- MG overrides `full_cell_type_figsize` to `[13.8, 8.925]` but has no point-size
  override.
- `run_sample()` currently reuses one `point_size` for cluster, cell-type and
  marker-spatial plots.

Preferred fix: introduce `cell_type_point_size`, falling back to `point_size`,
and use it for the full cell-type map, optional on-data map and cell-type axis
scout. Set the larger value only for MG so marker maps and other samples remain
unchanged. Compare sizes 3, 4 and 5 at poster scale and choose the smallest that
closes the washed-out gaps without obscuring local structure or rare cell types.

The helper supports a sample-restricted run:

```bash
python helper_scripts/conference_figures/plot_spatial_conference_figures_from_config.py \
  --config config/conference_figures/vbct_min15_area95_selected_resolutions.json \
  --only-sample MG_gastric_non_res
```

That command still regenerates every plot family enabled for MG. For a focused
comparison, use a temporary config or sample-level overrides that disable marker
spatial plots, dotplots and crops. Preserve the exact script 08/09 palette and
legend order. Visually inspect the PNG, verify PDF/SVG integrity, and confirm the
SVG parses before replacing the poster panel. The mounted H5AD is read-only; do
not modify it or attempt to reconnect the mount if unavailable.

## Methods-poster pivot, cohort figures and final proofreading (2026-10-28, later session)

This section supersedes the "Spatial analysis selected for the poster" panel,
the "Recommended future directions" list and the "Immediate next task" above.
The MG point-size task was not done in this session.

### Decision: methods poster, no new spatial results

Nathalie reviewed candidate neighbourhood panels (focus pairs, graph
sensitivity, 500 µm tiles, immune infiltration) and decided they were too new
to verify before the conference. The poster is now a **methods poster** built
only from reviewed outputs: annotation, composition, marker expression and
endothelial morphology. Her reasoning: the poster is a networking tool, so
every panel must be one she can explain and defend in person. Neighbourhood
results are kept as talking points and future work, not on the poster.

Exploratory observations kept for conversation only (unverified, confounded by
tissue site and n = 1 per group):

- Endothelial cells self-associate in all 8 lesions (log2 obs/exp +2.1 to +3.8),
  are depleted next to melanoma cells, and are near-neutral with T cells.
- Melanoma–T cell segregation varies widely and does not track response.
- In each tissue-matched pair, endothelial–melanoma segregation is stronger in
  the non-responder (skin −1.0/−1.8 vs −2.7; GI −1.4 vs −2.5; brain −1.8 vs −3.6).
- The old poster Fig 4 was labelled "GI responder" but used MG_gastric_non_res
  values; that panel has been removed.

### Neighbourhood focus pairs changed

`config/spatial_analysis/vbct/nhood_enrichment_comparison.json`
`focus_pairs` now centre on endothelial cells: Endo–Endo, Endo–Melanoma,
Endo–T, Endo–Macrophage, Melanoma–T. `compare_nhood_enrichment_from_config.py`
was rerun locally. The previous focus-pair, graph-sensitivity and tile figures
and CSVs were backed up to `previous_focus_pairs_2026-10-26/` subfolders in the
`cross_sample_nhood_enrichment/` figure and table directories.

### New poster helper scripts (read-only on script 08 objects)

Both read the local script 08 annotated H5ADs via h5py and never modify them.
Script 08 was **not** rerun (it rewrites the annotated objects).

1. Cohort UMAP grid
   - Script: `helper_scripts/poster/plot_cohort_umap_grid_from_config.py`
   - Config: `config/poster/vbct/poster_cohort_umap_grid.json`
   - Outputs: `figures/conference_min15_area95_selected_resolutions_poster/cohort_umap_grid/`
     `cohort_cell_type_umap_grid_2x4` (300 dpi) and `..._4x2_150dpi`, PNG/SVG/PDF.
   - Per-sample `point_size` overrides enlarge CK_skin_res (4) and
     GR_lung_non_res_roi (2). Total 2,083,396 cells. Not used on the final
     poster (no space), but available.

2. Cohort canonical marker dot plot (poster Fig 4)
   - Script: `helper_scripts/poster/plot_cohort_marker_dotplot_from_config.py`
   - Original config (Reds): `config/poster/vbct/poster_cohort_marker_dotplot.json`
     → `cohort_canonical_marker_dotplot_150dpi.{png,svg,pdf}`
   - Poster config (RdBu_r, larger fonts):
     `config/poster/vbct/poster_cohort_marker_dotplot_rdbu_large_font.json`
     → `cohort_canonical_marker_dotplot_rdbu_large_font_150dpi.{png,svg,pdf}`
     (fontsize 16, gene_fontsize 15, figsize 11 × 4.8 in). `plot.cmap` and
     `plot.colorbar_label` are optional config keys; defaults reproduce the
     original Reds figure. RdBu_r matches poster Fig 2 (script conference
     helper, Scanpy `standard_scale="var"`).
   - Figures: `figures/conference_min15_area95_selected_resolutions_poster/cohort_marker_dotplot/`
   - Tables: `data/xenium/processed/vbct/cell_type_annotation/min15_area95_selected_resolutions/poster_cohort_marker_dotplot/`
     (`marker_expression_per_sample.csv`, `marker_expression_cohort.csv`).
   - Method: log-normalised `X`; per sample, mean expression and fraction > 0
     per cell type; cohort value = equal-weight mean over samples containing
     the type (`sample_weighting: "equal"`); min-max scaled per gene.
   - Markers: Melanoma MITF/MLANA/PMEL; Macrophage C1QA/CD163/CD68; T CD3D/CD3E/CD8A;
     B CD79A/MZB1/TNFRSF17 (MZB1/TNFRSF17 are plasma-cell genes); Fibroblast
     COL6A1/LUM/PDGFRA; Endothelial CDH5/PLVAP/VWF; Epithelial KRT18/KRT5;
     Pericyte ACTA2/RGS5/NOTCH3. Not on the panel: CD2, MS4A1, KRT8, KRT19,
     EPCAM, PECAM1, COL1A1, COL1A2, PDGFRB, MCAM, C1QB, CD14.
   - Samples per type: 8 for most, B cell 6, Epithelial 5, Pericyte 1 (GR lung,
     207 cells). RGS5 is higher in endothelial than pericyte cells; the
     pericyte row is caveated in the caption rather than dropped.

### Workflow icon

New `marker_review_icon.svg` in
`/home/nnataren/Documents/PhD/Conferences/HCA_2026/Poster/Image_editing_for_poster/workflow/`
replaces "Spatial organisation" in the downstream row (Cohort landscape ·
Vascular evidence · Marker review). Same card, CSS and layer offset as
`spatial_organisation_icon.svg`; a mini dot plot with palette row swatches and
a shaded diagonal. `spatial_organisation_icon.svg` and the combined
`downstream_interpretation_icons.svg` were not changed. The workflow's Squidpy
node is now labelled "(in development)".

### Final poster state

Latest reviewed version:
`Copy of HCA_research_poster_2026_Nathalie_Nataren_v2.pdf` (Canva export).

- Title: **A QC-aware spatial workflow for mapping cellular diversity in metastatic melanoma**.
- Banner: 8 metastatic samples | 304-gene Xenium panel | 2.08 million post-QC
  cells | 8 broad cell types.
- Fig 1 GI UMAPs/composition/spatial maps; Fig 2 endothelial marker dot plots;
  Fig 3 CDH5, endothelial map, multimodal morphology; Fig 4 cohort canonical
  marker dot plot.
- Take-home: QC-aware BANKSY workflow → consistent annotations across 2.08
  million cells, foundation for vascular–immune analysis and ICI response.

Agreed future directions (each bullet one step: map → measure → test):

- Extend the workflow with Squidpy spatial analysis (neighbourhood enrichment,
  co-occurrence, centrality and Moran's I) to map endothelial–immune–tumour
  organisation across all eight lesions.
- Quantify immune infiltration by measuring T-cell and macrophage distances to
  melanoma cells and vessels, classifying lesions as infiltrated or excluded.
- Apply the workflow to a larger tissue-matched cohort to test whether these
  spatial patterns are associated with ICI response independent of tissue site.

Outstanding proofreading flagged on v2 (for Nathalie in Canva):

- "Quatify" → "Quantify" (superseded if the new future directions are pasted).
- Fig 1 caption: full stop before "(B)"; rephrase "show that spatially coherent
  clusters, mapping onto…" → "show spatially coherent clusters that map onto…";
  last caption line collides with the Fig 2 headings.
- Fig 3B: remove "spatial1/spatial2" axis labels.
- Fig 4 caption: "0–1" with en dash; update colour description for RdBu_r.
- Reference [1]: check against PubMed (likely N Engl J Med 2025;392(1):11–22,
  online Sep 2024); use en dashes in all page ranges.
- Open: response-criteria definition; verify unique patient count (CK appears
  in skin and bowel).
