# Annotation sets

An annotation set is one reviewed labelling of one cohort's clusters. Script 08
applies it; scripts 09, 10, 10a, and the spatial and poster helpers read the
objects and tables it produced. Added 2026-11-22. Design and rationale:
`context/tasks/cell_type_annotation.md` ("Design: Versioned Annotation Sets").

```text
config/annotations/<project>/<cohort>/<version>/
  annotation.json   version, cohort, parent, status (draft | released), release hash
  mapping.csv       one row per source cluster
  palette.json      per tier: label order, colours, and (tier 2+) parent labels
  CHANGES.md        what this version represents and why (written by hand)
```

`<cohort>` is the stem of the cohort config in `config/cohorts/<project>/`;
the cohort holds the clustering facts (samples, branch, resolution, script 06
objects, cluster and UMAP keys). Versions are `v1`, `v2`, ...

## Tiers

- `cell_type_tier_1`: broad cell type. Required for every cluster.
- `cell_type_tier_2` (optional): a finer label. Each tier-2 label is declared in
  `palette.json` with its tier-1 parent, and script 08 rejects a row whose
  tier-2 label sits under a different tier-1 label.
- A blank tier-2 cell means "no finer split": the cluster inherits its tier-1
  label and colour, and `obs["cell_type_tier_2_inherited"]` is True.
- Merging clusters means giving them the same label in that tier.

Script 08 writes one obs column per tier with matching `uns["<column>_colors"]`.
Later steps choose a tier with `label_level` (`tier_1` by default) or
`--label-level`; the tier is part of their output path.

`mapping.csv` columns: `sample, branch, resolution, groupby, cluster_id,
cell_type_tier_1[, cell_type_tier_2], annotation_status, notes`.
`annotation_status` is `reviewed` or `proposed`; release needs all `reviewed`.

`palette.json`:

```json
{
  "tier_1": {"order": ["Melanoma cell", "T cell"], "colors": {"Melanoma cell": "#AA4499", "T cell": "#4477AA"}},
  "tier_2": {"order": ["CD8 T cell", "CD4 T cell"],
             "colors": {"CD8 T cell": "#33528A", "CD4 T cell": "#7FA6D6"},
             "parent": {"CD8 T cell": "T cell", "CD4 T cell": "T cell"}}
}
```

Keep tier-1 colours unchanged between versions; give tier-2 labels shades of
their parent's colour.

## Making a new version

```bash
A=config/annotations/vbct/min15_area95_selected_resolutions
python -m xenium_workflow.annotation_sets new $A/v1          # creates v2 as a draft copy (parent must be released)
# edit $A/v2/mapping.csv, palette.json, and CHANGES.md
python -m xenium_workflow.annotation_sets check $A/v2        # validate (read-only)
python -m xenium_workflow.annotation_sets diff  $A/v2        # clusters whose labels changed since the parent
```

Run script 08 on the draft (locally with `--input-root` on the mount, or on
the HPC with `ANNOTATION=$A/v2`); review UMAPs, counts and
`changes_vs_v1.csv`; rerun with `--overwrite` while it is a draft. When it is
final:

```bash
python -m xenium_workflow.annotation_sets release $A/v2      # locks mapping.csv + palette.json by SHA-256
```

A released set cannot be edited: scripts refuse to load it if its files no
longer match the hash. Further changes mean another version. Script 08 objects
written from a draft remain valid after release only if the files did not
change. Scripts 09 and 10 and the helpers need a released set; `--allow-draft`
is for local tests.

All step configs name a set with `"annotation"`; `--annotation <set dir>`
points the same config at another version, so new versions need no new step
configs.

## Outputs

```text
annotated H5AD  data/xenium/processed/<project>/cell_type_annotation/<cohort>/<version>/<sample>/
                  adata_expression_clean_<sample>_<branch>_cell_type_annotated_<version>_r<res>.h5ad
08 tables       data/xenium/processed/<project>/cell_type_annotation/<cohort>/<version>/
                  cluster_annotation_counts.csv, cell_type_counts_<tier>.csv,
                  changes_vs_<parent>.csv, run manifest (JSON)
step outputs    data/xenium/{output,processed}/<project>/<step>/<cohort>/<version>/<tier>/[<output_subdir>/]
```

Each annotated H5AD carries `uns["cell_type_annotation"]` (version, parent,
status, hash, cohort, tiers, grouping, resolution, UMAP key); later steps check
it against the set they were asked to use. `tools/pipeline_status.py` checks
the highest released version by default (`--annotation`, `--label-level` to
choose).

## Sets

| Set | Status | Represents |
|---|---|---|
| `vbct/min15_area95_selected_resolutions/v1` | released 2026-11-22 | first reviewed labelling, tier 1 only (behind the preliminary results and HCA poster) |
