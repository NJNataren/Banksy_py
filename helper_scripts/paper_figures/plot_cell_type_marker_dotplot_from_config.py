#!/usr/bin/env python
# coding: utf-8

"""
Title: Cell-type marker dot plot by sample (paper figure) and marker specificity tables
Created: 2026-11-24
Last modified: 2026-11-24
Summary: Label every cluster in the script 03 per-cluster expression exports
with an annotation set, merge the clusters of each cell type within each sample
(cell-count weighted, which is exact because script 03 averages over all cells
including zeros), and plot a short marker list either with one row per sample x
cell type or with cohort rows per cell type (samples weighted equally). Colour
is a per-gene z-score (Seurat DotPlot style, clipped) or min-max of mean
expression, within each sample or across all rows, as set in the config.
Also writes marker specificity tables (each cell type against all other cells
in the same sample, every panel gene) so markers can be chosen from the data.
Issue #64; design notes in context/tasks/dotplot/dotplot.md.
Inputs: JSON config (config/paper_figures/vbct/); annotation set and its cohort
(config/annotations/); script 03 dotplot summary CSVs listed in the config;
canonical gene list; script 09 cell-type counts (optional cross-check).
Outputs: figure (png/svg/pdf) in
data/xenium/output/<project>/paper_figures/<cohort>/<version>/<tier>/<output_subdir>/;
tables in the same path under data/xenium/processed/:
cell_type_expression_by_sample.csv, marker_specificity_per_sample.csv,
marker_specificity_summary.csv, <output_stem>_values.csv.
Usage: python helper_scripts/paper_figures/plot_cell_type_marker_dotplot_from_config.py \
    --config config/paper_figures/vbct/cell_type_marker_dotplot.json [--validate-only]
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
import pandas as pd
from matplotlib import rc_context
from matplotlib.lines import Line2D

# The shared xenium_workflow package lives at the repository root, two levels up.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from xenium_workflow.annotation_sets import (  # noqa: E402
    add_annotation_arguments, label_column, level_vocabulary, load_mapping, require_released,
    resolve_step, step_dir,
)
from xenium_workflow.cohorts import add_root_arguments, resolve_path, select_samples  # noqa: E402


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]

EXPORT_COLUMNS = ["cluster_id", "gene", "mean_expression", "percent_expressing", "n_cells",
                  "groupby", "resolution", "expression_source"]


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog="plot cell-type marker dot plot by sample")
    parser.add_argument("--config", required=True, help="Paper figure JSON config")
    parser.add_argument("--validate-only", action="store_true",
                        help="Load, label and merge the exports and check markers, without writing anything")
    parser.add_argument("--allow-draft", action="store_true", help="Use a draft annotation set")
    add_annotation_arguments(parser)
    add_root_arguments(parser)
    return parser.parse_args()


def load_inputs(cfg, args):
    """Resolve the annotation set, sample records (config order) and output folders.

    Returns:
        context dict (annotation, level, label column, vocabulary order, folders)
        and the list of sample records with their export path.
    """
    annotation, level = resolve_step(cfg, args)
    require_released(annotation, args.allow_draft, "The paper marker dot plot")
    order, _ = level_vocabulary(annotation, level)
    exports = cfg["exports"]
    records = []
    for record in select_samples(annotation["cohort"], cfg.get("samples")):
        if record["sample"] not in exports:
            raise KeyError(f"No export listed in the config for {record['sample']}")
        records.append({**record, "export": resolve_path(args.input_root, exports[record["sample"]])})
    subdir = cfg.get("output_subdir")
    context = {
        "annotation": annotation, "level": level, "column": label_column(level), "order": order,
        "figure_dir": resolve_path(args.output_root,
                                   step_dir(annotation, "output", "paper_figures", level, subdir)),
        "table_dir": resolve_path(args.output_root,
                                  step_dir(annotation, "processed", "paper_figures", level, subdir)),
        "counts_csv": resolve_path(args.output_root,
                                   step_dir(annotation, "processed", "cell_type_proportions", level)
                                   / "cell_type_proportions.csv"),
    }
    print(f"Annotation set {annotation['version']} ({annotation['status']}), level {level}")
    print(f"Samples (row order): {[r['sample'] for r in records]}")
    return context, records


def read_sample_export(record, source_rows):
    """Read one script 03 export and keep the cohort's clustering for this sample.

    Rows are kept only for the cohort's `groupby` and selected resolution, and
    for `source_rows` (e.g. the min15_area95 pseudo-sample of the
    filter-sensitivity exports, which also hold other QC branches).

    Returns:
        DataFrame with one row per cluster x gene: cluster_id (str), gene,
        mean_expression (log-normalised, averaged over all cells), percent_expressing (0-100), n_cells.
    """
    path = record["export"]
    if not path.is_file():
        raise FileNotFoundError(f"{record['sample']}: export not found: {path}")
    df = pd.read_csv(path)
    missing = sorted(set(EXPORT_COLUMNS) - set(df.columns))
    if missing:
        raise KeyError(f"{path}: missing columns {missing}")
    keep = (df["groupby"].astype(str) == record["groupby"]) & np.isclose(
        df["resolution"].astype(float), float(record["resolution"]))
    for column, value in (source_rows or {}).items():
        keep &= df[column].astype(str) == str(value)
    df = df[keep].copy()
    if df.empty:
        raise ValueError(f"{record['sample']}: no export rows for {record['groupby']} r{record['resolution']}")

    # Biological summaries must come from clean expression, never BANKSY-expanded features.
    sources = sorted(df["expression_source"].astype(str).unique())
    if sources != ["raw"]:
        raise ValueError(f"{path}: expression_source must be 'raw' (clean expression); found {sources}")
    df["cluster_id"] = df["cluster_id"].astype(str)
    if df.duplicated(["cluster_id", "gene"]).any():
        raise ValueError(f"{path}: duplicate cluster x gene rows after filtering")
    if (df.groupby("cluster_id")["n_cells"].nunique() > 1).any():
        raise ValueError(f"{path}: a cluster has more than one n_cells value")
    # Every cluster must carry every gene, so cell counts summed per gene are cell-type totals.
    if df.groupby("cluster_id")["gene"].nunique().nunique() != 1:
        raise ValueError(f"{path}: clusters do not all have the same genes")
    return df[["cluster_id", "gene", "mean_expression", "percent_expressing", "n_cells"]]


def label_clusters(df, mapping, column, sample):
    """Add the tier label to every cluster; fail if export and mapping disagree on clusters."""
    rows = mapping[mapping["sample"] == sample].set_index("cluster_id")
    exported = set(df["cluster_id"])
    unmapped, absent = sorted(exported - set(rows.index)), sorted(set(rows.index) - exported)
    if unmapped or absent:
        raise ValueError(f"{sample}: export and annotation set disagree; unmapped clusters {unmapped}, "
                         f"mapped clusters absent from export {absent}")
    df = df.copy()
    df["cell_type"] = df["cluster_id"].map(rows[column])
    return df


def merge_to_cell_types(labelled):
    """Merge clusters into sample x cell type x gene, weighting each cluster by its cell count.

    Script 03's mean and percent are taken over all cells of a cluster
    (zeros included), so sum(n * value) / sum(n) equals the value computed
    on the pooled cells of the cell type.
    """
    weighted = labelled.assign(
        w_mean=labelled["mean_expression"] * labelled["n_cells"],
        w_pct=labelled["percent_expressing"] * labelled["n_cells"],
    )
    grouped = weighted.groupby(["sample", "cell_type", "gene"], sort=False)
    merged = grouped.agg(w_mean=("w_mean", "sum"), w_pct=("w_pct", "sum"),
                         n_cells=("n_cells", "sum"), n_clusters=("cluster_id", "nunique")).reset_index()
    merged["mean_expression"] = merged["w_mean"] / merged["n_cells"]
    merged["percent_expressing"] = merged["w_pct"] / merged["n_cells"]
    return merged.drop(columns=["w_mean", "w_pct"])


def check_counts(merged, counts_csv):
    """Cross-check merged cell counts against script 09's counts for the same set, when present."""
    if not counts_csv.is_file():
        print(f"WARNING: script 09 counts not found, skipping count cross-check: {counts_csv}")
        return
    counts = pd.read_csv(counts_csv)
    counts = counts[counts["n_cells"] > 0].set_index(["sample", "cell_type"])["n_cells"]
    ours = merged.drop_duplicates(["sample", "cell_type"]).set_index(["sample", "cell_type"])["n_cells"]
    ours = ours[ours.index.get_level_values("sample").isin(counts.index.get_level_values("sample"))]
    joined = pd.concat([ours.rename("export"), counts.rename("script09")], axis=1)
    bad = joined[joined["export"] != joined["script09"]]
    if not bad.empty:
        raise ValueError(f"Merged cell counts differ from script 09 ({counts_csv}):\n{bad}")
    print(f"Cell counts match script 09 for {len(joined)} sample x cell type rows")


def specificity_tables(merged, canonical, chosen):
    """Each cell type against all other cells of the same sample, for every gene.

    Descriptive only: no thresholds are applied and nothing is selected.
    `diff_vs_rest` is the difference in mean log-normalised expression
    (target minus the pooled rest of the sample), i.e. roughly a log fold change.

    Returns:
        per_sample: one row per sample x cell type x gene.
        summary: one row per cell type x gene across the samples containing that type,
            ranked within cell type by samples where it is top, then median diff_vs_rest.
    """
    df = merged.copy()
    keys = ["sample", "gene"]
    df["w_mean"] = df["mean_expression"] * df["n_cells"]
    df["w_pct"] = df["percent_expressing"] * df["n_cells"]
    totals = df.groupby(keys)[["w_mean", "w_pct", "n_cells"]].transform("sum")
    rest_n = totals["n_cells"] - df["n_cells"]
    # A sample with a single cell type has no "rest"; leave those values empty.
    rest_n = rest_n.where(rest_n > 0)
    df["n_cells_rest"] = rest_n
    df["mean_rest"] = (totals["w_mean"] - df["w_mean"]) / rest_n
    df["percent_rest"] = (totals["w_pct"] - df["w_pct"]) / rest_n
    df["diff_vs_rest"] = df["mean_expression"] - df["mean_rest"]

    # Highest and runner-up cell type per sample x gene, to show what a marker competes with.
    ranked = df.sort_values(keys + ["mean_expression"], ascending=[True, True, False])
    ranked["rank_in_sample"] = ranked.groupby(keys).cumcount() + 1
    df["rank_in_sample"] = ranked["rank_in_sample"]
    df["is_top"] = df["rank_in_sample"] == 1
    best = ranked[ranked["rank_in_sample"] <= 2].set_index(keys + ["rank_in_sample"])
    first = best.xs(1, level="rank_in_sample")[["cell_type", "mean_expression"]]
    second = best.xs(2, level="rank_in_sample")[["cell_type", "mean_expression"]]
    joined = df.join(first.rename(columns=lambda c: f"top_{c}"), on=keys).join(
        second.rename(columns=lambda c: f"second_{c}"), on=keys)
    # Competitor: the best other cell type (the runner-up when this row is itself top).
    own_top = joined["top_cell_type"] == joined["cell_type"]
    joined["best_other_cell_type"] = joined["second_cell_type"].where(own_top, joined["top_cell_type"])
    joined["best_other_mean"] = joined["second_mean_expression"].where(own_top, joined["top_mean_expression"])
    per_sample = joined[["sample", "cell_type", "gene", "n_cells", "n_clusters", "mean_expression",
                         "percent_expressing", "n_cells_rest", "mean_rest", "percent_rest",
                         "diff_vs_rest", "rank_in_sample", "is_top", "best_other_cell_type",
                         "best_other_mean"]].copy()

    grouped = per_sample.groupby(["cell_type", "gene"])
    summary = grouped.agg(
        n_samples_present=("sample", "nunique"),
        n_samples_top=("is_top", "sum"),
        median_diff_vs_rest=("diff_vs_rest", "median"),
        min_diff_vs_rest=("diff_vs_rest", "min"),
        median_mean_target=("mean_expression", "median"),
        median_percent_target=("percent_expressing", "median"),
        median_percent_rest=("percent_rest", "median"),
    ).reset_index()
    top_samples = per_sample[per_sample["is_top"]].groupby(["cell_type", "gene"])["sample"].agg(
        lambda s: ";".join(sorted(s)))
    summary = summary.join(top_samples.rename("samples_top"), on=["cell_type", "gene"])
    summary["samples_top"] = summary["samples_top"].fillna("")
    summary["canonical"] = summary["gene"].isin(canonical)
    summary["chosen_marker"] = [gene in chosen.get(cell_type, []) for cell_type, gene
                                in zip(summary["cell_type"], summary["gene"])]
    summary = summary.sort_values(["cell_type", "n_samples_top", "median_diff_vs_rest"],
                                  ascending=[True, False, False])
    summary.insert(2, "rank_in_cell_type", summary.groupby("cell_type").cumcount() + 1)
    return per_sample, summary.reset_index(drop=True)


POOLED_LABEL = "Cohort"


def pool_samples(merged):
    """Cohort rows: each cell type averaged over the samples containing it, every sample weighted equally.

    Equal weights stop the largest sections (e.g. MF_skin_non_res, 656k cells)
    from dominating; within each sample, clusters were already merged by cell count.
    """
    grouped = merged.groupby(["cell_type", "gene"], sort=False)
    pooled = grouped.agg(mean_expression=("mean_expression", "mean"),
                         percent_expressing=("percent_expressing", "mean"),
                         n_cells=("n_cells", "sum"), n_samples=("sample", "nunique")).reset_index()
    pooled.insert(0, "sample", POOLED_LABEL)
    return pooled


def scale_expression(values, method, within_sample, clip):
    """Add `scaled_mean`: per-gene z-score (Seurat DotPlot style) or min-max of mean expression.

    Scaling runs per gene across the plotted rows, or per gene within each
    sample when `within_sample` is true. The z-score uses the sample standard
    deviation (ddof=1, as R's scale()) and is clipped at +/-`clip`, like
    Seurat's col.min/col.max defaults (-2.5/2.5).
    """
    keys = (["sample"] if within_sample else []) + ["gene"]
    grouped = values.groupby(keys)["mean_expression"]
    if method == "zscore":
        centre, spread = grouped.transform("mean"), grouped.transform(lambda s: s.std(ddof=1))
        scaled = ((values["mean_expression"] - centre) / spread.replace(0, np.nan)).clip(-clip, clip)
    elif method == "minmax":
        lo, hi = grouped.transform("min"), grouped.transform("max")
        scaled = (values["mean_expression"] - lo) / (hi - lo).replace(0, np.nan)
    else:
        raise ValueError(f"scaling.method must be 'zscore' or 'minmax', not {method!r}")
    return values.assign(scaled_mean=scaled)


def dotplot_values(merged, records, markers, order, cfg):
    """Plotted rows and marker columns, with the configured row layout and colour scaling.

    `rows`: "sample_cell_type" (one row per sample x cell type, sample order then
    vocabulary order) or "cell_type" (cohort rows, equal sample weight).
    `scaling`: {"method": "zscore" | "minmax", "within_sample": bool, "clip": float}.
    """
    present = set(merged["cell_type"])
    types = [cell_type for cell_type in order if cell_type in present and cell_type in markers]
    genes = [gene for cell_type in types for gene in markers[cell_type]]
    layout = cfg.get("rows", "sample_cell_type")
    if layout == "cell_type":
        subset = pool_samples(merged[merged["gene"].isin(genes)])
    elif layout == "sample_cell_type":
        subset = merged[merged["gene"].isin(genes)].copy()
    else:
        raise ValueError(f"rows must be 'sample_cell_type' or 'cell_type', not {layout!r}")
    scaling = cfg.get("scaling", {"method": "minmax", "within_sample": False})
    within = bool(scaling.get("within_sample", False)) and layout == "sample_cell_type"
    subset = scale_expression(subset, scaling.get("method", "minmax"), within, float(scaling.get("clip", 2.5)))

    rank = {cell_type: i for i, cell_type in enumerate(order)}
    sample_rank = {record["sample"]: i for i, record in enumerate(records)}
    sample_rank[POOLED_LABEL] = -1
    subset = subset.sort_values(by=["sample", "cell_type"],
                                key=lambda s: s.map(sample_rank if s.name == "sample" else rank))
    keep = ["sample", "cell_type", "n_cells"] + (["n_samples"] if layout == "cell_type" else [])
    rows = subset.drop_duplicates(["sample", "cell_type"])[keep]
    return subset, rows.reset_index(drop=True), types, genes


def plot_dotplot(values, rows, types, genes, markers, records, context, cfg):
    """Draw the dot plot: rows (grouped by sample unless pooled), marker columns grouped by cell type."""
    plot_cfg = cfg["plot"]
    scaling = cfg.get("scaling", {"method": "minmax"})
    pooled = cfg.get("rows", "sample_cell_type") == "cell_type"
    if scaling.get("method", "minmax") == "zscore":
        clip = float(scaling.get("clip", 2.5))
        vmin, vmax = -clip, clip
    else:
        vmin, vmax = 0, 1
    fontsize = float(plot_cfg.get("fontsize", 9))
    max_size = float(plot_cfg.get("max_dot_size", 110))
    n_rows, n_cols = len(rows), len(genes)
    fig, ax = plt.subplots(figsize=(
        n_cols * float(plot_cfg.get("column_width_in", 0.34)) + float(plot_cfg.get("extra_width_in", 4.5)),
        n_rows * float(plot_cfg.get("row_height_in", 0.24)) + float(plot_cfg.get("extra_height_in", 2.2)),
    ))
    row_index = {(r.sample, r.cell_type): i for i, r in enumerate(rows.itertuples())}
    col_index = {gene: j for j, gene in enumerate(genes)}
    x = values["gene"].map(col_index).to_numpy()
    y = [row_index[key] for key in zip(values["sample"], values["cell_type"])]
    dots = ax.scatter(x, y, s=values["percent_expressing"].to_numpy() / 100 * max_size,
                      c=values["scaled_mean"].to_numpy(), cmap=plot_cfg.get("cmap", "RdBu_r"),
                      vmin=vmin, vmax=vmax, edgecolors="#9a9a9a", linewidths=0.4, zorder=2)

    # Shade each row's own marker block so the expected on-target dots are easy to find.
    starts, start = {}, 0
    for cell_type in types:
        starts[cell_type] = start
        start += len(markers[cell_type])
    for i, r in enumerate(rows.itertuples()):
        if r.cell_type in starts:
            ax.add_patch(plt.Rectangle((starts[r.cell_type] - 0.5, i - 0.5), len(markers[r.cell_type]), 1,
                                       facecolor="#f1efe8", edgecolor="none", zorder=0))
    for cell_type in types[1:]:
        ax.axvline(starts[cell_type] - 0.5, color="#d8d6cf", linewidth=0.8, zorder=1)

    if pooled:
        labels = [f"{r.cell_type} ({r.n_samples} sample{'s' if r.n_samples != 1 else ''}, n={r.n_cells:,})"
                  for r in rows.itertuples()]
    else:
        labels = [f"{r.cell_type} (n={r.n_cells:,})" for r in rows.itertuples()]
    ax.set_yticks(range(n_rows), labels, fontsize=fontsize)
    ax.set_xticks(range(n_cols), genes, rotation=90, fontsize=float(plot_cfg.get("gene_fontsize", fontsize)))
    ax.set_xlim(-0.6, n_cols - 0.4)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Cell-type names above their marker columns.
    for cell_type in types:
        centre = starts[cell_type] + (len(markers[cell_type]) - 1) / 2
        ax.text(centre, -0.9, cell_type, ha="left", va="bottom", rotation=40, fontsize=fontsize,
                rotation_mode="anchor")

    # Sample blocks: thin lines between samples, thicker between tissue groups; sample names
    # are placed left of the row labels once their width is known.
    by_sample = {record["sample"]: record for record in records}
    blocks = {}
    for i, sample in enumerate(rows["sample"]):
        if sample in by_sample:  # pooled cohort rows have no sample block
            blocks[sample] = (blocks.get(sample, (i, i))[0], i)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    left = min(label.get_window_extent(renderer).x0 for label in ax.get_yticklabels())
    left_axes = ax.transAxes.inverted().transform((left, 0))[0]
    blend = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
    previous_group = None
    for sample, (first, last) in blocks.items():
        record = by_sample[sample]
        group = record.get("tissue_group", record.get("tissue"))
        if first > 0:
            thick = group != previous_group
            ax.axhline(first - 0.5, color="#555555" if thick else "#bbbbbb",
                       linewidth=1.4 if thick else 0.6, zorder=1, xmin=-0.0, xmax=1.0, clip_on=False)
        previous_group = group
        ax.text(left_axes - 0.02, (first + last) / 2, f"{sample}\n{record['tissue']}, {record['response']}",
                transform=blend, ha="right", va="center", fontsize=fontsize, fontweight="bold")

    annotation = context["annotation"]
    what = "cohort, equal sample weight" if pooled else "by sample"
    ax.set_title(plot_cfg.get("title", f"Cell-type marker expression ({what}) · annotation "
                                       f"{annotation['version']} ({annotation['status']}), "
                                       f"{context['level'].replace('_', ' ')}"),
                 fontsize=fontsize + 2, pad=float(plot_cfg.get("title_pad", 90)))

    # Short plots need a taller colour bar and a lower size legend to stay readable.
    bar_height = 0.16 if n_rows > 20 else 0.45
    cax = ax.inset_axes([1.04, 0.98 - bar_height, 0.04, bar_height])
    colorbar = fig.colorbar(dots, cax=cax)
    colorbar.set_label(plot_cfg.get("colorbar_label", "Scaled mean\nexpression"), fontsize=fontsize)
    colorbar.ax.tick_params(labelsize=fontsize - 1)
    handles = [Line2D([], [], marker="o", linestyle="", markersize=np.sqrt(f * max_size),
                      markerfacecolor="#9a9a9a", markeredgecolor="none", label=f"{int(f * 100)}%")
               for f in (0.25, 0.5, 0.75, 1.0)]
    ax.legend(handles=handles, title="Cells expressing", loc="upper left",
              bbox_to_anchor=(1.02, 0.96 - bar_height - 0.08 if n_rows > 20 else 0.4),
              frameon=False, fontsize=fontsize - 1, title_fontsize=fontsize, labelspacing=1.0)
    return fig


def main():
    """Load and merge exports, write specificity tables, and render the dot plot."""
    args = parse_args()
    cfg = json.loads(Path(args.config).resolve().read_text())
    context, records = load_inputs(cfg, args)
    order, column = context["order"], context["column"]
    markers = cfg["markers"]
    unknown = sorted(set(markers) - set(order))
    if unknown:
        raise ValueError(f"Marker cell types are not labels of {context['level']}: {unknown}")

    mapping = load_mapping(context["annotation"])
    labelled = []
    for record in records:
        df = label_clusters(read_sample_export(record, cfg.get("source_rows")), mapping, column, record["sample"])
        df.insert(0, "sample", record["sample"])
        labelled.append(df)
        print(f"{record['sample']} (r{record['resolution']}): {df['cluster_id'].nunique()} clusters, "
              f"{df.drop_duplicates('cluster_id')['n_cells'].sum():,} cells, "
              f"{df['cell_type'].nunique()} cell types, {df['gene'].nunique()} genes")
    labelled = pd.concat(labelled, ignore_index=True)
    panel = set(labelled["gene"])
    absent = sorted({gene for genes in markers.values() for gene in genes} - panel)
    if absent:
        raise KeyError(f"Marker genes not in the exports: {absent}")

    merged = merge_to_cell_types(labelled)
    check_counts(merged, context["counts_csv"])
    print(f"Merged {labelled['cluster_id'].groupby(labelled['sample']).nunique().sum()} clusters into "
          f"{len(merged.drop_duplicates(['sample', 'cell_type']))} sample x cell type rows")
    if args.validate_only:
        print("Validation passed.")
        return

    canonical = set()
    if cfg.get("canonical_gene_file"):
        path = resolve_path(args.input_root, cfg["canonical_gene_file"])
        canonical = set(pd.read_csv(path)[cfg.get("canonical_gene_column", "Gene")].astype(str).str.strip())
    per_sample, summary = specificity_tables(merged, canonical, markers)
    values, rows, types, genes = dotplot_values(merged, records, markers, order, cfg)
    print(f"Rows: {cfg.get('rows', 'sample_cell_type')}; scaling: {cfg.get('scaling')}")

    table_dir, figure_dir = context["table_dir"], context["figure_dir"]
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    merged.to_csv(table_dir / "cell_type_expression_by_sample.csv", index=False)
    per_sample.to_csv(table_dir / "marker_specificity_per_sample.csv", index=False)
    summary.to_csv(table_dir / "marker_specificity_summary.csv", index=False)
    # Plotted values are named by figure, since several configs share these tables' folder.
    values.to_csv(table_dir / f"{cfg['output_stem']}_values.csv", index=False)
    print(f"Saved tables: {table_dir}")

    fig = plot_dotplot(values, rows, types, genes, markers, records, context, cfg)
    try:
        with rc_context({"svg.fonttype": "none", "pdf.fonttype": 42}):
            for fmt in cfg.get("figure_formats", ["png"]):
                path = figure_dir / f"{cfg['output_stem']}.{fmt}"
                fig.savefig(path, dpi=int(cfg.get("figure_dpi", 200)), bbox_inches="tight", facecolor="white")
                print(f"Saved figure: {path}")
    finally:
        plt.close(fig)


if __name__ == "__main__":
    main()
