#!/usr/bin/env python
# coding: utf-8

"""
Title: Poster Dot Plot Of Canonical Markers Across The Cohort
Date: 2026-10-28
Updated: 2026-11-19 - samples, palette, and mapping come from a cohort config
(xenium_workflow/cohorts.py); --input-root/--output-root default to the repository root.
Updated: 2026-11-22 - names an annotation set and label level (xenium_workflow/annotation_sets.py);
reads `obs/cell_type_<tier>`, checks the annotation stamp, derives output folders.
Summary: Read reviewed cell-type labels and log-normalised expression for a
short canonical marker list from each script 08 annotated H5AD, and summarise
every cell type by mean expression and fraction of expressing cells. Samples
are weighted equally by default (the cohort value is the mean of per-sample
values over samples containing that type), so the largest sections do not
dominate. Mean expression is min-max scaled per gene across cell types, as in
Scanpy's `standard_scale="var"`. Inputs are never modified.
"""

import argparse
import json
import sys
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rc_context
from matplotlib.lines import Line2D
from scipy import sparse

try:
    from anndata.io import read_elem
except ImportError:
    from anndata.experimental import read_elem

# The shared xenium_workflow package lives at the repository root, two levels up.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from xenium_workflow.annotation_sets import (  # noqa: E402
    add_annotation_arguments, annotated_h5ad_path, check_stamp, label_column, level_vocabulary,
    require_released, resolve_step, step_dir,
)
from xenium_workflow.cohorts import add_root_arguments, resolve_path, select_samples  # noqa: E402


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]


def parse_args():
    """Parse command-line arguments for the poster marker dot plot."""
    parser = argparse.ArgumentParser(prog="plot cohort canonical marker dot plot for the poster")
    parser.add_argument("--config", required=True, help="Poster dot plot JSON config")
    parser.add_argument("--validate-only", action="store_true", help="Check inputs and genes without plotting")
    parser.add_argument("--allow-draft", action="store_true", help="Use a draft annotation set")
    add_annotation_arguments(parser)
    add_root_arguments(parser)
    return parser.parse_args()


def load_inputs(cfg, args):
    """Read the annotation set; return context (set, level, column, folders), records, and label order."""
    annotation, level = resolve_step(cfg, args)
    require_released(annotation, args.allow_draft, "The poster dot plot")
    order, _ = level_vocabulary(annotation, level)
    records = [
        {"sample": s["sample"], "resolution": s["resolution"],
         "h5ad": resolve_path(args.input_root, annotated_h5ad_path(annotation, s))}
        for s in select_samples(annotation["cohort"], cfg.get("samples"))
    ]
    subdir = cfg.get("output_subdir")
    context = {
        "annotation": annotation, "level": level, "column": label_column(level),
        "figure_dir": resolve_path(args.output_root, step_dir(annotation, "output", "poster", level, subdir)),
        "table_dir": resolve_path(args.output_root, step_dir(annotation, "processed", "poster", level, subdir)),
    }
    print(f"Annotation set {annotation['version']} ({annotation['status']}), level {level}")
    return context, records, order

def read_checked_labels(handle, column, annotation, h5ad_path):
    """Check the object's annotation stamp, then return the level's labels as strings."""
    found = read_elem(handle["uns/cell_type_annotation"]) if "uns/cell_type_annotation" in handle else {}
    check_stamp(found, annotation, str(h5ad_path))
    return pd.Series(read_elem(handle[f"obs/{column}"])).astype(str).to_numpy()



def read_sample(h5ad_path, genes, context):
    """Read the level's labels and the log-normalised `X` columns for `genes` only."""
    with h5py.File(h5ad_path, "r") as handle:
        labels = read_checked_labels(handle, context["column"], context["annotation"], h5ad_path)
        var_names = pd.Index(read_elem(handle["var"]).index.astype(str))
        group = handle["X"]
        if group.attrs.get("encoding-type") != "csr_matrix":
            raise ValueError(f"{h5ad_path}: expected CSR X, found {group.attrs.get('encoding-type')}")
        matrix = sparse.csr_matrix(
            (group["data"][:], group["indices"][:], group["indptr"][:]), shape=tuple(group.attrs["shape"]),
        )
    missing = [gene for gene in genes if gene not in var_names]
    if missing:
        raise KeyError(f"{h5ad_path}: marker genes absent from panel: {missing}")
    return labels, matrix[:, var_names.get_indexer(genes)].tocsc()


def summarise_sample(labels, expression, genes, order):
    """Per cell type: mean log-normalised expression and fraction of cells with expression > 0."""
    rows = []
    for cell_type in order:
        mask = labels == cell_type
        if not mask.any():
            continue
        subset = expression[mask]
        means = np.asarray(subset.mean(axis=0)).ravel()
        fractions = np.asarray((subset > 0).mean(axis=0)).ravel()
        for gene, mean, fraction in zip(genes, means, fractions):
            rows.append({"cell_type": cell_type, "gene": gene, "n_cells": int(mask.sum()),
                         "mean_expression": mean, "fraction_expressing": fraction})
    return rows


def cohort_summary(per_sample, weighting):
    """Combine per-sample summaries: equal sample weights, or pooled by cell number."""
    if weighting == "equal":
        grouped = per_sample.groupby(["cell_type", "gene"], sort=False)
        summary = grouped[["mean_expression", "fraction_expressing"]].mean()
    elif weighting == "cells":
        weighted = per_sample.assign(
            mean_expression=per_sample["mean_expression"] * per_sample["n_cells"],
            fraction_expressing=per_sample["fraction_expressing"] * per_sample["n_cells"],
        )
        grouped = weighted.groupby(["cell_type", "gene"], sort=False)
        summary = grouped[["mean_expression", "fraction_expressing"]].sum()
        summary = summary.div(grouped["n_cells"].sum(), axis=0)
    else:
        raise ValueError(f"sample_weighting must be 'equal' or 'cells', not {weighting!r}")
    summary["n_samples"] = per_sample.groupby(["cell_type", "gene"], sort=False)["sample"].nunique()
    summary["n_cells"] = per_sample.groupby(["cell_type", "gene"], sort=False)["n_cells"].sum()
    return summary.reset_index()


def plot_dotplot(summary, markers, order, plot_cfg):
    """Cell types as rows, marker genes grouped by their cell type as columns."""
    types = [cell_type for cell_type in order if cell_type in set(summary["cell_type"])]
    genes = [gene for cell_type in types for gene in markers.get(cell_type, [])]
    means = summary.pivot(index="cell_type", columns="gene", values="mean_expression").reindex(index=types, columns=genes)
    fractions = summary.pivot(index="cell_type", columns="gene", values="fraction_expressing").reindex(index=types, columns=genes)
    scaled = (means - means.min()) / (means.max() - means.min()).replace(0, np.nan)

    fontsize = float(plot_cfg["fontsize"])
    max_size = float(plot_cfg["max_dot_size"])
    fig, ax = plt.subplots(figsize=tuple(plot_cfg["figsize_in"]))
    x, y = np.meshgrid(np.arange(len(genes)), np.arange(len(types)))
    dots = ax.scatter(
        x.ravel(), y.ravel(), s=fractions.to_numpy().ravel() * max_size, c=scaled.to_numpy().ravel(),
        cmap=plot_cfg.get("cmap", "Reds"), vmin=0, vmax=1, edgecolors="#9a9a9a", linewidths=0.4,
    )
    ax.set_xticks(range(len(genes)), genes, rotation=90, fontsize=float(plot_cfg["gene_fontsize"]))
    ax.set_yticks(range(len(types)), types, fontsize=fontsize)
    ax.set_xlim(-0.6, len(genes) - 0.4)
    ax.set_ylim(len(types) - 0.5, -0.5)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Light boxes mark each cell type's own markers, so the expected diagonal is easy to read.
    start = 0
    for row, cell_type in enumerate(types):
        width = len(markers.get(cell_type, []))
        ax.add_patch(plt.Rectangle((start - 0.5, row - 0.5), width, 1, facecolor="#f1efe8",
                                   edgecolor="none", zorder=0))
        if start:
            ax.axvline(start - 0.5, color="#d8d6cf", linewidth=0.8, zorder=0)
        start += width

    colorbar = fig.colorbar(dots, ax=ax, fraction=0.025, pad=0.02)
    colorbar.set_label(plot_cfg.get("colorbar_label", "Scaled mean\nexpression"), fontsize=fontsize - 1)
    colorbar.ax.tick_params(labelsize=fontsize - 2)
    sizes = [0.25, 0.5, 0.75, 1.0]
    handles = [Line2D([], [], marker="o", linestyle="", markersize=np.sqrt(f * max_size),
                      markerfacecolor="#9a9a9a", markeredgecolor="none", label=f"{int(f * 100)}%") for f in sizes]
    ax.legend(handles=handles, title="Cells expressing", loc="upper left", bbox_to_anchor=(1.12, 1.02),
              frameon=False, fontsize=fontsize - 2, title_fontsize=fontsize - 1, labelspacing=1.0)
    fig.tight_layout()
    return fig


def main():
    """Validate inputs, summarise marker expression per cell type, and render the dot plot."""
    args = parse_args()
    cfg = json.loads(Path(args.config).resolve().read_text())
    context, records, order = load_inputs(cfg, args)
    markers = cfg["markers"]
    if set(markers) - set(order):
        raise ValueError(f"Marker cell types absent from palette: {sorted(set(markers) - set(order))}")
    genes = list(dict.fromkeys(gene for cell_type in order for gene in markers.get(cell_type, [])))

    rows = []
    for record in records:
        if not record["h5ad"].is_file():
            raise FileNotFoundError(f"{record['sample']}: annotated object not found: {record['h5ad']}")
        labels, expression = read_sample(record["h5ad"], genes, context)
        print(f"{record['sample']} (r{record['resolution']}): {len(labels):,} cells, {len(set(labels))} cell types")
        if args.validate_only:
            continue
        for row in summarise_sample(labels, expression, genes, order):
            rows.append({"sample": record["sample"], "resolution": record["resolution"], **row})
    if args.validate_only:
        print(f"Validation passed: {len(genes)} marker genes present in all samples.")
        return

    per_sample = pd.DataFrame(rows)
    summary = cohort_summary(per_sample, cfg.get("sample_weighting", "equal"))
    table_dir, figure_dir = context["table_dir"], context["figure_dir"]
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    per_sample.to_csv(table_dir / "marker_expression_per_sample.csv", index=False)
    summary.to_csv(table_dir / "marker_expression_cohort.csv", index=False)
    print(f"Saved tables: {table_dir}")

    fig = plot_dotplot(summary, markers, order, cfg["plot"])
    try:
        with rc_context({"svg.fonttype": "none", "pdf.fonttype": 42}):
            for fmt in cfg.get("figure_formats", ["png"]):
                path = figure_dir / f"{cfg['output_stem']}.{fmt}"
                fig.savefig(path, dpi=int(cfg.get("figure_dpi", 150)), bbox_inches="tight", facecolor="white")
                print(f"Saved figure: {path}")
    finally:
        plt.close(fig)


if __name__ == "__main__":
    main()
