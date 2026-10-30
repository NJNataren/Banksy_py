#!/usr/bin/env python
# coding: utf-8

"""
Title: Poster Dot Plot Of Canonical Markers Across The Cohort
Date: 2026-10-28
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


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]


def parse_args():
    """Parse command-line arguments for the poster marker dot plot."""
    parser = argparse.ArgumentParser(prog="plot cohort canonical marker dot plot for the poster")
    parser.add_argument("--config", required=True, help="Poster dot plot JSON config")
    parser.add_argument("--validate-only", action="store_true", help="Check inputs and genes without plotting")
    return parser.parse_args()


def resolve_path(root, value):
    """Resolve a configured absolute path or a path relative to its project root."""
    path = Path(value)
    return path if path.is_absolute() else root / path


def load_inputs(cfg, config_path):
    """Follow the poster config to script 08's config and palette; return sample records."""
    root = config_path.parents[int(cfg.get("project_root_parent_depth", 3))]
    annotation_path = resolve_path(root, cfg["annotation_config"])
    annotation = json.loads(annotation_path.read_text())
    annotation_project = annotation_path.parents[int(annotation.get("project_root_parent_depth", 3))]
    annotation_root = Path(annotation.get("output_root", annotation_project))
    palette = json.loads(resolve_path(annotation_project, annotation["palette_json"]).read_text())
    records = [
        {"sample": s["sample"], "resolution": s["resolution"],
         "h5ad": resolve_path(annotation_root, s["output_h5ad"])}
        for s in annotation["samples"]
    ]
    return root, records, palette["cell_type_order"]


def read_sample(h5ad_path, genes):
    """Read cell-type labels and the log-normalised `X` columns for `genes` only."""
    with h5py.File(h5ad_path, "r") as handle:
        labels = pd.Series(read_elem(handle["obs/cell_type"])).astype(str).to_numpy()
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
    config_path = Path(args.config).resolve()
    cfg = json.loads(config_path.read_text())
    root, records, order = load_inputs(cfg, config_path)
    markers = cfg["markers"]
    if set(markers) - set(order):
        raise ValueError(f"Marker cell types absent from palette: {sorted(set(markers) - set(order))}")
    genes = list(dict.fromkeys(gene for cell_type in order for gene in markers.get(cell_type, [])))

    rows = []
    for record in records:
        if not record["h5ad"].is_file():
            raise FileNotFoundError(f"{record['sample']}: annotated object not found: {record['h5ad']}")
        labels, expression = read_sample(record["h5ad"], genes)
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
    table_dir = resolve_path(root, cfg["table_dir"])
    figure_dir = resolve_path(root, cfg["figure_dir"])
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
