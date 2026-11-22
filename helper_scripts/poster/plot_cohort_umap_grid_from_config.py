#!/usr/bin/env python
# coding: utf-8

"""
Title: Poster Grid Of Reviewed Cell-Type UMAPs Across The Cohort
Date: 2026-10-28
Updated: 2026-11-19 - samples, palette, and mapping come from a cohort config
(xenium_workflow/cohorts.py); --input-root/--output-root default to the repository root.
Updated: 2026-11-22 - names an annotation set and label level (xenium_workflow/annotation_sets.py);
reads `obs/cell_type_<tier>`, checks the annotation stamp, derives output folders.
Summary: Read only the reviewed cell-type labels and saved UMAP embedding from
each script 08 annotated H5AD and draw one compact multi-sample grid for the
poster, using the shared cell-type palette. Panel order, titles, and styling
come from the poster config; sample IDs are replaced by tissue and response
titles with the post-QC cell count. Points are drawn in a fixed random order so
no cell type is systematically hidden under another. Inputs are never modified.
"""

import argparse
import json
import sys
import math
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rc_context
from matplotlib.lines import Line2D

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
    """Parse command-line arguments for the poster UMAP grid."""
    parser = argparse.ArgumentParser(prog="plot cohort cell-type UMAP grid for the poster")
    parser.add_argument("--config", required=True, help="Poster grid JSON config")
    parser.add_argument("--validate-only", action="store_true", help="Check inputs without plotting")
    parser.add_argument("--allow-draft", action="store_true", help="Use a draft annotation set")
    add_annotation_arguments(parser)
    add_root_arguments(parser)
    return parser.parse_args()


def load_inputs(cfg, args):
    """Read the annotation set; return context, records in poster order, label order, and colours."""
    annotation, level = resolve_step(cfg, args)
    require_released(annotation, args.allow_draft, "The poster UMAP grid")
    order, colors = level_vocabulary(annotation, level)
    records = []
    # Poster titles and point sizes come from this config; paths and keys from the set's cohort.
    for sample in select_samples(annotation["cohort"], cfg["samples"]):
        records.append({
            "sample": sample["sample"],
            "title": sample["title"],
            "h5ad": resolve_path(args.input_root, annotated_h5ad_path(annotation, sample)),
            "umap_key": sample["umap_key"],
            "point_size": sample.get("point_size"),
        })
    context = {
        "annotation": annotation, "level": level, "column": label_column(level),
        "figure_dir": resolve_path(args.output_root, step_dir(annotation, "output", "poster", level,
                                                              cfg.get("output_subdir"))),
    }
    print(f"Annotation set {annotation['version']} ({annotation['status']}), level {level}")
    return context, records, order, colors


def read_checked_labels(handle, column, annotation, h5ad_path):
    """Check the object's annotation stamp, then return the level's labels as strings."""
    found = read_elem(handle["uns/cell_type_annotation"]) if "uns/cell_type_annotation" in handle else {}
    check_stamp(found, annotation, str(h5ad_path))
    return pd.Series(read_elem(handle[f"obs/{column}"])).astype(str).to_numpy()


def read_labels_and_umap(h5ad_path, umap_key, order, context):
    """Read only the level's labels and the saved UMAP embedding from an annotated H5AD."""
    with h5py.File(h5ad_path, "r") as handle:
        labels = read_checked_labels(handle, context["column"], context["annotation"], h5ad_path)
        umap = np.asarray(handle[f"obsm/{umap_key}"][:])
    if umap.shape != (len(labels), 2) or not np.isfinite(umap).all():
        raise ValueError(f"{h5ad_path}: UMAP shape or values invalid ({umap.shape})")
    if set(labels) - set(order):
        raise ValueError(f"{h5ad_path}: unknown cell types {sorted(set(labels) - set(order))}")
    return labels, umap


def add_centroid_labels(ax, umap, labels, present, fontsize):
    """Place one haloed label at each population's median UMAP location."""
    for cell_type in present:
        x, y = np.median(umap[labels == cell_type], axis=0)
        text = ax.text(x, y, cell_type, fontsize=fontsize, fontweight="bold", ha="center", va="center")
        text.set_path_effects([path_effects.withStroke(linewidth=2.5, foreground="white")])


def plot_grid(panels, order, colors, plot_cfg, layout):
    """Draw one UMAP per sample with a shared bottom legend of the cell types present."""
    ncols = int(layout["ncols"])
    nrows = math.ceil(len(panels) / ncols)
    width, height = plot_cfg["panel_size_in"]
    legend_ncol = int(layout.get("legend_ncol", len(order)))
    legend_rows = math.ceil(len(order) / legend_ncol)
    fig_height = height * nrows + 0.45 * legend_rows + 0.35
    fig, axes = plt.subplots(nrows, ncols, figsize=(width * ncols, fig_height), squeeze=False)
    rng = np.random.default_rng(int(plot_cfg.get("shuffle_seed", 1234)))
    present_any = set()
    for ax, (record, labels, umap) in zip(axes.flat, panels):
        draw = rng.permutation(len(labels))
        # Small samples can override the shared point size so sparse panels stay readable.
        point_size = record["point_size"] or plot_cfg["point_size"]
        ax.scatter(
            umap[draw, 0], umap[draw, 1], c=pd.Series(labels[draw]).map(colors).to_numpy(),
            s=float(point_size), alpha=float(plot_cfg["alpha"]), linewidths=0, rasterized=True,
        )
        present = [label for label in order if label in set(labels)]
        present_any.update(present)
        if plot_cfg.get("label_on_data", False):
            add_centroid_labels(ax, umap, labels, present, float(plot_cfg["label_fontsize"]))
        ax.set_title(record["title"], fontsize=float(plot_cfg["title_fontsize"]), fontweight="bold", pad=2)
        ax.text(0.5, -0.02, f"{len(labels):,} cells", transform=ax.transAxes, ha="center", va="top",
                fontsize=float(plot_cfg["subtitle_fontsize"]), color="#52514e")
        ax.set_aspect("equal", adjustable="datalim")
        ax.set_axis_off()
    for ax in axes.flat[len(panels):]:
        ax.set_visible(False)

    handles = [
        Line2D([], [], marker="o", linestyle="", markersize=9, markerfacecolor=colors[label],
               markeredgecolor="none", label=label)
        for label in order if label in present_any
    ]
    fig.legend(handles=handles, loc="lower center", ncol=min(legend_ncol, len(handles)), frameon=False,
               fontsize=float(plot_cfg["legend_fontsize"]), handletextpad=0.3, columnspacing=1.2)
    # Reserve a fixed band in inches for the legend so it never overlaps the bottom row.
    bottom = (0.45 * legend_rows + 0.25) / fig_height
    fig.subplots_adjust(left=0.01, right=0.99, top=1 - 0.3 / fig_height, bottom=bottom,
                        wspace=0.05, hspace=0.18)
    return fig


def main():
    """Validate inputs and render the poster UMAP grid."""
    args = parse_args()
    cfg = json.loads(Path(args.config).resolve().read_text())
    context, records, order, colors = load_inputs(cfg, args)

    panels = []
    for record in records:
        if not record["h5ad"].is_file():
            raise FileNotFoundError(f"{record['sample']}: annotated object not found: {record['h5ad']}")
        labels, umap = read_labels_and_umap(record["h5ad"], record["umap_key"], order, context)
        print(f"{record['sample']}: {len(labels):,} cells, {len(set(labels))} cell types")
        panels.append((record, labels, umap))
    print(f"Total cells: {sum(len(labels) for _, labels, _ in panels):,}")
    if args.validate_only:
        return

    figure_dir = context["figure_dir"]
    figure_dir.mkdir(parents=True, exist_ok=True)
    for layout in cfg["layouts"]:
        fig = plot_grid(panels, order, colors, cfg["plot"], layout)
        try:
            with rc_context({"svg.fonttype": "none", "pdf.fonttype": 42}):
                for fmt in cfg.get("figure_formats", ["png"]):
                    path = figure_dir / f"{layout['output_stem']}.{fmt}"
                    fig.savefig(path, dpi=int(layout.get("dpi", 300)), bbox_inches="tight", facecolor="white")
                    print(f"Saved figure: {path}")
        finally:
            plt.close(fig)


if __name__ == "__main__":
    main()
