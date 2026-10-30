#!/usr/bin/env python
# coding: utf-8

"""
Title: Plot Immune Infiltration Distances From Config
Date: 2026-10-26
Summary: For each sample, measure the distance from every immune cell of the
configured types (e.g. T cells) to its nearest tumour cell, then compare how
much of each population lies in contact with tumour, at the tumour margin, or
distant from it. All non-tumour cells are summarised as a reference, so a
population can be judged against what tumour geometry alone would give.
Samples, response groups, and tissues come from the neighbourhood enrichment
comparison config, which follows the script 10 -> script 08 config chain.
Only labels and coordinates are read; inputs are never modified.
"""

import argparse
import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.spatial import cKDTree

from spatial_common import (
    RESPONSE_COLORS, RESPONSE_ORDER, SURFACE, TEXT_DARK, TEXT_MUTED,
    load_sample_records, read_labels_and_coords, resolve_path, save_figure,
)


REFERENCE_LABEL = "All non-tumour cells"
# Ordinal blue ramp: darkest = closest to tumour. Lighter steps fall below 2:1
# contrast on the surface, so at most three distance bins are supported.
BIN_COLORS = ["#184f95", "#2a78d6", "#86b6ef"]
LINE_STYLES = ["-", (0, (5, 2)), (0, (1.5, 1.5)), (0, (6, 2, 1.5, 2))]
BACKGROUND_POINT = "#dddcd6"
TUMOUR_POINT = "#9f9e98"


def parse_args():
    """Parse command-line arguments for the infiltration-distance analysis."""
    parser = argparse.ArgumentParser(prog="plot immune-cell distance to tumour across samples")
    parser.add_argument("--config", required=True, help="Infiltration JSON config")
    parser.add_argument("--sample", action="append", help="Restrict to a sample; repeat as needed")
    parser.add_argument("--validate-only", action="store_true", help="Check inputs without computing distances")
    parser.add_argument("--skip-maps", action="store_true", help="Skip per-sample spatial maps")
    return parser.parse_args()


def bin_labels(edges):
    """Readable distance-bin labels, e.g. ['≤20 µm', '20–100 µm', '>100 µm']."""
    labels = [f"≤{edges[0]:g} µm"]
    labels += [f"{low:g}–{high:g} µm" for low, high in zip(edges[:-1], edges[1:])]
    labels.append(f">{edges[-1]:g} µm")
    return labels


def validate_config(cfg, order):
    """Check cell types and distance bins before any data are read."""
    tumour = cfg["tumour_cell_type"]
    queries = cfg["query_cell_types"]
    if tumour not in order or set(queries) - set(order) or tumour in queries:
        raise ValueError("tumour_cell_type and query_cell_types must be distinct palette cell types")
    edges = cfg.get("distance_bins_um", [20, 100])
    if not edges or any(b <= a for a, b in zip(edges[:-1], edges[1:])) or edges[0] <= 0:
        raise ValueError("distance_bins_um must be positive and strictly increasing")
    if len(edges) + 1 > len(BIN_COLORS):
        raise ValueError(f"At most {len(BIN_COLORS) - 1} distance-bin edges are supported")
    map_type = cfg.get("map_cell_type", queries[0])
    if map_type not in queries:
        raise ValueError("map_cell_type must be one of query_cell_types")
    return tumour, queries, edges, map_type


def tumour_distances(labels, coords, tumour, groups):
    """Distance (µm) from each cell in every group to its nearest tumour cell.

    `groups` maps a group name to a boolean cell mask. Tumour cells themselves
    are never included, so distances are always to a different cell.
    """
    tree = cKDTree(coords[labels == tumour])
    return {name: tree.query(coords[mask], k=1, workers=-1)[0] for name, mask in groups.items()}


def summarise(record, distances, edges, histogram_bin, max_distance):
    """Per-group summary rows and binned distance histograms for one sample."""
    labels = bin_labels(edges)
    reference_close = np.mean(distances[REFERENCE_LABEL] <= edges[0])
    summary, histogram = [], []
    hist_edges = np.arange(0, max_distance + histogram_bin, histogram_bin)
    for group, values in distances.items():
        if len(values) == 0:
            continue
        fractions = np.bincount(np.digitize(values, edges, right=True), minlength=len(labels)) / len(values)
        row = {
            "sample": record["sample"], "response": record["response"], "tissue": record["tissue"],
            "resolution": record["resolution"], "cell_type": group, "n_cells": len(values),
            "median_distance_um": np.median(values),
            "q25_distance_um": np.quantile(values, 0.25), "q75_distance_um": np.quantile(values, 0.75),
        }
        row.update({f"fraction_{label}": fraction for label, fraction in zip(labels, fractions)})
        # >1 means this population is closer to tumour than non-tumour cells overall.
        row[f"contact_ratio_vs_reference_{labels[0]}"] = fractions[0] / reference_close if reference_close else np.nan
        summary.append(row)
        counts, _ = np.histogram(np.minimum(values, max_distance), bins=hist_edges)
        for low, high, count in zip(hist_edges[:-1], hist_edges[1:], counts):
            histogram.append({
                "sample": record["sample"], "response": record["response"], "cell_type": group,
                "bin_start_um": low, "bin_end_um": high, "n_cells": int(count),
                "includes_overflow": bool(high >= max_distance),
            })
    return summary, histogram


def ordered_records(records):
    """Responders first, then non-responders, keeping config order within each group."""
    return sorted(records, key=lambda r: RESPONSE_ORDER.index(r["response"]))


def plot_bin_fractions(summary_df, records, groups, edges, figure_dir, cfg):
    """Stacked horizontal bars of distance-bin fractions, one panel per cell group."""
    labels = bin_labels(edges)
    samples = [r for r in ordered_records(records) if r["sample"] in set(summary_df["sample"])]
    fig, axes = plt.subplots(
        1, len(groups), figsize=(3.3 * len(groups) + 1.8, 0.42 * len(samples) + 1.6),
        sharey=True, squeeze=False, layout="constrained",
    )
    y = np.arange(len(samples))[::-1]
    for ax, group in zip(axes[0], groups):
        subset = summary_df[summary_df["cell_type"] == group].set_index("sample")
        left = np.zeros(len(samples))
        for label, color in zip(labels, BIN_COLORS):
            widths = np.array([subset.at[s["sample"], f"fraction_{label}"] if s["sample"] in subset.index else 0
                               for s in samples])
            ax.barh(y, widths, left=left, color=color, edgecolor=SURFACE, linewidth=1.5, height=0.72)
            for yi, li, wi in zip(y, left, widths):
                if wi >= 0.1:
                    ink = "#ffffff" if color in BIN_COLORS[:2] else TEXT_DARK
                    ax.text(li + wi / 2, yi, f"{wi:.0%}", ha="center", va="center", fontsize=6.5, color=ink)
            left += widths
        ax.set_xlim(0, 1)
        ax.set_xticks([0, 0.5, 1], ["0%", "50%", "100%"], fontsize=8)
        ax.set_title(group, fontsize=10)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0, 0].set_yticks(y, [f"{s['sample']} · {s['response']}" for s in samples], fontsize=8)
    # Separate the response groups with a thin rule.
    n_res = sum(s["response"] == "res" for s in samples)
    if 0 < n_res < len(samples):
        for ax in axes[0]:
            ax.axhline(y[n_res - 1] - 0.5, color=TEXT_MUTED, linewidth=0.8)
    handles = [Patch(facecolor=c, label=l) for l, c in zip(labels, BIN_COLORS)]
    fig.legend(handles=handles, title="distance to nearest tumour cell", fontsize=8, title_fontsize=8,
               frameon=False, loc="outside lower center", ncol=len(labels))
    fig.suptitle("Share of each cell population by distance to tumour", fontsize=12)
    save_figure(fig, figure_dir / "infiltration_distance_bin_fractions", cfg)


def plot_distance_ecdf(distances_by_sample, records, queries, edges, max_distance, figure_dir, cfg):
    """Cumulative distance-to-tumour curves per sample, one panel per immune cell type.

    Curves converge near 100%, so end-of-line labels collide; samples are told
    apart by colour (response) plus a line style that is fixed per sample.
    """
    samples = ordered_records(records)
    styles = {}
    for group in RESPONSE_ORDER:
        members = [r["sample"] for r in samples if r["response"] == group]
        styles.update({name: LINE_STYLES[i % len(LINE_STYLES)] for i, name in enumerate(members)})
    fig, axes = plt.subplots(
        1, len(queries), figsize=(4.4 * len(queries) + 2.5, 4.2), squeeze=False, layout="constrained",
    )
    for ax, group in zip(axes[0], queries):
        for edge in edges:
            ax.axvline(edge, color="#c3c2b7", linewidth=1, linestyle="--", zorder=0)
        for record in samples:
            values = distances_by_sample.get(record["sample"], {}).get(group)
            if values is None or len(values) == 0:
                continue
            ordered = np.sort(values)
            cumulative = np.arange(1, len(ordered) + 1) / len(ordered)
            ax.plot(ordered, cumulative, color=RESPONSE_COLORS[record["response"]],
                    linestyle=styles[record["sample"]], linewidth=1.8, alpha=0.9)
        ax.set_xlim(0, max_distance)
        ax.set_ylim(0, 1.02)
        ax.set_xticks(np.linspace(0, max_distance, 6))
        ax.tick_params(labelsize=8)
        ax.set_xlabel("distance to nearest tumour cell (µm)", fontsize=9)
        ax.set_title(group, fontsize=10)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0, 0].set_ylabel("cumulative fraction of cells", fontsize=9)
    handles = [Line2D([], [], color=RESPONSE_COLORS[r["response"]], linestyle=styles[r["sample"]], linewidth=1.8,
                      label=f"{r['sample']} · {r['response']}") for r in samples]
    fig.legend(handles=handles, frameon=False, fontsize=7.5, loc="outside right center", handlelength=3)
    fig.suptitle("Distance from immune cells to tumour (steeper rise = closer to tumour)", fontsize=12)
    save_figure(fig, figure_dir / "infiltration_distance_ecdf", cfg)


def plot_spatial_map(record, labels, coords, tumour, map_type, distances, edges, figure_dir, cfg):
    """Spatial map: tumour in grey, `map_type` cells coloured by distance bin."""
    bins = np.digitize(distances, edges, right=True)
    point_size = float(cfg.get("map_point_size", 0.6))
    width, height = np.ptp(coords, axis=0)
    scale = 9 / max(width, height)
    fig, ax = plt.subplots(figsize=(width * scale + 2.2, height * scale + 0.8), layout="constrained")
    other = (labels != tumour) & (labels != map_type)
    ax.scatter(*coords[other].T, s=point_size, color=BACKGROUND_POINT, linewidth=0, rasterized=True)
    ax.scatter(*coords[labels == tumour].T, s=point_size, color=TUMOUR_POINT, linewidth=0, rasterized=True)
    query_coords = coords[labels == map_type]
    # Draw distant cells first so contact cells stay visible on top.
    for bin_index in reversed(range(len(edges) + 1)):
        mask = bins == bin_index
        ax.scatter(*query_coords[mask].T, s=point_size * 4, color=BIN_COLORS[bin_index], linewidth=0,
                   rasterized=True)
    ax.set_aspect("equal")
    ax.invert_yaxis()  # match image orientation used by Squidpy spatial plots
    ax.axis("off")
    handles = [Line2D([], [], marker="o", linestyle="", color=TUMOUR_POINT, label=tumour),
               Line2D([], [], marker="o", linestyle="", color=BACKGROUND_POINT, label="other cells")]
    handles += [Line2D([], [], marker="o", linestyle="", color=c, label=f"{map_type} {l}")
                for l, c in zip(bin_labels(edges), BIN_COLORS)]
    fig.legend(handles=handles, frameon=False, fontsize=8, loc="outside right center")
    ax.set_title(f"{record['sample']} ({record['response']}): {map_type} by distance to tumour", fontsize=11)
    save_figure(fig, figure_dir / "maps" / f"infiltration_map_{record['sample']}_{map_type.replace(' ', '_')}", cfg)


def main():
    """Compute distances to tumour, write summary tables, and render figures."""
    args = parse_args()
    config_path = Path(args.config).resolve()
    cfg = json.loads(config_path.read_text())
    root = config_path.parents[int(cfg.get("project_root_parent_depth", 4))]
    comparison_path = resolve_path(root, cfg["comparison_config"]).resolve()
    comparison = json.loads(comparison_path.read_text())
    context, records = load_sample_records(comparison, comparison_path, args.sample)
    tumour, queries, edges, map_type = validate_config(cfg, context["order"])
    max_distance = float(cfg.get("max_plot_distance_um", 500))
    histogram_bin = float(cfg.get("histogram_bin_um", 10))

    if args.validate_only:
        for record in records:
            if not record["annotated_h5ad"].is_file():
                raise FileNotFoundError(f"{record['sample']}: missing {record['annotated_h5ad']}")
            print(f"{record['sample']} ({record['response']}): {record['annotated_h5ad']}")
        print(f"Tumour: {tumour}; query types: {queries}; bins: {bin_labels(edges)}")
        return

    table_dir = resolve_path(root, cfg["table_dir"])
    figure_dir = resolve_path(root, cfg["figure_dir"])
    os.makedirs(table_dir, exist_ok=True)
    os.makedirs(figure_dir / "maps", exist_ok=True)
    summary_rows, histogram_rows, distances_by_sample = [], [], {}
    for record in records:
        labels, coords, present = read_labels_and_coords(
            record["annotated_h5ad"], context["order"], context["coord_keys"],
        )
        if tumour not in present:
            print(f"Skipping {record['sample']}: no {tumour} cells")
            continue
        groups = {group: labels == group for group in queries if group in present}
        groups[REFERENCE_LABEL] = labels != tumour
        distances = tumour_distances(labels, coords, tumour, groups)
        distances_by_sample[record["sample"]] = distances
        summary, histogram = summarise(record, distances, edges, histogram_bin, max_distance)
        summary_rows.extend(summary)
        histogram_rows.extend(histogram)
        print(f"Computed distances for {record['sample']}: " + ", ".join(
            f"{g} n={len(d):,}" for g, d in distances.items()))
        if cfg.get("spatial_maps", True) and not args.skip_maps and map_type in groups:
            plot_spatial_map(record, labels, coords, tumour, map_type, distances[map_type], edges, figure_dir, cfg)

    summary_df = pd.DataFrame(summary_rows)
    summary_csv = table_dir / "infiltration_distance_summary.csv"
    summary_df.to_csv(summary_csv, index=False)
    histogram_csv = table_dir / "infiltration_distance_histogram.csv"
    pd.DataFrame(histogram_rows).to_csv(histogram_csv, index=False)
    print(f"Saved summary table: {summary_csv}")
    print(f"Saved histogram table: {histogram_csv}")

    plot_bin_fractions(summary_df, records, queries + [REFERENCE_LABEL], edges, figure_dir, cfg)
    plot_distance_ecdf(distances_by_sample, records, queries, edges, max_distance, figure_dir, cfg)


if __name__ == "__main__":
    main()
