#!/usr/bin/env python
# coding: utf-8

"""
Title: Compare Squidpy Neighbourhood Enrichment Across Samples
Date: 2026-10-26
Summary: Collect script 10 cell-type neighbourhood enrichment matrices for all
samples, convert them to a sample-size-independent metric, and plot responder
versus non-responder comparisons. Squidpy z-scores grow with cell number, so
the comparable metric is log2(observed / expected) neighbour edges from the
edge-count matrix. Samples without a script 10 result H5AD can have enrichment
recomputed from their script 08 annotated object with the same graph settings
as script 10 and a fixed seed. Optional sensitivity analyses repeat the edge
counts for wider neighbourhoods (more neighbours or fixed radii in µm) and
within spatial tiles, to show whether a pattern depends on graph scale or on
one region of the section. Inputs are never modified.
"""

import argparse
import json
import os
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree

from spatial_common import (
    RESPONSE_COLORS, RESPONSE_ORDER, SURFACE, TEXT_DARK, TEXT_MUTED,
    load_sample_records, read_elem, read_labels_and_coords, resolve_path,
    sample_title, save_figure, spread_label_positions, tissue_markers,
)


# Blue = depleted, red = enriched, neutral grey at zero (no hue at the midpoint).
DIVERGING_CMAP = LinearSegmentedColormap.from_list(
    "depleted_enriched", ["#184f95", "#2a78d6", "#f0efec", "#e34948", "#a32a2a"]
).with_extremes(bad="#ffffff")
METRIC_LABEL = "log2 observed / expected neighbour edges"
SCRIPT10_GRAPH = {"n_neighs": 6}


def parse_args():
    """Parse command-line arguments for the cross-sample comparison."""
    parser = argparse.ArgumentParser(
        prog="compare script 10 neighbourhood enrichment across samples"
    )
    parser.add_argument("--config", required=True, help="Comparison JSON config")
    parser.add_argument("--sample", action="append", help="Restrict to a sample; repeat as needed")
    parser.add_argument("--validate-only", action="store_true", help="Check inputs and report sources without plotting")
    parser.add_argument(
        "--recompute-missing", action="store_true",
        help="Recompute enrichment from script 08 objects for samples without a script 10 result",
    )
    parser.add_argument(
        "--force-recompute", action="store_true",
        help="Ignore cached recomputed matrices (only affects samples without a script 10 result)",
    )
    return parser.parse_args()


def graph_settings(cfg):
    """Validate neighbourhood graph settings and return them keyed by name.

    Each setting has exactly one of `n_neighs` (k nearest neighbours) or
    `radius_um`. The primary graph must match script 10's 6-nearest-neighbour
    graph because its z-scores and comparison figures come from script 10.
    """
    settings = cfg.get("graph_settings", [{"name": "knn6", **SCRIPT10_GRAPH}])
    by_name = {}
    for setting in settings:
        if ("n_neighs" in setting) == ("radius_um" in setting):
            raise ValueError(f"Graph setting needs exactly one of n_neighs or radius_um: {setting}")
        if setting["name"] in by_name:
            raise ValueError(f"Duplicate graph setting name: {setting['name']}")
        by_name[setting["name"]] = setting
    primary = cfg.get("primary_graph", settings[0]["name"])
    if primary not in by_name or by_name[primary].get("n_neighs") != SCRIPT10_GRAPH["n_neighs"]:
        raise ValueError("primary_graph must name a 6-nearest-neighbour setting (script 10's graph)")
    tiling = cfg.get("tiling")
    if tiling and tiling.get("graph", primary) not in by_name:
        raise ValueError("tiling.graph must name a configured graph setting")
    return by_name, primary


def validate_comparison_config(cfg, order):
    """Check focus pairs and comparisons against the palette and configured samples."""
    for pair in cfg.get("focus_pairs", []):
        if set(pair) - set(order):
            raise ValueError(f"Focus pair uses unknown cell types: {pair}")
    for comparison in cfg.get("comparisons", []):
        if set(comparison["samples"]) - set(cfg["samples"]):
            raise ValueError(f"Comparison {comparison['name']!r} uses unconfigured samples")


def read_stored_enrichment(result_h5ad):
    """Read script 10's enrichment matrices and cell-type counts without loading `.X`."""
    with h5py.File(result_h5ad, "r") as handle:
        enrichment = read_elem(handle["uns/cell_type_nhood_enrichment"])
        labels = read_elem(handle["obs/cell_type"])
    categories = list(labels.categories)
    n_cells = pd.Series(labels).value_counts().reindex(categories, fill_value=0)
    zscore, count = np.asarray(enrichment["zscore"]), np.asarray(enrichment["count"])
    if zscore.shape != (len(categories),) * 2 or count.shape != zscore.shape:
        raise ValueError(f"{result_h5ad}: enrichment shape does not match cell_type categories")
    return zscore, count, categories, n_cells.to_numpy()


def recompute_enrichment(annotated_h5ad, order, coord_keys, seed, n_perms):
    """Recompute enrichment exactly as script 10 does, from labels and coordinates only.

    Absent cell types are dropped and present types keep palette order, matching
    script 10's categories. Only `.obs` columns are read, so memory stays small.
    """
    import anndata as ad
    import squidpy as sq

    labels, coords, present = read_labels_and_coords(annotated_h5ad, order, coord_keys)
    adata = ad.AnnData(
        obs=pd.DataFrame({"cell_type": pd.Categorical(labels, categories=present)}),
        obsm={"spatial": coords},
    )
    sq.gr.spatial_neighbors(adata, coord_type="generic")
    sq.gr.nhood_enrichment(adata, cluster_key="cell_type", seed=seed, n_perms=n_perms)
    enrichment = adata.uns["cell_type_nhood_enrichment"]
    n_cells = pd.Series(labels).value_counts().reindex(present, fill_value=0).to_numpy()
    return np.asarray(enrichment["zscore"]), np.asarray(enrichment["count"]), present, n_cells


def load_enrichment(record, cfg, args, table_dir, context):
    """Return primary-graph enrichment for one sample from script 10, the cache, or a recompute."""
    if record["result_h5ad"].is_file():
        return (*read_stored_enrichment(record["result_h5ad"]), "stored_script10")

    seed, n_perms = int(cfg.get("recompute_seed", 1234)), int(cfg.get("recompute_n_perms", 1000))
    source = f"recomputed_seed{seed}_perms{n_perms}"
    cache = table_dir / "recomputed" / f"{record['sample']}_r{record['resolution']}_{source}.npz"
    if cache.is_file() and not args.force_recompute:
        cached = np.load(cache, allow_pickle=False)
        return cached["zscore"], cached["count"], cached["categories"].tolist(), cached["n_cells"], source
    if not args.recompute_missing:
        raise FileNotFoundError(
            f"{record['sample']}: no script 10 result at {record['result_h5ad']}; "
            "rerun with --recompute-missing"
        )
    print(f"Recomputing enrichment for {record['sample']} from {record['annotated_h5ad']}")
    zscore, count, categories, n_cells = recompute_enrichment(
        record["annotated_h5ad"], context["order"], context["coord_keys"], seed, n_perms,
    )
    os.makedirs(cache.parent, exist_ok=True)
    np.savez(cache, zscore=zscore, count=count, categories=np.array(categories), n_cells=n_cells)
    print(f"Cached recomputed enrichment: {cache}")
    return zscore, count, categories, n_cells, source


def neighbour_type_counts(coords, codes, n_types, setting):
    """Count each cell's neighbours of every type under one graph setting.

    Returns an (n_cells, n_types) integer array. The kNN branch reproduces
    Squidpy's generic kNN graph (verified edge-for-edge against script 10);
    the radius branch uses one KD-tree per type so no edge list is stored.
    """
    counts = np.zeros((len(coords), n_types), dtype=np.int64)
    if "n_neighs" in setting:
        n_neighs = int(setting["n_neighs"])
        _, index = cKDTree(coords).query(coords, k=n_neighs + 1, workers=-1)
        is_self = index == np.arange(len(coords))[:, None]
        # Coincident points can push self out of the result; drop the farthest instead.
        is_self[~is_self.any(axis=1), -1] = True
        neighbour_codes = codes[index[~is_self].reshape(len(coords), n_neighs)]
        for type_code in range(n_types):
            counts[:, type_code] = (neighbour_codes == type_code).sum(axis=1)
    else:
        radius = float(setting["radius_um"])
        for type_code in range(n_types):
            tree = cKDTree(coords[codes == type_code])
            within = tree.query_ball_point(coords, radius, return_length=True, workers=-1)
            counts[:, type_code] = within - (codes == type_code)  # remove each cell itself
    return counts


def type_edge_matrix(per_cell, codes, n_types):
    """Sum per-cell neighbour counts by source cell type into a type-by-type edge matrix."""
    return np.vstack([per_cell[codes == source].sum(axis=0) for source in range(n_types)])


def log2_obs_exp(count, pseudocount):
    """Symmetrise directed edge counts and compare with edges expected from type edge shares."""
    symmetric = count + count.T
    edge_totals = symmetric.sum(axis=1)
    total = symmetric.sum()
    expected = np.outer(edge_totals, edge_totals) / total if total else np.zeros_like(symmetric, dtype=float)
    return symmetric, expected, np.log2((symmetric + pseudocount) / (expected + pseudocount))


def enrichment_to_long(record, graph, zscore, count, categories, n_cells, source, pseudocount):
    """Convert one sample's matrices for one graph to long format with log2 O/E."""
    symmetric, expected, log2_oe = log2_obs_exp(count, pseudocount)
    cells = dict(zip(categories, n_cells))
    rows = []
    for i, cell_a in enumerate(categories):
        for j, cell_b in enumerate(categories):
            rows.append({
                "sample": record["sample"], "response": record["response"],
                "tissue": record["tissue"], "branch": record["branch"],
                "resolution": record["resolution"], "graph": graph,
                "cell_type_a": cell_a, "cell_type_b": cell_b,
                "n_cells_a": int(cells[cell_a]), "n_cells_b": int(cells[cell_b]),
                "zscore": np.nan if zscore is None else zscore[i, j],
                "edge_count": int(count[i, j]), "symmetric_edge_count": int(symmetric[i, j]),
                "expected_symmetric_edges": expected[i, j], "log2_obs_exp": log2_oe[i, j],
                "source": source,
            })
    return rows


def tile_focus_pairs(record, graph, per_cell, codes, coords, present, pairs, tiling, pseudocount):
    """Compute focus-pair log2 O/E inside square spatial tiles.

    Edges are taken from the whole-section graph and assigned to the tile of
    their source cell, so cells at tile borders keep their real neighbours.
    Expected edges use each tile's own composition; a pair is reported only
    when both types have at least `min_cells_per_type` cells in the tile.
    """
    size = float(tiling.get("tile_size_um", 500))
    min_cells = int(tiling.get("min_cells_per_type", 20))
    n_types = len(present)
    tile_xy = np.floor((coords - coords.min(axis=0)) / size).astype(np.int64)
    tile_ids, tile_index = np.unique(tile_xy, axis=0, return_inverse=True)
    tile_index = tile_index.ravel()
    flat = tile_index * n_types + codes
    n_bins = len(tile_ids) * n_types
    # tile_counts[t, a, b]: edges from type-a cells in tile t to type-b neighbours.
    tile_counts = np.stack(
        [np.bincount(flat, weights=per_cell[:, b], minlength=n_bins) for b in range(n_types)], axis=1,
    ).reshape(len(tile_ids), n_types, n_types)
    tile_cells = np.bincount(flat, minlength=n_bins).reshape(len(tile_ids), n_types)

    rows = []
    usable = [(a, b) for a, b in pairs if a in present and b in present]
    for tile, (tile_x, tile_y) in enumerate(tile_ids):
        _, _, log2_oe = log2_obs_exp(tile_counts[tile], pseudocount)
        for cell_a, cell_b in usable:
            i, j = present.index(cell_a), present.index(cell_b)
            if min(tile_cells[tile, i], tile_cells[tile, j]) < min_cells:
                continue
            rows.append({
                "sample": record["sample"], "response": record["response"],
                "tissue": record["tissue"], "graph": graph, "tile_size_um": size,
                "tile_x": int(tile_x), "tile_y": int(tile_y),
                "n_cells_tile": int(tile_cells[tile].sum()),
                "cell_type_a": cell_a, "cell_type_b": cell_b,
                "n_cells_a": int(tile_cells[tile, i]), "n_cells_b": int(tile_cells[tile, j]),
                "log2_obs_exp": log2_oe[i, j],
            })
    return rows


def present_types(long_df, samples, order):
    """Palette-ordered cell types found in any of `samples`, so panels share one layout."""
    found = set(long_df.loc[long_df["sample"].isin(samples), "cell_type_a"])
    return [label for label in order if label in found]


def metric_matrix(long_df, sample, present, metric="log2_obs_exp"):
    """Return a square matrix for one sample over `present`; absent types are NaN."""
    subset = long_df[long_df["sample"] == sample]
    matrix = subset.pivot(index="cell_type_a", columns="cell_type_b", values=metric)
    return matrix.reindex(index=present, columns=present)


def symmetric_limit(matrices, off_diagonal=True):
    """Pick one symmetric colour limit; the diagonal is excluded so it cannot dominate."""
    values = []
    for matrix in matrices:
        array = matrix.to_numpy(dtype=float).copy()
        if off_diagonal:
            np.fill_diagonal(array, np.nan)
        values.append(np.abs(array[np.isfinite(array)]))
    limit = np.concatenate(values).max() if values else 1.0
    return max(np.ceil(limit * 4) / 4, 0.25)


def draw_heatmap(ax, matrix, limit, title, show_y=True):
    """Draw one annotated diverging heatmap; NaN cells (absent types) stay blank."""
    image = ax.imshow(matrix.to_numpy(dtype=float), cmap=DIVERGING_CMAP, vmin=-limit, vmax=limit)
    labels = list(matrix.index)
    ax.set_xticks(range(len(labels)), labels, rotation=90, fontsize=8)
    ax.set_yticks(range(len(labels)), labels if show_y else [], fontsize=8)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(title, fontsize=10, color=TEXT_DARK)
    for i in range(len(labels)):
        for j in range(len(labels)):
            value = matrix.iat[i, j]
            if np.isfinite(value):
                ink = "#ffffff" if abs(value) > 0.6 * limit else TEXT_DARK
                ax.text(j, i, f"{value:.1f}", ha="center", va="center", fontsize=6.5, color=ink)
    return image


def response_tissue_legend(ax, records, markers, **kwargs):
    """Legend with response colours and tissue markers (identity never colour-only)."""
    handles = [Line2D([], [], marker="o", linestyle="", color=RESPONSE_COLORS[g], label=g) for g in RESPONSE_ORDER]
    handles += [Line2D([], [], marker=m, linestyle="", color=TEXT_MUTED, label=t) for t, m in markers.items()]
    ax.legend(handles=handles, frameon=False, fontsize=8, **kwargs)


def plot_sample_grid(long_df, records, order, figure_dir, cfg):
    """Plot every sample on one shared colour scale: responders top, non-responders bottom."""
    rows = [[r for r in records if r["response"] == group] for group in RESPONSE_ORDER]
    rows = [row for row in rows if row]
    ncols = max(len(row) for row in rows)
    present = present_types(long_df, [r["sample"] for r in records], order)
    matrices = {r["sample"]: metric_matrix(long_df, r["sample"], present) for r in records}
    limit = symmetric_limit(matrices.values())

    fig, axes = plt.subplots(
        len(rows), ncols, figsize=(3.4 * ncols, 4.2 * len(rows)), squeeze=False, layout="constrained",
    )
    image = None
    for row_index, row in enumerate(rows):
        for col_index in range(ncols):
            ax = axes[row_index, col_index]
            if col_index >= len(row):
                ax.axis("off")
                continue
            record = row[col_index]
            image = draw_heatmap(ax, matrices[record["sample"]], limit, sample_title(record), show_y=col_index == 0)
    colorbar = fig.colorbar(image, ax=axes, shrink=0.6, extend="both")
    colorbar.set_label(f"{METRIC_LABEL}\n(scale set by off-diagonal pairs; diagonal may exceed it)", fontsize=8)
    fig.suptitle("Neighbourhood enrichment by response group", fontsize=12)
    save_figure(fig, figure_dir / "nhood_enrichment_log2_obs_exp_all_samples", cfg)


def plot_focus_pairs(long_df, records, pairs, figure_dir, cfg):
    """Strip plot of selected cell-type pairs: colour = response, marker = tissue."""
    markers = tissue_markers(records)
    offsets = {"res": -0.17, "non_res": 0.17}

    fig, ax = plt.subplots(figsize=(1.6 * len(pairs) + 2.5, 4.5))
    ax.axhline(0, color="#c3c2b7", linewidth=1, zorder=0)
    for pair_index, (cell_a, cell_b) in enumerate(pairs):
        for record in records:
            hit = long_df[(long_df["sample"] == record["sample"]) & (long_df["cell_type_a"] == cell_a)
                          & (long_df["cell_type_b"] == cell_b)]
            if hit.empty:
                continue
            group = [r for r in records if r["response"] == record["response"]]
            # Small deterministic spread so overlapping samples stay visible.
            spread = (group.index(record) - (len(group) - 1) / 2) * 0.06
            ax.scatter(
                pair_index + offsets[record["response"]] + spread, hit["log2_obs_exp"].iat[0],
                s=55, marker=markers[record["tissue"]], color=RESPONSE_COLORS[record["response"]],
                edgecolor=SURFACE, linewidth=1.5, zorder=3,
            )
    ax.set_xticks(range(len(pairs)), [f"{a}\n– {b}" for a, b in pairs], fontsize=8)
    ax.set_ylabel(METRIC_LABEL, fontsize=9)
    ax.tick_params(axis="y", labelsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    response_tissue_legend(ax, records, markers, loc="center left", bbox_to_anchor=(1.01, 0.5))
    ax.set_title("Focus cell-type pairs (above 0 = co-localised, below 0 = segregated)", fontsize=10)
    save_figure(fig, figure_dir / "nhood_enrichment_focus_pairs", cfg)


def plot_graph_sensitivity(long_df, records, pairs, graph_names, figure_dir, cfg):
    """One panel per focus pair: each sample's log2 O/E across neighbourhood sizes."""
    markers = tissue_markers(records)
    ncols = min(len(pairs), 3)
    nrows = int(np.ceil(len(pairs) / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(4.6 * ncols, 3.9 * nrows), squeeze=False, layout="constrained",
    )
    flat = axes.ravel()
    x = np.arange(len(graph_names))
    for ax, (cell_a, cell_b) in zip(flat, pairs):
        ax.axhline(0, color="#c3c2b7", linewidth=1, zorder=0)
        ends, end_labels = [], []
        for record in records:
            hit = long_df[(long_df["sample"] == record["sample"]) & (long_df["cell_type_a"] == cell_a)
                          & (long_df["cell_type_b"] == cell_b)].set_index("graph")
            if hit.empty:
                continue
            values = hit["log2_obs_exp"].reindex(graph_names).to_numpy(dtype=float)
            color = RESPONSE_COLORS[record["response"]]
            ax.plot(x, values, color=color, linewidth=2, alpha=0.85, zorder=2)
            ax.scatter(x, values, marker=markers[record["tissue"]], s=40, color=color,
                       edgecolor=SURFACE, linewidth=1.2, zorder=3)
            ends.append(values[-1])
            end_labels.append((record["sample"], color))
        if ends:
            span = np.ptp(ax.get_ylim())
            for y, (name, color) in zip(spread_label_positions(ends, span * 0.055), end_labels):
                ax.text(x[-1] + 0.12, y, name, fontsize=6.5, va="center", color=TEXT_DARK)
        ax.set_xticks(x, graph_names, fontsize=8, rotation=30, ha="right")
        ax.set_xlim(-0.3, len(graph_names) + 0.9)
        ax.tick_params(axis="y", labelsize=8)
        ax.set_title(f"{cell_a} – {cell_b}", fontsize=10)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    for ax in flat[len(pairs):]:
        ax.axis("off")
    for row in axes:
        row[0].set_ylabel("log2 O/E neighbour edges", fontsize=9)
    response_tissue_legend(fig, records, markers, loc="outside right center")
    fig.suptitle("Focus pairs across neighbourhood sizes (left = contact, right = wider area)", fontsize=12)
    save_figure(fig, figure_dir / "nhood_enrichment_graph_sensitivity", cfg)


def plot_tile_distributions(tile_df, sample_df, records, pairs, graph, tiling, figure_dir, cfg):
    """Per focus pair, each sample's tile-level log2 O/E with its whole-section value."""
    ncols = min(len(pairs), 3)
    nrows = int(np.ceil(len(pairs) / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(4.8 * ncols, 4.4 * nrows), squeeze=False, layout="constrained",
    )
    flat = axes.ravel()
    ordered = sorted(records, key=lambda r: RESPONSE_ORDER.index(r["response"]))
    rng = np.random.default_rng(1234)
    for ax, (cell_a, cell_b) in zip(flat, pairs):
        ax.axhline(0, color="#c3c2b7", linewidth=1, zorder=0)
        tick_labels = []
        for position, record in enumerate(ordered):
            tiles = tile_df[(tile_df["sample"] == record["sample"]) & (tile_df["cell_type_a"] == cell_a)
                            & (tile_df["cell_type_b"] == cell_b)]["log2_obs_exp"].to_numpy()
            color = RESPONSE_COLORS[record["response"]]
            if len(tiles):
                jitter = rng.uniform(-0.22, 0.22, len(tiles))
                ax.scatter(position + jitter, tiles, s=7, color=color, alpha=0.35, linewidth=0, zorder=2)
                ax.plot([position - 0.28, position + 0.28], [np.median(tiles)] * 2,
                        color=TEXT_DARK, linewidth=2, zorder=4)
            whole = sample_df[(sample_df["sample"] == record["sample"]) & (sample_df["cell_type_a"] == cell_a)
                              & (sample_df["cell_type_b"] == cell_b)]["log2_obs_exp"]
            if not whole.empty:
                ax.scatter(position, whole.iat[0], marker="D", s=36, facecolor=SURFACE,
                           edgecolor=TEXT_DARK, linewidth=1.5, zorder=5)
            tick_labels.append(f"{record['sample']}\n(n={len(tiles)})")
        ax.set_xticks(range(len(ordered)), tick_labels, fontsize=6.5, rotation=45, ha="right")
        ax.tick_params(axis="y", labelsize=8)
        ax.set_title(f"{cell_a} – {cell_b}", fontsize=10)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    for ax in flat[len(pairs):]:
        ax.axis("off")
    for row in axes:
        row[0].set_ylabel("log2 O/E neighbour edges (per tile)", fontsize=9)
    handles = [Line2D([], [], marker="o", linestyle="", color=RESPONSE_COLORS[g], label=f"{g} tiles") for g in RESPONSE_ORDER]
    handles += [Line2D([], [], color=TEXT_DARK, linewidth=2, label="tile median"),
                Line2D([], [], marker="D", linestyle="", markerfacecolor=SURFACE, markeredgecolor=TEXT_DARK,
                       label="whole section")]
    fig.legend(handles=handles, frameon=False, fontsize=8, loc="outside right center")
    size = tiling.get("tile_size_um", 500)
    min_cells = tiling.get("min_cells_per_type", 20)
    fig.suptitle(
        f"Focus pairs within {size:g} µm tiles ({graph}; tiles with ≥{min_cells} cells of both types)",
        fontsize=12,
    )
    save_figure(fig, figure_dir / f"nhood_enrichment_tiles_{graph}_{size:g}um", cfg)


def comparison_difference(long_df, comparison_samples, present):
    """Mean responder minus mean non-responder log2 O/E over types present in every sample."""
    group_means = {}
    for group in RESPONSE_ORDER:
        samples = [s for s in comparison_samples if s["response"] == group]
        stacked = [metric_matrix(long_df, s["sample"], present) for s in samples]
        group_means[group] = sum(stacked) / len(stacked)
    # Pandas aligns on labels; a type missing from any sample becomes NaN.
    return group_means["res"] - group_means["non_res"]


def comparison_samples(comparison, records):
    """Selected samples of one comparison, responders first; None if a group is missing."""
    by_name = {r["sample"]: r for r in records}
    samples = [by_name[name] for name in comparison["samples"] if name in by_name]
    if {s["response"] for s in samples} != set(RESPONSE_ORDER):
        return None
    return sorted(samples, key=lambda s: RESPONSE_ORDER.index(s["response"]))


def difference_table(long_df, comparison, samples, order, graph_names):
    """Long table of responder-minus-non-responder log2 O/E for every graph setting."""
    tables = []
    for graph in graph_names:
        graph_df = long_df[long_df["graph"] == graph]
        present = present_types(graph_df, [s["sample"] for s in samples], order)
        stacked = comparison_difference(graph_df, samples, present).stack().reset_index()
        stacked.columns = ["cell_type_a", "cell_type_b", "res_minus_non_res_log2_obs_exp"]
        stacked.insert(0, "graph", graph)
        tables.append(stacked)
    table = pd.concat(tables, ignore_index=True)
    table.insert(0, "comparison", comparison["name"])
    table.insert(1, "n_res", sum(s["response"] == "res" for s in samples))
    table.insert(2, "n_non_res", sum(s["response"] == "non_res" for s in samples))
    return table


def plot_comparison(primary_df, comparison, samples, order, figure_dir, cfg):
    """Per-sample panels plus a responder-minus-non-responder difference panel (primary graph)."""
    present = present_types(primary_df, [s["sample"] for s in samples], order)
    matrices = [metric_matrix(primary_df, s["sample"], present) for s in samples]
    difference = comparison_difference(primary_df, samples, present)
    shared = [label for label in present if difference.loc[label].notna().any()]
    difference = difference.reindex(index=shared, columns=shared)
    n_res = sum(s["response"] == "res" for s in samples)

    n_panels = len(samples) + 1
    ncols = min(n_panels, 5)
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(3.6 * ncols + 1.5, 4.2 * nrows), squeeze=False, layout="constrained",
    )
    flat = axes.ravel()
    limit = symmetric_limit(matrices)
    image = None
    for index, (record, matrix) in enumerate(zip(samples, matrices)):
        image = draw_heatmap(flat[index], matrix, limit, sample_title(record), show_y=index % ncols == 0)
    fig.colorbar(image, ax=flat[: len(samples)].tolist(), shrink=0.6, extend="both").set_label("log2 O/E", fontsize=8)

    diff_ax = flat[len(samples)]
    diff_limit = symmetric_limit([difference], off_diagonal=False)
    diff_image = draw_heatmap(
        diff_ax, difference, diff_limit,
        f"res − non_res (mean)\n{n_res} vs {len(samples) - n_res}, descriptive only",
        show_y=True,
    )
    fig.colorbar(diff_image, ax=diff_ax, shrink=0.6).set_label("difference in log2 O/E", fontsize=8)
    for ax in flat[n_panels:]:
        ax.axis("off")
    fig.suptitle(f"Neighbourhood enrichment comparison: {comparison['name']}", fontsize=12)
    save_figure(fig, figure_dir / f"nhood_enrichment_comparison_{comparison['name']}", cfg)


def main():
    """Collect enrichment matrices, write tables, and render comparison figures."""
    args = parse_args()
    config_path = Path(args.config).resolve()
    cfg = json.loads(config_path.read_text())
    context, records = load_sample_records(cfg, config_path, args.sample)
    order, root = context["order"], context["root"]
    validate_comparison_config(cfg, order)
    graphs, primary = graph_settings(cfg)
    tiling = cfg.get("tiling")
    table_dir = resolve_path(root, cfg["table_dir"])
    figure_dir = resolve_path(root, cfg["figure_dir"])
    pseudocount = float(cfg.get("pseudocount", 1))
    pairs = [tuple(pair) for pair in cfg.get("focus_pairs", [])]

    if args.validate_only:
        for record in records:
            if record["result_h5ad"].is_file():
                status = f"stored script 10 result: {record['result_h5ad']}"
            elif record["annotated_h5ad"].is_file():
                status = f"no script 10 result; recompute from {record['annotated_h5ad']}"
            else:
                raise FileNotFoundError(f"{record['sample']}: neither script 10 result nor script 08 object found")
            if (len(graphs) > 1 or tiling) and not record["annotated_h5ad"].is_file():
                raise FileNotFoundError(f"{record['sample']}: sensitivity analyses need {record['annotated_h5ad']}")
            print(f"{record['sample']} ({record['response']}, r{record['resolution']}): {status}")
        print(f"Graphs: {list(graphs)} (primary {primary}); tiling: {tiling or 'off'}")
        return

    os.makedirs(table_dir, exist_ok=True)
    os.makedirs(figure_dir, exist_ok=True)
    rows, tile_rows = [], []
    for record in records:
        zscore, count, categories, n_cells, source = load_enrichment(record, cfg, args, table_dir, context)
        print(f"Loaded {record['sample']}: {len(categories)} cell types ({source})")
        rows.extend(enrichment_to_long(record, primary, zscore, count, categories, n_cells, source, pseudocount))
        if len(graphs) == 1 and not tiling:
            continue

        labels, coords, present = read_labels_and_coords(record["annotated_h5ad"], order, context["coord_keys"])
        if present != categories:
            raise ValueError(f"{record['sample']}: annotated cell types differ from script 10 categories")
        codes = pd.Categorical(labels, categories=present).codes.astype(np.int64)
        per_cell_by_graph = {}
        for name, setting in graphs.items():
            needed = name != primary or (tiling and tiling.get("graph", primary) == name)
            if not needed:
                continue
            per_cell_by_graph[name] = neighbour_type_counts(coords, codes, len(present), setting)
            if name != primary:
                graph_count = type_edge_matrix(per_cell_by_graph[name], codes, len(present))
                rows.extend(enrichment_to_long(
                    record, name, None, graph_count, categories, n_cells, "kdtree_edge_counts", pseudocount,
                ))
        if tiling:
            tile_graph = tiling.get("graph", primary)
            tile_rows.extend(tile_focus_pairs(
                record, tile_graph, per_cell_by_graph[tile_graph], codes, coords, present, pairs, tiling, pseudocount,
            ))
        print(f"Computed sensitivity analyses for {record['sample']}")

    long_df = pd.DataFrame(rows)
    long_df.to_csv(table_dir / "nhood_enrichment_long.csv", index=False)
    print(f"Saved long-format enrichment table: {table_dir / 'nhood_enrichment_long.csv'}")
    primary_df = long_df[long_df["graph"] == primary]
    graph_names = [name for name in graphs if name in set(long_df["graph"])]

    if pairs:
        is_focus = pd.Series(list(zip(long_df["cell_type_a"], long_df["cell_type_b"]))).isin(pairs)
        long_df[is_focus.to_numpy()].to_csv(table_dir / "nhood_enrichment_focus_pairs.csv", index=False)
        plot_focus_pairs(primary_df, records, pairs, figure_dir, cfg)
        if len(graph_names) > 1:
            plot_graph_sensitivity(long_df, records, pairs, graph_names, figure_dir, cfg)
    if tile_rows:
        tile_df = pd.DataFrame(tile_rows)
        tile_graph = tiling.get("graph", primary)
        tile_csv = table_dir / f"nhood_enrichment_tiles_{tile_graph}_{tiling.get('tile_size_um', 500):g}um.csv"
        tile_df.to_csv(tile_csv, index=False)
        print(f"Saved tile table: {tile_csv}")
        plot_tile_distributions(
            tile_df, long_df[long_df["graph"] == tile_graph], records, pairs, tile_graph, tiling, figure_dir, cfg,
        )

    plot_sample_grid(primary_df, records, order, figure_dir, cfg)
    for comparison in cfg.get("comparisons", []):
        samples = comparison_samples(comparison, records)
        if samples is None:
            print(f"Skipping comparison {comparison['name']!r}: needs both response groups among selected samples")
            continue
        table = difference_table(long_df, comparison, samples, order, graph_names)
        csv_path = table_dir / f"nhood_enrichment_difference_{comparison['name']}.csv"
        table.to_csv(csv_path, index=False)
        print(f"Saved comparison table: {csv_path}")
        plot_comparison(primary_df, comparison, samples, order, figure_dir, cfg)


if __name__ == "__main__":
    main()
