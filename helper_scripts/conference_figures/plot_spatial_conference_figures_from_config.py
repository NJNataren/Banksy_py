#!/usr/bin/env python
# coding: utf-8

"""
Title: Plot Spatial Conference Figures From Config
Date: 2026-08-27
Summary: Generate presentation-oriented spatial, UMAP, marker, dotplot, and
cluster abundance figures for selected Xenium AnnData samples. Sample paths,
cluster columns, marker sets, and optional cluster-to-cell-type labels are
provided in a JSON config so the same helper can run locally or on the HPC.
"""

import argparse
import json
from pathlib import Path


DEFAULT_UMAP_KEYS = [
    "X_umap_scaled_gaussian_pc30_nc0.20",
    "X_umap_scaled_gaussian_pc35_nc0.20",
    "X_umap_scaled_gaussian_pc55_nc0.20",
    "X_umap",
]


pd = None
plt = None
sc = None


def load_plotting_stack():
    """Import plotting dependencies only when figure generation runs."""
    global pd, plt, sc

    if sc is not None:
        return

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as pyplot
    import pandas as pandas
    import scanpy as scanpy

    pd = pandas
    plt = pyplot
    sc = scanpy


def parse_args():
    """Parse command-line arguments for config-driven figure generation."""
    parser = argparse.ArgumentParser(
        description="Generate spatial conference figures from a JSON config."
    )
    parser.add_argument("--config", required=True, help="Path to JSON figure config.")
    return parser.parse_args()


def read_config(config_path):
    """Read a JSON config file."""
    with open(config_path) as handle:
        return json.load(handle)


def ensure_output_dir(path):
    """Create the output directory and return it as a Path."""
    output_dir = Path(path)
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def save_current_figure(output_dir, stem, dpi, save_pdf=True):
    """Save the current Matplotlib figure as PNG and optionally PDF."""
    png_path = output_dir / f"{stem}.png"
    plt.savefig(png_path, dpi=dpi, bbox_inches="tight")
    print(f"Wrote {png_path}")

    if save_pdf:
        pdf_path = output_dir / f"{stem}.pdf"
        plt.savefig(pdf_path, bbox_inches="tight")
        print(f"Wrote {pdf_path}")

    plt.close()


def get_sample_value(sample_cfg, cfg, key, default=None):
    """Resolve a sample-level config value with a project-level fallback."""
    return sample_cfg.get(key, cfg.get(key, default))


def prepare_adata(adata_path, cluster_col):
    """Load AnnData and add Scanpy-compatible spatial coordinates."""
    print(f"Reading {adata_path}")
    adata = sc.read_h5ad(adata_path)

    if cluster_col not in adata.obs.columns:
        raise KeyError(
            f"Cluster column {cluster_col!r} was not found in {adata_path}. "
            f"Available obs columns include: {list(adata.obs.columns[:25])}"
        )

    # Scanpy's generic embedding plot expects X_<basis>. Script 00 clean
    # expression objects usually store coordinates as obsm['spatial'] and 'xy'.
    if "X_spatial" not in adata.obsm:
        if "spatial" in adata.obsm:
            adata.obsm["X_spatial"] = adata.obsm["spatial"]
        elif "xy" in adata.obsm:
            adata.obsm["X_spatial"] = adata.obsm["xy"]
        else:
            raise KeyError("No spatial coordinates found in obsm['spatial'] or obsm['xy'].")

    adata.obs["conference_cluster"] = (
        adata.obs[cluster_col].astype(str).astype("category")
    )
    return adata


def read_label_map_from_csv(label_file, cluster_column, label_column):
    """Read a cluster-to-label mapping from CSV."""
    label_df = pd.read_csv(label_file)
    missing = [
        col for col in [cluster_column, label_column] if col not in label_df.columns
    ]
    if missing:
        raise ValueError(
            f"Label file {label_file} is missing columns {missing}. "
            f"Available columns: {list(label_df.columns)}"
        )

    label_df = label_df.dropna(subset=[cluster_column, label_column]).copy()
    label_df[cluster_column] = label_df[cluster_column].astype(str)
    label_df[label_column] = label_df[label_column].astype(str)
    return dict(zip(label_df[cluster_column], label_df[label_column]))


def read_label_map_from_json(label_config, map_key):
    """Read a cluster-to-label mapping from a JSON config key."""
    with open(label_config) as handle:
        cfg = json.load(handle)

    if map_key not in cfg:
        raise KeyError(
            f"Label config {label_config} does not contain key {map_key!r}. "
            f"Available keys: {list(cfg.keys())}"
        )

    return {str(key): str(value) for key, value in cfg[map_key].items()}


def resolve_label_map(sample_cfg):
    """Resolve optional cell-type labels from an inline map, JSON config, or CSV file."""
    label_map = sample_cfg.get("label_map")
    if label_map:
        return {str(key): str(value) for key, value in label_map.items()}

    label_config = sample_cfg.get("label_config")
    if label_config:
        return read_label_map_from_json(
            label_config,
            sample_cfg.get("label_config_map_key", "new_labels"),
        )

    label_file = sample_cfg.get("label_file")
    if not label_file:
        return None

    return read_label_map_from_csv(
        label_file,
        sample_cfg.get("label_cluster_column", "cluster_id"),
        sample_cfg.get("label_column", "cell_type_label"),
    )


def add_cell_type_labels(adata, sample_cfg):
    """Add a conference cell-type label column when labels are provided."""
    label_obs_col = sample_cfg.get("label_obs_col")
    if label_obs_col:
        if label_obs_col not in adata.obs.columns:
            raise KeyError(f"Configured label_obs_col {label_obs_col!r} not in adata.obs.")
        adata.obs["conference_cell_type"] = (
            adata.obs[label_obs_col].astype(str).astype("category")
        )
        return True

    label_map = resolve_label_map(sample_cfg)
    if not label_map:
        return False

    labels = adata.obs["conference_cluster"].astype(str).map(label_map)
    labels = labels.fillna("cluster_" + adata.obs["conference_cluster"].astype(str))
    adata.obs["conference_cell_type"] = labels.astype("category")
    return True


def available_genes(adata, requested_genes):
    """Return requested genes present in `adata`, preserving requested order."""
    present = [gene for gene in requested_genes if gene in adata.var_names]
    missing = [gene for gene in requested_genes if gene not in adata.var_names]

    if missing:
        print(f"Skipping missing marker genes: {', '.join(missing)}")
    if not present:
        raise ValueError("None of the requested marker genes were found in adata.var_names.")

    return present


def plot_spatial_category(adata, color_col, title, output_dir, stem, point_size, dpi):
    """Plot spatial coordinates colored by a categorical obs column."""
    sc.pl.embedding(
        adata,
        basis="spatial",
        color=color_col,
        size=point_size,
        frameon=False,
        legend_loc="right margin",
        title=title,
        show=False,
    )
    save_current_figure(output_dir, stem, dpi)


def plot_spatial_on_data(adata, color_col, title, output_dir, stem, point_size, dpi):
    """Plot spatial coordinates with category labels placed on the tissue."""
    sc.pl.embedding(
        adata,
        basis="spatial",
        color=color_col,
        size=point_size,
        frameon=False,
        legend_loc="on data",
        title=title,
        show=False,
    )
    save_current_figure(output_dir, stem, dpi)


def plot_umap_category(adata, color_col, title, output_dir, stem, point_size, dpi, umap_keys):
    """Plot a UMAP embedding when one of the configured keys is available."""
    for obsm_key in umap_keys:
        if obsm_key in adata.obsm:
            adata.obsm["X_conference_umap"] = adata.obsm[obsm_key]
            print(f"Using {obsm_key} for UMAP plot.")
            break
    else:
        print("Skipping UMAP plot: no configured UMAP key found.")
        return

    sc.pl.embedding(
        adata,
        basis="conference_umap",
        color=color_col,
        size=point_size,
        frameon=False,
        legend_loc="right margin",
        title=title,
        show=False,
    )
    save_current_figure(output_dir, stem, dpi)


def plot_marker_spatial_panels(adata, marker_set, output_dir, sample_name, point_size, dpi):
    """Plot spatial expression panels for one marker set."""
    markers = available_genes(adata, marker_set["genes"])
    marker_name = marker_set["name"]
    safe_name = marker_name.replace(" ", "_")

    sc.pl.embedding(
        adata,
        basis="spatial",
        color=markers,
        size=point_size,
        frameon=False,
        cmap=marker_set.get("cmap", "viridis"),
        ncols=marker_set.get("ncols", 4),
        title=markers,
        show=False,
    )
    save_current_figure(
        output_dir,
        f"{sample_name}_spatial_{safe_name}_markers",
        dpi,
    )


def plot_marker_dotplot(adata, marker_set, groupby, output_dir, sample_name, resolution, dpi):
    """Plot marker expression summarized across clusters or cell-type labels."""
    markers = available_genes(adata, marker_set["genes"])
    marker_name = marker_set["name"]
    safe_name = marker_name.replace(" ", "_")

    dotplot = sc.pl.dotplot(
        adata,
        var_names=markers,
        groupby=groupby,
        standard_scale="var",
        cmap=marker_set.get("dotplot_cmap", "RdBu_r"),
        dendrogram=False,
        show=False,
        return_fig=True,
    )
    output_path = output_dir / f"{sample_name}_{safe_name}_dotplot_{resolution}.png"
    dotplot.savefig(output_path, dpi=dpi, bbox_inches="tight")
    print(f"Wrote {output_path}")


def plot_cluster_abundance(adata, output_dir, sample_name, resolution, dpi):
    """Plot the number of cells assigned to each configured cluster."""
    counts = (
        adata.obs["conference_cluster"]
        .value_counts()
        .rename_axis("cluster")
        .reset_index(name="n_cells")
    )
    counts["cluster_sort"] = pd.to_numeric(counts["cluster"], errors="coerce")
    counts = counts.sort_values(["cluster_sort", "cluster"], na_position="last")

    fig_width = max(8, 0.35 * len(counts))
    plt.figure(figsize=(fig_width, 4.5))
    plt.bar(counts["cluster"], counts["n_cells"], color="#4c78a8")
    plt.xlabel(f"BANKSY cluster, {resolution}")
    plt.ylabel("Cells")
    plt.title(f"{sample_name} {resolution} cluster abundance")
    plt.xticks(rotation=90)
    save_current_figure(output_dir, f"{sample_name}_cluster_abundance_{resolution}", dpi)


def run_sample(sample_cfg, cfg):
    """Generate all requested conference plots for one sample."""
    sample_name = sample_cfg["sample"]
    resolution = sample_cfg.get("resolution", "selected_resolution")
    output_root = get_sample_value(sample_cfg, cfg, "output_dir", "figures/conference")
    output_subdir = sample_cfg.get("output_subdir", f"{sample_name}_{resolution}")
    output_dir = ensure_output_dir(Path(output_root) / output_subdir)
    point_size = float(get_sample_value(sample_cfg, cfg, "point_size", 3.0))
    dpi = int(get_sample_value(sample_cfg, cfg, "dpi", 300))
    umap_keys = get_sample_value(sample_cfg, cfg, "umap_keys", DEFAULT_UMAP_KEYS)
    marker_sets = get_sample_value(sample_cfg, cfg, "marker_sets", [])

    adata = prepare_adata(sample_cfg["adata_path"], sample_cfg["cluster_col"])
    has_cell_type_labels = add_cell_type_labels(adata, sample_cfg)

    print(f"Running {sample_name}: {adata.n_obs:,} cells, {adata.n_vars:,} genes.")

    plot_spatial_category(
        adata,
        "conference_cluster",
        f"{sample_name} BANKSY clusters, {resolution}",
        output_dir,
        f"{sample_name}_spatial_clusters_{resolution}",
        point_size,
        dpi,
    )
    plot_spatial_on_data(
        adata,
        "conference_cluster",
        f"{sample_name} BANKSY clusters, {resolution}",
        output_dir,
        f"{sample_name}_spatial_clusters_{resolution}_on_data",
        point_size,
        dpi,
    )
    plot_umap_category(
        adata,
        "conference_cluster",
        f"{sample_name} BANKSY UMAP, {resolution} clusters",
        output_dir,
        f"{sample_name}_umap_clusters_{resolution}",
        point_size,
        dpi,
        umap_keys,
    )

    if has_cell_type_labels:
        plot_spatial_category(
            adata,
            "conference_cell_type",
            f"{sample_name} cell-type labels, {resolution}",
            output_dir,
            f"{sample_name}_spatial_cell_type_labels_{resolution}",
            point_size,
            dpi,
        )
        plot_spatial_on_data(
            adata,
            "conference_cell_type",
            f"{sample_name} cell-type labels, {resolution}",
            output_dir,
            f"{sample_name}_spatial_cell_type_labels_{resolution}_on_data",
            point_size,
            dpi,
        )

    dotplot_groupby = "conference_cell_type" if has_cell_type_labels else "conference_cluster"
    for marker_set in marker_sets:
        plot_marker_spatial_panels(adata, marker_set, output_dir, sample_name, point_size, dpi)
        plot_marker_dotplot(adata, marker_set, dotplot_groupby, output_dir, sample_name, resolution, dpi)

    plot_cluster_abundance(adata, output_dir, sample_name, resolution, dpi)


def main():
    """Run configured conference figure generation."""
    args = parse_args()
    load_plotting_stack()
    cfg = read_config(args.config)

    for sample_cfg in cfg["samples"]:
        run_sample(sample_cfg, cfg)


if __name__ == "__main__":
    main()
