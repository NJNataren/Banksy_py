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
import re
from pathlib import Path


DEFAULT_UMAP_KEYS = [
    "X_umap_scaled_gaussian_pc30_nc0.20",
    "X_umap_scaled_gaussian_pc35_nc0.20",
    "X_umap_scaled_gaussian_pc55_nc0.20",
    "X_umap",
]


np = None
pd = None
plt = None
sc = None
# DAPI/OME-TIFF overlay support is currently disabled for conference figures.
# tifffile = None


def load_plotting_stack():
    """Import plotting dependencies only when figure generation runs."""
    global np, pd, plt, sc

    if sc is not None:
        return

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as pyplot
    import numpy as numpy
    import pandas as pandas
    import scanpy as scanpy
    # import tifffile as tiff

    np = numpy
    pd = pandas
    plt = pyplot
    sc = scanpy
    # tifffile = tiff


def parse_args():
    """Parse command-line arguments for config-driven figure generation."""
    parser = argparse.ArgumentParser(
        description="Generate spatial conference figures from a JSON config."
    )
    parser.add_argument("--config", required=True, help="Path to JSON figure config.")
    parser.add_argument(
        "--only-sample",
        action="append",
        default=[],
        help="Sample name to run. Repeat to run several samples from one config.",
    )
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


def figure_type_dir(output_dir, cfg, dirname):
    """Return the output directory for one plot family."""
    if not bool(cfg.get("organize_output_subdirs", True)):
        return output_dir
    return ensure_output_dir(Path(output_dir) / dirname)


def save_current_figure(output_dir, stem, dpi, save_pdf=True):
    """Save the current Matplotlib figure as PNG and optionally PDF."""
    output_dir = ensure_output_dir(output_dir)
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


def natural_sort_key(value):
    """Return a key that sorts embedded numbers by numeric value."""
    parts = re.split(r"(\d+)", str(value))
    return [int(part) if part.isdigit() else part.lower() for part in parts]


def sorted_category_values(values):
    """Return unique category values in natural, human-readable order."""
    categories = list(pd.Categorical(values.astype(str)).categories)
    return sorted(categories, key=natural_sort_key)


def reorder_obs_category(adata, obs_col):
    """Reorder an observation column's categories using natural sorting."""
    if obs_col in adata.obs.columns:
        categories = sorted_category_values(adata.obs[obs_col])
        adata.obs[obs_col] = adata.obs[obs_col].astype("category").cat.reorder_categories(categories)


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

    adata.obs["conference_cluster"] = adata.obs[cluster_col].astype(str).astype("category")
    reorder_obs_category(adata, "conference_cluster")
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
        labels = adata.obs[label_obs_col].astype(str)
        adata.obs["conference_cell_type_detail"] = labels.astype("category")
        adata.obs["conference_cell_type"] = labels.astype("category")
        reorder_obs_category(adata, "conference_cell_type_detail")
        reorder_obs_category(adata, "conference_cell_type")
        return True

    label_map = resolve_label_map(sample_cfg)
    if not label_map:
        return False

    labels = adata.obs["conference_cluster"].astype(str).map(label_map)
    labels = labels.fillna("cluster_" + adata.obs["conference_cluster"].astype(str))
    adata.obs["conference_cell_type_detail"] = labels.astype("category")
    adata.obs["conference_cell_type"] = labels.astype("category")
    reorder_obs_category(adata, "conference_cell_type_detail")
    reorder_obs_category(adata, "conference_cell_type")
    return True


def collapse_cell_type_labels(adata, collapse_rules):
    """Collapse detailed cell-type labels into simpler presentation categories."""
    if not collapse_rules or "conference_cell_type" not in adata.obs.columns:
        return

    labels = adata.obs["conference_cell_type"].astype(str).copy()
    for rule in collapse_rules:
        collapsed_label = rule["label"]
        contains = rule.get("contains", [])
        exact = {str(value).lower() for value in rule.get("exact", [])}

        mask = pd.Series(False, index=labels.index)
        if exact:
            mask = labels.str.lower().isin(exact)
        for pattern in contains:
            mask = mask | labels.str.contains(str(pattern), case=False, regex=False)

        labels.loc[mask] = collapsed_label

    adata.obs["conference_cell_type"] = labels.astype("category")
    reorder_obs_category(adata, "conference_cell_type")


def resolve_cell_type_palette(sample_cfg, cfg):
    """Return the configured cell-type palette with sample-level overrides."""
    palette = dict(cfg.get("cell_type_palette", {}))
    palette.update(sample_cfg.get("cell_type_palette", {}))
    return palette


def color_for_cell_type(label, palette, fallback):
    """Return a configured color for a cell-type label when available."""
    label = str(label)
    if label in palette:
        return palette[label]
    for key, color in palette.items():
        if key and key.lower() in label.lower():
            return color
    return fallback


def apply_conference_palette(adata, sample_cfg, cfg):
    """Apply configured categorical palettes to conference plotting columns."""
    palette = resolve_cell_type_palette(sample_cfg, cfg)
    if not palette:
        return

    fallback = plt.get_cmap("tab20")
    if "conference_cell_type" in adata.obs.columns:
        categories = sorted(adata.obs["conference_cell_type"].cat.categories, key=natural_sort_key)
        adata.uns["conference_cell_type_colors"] = [
            color_for_cell_type(category, palette, fallback(i % fallback.N))
            for i, category in enumerate(categories)
        ]

    if "conference_cluster" in adata.obs.columns and "conference_cell_type" in adata.obs.columns:
        cluster_categories = sorted(adata.obs["conference_cluster"].cat.categories, key=natural_sort_key)
        labels = adata.obs["conference_cell_type"].astype(str)
        clusters = adata.obs["conference_cluster"].astype(str)
        colors = []
        for i, cluster in enumerate(cluster_categories):
            cluster_labels = labels.loc[clusters == str(cluster)]
            label = cluster_labels.mode().iat[0] if not cluster_labels.empty else cluster
            colors.append(color_for_cell_type(label, palette, fallback(i % fallback.N)))
        adata.uns["conference_cluster_colors"] = colors


def available_genes(adata, requested_genes):
    """Return requested genes present in `adata`, preserving requested order."""
    present = [gene for gene in requested_genes if gene in adata.var_names]
    missing = [gene for gene in requested_genes if gene not in adata.var_names]

    if missing:
        print(f"Skipping missing marker genes: {', '.join(missing)}")
    if not present:
        raise ValueError("None of the requested marker genes were found in adata.var_names.")

    return present


def add_coordinate_axis_guides(ax):
    """Make spatial coordinate axes readable for manual crop selection."""
    ax.set_xlabel("x coordinate")
    ax.set_ylabel("y coordinate")
    ax.ticklabel_format(axis="both", style="plain", useOffset=False)
    ax.tick_params(
        axis="both",
        which="major",
        labelbottom=True,
        labelleft=True,
        bottom=True,
        left=True,
        labelsize=8,
    )
    ax.grid(True, color="0.85", linewidth=0.5, alpha=0.7)


def plot_spatial_axis_scout(adata, color_col, title, output_dir, stem, point_size, dpi):
    """Plot spatial coordinates with explicit numeric axes for crop selection."""
    coords = spatial_coordinate_frame(adata)
    values = adata.obs[color_col].astype(str)
    colors = categorical_color_map(values)

    fig, ax = plt.subplots(figsize=(10, 8))
    for category, color in colors.items():
        mask = values == category
        ax.scatter(
            coords.loc[mask, "x"],
            coords.loc[mask, "y"],
            s=point_size,
            c=[color],
            alpha=0.8,
            linewidths=0,
            label=category,
        )

    ax.set_title(title)
    add_coordinate_axis_guides(ax)
    ax.set_aspect("equal", adjustable="box")
    ax.margins(0.02)
    ax.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        markerscale=5,
        fontsize=9,
        ncol=2 if len(colors) > 12 else 1,
    )
    save_current_figure(output_dir, stem, dpi)


def highlight_obs_mask(adata, highlight_cfg):
    """Return cells matching configured label-highlight rules."""
    if not highlight_cfg:
        return None

    obs_col = highlight_cfg.get("obs_col", "conference_cell_type")
    if obs_col not in adata.obs.columns:
        return None

    labels = adata.obs[obs_col].astype(str)
    mask = pd.Series(False, index=labels.index)
    exact = {str(value).lower() for value in highlight_cfg.get("exact", [])}
    if exact:
        mask = labels.str.lower().isin(exact)
    for pattern in highlight_cfg.get("contains", []):
        mask = mask | labels.str.contains(str(pattern), case=False, regex=False)

    return mask


def color_values_for_cells(adata, color_col, mask):
    """Return colors matching Scanpy's palette for highlighted cells."""
    values = adata.obs.loc[mask, color_col]
    categories = list(adata.obs[color_col].cat.categories)
    palette = adata.uns.get(f"{color_col}_colors")
    if palette is None or len(palette) == 0:
        colors = categorical_color_map(adata.obs[color_col].astype(str))
        return [colors[str(value)] for value in values.astype(str)]

    color_by_category = dict(zip(categories, list(palette)))
    return [color_by_category[value] for value in values]


def add_spatial_highlight_overlay(adata, color_col, highlight_cfg, default_size):
    """Redraw selected spatial cells on top of a Scanpy categorical plot."""
    mask = highlight_obs_mask(adata, highlight_cfg)
    if mask is None or not bool(mask.any()):
        return

    coords = spatial_coordinate_frame(adata)
    size = float(highlight_cfg.get("point_size", default_size))
    edgecolor = highlight_cfg.get("edgecolor", "none")
    linewidth = float(highlight_cfg.get("linewidth", 0.0))
    alpha = float(highlight_cfg.get("alpha", 1.0))
    colors = color_values_for_cells(adata, color_col, mask)

    plt.gca().scatter(
        coords.loc[mask, "x"],
        coords.loc[mask, "y"],
        s=size,
        c=colors,
        alpha=alpha,
        edgecolors=edgecolor,
        linewidths=linewidth,
        zorder=20,
    )
    print(f"Highlighted {int(mask.sum()):,} cells on {color_col} plot.")


def force_categorical_legend(adata, color_col, point_size, label_col=None):
    """Rebuild a right-margin legend for categorical spatial plots."""
    from matplotlib.lines import Line2D

    ax = plt.gca()
    existing_legend = ax.get_legend()
    if existing_legend is not None:
        existing_legend.remove()

    values = adata.obs[color_col]
    categories = list(values.cat.categories) if hasattr(values, "cat") else list(pd.Categorical(values).categories)
    categories = sorted(categories, key=natural_sort_key)
    palette = adata.uns.get(f"{color_col}_colors")
    if palette is None or len(palette) == 0:
        color_map = categorical_color_map(values.astype(str))
        colors = [color_map[str(category)] for category in categories]
    else:
        colors = list(palette)[:len(categories)]

    legend_entries = list(zip([str(category) for category in categories], colors))
    if label_col and label_col in adata.obs.columns:
        label_values = adata.obs[label_col].astype(str)
        label_to_color = {}
        for category, color in zip(categories, colors):
            category_mask = values.astype(str) == str(category)
            category_labels = label_values.loc[category_mask]
            display_label = category_labels.mode().iat[0] if not category_labels.empty else str(category)
            label_to_color.setdefault(str(display_label), color)
        legend_entries = sorted(label_to_color.items(), key=lambda item: item[0].lower())

    marker_size = max(4.0, min(10.0, float(point_size) ** 0.5 * 2.0))
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=color,
            markeredgecolor="none",
            markersize=marker_size,
            label=str(label),
        )
        for label, color in legend_entries
    ]
    ax.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=8,
        ncol=2 if len(handles) > 12 else 1,
    )


def plot_spatial_category(
    adata,
    color_col,
    title,
    output_dir,
    stem,
    point_size,
    dpi,
    frameon=False,
    coordinate_axes=False,
    highlight_cfg=None,
    force_legend=False,
    legend_label_col=None,
):
    """Plot spatial coordinates colored by a categorical obs column."""
    sc.pl.embedding(
        adata,
        basis="spatial",
        color=color_col,
        size=point_size,
        frameon=frameon,
        legend_loc="right margin",
        title=title,
        show=False,
    )
    if highlight_cfg:
        add_spatial_highlight_overlay(adata, color_col, highlight_cfg, point_size)
    if force_legend:
        force_categorical_legend(adata, color_col, point_size, legend_label_col)
    if coordinate_axes:
        add_coordinate_axis_guides(plt.gca())
    save_current_figure(output_dir, stem, dpi)


def plot_spatial_on_data(
    adata, color_col, title, output_dir, stem, point_size, dpi, frameon=False
):
    """Plot spatial coordinates with category labels placed on the tissue."""
    sc.pl.embedding(
        adata,
        basis="spatial",
        color=color_col,
        size=point_size,
        frameon=frameon,
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


def plot_marker_spatial_panels(
    adata, marker_set, output_dir, sample_name, point_size, dpi, stem_suffix=""
):
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
        f"{sample_name}_spatial_{safe_name}_markers{stem_suffix}",
        dpi,
    )


def spatial_coordinate_frame(adata):
    """Return spatial x/y coordinates as a DataFrame aligned to observations."""
    coords = adata.obsm["X_spatial"]
    return pd.DataFrame(coords[:, :2], index=adata.obs_names, columns=["x", "y"])


def write_spatial_coordinate_summary(adata, summary_dir, sample_name):
    """Write min/max spatial coordinates to help choose manual crop windows."""
    coords = spatial_coordinate_frame(adata)
    summary = pd.DataFrame(
        {
            "axis": ["x", "y"],
            "min": [coords["x"].min(), coords["y"].min()],
            "max": [coords["x"].max(), coords["y"].max()],
            "p01": [coords["x"].quantile(0.01), coords["y"].quantile(0.01)],
            "p50": [coords["x"].quantile(0.50), coords["y"].quantile(0.50)],
            "p99": [coords["x"].quantile(0.99), coords["y"].quantile(0.99)],
        }
    )
    output_path = summary_dir / f"{sample_name}_spatial_coordinate_summary.csv"
    summary.to_csv(output_path, index=False)
    print(f"Wrote {output_path}")


def resolve_crop_limits(window):
    """Resolve crop x/y limits from xlim/ylim or center plus width/height."""
    if "xlim" in window and "ylim" in window:
        return tuple(window["xlim"]), tuple(window["ylim"])

    if "center" in window and "width" in window and "height" in window:
        center_x, center_y = window["center"]
        half_width = float(window["width"]) / 2
        half_height = float(window["height"]) / 2
        return (center_x - half_width, center_x + half_width), (center_y - half_height, center_y + half_height)

    raise ValueError(
        "Each crop window must define either xlim/ylim or center/width/height."
    )


def crop_adata_to_window(adata, window):
    """Return an AnnData view copied to cells inside a manual spatial crop window."""
    xlim, ylim = resolve_crop_limits(window)
    coords = spatial_coordinate_frame(adata)
    mask = (
        coords["x"].between(min(xlim), max(xlim))
        & coords["y"].between(min(ylim), max(ylim))
    )
    n_cells = int(mask.sum())
    if n_cells == 0:
        raise ValueError(f"Crop window {window.get('name', 'unnamed')!r} contains 0 cells.")

    print(f"Crop {window.get('name', 'unnamed')}: {n_cells:,} cells")
    return adata[mask].copy()


# DAPI/OME-TIFF overlay support disabled for now.
# def resolve_image_overlay(sample_cfg, cfg):
#     """Return sample-level image overlay settings with top-level fallback."""
#     return sample_cfg.get("image_overlay", cfg.get("image_overlay"))
#
#
# def crop_limits_to_pixel_bounds(window, image_cfg):
#     """Convert micron crop limits to integer image pixel bounds."""
#     xlim, ylim = resolve_crop_limits(window)
#     pixel_size_um = float(image_cfg["pixel_size_um"])
#     origin_x, origin_y = image_cfg.get("origin_um", [0.0, 0.0])
#
#     x0 = int(np.floor((min(xlim) - origin_x) / pixel_size_um))
#     x1 = int(np.ceil((max(xlim) - origin_x) / pixel_size_um))
#     y0 = int(np.floor((min(ylim) - origin_y) / pixel_size_um))
#     y1 = int(np.ceil((max(ylim) - origin_y) / pixel_size_um))
#
#     return max(0, x0), max(0, x1), max(0, y0), max(0, y1)
#
#
# def read_ome_crop_zarr(image_path, x0, x1, y0, y1, pyramid_level=0):
#     """Read a cropped OME-TIFF region through tifffile's zarr interface."""
#     import zarr
#
#     with tifffile.TiffFile(image_path) as tif:
#         store = tif.series[0].aszarr()
#         try:
#             root = zarr.open(store, mode="r")
#             arr = root
#             if not hasattr(arr, "shape") or arr.shape is None:
#                 arr = root[str(pyramid_level)]
#             return np.asarray(arr[:, y0:y1, x0:x1])
#         finally:
#             store.close()
#
#
# def read_ome_crop_memmap(image_path, x0, x1, y0, y1):
#     """Read a cropped OME-TIFF region through memory mapping when possible."""
#     arr = tifffile.memmap(image_path)
#     return np.asarray(arr[:, y0:y1, x0:x1])
#
#
# def read_ome_crop(image_cfg, window):
#     """Read and project a DAPI OME-TIFF crop for a spatial window."""
#     image_path = image_cfg["image_path"]
#     x0, x1, y0, y1 = crop_limits_to_pixel_bounds(window, image_cfg)
#     if x1 <= x0 or y1 <= y0:
#         raise ValueError(f"Invalid image crop bounds for {window.get('name', 'unnamed')!r}.")
#
#     pyramid_level = int(image_cfg.get("pyramid_level", 0))
#     try:
#         crop_stack = read_ome_crop_zarr(image_path, x0, x1, y0, y1, pyramid_level)
#     except Exception as exc:
#         print(f"zarr crop read failed ({exc}); trying tifffile.memmap fallback.")
#         crop_stack = read_ome_crop_memmap(image_path, x0, x1, y0, y1)
#
#     projection = image_cfg.get("z_projection", "max")
#     if crop_stack.ndim == 3:
#         if projection == "max":
#             image = crop_stack.max(axis=0)
#         elif projection == "mean":
#             image = crop_stack.mean(axis=0)
#         elif isinstance(projection, int):
#             image = crop_stack[projection]
#         else:
#             raise ValueError(f"Unsupported z_projection: {projection!r}")
#     elif crop_stack.ndim == 2:
#         image = crop_stack
#     else:
#         raise ValueError(f"Unsupported image crop shape: {crop_stack.shape}")
#
#     return image, (x0, x1, y0, y1)
#
#
# def scale_background_image(image, percentile_clip):
#     """Contrast-scale an image crop to 0-1 for display."""
#     low_pct, high_pct = percentile_clip
#     low, high = np.percentile(image, [low_pct, high_pct])
#     if high <= low:
#         high = image.max()
#         low = image.min()
#     if high <= low:
#         return np.zeros_like(image, dtype=float)
#     return np.clip((image.astype(float) - low) / (high - low), 0, 1)
#
#
def categorical_color_map(values):
    """Build a stable categorical color map for overlay scatter points."""
    categories = sorted(pd.Categorical(values).categories, key=natural_sort_key)
    cmap = plt.get_cmap("tab20")
    return {category: cmap(i % cmap.N) for i, category in enumerate(categories)}


# DAPI/OME-TIFF overlay plotting disabled for now.
# def prepare_background_crop(image_cfg, window):
#     """Read, contrast-scale, and return a DAPI crop for plotting."""
#     image, bounds = read_ome_crop(image_cfg, window)
#     clip = image_cfg.get("background_percentile_clip", [1, 99.8])
#     return scale_background_image(image, clip), bounds
#
#
# def plot_background_crop(image_cfg, window, output_dir, stem, title, dpi):
#     """Write a DAPI crop without cell or cluster overlays."""
#     background, _ = prepare_background_crop(image_cfg, window)
#     fig_width = float(image_cfg.get("fig_width", 8.0))
#     fig_height = fig_width * background.shape[0] / max(background.shape[1], 1)
#
#     plt.figure(figsize=(fig_width, fig_height))
#     plt.imshow(
#         background,
#         cmap=image_cfg.get("background_cmap", "gray"),
#         origin="upper",
#         alpha=float(image_cfg.get("background_alpha", 1.0)),
#     )
#     plt.title(title)
#     plt.axis("off")
#     save_current_figure(output_dir, stem, dpi)
#
#
# def plot_overlay_category(
#     adata,
#     image_cfg,
#     window,
#     color_col,
#     output_dir,
#     stem,
#     title,
#     point_size,
#     dpi,
# ):
#     """Plot crop-window cell categories over a DAPI OME-TIFF background."""
#     background, bounds = prepare_background_crop(image_cfg, window)
#     x0, x1, y0, y1 = bounds
#     pixel_size_um = float(image_cfg["pixel_size_um"])
#     origin_x, origin_y = image_cfg.get("origin_um", [0.0, 0.0])
#
#     coords = spatial_coordinate_frame(adata)
#     xs = ((coords["x"] - origin_x) / pixel_size_um) - x0
#     ys = ((coords["y"] - origin_y) / pixel_size_um) - y0
#     if image_cfg.get("invert_y", False):
#         ys = background.shape[0] - ys
#
#     values = adata.obs[color_col].astype(str)
#     colors = categorical_color_map(values)
#
#     fig_width = float(image_cfg.get("fig_width", 8.0))
#     fig_height = fig_width * background.shape[0] / max(background.shape[1], 1)
#     plt.figure(figsize=(fig_width, fig_height))
#     plt.imshow(
#         background,
#         cmap=image_cfg.get("background_cmap", "gray"),
#         origin="upper",
#         alpha=float(image_cfg.get("background_alpha", 1.0)),
#     )
#
#     for category, color in colors.items():
#         mask = values == category
#         plt.scatter(
#             xs.loc[mask],
#             ys.loc[mask],
#             s=point_size,
#             c=[color],
#             alpha=float(image_cfg.get("cell_alpha", 0.85)),
#             linewidths=0,
#             label=category,
#         )
#
#     plt.xlim(0, background.shape[1])
#     plt.ylim(background.shape[0], 0)
#     plt.title(title)
#     plt.axis("off")
#     if bool(image_cfg.get("show_legend", True)):
#         plt.legend(
#             loc="center left",
#             bbox_to_anchor=(1.02, 0.5),
#             frameon=False,
#             markerscale=float(image_cfg.get("legend_markerscale", 3.0)),
#             fontsize=float(image_cfg.get("legend_fontsize", 8.0)),
#         )
#     save_current_figure(output_dir, stem, dpi)
#
#
# def plot_image_overlays_for_crop(
#     cropped, image_cfg, window, output_dir, sample_name, resolution, name, point_size, dpi
# ):
#     """Write DAPI-only and labeled overlay plots for one crop window."""
#     plot_background_crop(
#         image_cfg,
#         window,
#         output_dir,
#         f"{sample_name}_dapi_{resolution}_{name}",
#         f"{sample_name} DAPI, {resolution}, {window['name']}",
#         dpi,
#     )
#     plot_overlay_category(
#         cropped,
#         image_cfg,
#         window,
#         "conference_cluster",
#         output_dir,
#         f"{sample_name}_spatial_clusters_{resolution}_{name}_overlay_dapi",
#         f"{sample_name} clusters on DAPI, {resolution}, {window['name']}",
#         point_size,
#         dpi,
#     )
#     if "conference_cell_type" in cropped.obs.columns:
#         plot_overlay_category(
#             cropped,
#             image_cfg,
#             window,
#             "conference_cell_type",
#             output_dir,
#             f"{sample_name}_spatial_cell_type_labels_{resolution}_{name}_overlay_dapi",
#             f"{sample_name} cell types on DAPI, {resolution}, {window['name']}",
#             point_size,
#             dpi,
#         )
#
#
def plot_crop_windows(adata, sample_cfg, cfg, crop_dir, sample_name, resolution, marker_sets):
    """Generate manual close-up plots for configured spatial crop windows."""
    crop_windows = get_sample_value(sample_cfg, cfg, "crop_windows", [])
    if not crop_windows:
        return

    default_point_size = float(get_sample_value(sample_cfg, cfg, "crop_point_size", 8.0))
    dpi = int(get_sample_value(sample_cfg, cfg, "dpi", 300))
    frameon = bool(get_sample_value(sample_cfg, cfg, "crop_show_axes", True))
    # image_cfg = resolve_image_overlay(sample_cfg, cfg)

    for window in crop_windows:
        name = window["name"].replace(" ", "_")
        cropped = crop_adata_to_window(adata, window)
        point_size = float(window.get("point_size", default_point_size))
        cluster_point_size = float(window.get("cluster_point_size", point_size))
        stem_suffix = f"_{name}"

        plot_spatial_category(
            cropped,
            "conference_cluster",
            f"{sample_name} BANKSY clusters, {resolution}, {window['name']}",
            crop_dir,
            f"{sample_name}_spatial_clusters_{resolution}{stem_suffix}",
            cluster_point_size,
            dpi,
            frameon=frameon,
            highlight_cfg=window.get("cluster_highlight"),
            force_legend=True,
            legend_label_col="conference_cell_type",
        )
        if "conference_cell_type" in cropped.obs.columns:
            plot_spatial_category(
                cropped,
                "conference_cell_type",
                f"{sample_name} cell-type labels, {resolution}, {window['name']}",
                crop_dir,
                f"{sample_name}_spatial_cell_type_labels_{resolution}{stem_suffix}",
                point_size,
                dpi,
                frameon=frameon,
                force_legend=True,
            )

        for marker_set in window.get("marker_sets", marker_sets):
            plot_marker_spatial_panels(
                cropped,
                marker_set,
                crop_dir,
                sample_name,
                point_size,
                dpi,
                stem_suffix=stem_suffix,
            )

        # DAPI/OME-TIFF overlay plotting is disabled for now.


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


def plot_cluster_abundance(adata, abundance_dir, sample_name, resolution, dpi):
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
    save_current_figure(abundance_dir, f"{sample_name}_cluster_abundance_{resolution}", dpi)


def run_sample(sample_cfg, cfg):
    """Generate all requested conference plots for one sample."""
    sample_name = sample_cfg["sample"]
    resolution = sample_cfg.get("resolution", "selected_resolution")
    output_root = get_sample_value(sample_cfg, cfg, "output_dir", "figures/conference")
    output_subdir = sample_cfg.get("output_subdir", f"{sample_name}_{resolution}")
    output_dir = ensure_output_dir(Path(output_root) / output_subdir)
    summary_dir = figure_type_dir(output_dir, cfg, "summary")
    cluster_dir = figure_type_dir(output_dir, cfg, "spatial_clusters")
    cell_type_dir = figure_type_dir(output_dir, cfg, "cell_type_labels")
    umap_dir = figure_type_dir(output_dir, cfg, "umap")
    marker_dir = figure_type_dir(output_dir, cfg, "marker_spatial")
    dotplot_dir = figure_type_dir(output_dir, cfg, "dotplots")
    crop_dir = figure_type_dir(output_dir, cfg, "crops")
    abundance_dir = figure_type_dir(output_dir, cfg, "cluster_abundance")
    point_size = float(get_sample_value(sample_cfg, cfg, "point_size", 3.0))
    dpi = int(get_sample_value(sample_cfg, cfg, "dpi", 300))
    umap_keys = get_sample_value(sample_cfg, cfg, "umap_keys", DEFAULT_UMAP_KEYS)
    marker_sets = get_sample_value(sample_cfg, cfg, "marker_sets", [])

    adata = prepare_adata(sample_cfg["adata_path"], sample_cfg["cluster_col"])
    has_cell_type_labels = add_cell_type_labels(adata, sample_cfg)
    collapse_rules = get_sample_value(sample_cfg, cfg, "label_collapse_rules", [])
    collapse_cell_type_labels(adata, collapse_rules)
    apply_conference_palette(adata, sample_cfg, cfg)

    print(f"Running {sample_name}: {adata.n_obs:,} cells, {adata.n_vars:,} genes.")
    if bool(get_sample_value(sample_cfg, cfg, "write_coordinate_summary", True)):
        write_spatial_coordinate_summary(adata, summary_dir, sample_name)

    plot_spatial_category(
        adata,
        "conference_cluster",
        f"{sample_name} BANKSY clusters, {resolution}",
        cluster_dir,
        f"{sample_name}_spatial_clusters_{resolution}",
        point_size,
        dpi,
    )
    plot_spatial_on_data(
        adata,
        "conference_cluster",
        f"{sample_name} BANKSY clusters, {resolution}",
        cluster_dir,
        f"{sample_name}_spatial_clusters_{resolution}_on_data",
        point_size,
        dpi,
    )
    if bool(get_sample_value(sample_cfg, cfg, "write_axis_scout_plots", True)):
        plot_spatial_axis_scout(
            adata,
            "conference_cluster",
            f"{sample_name} BANKSY clusters, {resolution}, axis scout",
            cluster_dir,
            f"{sample_name}_spatial_clusters_{resolution}_axis_scout",
            point_size,
            dpi,
        )
    plot_umap_category(
        adata,
        "conference_cluster",
        f"{sample_name} BANKSY UMAP, {resolution} clusters",
        umap_dir,
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
            cell_type_dir,
            f"{sample_name}_spatial_cell_type_labels_{resolution}",
            point_size,
            dpi,
            force_legend=True,
        )
        plot_spatial_on_data(
            adata,
            "conference_cell_type",
            f"{sample_name} cell-type labels, {resolution}",
            cell_type_dir,
            f"{sample_name}_spatial_cell_type_labels_{resolution}_on_data",
            point_size,
            dpi,
        )
        if bool(get_sample_value(sample_cfg, cfg, "write_axis_scout_plots", True)):
            plot_spatial_axis_scout(
                adata,
                "conference_cell_type",
                f"{sample_name} cell-type labels, {resolution}, axis scout",
                cell_type_dir,
                f"{sample_name}_spatial_cell_type_labels_{resolution}_axis_scout",
                point_size,
                dpi,
            )

    dotplot_groupby = "conference_cell_type" if has_cell_type_labels else "conference_cluster"
    for marker_set in marker_sets:
        plot_marker_spatial_panels(adata, marker_set, marker_dir, sample_name, point_size, dpi)
        plot_marker_dotplot(adata, marker_set, dotplot_groupby, dotplot_dir, sample_name, resolution, dpi)

    plot_crop_windows(adata, sample_cfg, cfg, crop_dir, sample_name, resolution, marker_sets)
    plot_cluster_abundance(adata, abundance_dir, sample_name, resolution, dpi)


def main():
    """Run configured conference figure generation."""
    args = parse_args()
    load_plotting_stack()
    cfg = read_config(args.config)

    requested_samples = set(args.only_sample)
    for sample_cfg in cfg["samples"]:
        if requested_samples and sample_cfg["sample"] not in requested_samples:
            continue
        run_sample(sample_cfg, cfg)

    missing_samples = requested_samples - {sample_cfg["sample"] for sample_cfg in cfg["samples"]}
    if missing_samples:
        raise ValueError(f"Requested samples were not found in config: {sorted(missing_samples)}")


if __name__ == "__main__":
    main()
