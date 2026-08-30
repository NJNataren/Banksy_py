#!/usr/bin/env python
# coding: utf-8

"""
Title: Plot Conference UMAP Figures From Config
Date: 2026-08-29
Summary: Generate UMAP-only conference figure variants from the shared JSON
config. The script reuses the main conference helper for AnnData loading,
cluster masking, label cleanup, and palette assignment, then writes lighter
right-margin, on-data labelled, and optional black-underlay UMAP plots for quick iteration. Dense point layers are rasterized in vector outputs by the shared saver so labels stay editable.
"""

import argparse
from pathlib import Path

import matplotlib.patheffects as path_effects

import plot_spatial_conference_figures_from_config as conference


def parse_args():
    """Parse command-line options for UMAP-only conference plotting."""
    parser = argparse.ArgumentParser(
        description="Generate UMAP-only conference figures from the shared JSON config."
    )
    parser.add_argument("--config", required=True, help="Path to JSON figure config.")
    parser.add_argument(
        "--only-sample",
        action="append",
        default=[],
        help="Sample name to run. Repeat to run several samples from one config.",
    )
    parser.add_argument(
        "--output-dirname",
        default="umap_rework",
        help="Plot-family subdirectory name under each sample output directory.",
    )
    parser.add_argument(
        "--point-size",
        type=float,
        default=None,
        help="Override UMAP point size. Defaults to config umap_point_size, then 1.0.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.35,
        help="Point alpha for less-solid dense UMAP clusters.",
    )
    parser.add_argument(
        "--figsize",
        type=float,
        nargs=2,
        default=None,
        metavar=("WIDTH", "HEIGHT"),
        help="Override UMAP figure size in inches.",
    )
    parser.add_argument(
        "--outlined",
        action="store_true",
        help="Also write cartoon-style cluster UMAPs with black point underlays and on-plot labels.",
    )
    parser.add_argument(
        "--skip-standard",
        action="store_true",
        help="Only write requested non-standard UMAP styles, such as --outlined.",
    )
    parser.add_argument(
        "--outline-point-size",
        type=float,
        default=3.0,
        help="Point size for outlined cluster UMAPs.",
    )
    parser.add_argument(
        "--outline-alpha",
        type=float,
        default=0.45,
        help="Point alpha for outlined cluster UMAPs.",
    )
    parser.add_argument(
        "--outline-width",
        type=float,
        default=2.0,
        help="Line width for optional convex hull outlines around clusters.",
    )
    parser.add_argument(
        "--hull-outlines",
        action="store_true",
        help="Draw convex hull lines around clusters. Off by default because dispersed clusters can make messy crossings.",
    )
    return parser.parse_args()


def prepare_sample_adata(sample_cfg, cfg):
    """Load one sample and apply the same presentation transforms as the main helper."""
    adata = conference.prepare_adata(sample_cfg["adata_path"], sample_cfg["cluster_col"])
    adata = conference.exclude_configured_clusters(adata, sample_cfg, cfg)
    has_cell_type_labels = conference.add_cell_type_labels(adata, sample_cfg)
    collapse_rules = conference.get_sample_value(sample_cfg, cfg, "label_collapse_rules", [])
    conference.collapse_cell_type_labels(adata, collapse_rules)
    conference.apply_conference_palette(adata, sample_cfg, cfg)
    return adata, has_cell_type_labels


def select_umap_basis(adata, sample_cfg, cfg):
    """Copy the first available configured UMAP embedding to X_conference_umap."""
    umap_keys = conference.get_sample_value(
        sample_cfg, cfg, "umap_keys", conference.DEFAULT_UMAP_KEYS
    )
    for obsm_key in umap_keys:
        if obsm_key in adata.obsm:
            adata.obsm["X_conference_umap"] = adata.obsm[obsm_key]
            print(f"Using {obsm_key} for UMAP plot.")
            return True

    print("Skipping UMAP plot: no configured UMAP key found.")
    return False


def save_umap_variant(
    adata,
    color_col,
    output_dir,
    stem,
    title,
    point_size,
    alpha,
    dpi,
    legend_loc,
    figsize,
):
    """Write one UMAP variant for a categorical observation column."""
    plot_kwargs = {
        "basis": "conference_umap",
        "color": color_col,
        "size": point_size,
        "alpha": alpha,
        "frameon": False,
        "legend_loc": legend_loc,
        "title": title,
        "show": False,
    }
    if legend_loc == "on data":
        plot_kwargs["legend_fontoutline"] = 2

    conference.sc.pl.embedding(adata, **plot_kwargs)
    if figsize:
        conference.plt.gcf().set_size_inches(float(figsize[0]), float(figsize[1]), forward=True)
    conference.save_current_figure(output_dir, stem, dpi)


def categorical_color_map(adata, obs_col):
    """Return a category-to-color mapping for a plotted observation column."""
    categories = list(adata.obs[obs_col].cat.categories)
    uns_key = f"{obs_col}_colors"
    colors = list(adata.uns.get(uns_key, []))
    if len(colors) >= len(categories):
        return dict(zip(categories, colors))

    fallback = conference.plt.get_cmap("tab20")
    return {category: fallback(i % fallback.N) for i, category in enumerate(categories)}


def convex_hull(points):
    """Return convex hull points using the monotonic chain algorithm."""
    if len(points) < 3:
        return points

    ordered = sorted({(float(x), float(y)) for x, y in points})
    if len(ordered) < 3:
        return conference.np.asarray(ordered)

    def cross(origin, first, second):
        return (first[0] - origin[0]) * (second[1] - origin[1]) - (
            first[1] - origin[1]
        ) * (second[0] - origin[0])

    lower = []
    for point in ordered:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)

    upper = []
    for point in reversed(ordered):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)

    return conference.np.asarray(lower[:-1] + upper[:-1])


def save_outlined_umap(
    adata,
    output_dir,
    stem,
    title,
    point_size,
    alpha,
    outline_width,
    draw_hulls,
    dpi,
    figsize,
):
    """Write a cluster UMAP with black dot underlays and on-plot labels."""
    coords = conference.np.asarray(adata.obsm["X_conference_umap"])
    labels = adata.obs["conference_cluster"]
    colors = categorical_color_map(adata, "conference_cluster")

    fig, ax = conference.plt.subplots(figsize=tuple(figsize or [7, 7]))
    ax.set_title(title)
    ax.set_axis_off()

    for category in labels.cat.categories:
        mask = (labels == category).to_numpy()
        cluster_coords = coords[mask]
        if cluster_coords.size == 0:
            continue

        # Draw a slightly larger black point cloud underneath the coloured points
        # so dense clusters have a clean cartoon-style edge without global hull artifacts.
        ax.scatter(
            cluster_coords[:, 0],
            cluster_coords[:, 1],
            s=point_size * 2.4,
            c="black",
            alpha=alpha,
            linewidths=0,
            rasterized=True,
        )
        ax.scatter(
            cluster_coords[:, 0],
            cluster_coords[:, 1],
            s=point_size,
            c=[colors[category]],
            alpha=alpha,
            linewidths=0,
            rasterized=True,
        )

        if draw_hulls:
            hull = convex_hull(cluster_coords)
            if len(hull) >= 3:
                closed_hull = conference.np.vstack([hull, hull[0]])
                ax.plot(
                    closed_hull[:, 0],
                    closed_hull[:, 1],
                    color="black",
                    linewidth=outline_width,
                    solid_capstyle="round",
                    solid_joinstyle="round",
                )

        center = conference.np.median(cluster_coords, axis=0)
        text = ax.text(
            center[0],
            center[1],
            str(category),
            color="white",
            fontsize=11,
            fontweight="bold",
            ha="center",
            va="center",
        )
        text.set_path_effects([path_effects.withStroke(linewidth=3.0, foreground="black")])

    ax.margins(0.06)
    conference.save_current_figure(output_dir, stem, dpi)


def run_sample(sample_cfg, cfg, args):
    """Generate UMAP variants for one configured sample."""
    sample_name = sample_cfg["sample"]
    resolution = sample_cfg.get("resolution", "selected_resolution")
    output_root = conference.get_sample_value(sample_cfg, cfg, "output_dir", "figures/conference")
    output_subdir = sample_cfg.get("output_subdir", f"{sample_name}_{resolution}")
    output_dir = conference.figure_type_dir(
        conference.ensure_output_dir(Path(output_root) / output_subdir),
        cfg,
        args.output_dirname,
    )
    dpi = int(conference.get_sample_value(sample_cfg, cfg, "dpi", 300))
    point_size = args.point_size
    if point_size is None:
        point_size = float(conference.get_sample_value(sample_cfg, cfg, "umap_point_size", 1.0))
    figsize = args.figsize or conference.get_sample_value(sample_cfg, cfg, "umap_figsize", [8, 7])

    adata, has_cell_type_labels = prepare_sample_adata(sample_cfg, cfg)
    if not select_umap_basis(adata, sample_cfg, cfg):
        return

    print(f"Running UMAP rework for {sample_name}: {adata.n_obs:,} cells.")
    if not args.skip_standard:
        for color_col, label in [("conference_cluster", "BANKSY clusters")]:
            save_umap_variant(
                adata,
                color_col,
                output_dir,
                f"{sample_name}_umap_{color_col}_{resolution}_right_margin",
                f"{sample_name} UMAP, {resolution} {label}",
                point_size,
                args.alpha,
                dpi,
                "right margin",
                figsize,
            )
            save_umap_variant(
                adata,
                color_col,
                output_dir,
                f"{sample_name}_umap_{color_col}_{resolution}_on_data",
                f"{sample_name} UMAP, {resolution} {label}",
                point_size,
                args.alpha,
                dpi,
                "on data",
                figsize,
            )

        if has_cell_type_labels:
            save_umap_variant(
                adata,
                "conference_cell_type",
                output_dir,
                f"{sample_name}_umap_cell_type_labels_{resolution}_right_margin",
                f"{sample_name} UMAP, {resolution} cell-type labels",
                point_size,
                args.alpha,
                dpi,
                "right margin",
                figsize,
            )
            save_umap_variant(
                adata,
                "conference_cell_type",
                output_dir,
                f"{sample_name}_umap_cell_type_labels_{resolution}_on_data",
                f"{sample_name} UMAP, {resolution} cell-type labels",
                point_size,
                args.alpha,
                dpi,
                "on data",
                figsize,
            )

    if args.outlined:
        save_outlined_umap(
            adata,
            output_dir,
            f"{sample_name}_umap_conference_cluster_{resolution}_outlined",
            f"{sample_name} UMAP, {resolution} outlined clusters",
            args.outline_point_size,
            args.outline_alpha,
            args.outline_width,
            args.hull_outlines,
            dpi,
            figsize,
        )


def main():
    """Run UMAP-only conference figure generation."""
    args = parse_args()
    conference.load_plotting_stack()
    cfg = conference.read_config(args.config)

    requested_samples = set(args.only_sample)
    for sample_cfg in cfg["samples"]:
        if requested_samples and sample_cfg["sample"] not in requested_samples:
            continue
        run_sample(sample_cfg, cfg, args)

    missing_samples = requested_samples - {sample_cfg["sample"] for sample_cfg in cfg["samples"]}
    if missing_samples:
        raise ValueError(f"Requested samples were not found in config: {sorted(missing_samples)}")


if __name__ == "__main__":
    main()
