#!/usr/bin/env python
# coding: utf-8

"""
Title: Plot Individual Cluster Spatial Figures From Config
Date: 2026-08-31
Summary: Generate one spatial highlight plot per configured reclustered BANKSY
cluster from the shared conference JSON config. Each plot shows all QC-passing
cells in grey and redraws the selected cluster in its conference palette colour,
then writes PNG, PDF, and SVG files under each sample's individual_clusters
output directory.
"""

import argparse
import re
from pathlib import Path

import plot_spatial_conference_figures_from_config as conference


def parse_args():
    """Parse command-line options for individual cluster plotting."""
    parser = argparse.ArgumentParser(
        description="Generate one spatial highlight plot per cluster from a JSON config."
    )
    parser.add_argument(
        "--config",
        default=(
            "config/conference_figures/"
            "vbct_exemplar_ck_bowel_mg_gastric_filtered_reclustered_qc_v1.json"
        ),
        help="Path to JSON figure config.",
    )
    parser.add_argument(
        "--only-sample",
        action="append",
        default=[],
        help="Sample name to run. Repeat to run several samples from one config.",
    )
    parser.add_argument(
        "--background-color",
        default="#eeeeee",
        help="Colour for non-selected cells.",
    )
    parser.add_argument(
        "--background-alpha",
        type=float,
        default=0.55,
        help="Alpha for non-selected cells.",
    )
    parser.add_argument(
        "--highlight-alpha",
        type=float,
        default=1.0,
        help="Alpha for selected-cluster cells.",
    )
    parser.add_argument(
        "--point-size",
        type=float,
        default=None,
        help="Override spatial point size. Defaults to the shared config point_size.",
    )
    parser.add_argument(
        "--highlight-point-size",
        type=float,
        default=None,
        help="Override selected-cluster point size. Defaults to background point size.",
    )
    return parser.parse_args()


def safe_filename_part(value):
    """Return a filesystem-safe label fragment."""
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    return re.sub(r"_+", "_", safe).strip("_") or "unlabeled"


def prepare_sample_adata(sample_cfg, cfg):
    """Load one sample and apply shared conference presentation transforms."""
    adata = conference.prepare_adata(sample_cfg["adata_path"], sample_cfg["cluster_col"])
    adata = conference.exclude_configured_clusters(adata, sample_cfg, cfg)
    conference.add_cell_type_labels(adata, sample_cfg)
    collapse_rules = conference.get_sample_value(sample_cfg, cfg, "label_collapse_rules", [])
    conference.collapse_cell_type_labels(adata, collapse_rules)
    conference.apply_conference_palette(adata, sample_cfg, cfg)
    return adata


def color_by_cluster(adata):
    """Return cluster IDs mapped to their configured conference colours."""
    categories = list(adata.obs["conference_cluster"].cat.categories)
    palette = list(adata.uns.get("conference_cluster_colors", []))
    if len(palette) >= len(categories):
        return dict(zip([str(category) for category in categories], palette))

    fallback = conference.plt.get_cmap("tab20")
    return {
        str(category): fallback(i % fallback.N)
        for i, category in enumerate(categories)
    }


def annotation_for_cluster(adata, cluster_id):
    """Return the presentation annotation for a cluster."""
    clusters = adata.obs["conference_cluster"].astype(str)
    mask = clusters == str(cluster_id)
    if "conference_cell_type" not in adata.obs.columns or not bool(mask.any()):
        return f"cluster {cluster_id}"

    labels = adata.obs.loc[mask, "conference_cell_type"].astype(str)
    return labels.mode().iat[0] if not labels.empty else f"cluster {cluster_id}"


def plot_individual_cluster(
    adata,
    cluster_id,
    annotation,
    color,
    output_dir,
    sample_name,
    resolution,
    point_size,
    highlight_point_size,
    background_color,
    background_alpha,
    highlight_alpha,
    dpi,
    figsize=None,
):
    """Write one grey-background spatial plot with a selected cluster highlighted."""
    coords = conference.spatial_coordinate_frame(adata)
    clusters = adata.obs["conference_cluster"].astype(str)
    mask = clusters == str(cluster_id)
    n_highlighted = int(mask.sum())
    if n_highlighted == 0:
        print(f"Skipping empty cluster {cluster_id} for {sample_name}.")
        return

    fig_size = tuple(float(value) for value in figsize) if figsize else None
    _, ax = conference.plt.subplots(figsize=fig_size)
    ax.scatter(
        coords["x"],
        coords["y"],
        s=point_size,
        c=background_color,
        alpha=background_alpha,
        linewidths=0,
        rasterized=True,
    )
    ax.scatter(
        coords.loc[mask, "x"],
        coords.loc[mask, "y"],
        s=highlight_point_size,
        c=[color],
        alpha=highlight_alpha,
        linewidths=0,
        rasterized=True,
    )

    title = (
        f"{sample_name}, {resolution}, cluster {cluster_id}: "
        f"{annotation} ({n_highlighted:,} cells)"
    )
    ax.set_title(title)
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()
    ax.margins(0.02)

    stem = (
        f"{sample_name}_spatial_individual_cluster_{resolution}_"
        f"cluster_{safe_filename_part(cluster_id)}_{safe_filename_part(annotation)}"
    )
    conference.save_current_figure(output_dir, stem, dpi)


def run_sample(sample_cfg, cfg, args):
    """Generate individual cluster spatial plots for one configured sample."""
    sample_name = sample_cfg["sample"]
    resolution = sample_cfg.get("resolution", "selected_resolution")
    output_root = conference.get_sample_value(
        sample_cfg, cfg, "output_dir", "figures/conference"
    )
    output_subdir = sample_cfg.get("output_subdir", f"{sample_name}_{resolution}")
    output_dir = conference.figure_type_dir(
        conference.ensure_output_dir(Path(output_root) / output_subdir),
        cfg,
        "individual_clusters",
    )

    point_size = args.point_size
    if point_size is None:
        point_size = float(conference.get_sample_value(sample_cfg, cfg, "point_size", 3.0))
    highlight_point_size = args.highlight_point_size or point_size
    dpi = int(conference.get_sample_value(sample_cfg, cfg, "dpi", 300))
    figsize = conference.get_sample_value(sample_cfg, cfg, "full_spatial_figsize")

    adata = prepare_sample_adata(sample_cfg, cfg)
    cluster_colors = color_by_cluster(adata)
    categories = list(adata.obs["conference_cluster"].cat.categories)

    print(
        f"Running individual cluster plots for {sample_name}: "
        f"{len(categories)} clusters, {adata.n_obs:,} cells."
    )
    for cluster_id in categories:
        annotation = annotation_for_cluster(adata, cluster_id)
        plot_individual_cluster(
            adata,
            str(cluster_id),
            annotation,
            cluster_colors[str(cluster_id)],
            output_dir,
            sample_name,
            resolution,
            point_size,
            highlight_point_size,
            args.background_color,
            args.background_alpha,
            args.highlight_alpha,
            dpi,
            figsize=figsize,
        )


def main():
    """Run configured individual cluster spatial figure generation."""
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
