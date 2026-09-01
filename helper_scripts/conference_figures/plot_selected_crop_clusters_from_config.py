#!/usr/bin/env python
# coding: utf-8

"""
Title: Plot Selected Crop Cluster Figures From Config
Date: 2026-08-31
Summary: Regenerate only selected spatial crop cluster plots from the filtered
VBCT conference config. The helper is intended for quick presentation-figure
iteration and currently redraws CK_bowel_res view 2 and MG_gastric_non_res
view 3 with configurable base and endothelial-cell dot sizes.
"""

import argparse
from pathlib import Path

import plot_spatial_conference_figures_from_config as conference


DEFAULT_SELECTED_CROPS = {
    "CK_bowel_res": {"view 2"},
    "MG_gastric_non_res": {"view 3"},
}


def parse_args():
    """Parse command-line options for selected crop cluster plotting."""
    parser = argparse.ArgumentParser(
        description="Regenerate selected crop cluster plots from a JSON config."
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
        "--cluster-point-size",
        type=float,
        default=17.0,
        help="Base point size for all cells in the selected cluster crops.",
    )
    parser.add_argument(
        "--endothelial-point-size",
        type=float,
        default=20.0,
        help="Point size for endothelial-labelled cells redrawn on top.",
    )
    return parser.parse_args()


def prepare_sample_adata(sample_cfg, cfg):
    """Load one sample and apply shared conference presentation transforms."""
    adata = conference.prepare_adata(sample_cfg["adata_path"], sample_cfg["cluster_col"])
    adata = conference.exclude_configured_clusters(adata, sample_cfg, cfg)
    conference.add_cell_type_labels(adata, sample_cfg)
    collapse_rules = conference.get_sample_value(sample_cfg, cfg, "label_collapse_rules", [])
    conference.collapse_cell_type_labels(adata, collapse_rules)
    conference.apply_conference_palette(adata, sample_cfg, cfg)
    return adata


def selected_windows(sample_cfg, cfg):
    """Yield configured crop windows selected for quick regeneration."""
    selected_names = DEFAULT_SELECTED_CROPS.get(sample_cfg["sample"], set())
    for window in conference.get_sample_value(sample_cfg, cfg, "crop_windows", []):
        if window.get("name") in selected_names:
            yield window


def endothelial_highlight_config(point_size):
    """Return a highlight rule for endothelial-labelled cells."""
    return {
        "obs_col": "conference_cell_type",
        "contains": ["Endothelial"],
        "point_size": point_size,
        "alpha": 1.0,
    }


def run_sample(sample_cfg, cfg, args):
    """Regenerate selected crop cluster plots for one sample."""
    sample_name = sample_cfg["sample"]
    if sample_name not in DEFAULT_SELECTED_CROPS:
        return

    resolution = sample_cfg.get("resolution", "selected_resolution")
    output_root = conference.get_sample_value(
        sample_cfg, cfg, "output_dir", "figures/conference"
    )
    output_subdir = sample_cfg.get("output_subdir", f"{sample_name}_{resolution}")
    output_dir = conference.ensure_output_dir(Path(output_root) / output_subdir)
    crop_dir = conference.figure_type_dir(output_dir, cfg, "crops")
    dpi = int(conference.get_sample_value(sample_cfg, cfg, "dpi", 300))
    crop_figsize = conference.get_sample_value(sample_cfg, cfg, "crop_figsize")
    frameon = bool(conference.get_sample_value(sample_cfg, cfg, "crop_show_axes", True))

    adata = prepare_sample_adata(sample_cfg, cfg)
    for window in selected_windows(sample_cfg, cfg):
        name = window["name"].replace(" ", "_")
        cropped = conference.crop_adata_to_window(adata, window)
        print(
            f"Regenerating {sample_name} {window['name']} cluster crop "
            f"with base point size {args.cluster_point_size:g} and "
            f"endothelial point size {args.endothelial_point_size:g}."
        )
        conference.plot_spatial_category(
            cropped,
            "conference_cluster",
            f"{sample_name} BANKSY clusters, {resolution}, {window['name']}",
            crop_dir,
            f"{sample_name}_spatial_clusters_{resolution}_{name}",
            args.cluster_point_size,
            dpi,
            frameon=frameon,
            highlight_cfg=endothelial_highlight_config(args.endothelial_point_size),
            force_legend=True,
            legend_label_col="conference_cell_type",
            figsize=window.get("figsize", crop_figsize),
        )


def main():
    """Run selected crop cluster figure generation."""
    args = parse_args()
    conference.load_plotting_stack()
    cfg = conference.read_config(args.config)

    for sample_cfg in cfg["samples"]:
        run_sample(sample_cfg, cfg, args)


if __name__ == "__main__":
    main()
