#!/usr/bin/env python
# coding: utf-8

"""
Title: Plot Zero-Count Xenium Cells From Raw AnnData
Date: 2026-08-23
Summary: Read raw Xenium AnnData from one config or a config directory, plot
the spatial distribution of zero-count cells before script 00 filtering, and
write a small per-sample summary CSV. This lightweight helper stops before
normalisation, PCA, BANKSY, or clustering.
"""

import argparse
import json
import os
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    """Parse command-line arguments for a single config or config directory."""
    parser = argparse.ArgumentParser(
        description=(
            "Plot pre-filter zero-count Xenium cells from raw AnnData without "
            "running script 00 clustering."
        )
    )
    config_group = parser.add_mutually_exclusive_group(required=True)
    config_group.add_argument(
        "--config",
        help="Path to one script 00-style JSON config.",
    )
    config_group.add_argument(
        "--config-dir",
        help="Directory of script 00-style JSON configs to run in sorted order.",
    )
    return parser.parse_args()


def load_config(config_path):
    """Load a JSON config file."""
    with open(config_path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def ensure_directory(path, label):
    """Create an output directory and print its absolute path."""
    os.makedirs(path, exist_ok=True)
    print(f"{label} directory ready: {os.path.abspath(path)}")


def resolve_config_paths(args):
    """Return config paths from either `--config` or sorted `--config-dir` JSONs."""
    if args.config:
        return [Path(args.config)]

    config_dir = Path(args.config_dir)
    config_paths = sorted(config_dir.glob("*.json"))
    if not config_paths:
        raise FileNotFoundError(f"No JSON configs found in: {config_dir}")
    return config_paths


def plot_zero_count_cells(config_path):
    """Plot and summarise cells with `nCount_Xenium <= 0` for one sample."""
    cfg = load_config(config_path)
    dataset_name = cfg["dataset_name"]
    project = cfg.get("project", "")
    raw_subdir = cfg.get("raw_subdir", "")
    base_dir = cfg.get("base_dir", "data/xenium")

    raw_path = os.path.join(base_dir, "raw_data", raw_subdir)
    raw_adata_path = os.path.join(raw_path, f"{dataset_name}_raw.h5ad")
    qc_path = os.path.join(base_dir, "output", project, "QC_testing", dataset_name)

    if not os.path.isfile(raw_adata_path):
        raise FileNotFoundError(
            f"Raw AnnData file '{raw_adata_path}' does not exist. "
            "Check dataset_name and raw_subdir in the config."
        )

    ensure_directory(qc_path, "Zero-count QC output")

    print(f"Reading raw AnnData from: {raw_adata_path}")
    adata = ad.read_h5ad(raw_adata_path)

    required_obs_columns = ["x", "y", "nCount_Xenium"]
    missing_columns = [
        column for column in required_obs_columns if column not in adata.obs.columns
    ]
    if missing_columns:
        raise KeyError(
            f"Missing required obs columns for {dataset_name}: {missing_columns}. "
            f"Available obs columns: {list(adata.obs.columns)}"
        )

    # Keep this before any zero-count filtering so empty cells are still visible.
    adata.obsm["xy"] = np.vstack([adata.obs["x"], adata.obs["y"]]).T

    zero_count_mask = adata.obs["nCount_Xenium"] <= 0
    adata.obs["zero_count_cell_cat"] = (
        zero_count_mask.map({True: "Fail", False: "Pass"}).astype("category")
    )

    zero_count_summary = (
        adata.obs["zero_count_cell_cat"]
        .value_counts()
        .rename_axis("qc_status")
        .reset_index(name="n_cells")
    )
    zero_count_summary["percent_cells"] = (
        zero_count_summary["n_cells"] / adata.n_obs * 100
    )

    zero_count_summary_path = os.path.join(
        qc_path,
        f"{dataset_name}_zero_count_cell_summary.csv",
    )
    zero_count_summary.to_csv(zero_count_summary_path, index=False)

    xy = adata.obsm["xy"]
    fig, ax = plt.subplots(figsize=(12, 8))

    for status, color, zorder, alpha in [
        ("Pass", "orange", 2, 0.45),
        ("Fail", "dodgerblue", 3, 0.9),
    ]:
        mask = np.asarray(adata.obs["zero_count_cell_cat"].astype(str) == status)
        if not mask.any():
            continue

        ax.scatter(
            xy[mask, 0],
            xy[mask, 1],
            s=2,
            c=color,
            label=status,
            linewidths=0,
            alpha=alpha,
            rasterized=True,
            zorder=zorder,
        )

    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"{dataset_name}: zero-count cells before filtering")
    ax.legend(fontsize=14, markerscale=4, frameon=False)

    zero_count_plot_path = os.path.join(
        qc_path,
        f"tissue_spatial_scatter_zero_count_cell_cat_{dataset_name}.png",
    )
    fig.savefig(zero_count_plot_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved zero-count cell QC summary to: {zero_count_summary_path}")
    print(f"Saved zero-count cell QC plot to: {zero_count_plot_path}")


def main():
    """Run zero-count QC plotting for all requested configs."""
    args = parse_args()
    config_paths = resolve_config_paths(args)

    for index, config_path in enumerate(config_paths, start=1):
        print("=" * 72)
        print(f"Running zero-count QC {index}/{len(config_paths)}: {config_path}")
        plot_zero_count_cells(config_path)


if __name__ == "__main__":
    main()
