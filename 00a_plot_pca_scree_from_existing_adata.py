#!/usr/bin/env python
# coding: utf-8

"""
Title: Plot PCA Scree From Existing Xenium AnnData
Date: 2026-08-04
Summary: Generate PCA scree plots and variance summary CSVs from existing
processed Xenium AnnData objects without running BANKSY clustering. The script
can run one or many clustering-style JSON configs and writes pca_qc outputs
under each sample's existing output directory.
"""

import argparse
import json
import os
from pathlib import Path


np = None
pd = None
plt = None
sc = None
seed = 1234


def load_scientific_stack():
    """Import plotting and single-cell dependencies only when scree plotting runs."""
    global np, pd, plt, sc

    if sc is not None:
        return

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as pyplot
    import numpy as numpy
    import pandas as pandas
    import scanpy as scanpy

    np = numpy
    pd = pandas
    plt = pyplot
    sc = scanpy




def ensure_directory(path, label):
    """Create `path` and any missing parent directories, then log its location."""
    os.makedirs(path, exist_ok=True)
    print(f"{label} directory ready: {os.path.abspath(path)}")


def read_json_config(config_path):
    """Read one JSON config file and add its source path for logging."""
    with open(config_path) as handle:
        cfg = json.load(handle)
    cfg["_config_path"] = str(config_path)
    return cfg


def collect_config_paths(configs, config_dirs, recursive):
    """Collect sorted JSON config paths from explicit files and directories."""
    config_paths = [Path(config) for config in configs]

    for config_dir in config_dirs:
        directory = Path(config_dir)
        pattern = "**/*.json" if recursive else "*.json"
        config_paths.extend(sorted(directory.glob(pattern)))

    unique_paths = []
    seen = set()
    for path in config_paths:
        resolved = path.resolve()
        if resolved in seen:
            continue
        if not path.is_file():
            raise FileNotFoundError(f"Config path is not a file: {path}")
        unique_paths.append(path)
        seen.add(resolved)

    if not unique_paths:
        raise ValueError("No config files supplied. Use --config or --config-dir.")

    return unique_paths


def resolve_scree_n_pcs(adata, selected_n_pcs, requested_n_pcs):
    """Choose a valid PCA depth for scree plotting beyond the selected PC count."""
    max_valid_n_pcs = min(adata.n_obs - 1, adata.n_vars - 1)
    scree_n_pcs = min(requested_n_pcs, max_valid_n_pcs)

    if scree_n_pcs < 2:
        print(
            "Skipping PCA scree plot: fewer than two valid components are "
            f"available for object shape {adata.shape}."
        )
        return None

    if scree_n_pcs < requested_n_pcs:
        print(
            f"PCA scree requested {requested_n_pcs} PCs, capped at "
            f"{scree_n_pcs} for object shape {adata.shape}."
        )

    return scree_n_pcs


def resolve_existing_adata_path(processed_path, dataset_name, input_h5ad=None):
    """Find the best existing AnnData object for scree plotting."""
    if input_h5ad is not None:
        return Path(input_h5ad), "custom"

    clean_expression_h5ad = (
        processed_path / f"adata_expression_clean_{dataset_name}_normalized_log1p.h5ad"
    )
    if clean_expression_h5ad.is_file():
        return clean_expression_h5ad, "normalized_log1p"

    float_h5ad = processed_path / f"{dataset_name}_float_32.h5ad"
    if float_h5ad.is_file():
        return float_h5ad, "float32_counts"

    raise FileNotFoundError(
        "Could not find an existing AnnData object for scree plotting. Checked: "
        f"{clean_expression_h5ad} and {float_h5ad}"
    )


def load_scree_adata(adata_path, input_kind):
    """Load AnnData and normalize/log1p only when using a count-matrix fallback."""
    print(f"Reading AnnData for PCA scree plot: {adata_path}")
    adata = sc.read_h5ad(adata_path)

    if input_kind == "float32_counts":
        print("Using float32 fallback object; applying normalize_total and log1p in memory.")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        print(f"Using existing {input_kind} object as stored.")

    return adata


def plot_pca_scree(adata, output_dir, dataset_name, selected_n_pcs, scree_n_pcs):
    """Run PCA and save scree plot plus variance summary for one AnnData object."""
    sc.tl.pca(
        adata,
        n_comps=scree_n_pcs,
        svd_solver="arpack",
        random_state=seed,
    )
    variance_ratio = np.asarray(adata.uns["pca"]["variance_ratio"]).copy()

    n_show = min(scree_n_pcs, len(variance_ratio))
    pcs = np.arange(1, n_show + 1)
    variance_percent = variance_ratio[:n_show] * 100
    cumulative_percent = np.cumsum(variance_ratio[:n_show]) * 100

    scree_dir = os.path.join(output_dir, "pca_qc")
    ensure_directory(scree_dir, "PCA QC")

    variance_df = pd.DataFrame(
        {
            "pc": pcs,
            "variance_ratio": variance_ratio[:n_show],
            "variance_percent": variance_percent,
            "cumulative_variance_percent": cumulative_percent,
            "selected_for_banksy": pcs <= selected_n_pcs,
        }
    )
    variance_csv = os.path.join(
        scree_dir,
        f"{dataset_name}_pca_scree_variance.csv",
    )
    variance_df.to_csv(variance_csv, index=False)

    fig, ax1 = plt.subplots(figsize=(9, 5))
    bar_container = ax1.bar(
        pcs,
        variance_percent,
        color="#4C78A8",
        alpha=0.85,
        label="Variance explained per PC",
    )
    ax1.set_xlabel("Principal component")
    ax1.set_ylabel("Variance explained (%)")
    ax1.set_xlim(0.5, n_show + 0.5)
    x_ticks = [1] + list(range(10, n_show + 1, 10))
    if selected_n_pcs <= n_show and selected_n_pcs not in x_ticks:
        x_ticks.append(selected_n_pcs)
    if n_show not in x_ticks:
        x_ticks.append(n_show)
    ax1.set_xticks(sorted(x_ticks))
    ax1.set_title(
        f"{dataset_name}: PCA scree plot before BANKSY clustering",
        fontsize=13,
    )

    ax2 = ax1.twinx()
    cumulative_line, = ax2.plot(
        pcs,
        cumulative_percent,
        color="#222222",
        marker="o",
        linewidth=1.4,
        label="Cumulative variance explained",
    )
    ax2.set_ylabel("Cumulative variance explained (%)")

    legend_handles = [bar_container, cumulative_line]
    if selected_n_pcs <= n_show:
        selected_line = ax1.axvline(
            selected_n_pcs,
            color="#D62728",
            linestyle="--",
            linewidth=1.2,
            label=f"Selected PC count ({selected_n_pcs})",
        )
        legend_handles.append(selected_line)
        ax1.text(
            selected_n_pcs,
            ax1.get_ylim()[1] * 0.95,
            f"selected PC {selected_n_pcs}",
            color="#D62728",
            ha="right",
            va="top",
            fontsize=9,
        )
    else:
        print(
            f"Selected PC count {selected_n_pcs} is beyond the scree plot range "
            f"of {n_show} PCs."
        )

    ax1.legend(handles=legend_handles, loc="lower right", frameon=False)
    fig.tight_layout()
    scree_png = os.path.join(
        scree_dir,
        f"{dataset_name}_pca_scree_plot.png",
    )
    fig.savefig(scree_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved PCA scree plot to: {scree_png}")
    print(f"Saved PCA variance summary to: {variance_csv}")


def run_config(cfg, args):
    """Generate a PCA scree plot for one clustering-style config."""
    dataset_name = cfg["dataset_name"]
    project = cfg.get("project", "")
    selected_n_pcs = int(cfg.get("pc_label", args.selected_n_pcs))
    requested_n_pcs = int(args.scree_n_pcs or cfg.get("scree_n_pcs", max(selected_n_pcs + 20, 50)))

    processed_path = Path(args.base_dir) / "processed" / project / dataset_name
    output_dir = Path(args.output_dir) if args.output_dir else Path(args.base_dir) / "output" / project / dataset_name
    ensure_directory(output_dir, "Output sample")

    adata_path, input_kind = resolve_existing_adata_path(
        processed_path=processed_path,
        dataset_name=dataset_name,
        input_h5ad=args.input_h5ad,
    )
    adata = load_scree_adata(adata_path, input_kind)
    scree_n_pcs = resolve_scree_n_pcs(adata, selected_n_pcs, requested_n_pcs)
    if scree_n_pcs is None:
        return

    plot_pca_scree(
        adata=adata,
        output_dir=output_dir,
        dataset_name=dataset_name,
        selected_n_pcs=selected_n_pcs,
        scree_n_pcs=scree_n_pcs,
    )


def parse_args():
    """Parse command-line arguments for single-sample or batch scree plotting."""
    parser = argparse.ArgumentParser(
        description=(
            "Generate PCA scree plots from existing processed Xenium AnnData "
            "objects without rerunning BANKSY clustering."
        )
    )
    parser.add_argument(
        "--config",
        action="append",
        default=[],
        help="Path to a clustering-style JSON config. May be supplied more than once.",
    )
    parser.add_argument(
        "--config-dir",
        action="append",
        default=[],
        help="Directory containing JSON configs to process.",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Find JSON configs recursively inside each --config-dir.",
    )
    parser.add_argument(
        "--base-dir",
        default="data/xenium",
        help="Base Xenium data directory containing processed/ and output/.",
    )
    parser.add_argument(
        "--input-h5ad",
        default=None,
        help="Optional explicit AnnData path. Intended for one --config run only.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Optional explicit output directory. Intended for one --config run only.",
    )
    parser.add_argument(
        "--scree-n-pcs",
        type=int,
        default=None,
        help="Number of PCs to calculate for the scree plot. Overrides config scree_n_pcs.",
    )
    parser.add_argument(
        "--selected-n-pcs",
        type=int,
        default=30,
        help="Fallback selected PC count if the config has no pc_label.",
    )
    return parser.parse_args()


def main():
    """Run PCA scree plotting for all supplied configs."""
    args = parse_args()
    config_paths = collect_config_paths(args.config, args.config_dir, args.recursive)

    if len(config_paths) > 1 and (args.input_h5ad or args.output_dir):
        raise ValueError("--input-h5ad and --output-dir are only supported for one config.")

    load_scientific_stack()

    for config_path in config_paths:
        print(f"\n--- PCA scree for config: {config_path} ---")
        cfg = read_json_config(config_path)
        run_config(cfg, args)


if __name__ == "__main__":
    main()
