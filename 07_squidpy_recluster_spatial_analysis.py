#!/usr/bin/env python
# coding: utf-8

"""
Title: Squidpy Spatial Analysis For BANKSY Reclusters
Date: 2026-08-18
Summary: Read the clean expression AnnData used for BANKSY by script 06,
then run Squidpy neighbourhood enrichment, centrality, co-occurrence, and
Moran's I analyses for selected BANKSY reclustering resolutions.
"""

import argparse
import json
import os
import random
from datetime import datetime
from matplotlib import rc_context

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import squidpy as sq


sc.logging.print_header()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 1
plt.rcParams["font.family"] = "Arial"
sns.set_style("white")

seed = 1234
np.random.seed(seed)
random.seed(seed)


def parse_args():
    """Parse command-line arguments for the Squidpy recluster analysis."""
    parser = argparse.ArgumentParser(
        prog="run Squidpy spatial analyses on script 06 BANKSY reclusters"
    )
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="JSON config with project, dataset_name, input_label, and reclustering settings.",
    )
    return parser.parse_args()


def log_time(step):
    """Print a timestamped workflow message."""
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {step}")


def ensure_directory(path, label):
    """Create `path` and log the resolved location."""
    os.makedirs(path, exist_ok=True)
    print(f"{label} directory ready: {os.path.abspath(path)}")




def make_cluster_palette(n_colors):
    """Create enough distinct categorical colours for cluster plots."""
    base_cmaps = ["tab20", "tab20b", "tab20c"]
    colors = []

    for cmap_name in base_cmaps:
        cmap = plt.get_cmap(cmap_name)
        colors.extend([mpl.colors.to_hex(cmap(i)) for i in range(cmap.N)])

    if n_colors > len(colors):
        extra_cmap = plt.get_cmap("nipy_spectral")
        extra_count = n_colors - len(colors)
        colors.extend(
            [
                mpl.colors.to_hex(extra_cmap(i / max(extra_count - 1, 1)))
                for i in range(extra_count)
            ]
        )

    return colors[:n_colors]


def make_banksy_label_col(nbr_weight_decay, pc_label, lambda_label, resolution):
    """Build the native BANKSY cluster-label column name for one resolution."""
    return (
        f"labels_{nbr_weight_decay}"
        f"_pc{pc_label}"
        f"_nc{lambda_label}"
        f"_r{resolution:.2f}"
    )


def make_recluster_label_col(
    nbr_weight_decay,
    pc_label,
    lambda_label,
    resolution,
    input_label,
):
    """Build the provenance-safe reclustered label column name from script 06."""
    banksy_col = make_banksy_label_col(
        nbr_weight_decay,
        pc_label,
        lambda_label,
        resolution,
    )
    return f"{banksy_col}_recluster_{input_label}"


def resolve_analysis_resolutions(cfg, res_label_list):
    """Return the reclustering resolutions that should receive Squidpy analyses."""
    analysis_res_label = cfg.get("analysis_res_label", res_label_list[-1])
    if analysis_res_label == "all":
        return [float(res) for res in res_label_list]
    if isinstance(analysis_res_label, list):
        return [float(res) for res in analysis_res_label]
    return [float(analysis_res_label)]


def prepare_spatial_coordinates(adata, coord_keys):
    """Store Xenium x/y coordinates in `.obsm['xy']` and `.obsm['spatial']`."""
    for col in [coord_keys[0], coord_keys[1]]:
        if col not in adata.obs.columns:
            raise KeyError(f"Missing coordinate column in adata.obs: {col}")

    coords = adata.obs[[coord_keys[0], coord_keys[1]]].to_numpy()
    adata.obsm["xy"] = coords
    adata.obsm["spatial"] = coords


def run_squidpy_for_resolution(
    adata,
    cluster_key,
    output_dir,
    dataset_name,
    pc_label,
    lambda_label,
    res_str_single,
    cfg,
):
    """Run Squidpy spatial analyses for one reclustering label column."""
    if cluster_key not in adata.obs.columns:
        raise KeyError(
            f"{cluster_key!r} was not found in adata.obs. "
            f"Available obs columns: {list(adata.obs.columns)}"
        )

    if not pd.api.types.is_categorical_dtype(adata.obs[cluster_key]):
        adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")

    # Squidpy plotting looks for this colour vector in `.uns`. Supplying it
    # explicitly avoids noisy palette warnings and keeps plots consistent.
    color_key = f"{cluster_key}_colors"
    if color_key not in adata.uns:
        adata.uns[color_key] = make_cluster_palette(
            len(adata.obs[cluster_key].cat.categories)
        )

    squidpy_path = os.path.join(output_dir, "squidpy")
    ensure_directory(squidpy_path, "Squidpy")
    sc.settings.figdir = squidpy_path

    # All Squidpy analyses use the clean expression object used for BANKSY from
    # script 06. In qc_pass_only mode, failed cells are absent here, so they do not contribute to any
    # spatial graph, neighbourhood summary, or expression autocorrelation.
    log_time(f"Running Squidpy spatial graph for {cluster_key}")
    sq.gr.spatial_neighbors(adata, coord_type="generic")

    log_time(f"Running neighbourhood enrichment for {cluster_key}")
    sq.gr.nhood_enrichment(adata, cluster_key=cluster_key)
    with rc_context({"figure.figsize": (10, 10)}):
        sq.pl.nhood_enrichment(
            adata,
            cluster_key=cluster_key,
            show=False,
            save=(
                f"neighbourhood_enrichment_{dataset_name}"
                f"_pc{pc_label}_nc{lambda_label}_r{res_str_single}.png"
            ),
        )
    plt.close("all")

    log_time(f"Running centrality scores for {cluster_key}")
    sq.gr.centrality_scores(adata, cluster_key=cluster_key)
    centrality_key = f"{cluster_key}_centrality_scores"
    if centrality_key in adata.uns:
        centrality_csv = os.path.join(
            squidpy_path,
            f"centrality_scores_{dataset_name}_pc{pc_label}_nc{lambda_label}_r{res_str_single}.csv",
        )
        adata.uns[centrality_key].to_csv(centrality_csv)
        print(f"Saved centrality scores to: {centrality_csv}")

    with rc_context({"figure.figsize": (22, 4)}):
        sq.pl.centrality_scores(adata, cluster_key=cluster_key)
        fig = plt.gcf()
        centrality_png = os.path.join(
            squidpy_path,
            f"centrality_scores_{dataset_name}_pc{pc_label}_nc{lambda_label}_r{res_str_single}.png",
        )
        fig.savefig(centrality_png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved centrality plot to: {centrality_png}")

    log_time(f"Running co-occurrence for {cluster_key}")
    sq.gr.co_occurrence(adata, cluster_key=cluster_key)
    cluster_categories = list(adata.obs[cluster_key].cat.categories)
    max_cooccurrence_clusters = int(
        cfg.get("max_cooccurrence_clusters", len(cluster_categories))
    )
    for cluster_number in cluster_categories[:max_cooccurrence_clusters]:
        with rc_context({"figure.figsize": (10, 10)}):
            sq.pl.co_occurrence(
                adata,
                cluster_key=cluster_key,
                clusters=cluster_number,
                save=(
                    f"co_occurrence_cluster_{cluster_number}_{dataset_name}"
                    f"_pc{pc_label}_nc{lambda_label}_r{res_str_single}.png"
                ),
            )
        plt.close("all")

    # Moran's I uses `.X`, so this must be the clean normalized/log1p expression
    # matrix from script 06 rather than any BANKSY-expanded spatial object.
    log_time("Running Moran's I on clean log-normalized expression")
    moran_fraction = float(cfg.get("moran_subsample_fraction", 0.5))
    moran_n_perms = int(cfg.get("moran_n_perms", 100))
    moran_n_jobs = int(cfg.get("moran_n_jobs", 1))
    moran_fdr_threshold = float(cfg.get("moran_fdr_threshold", 0.10))
    top_moran_n = int(cfg.get("top_moran_n", 10))

    if moran_fraction < 1:
        adata_moran = sc.pp.subsample(
            adata,
            fraction=moran_fraction,
            copy=True,
            random_state=seed,
        )
    else:
        adata_moran = adata.copy()

    sq.gr.spatial_neighbors(adata_moran, coord_type="generic", delaunay=True)
    sq.gr.spatial_autocorr(
        adata_moran,
        mode="moran",
        n_perms=moran_n_perms,
        n_jobs=moran_n_jobs,
    )

    moran_scores = adata_moran.uns["moranI"].reset_index(names=["Gene"])
    moran_scores_raw = moran_scores[
        ~moran_scores["Gene"].str.contains("_nbr_0", case=False, na=False)
    ]
    moran_scores_raw = moran_scores_raw[
        ~moran_scores_raw["Gene"].str.contains("_nbr_1", case=False, na=False)
    ]
    moran_scores_sig = moran_scores_raw[
        moran_scores_raw["pval_sim_fdr_bh"] <= moran_fdr_threshold
    ]

    moran_all_csv = os.path.join(
        squidpy_path,
        f"moran_scores_all_{dataset_name}_pc{pc_label}_nc{lambda_label}_r{res_str_single}.csv",
    )
    moran_sig_csv = os.path.join(
        squidpy_path,
        f"moran_scores_fdr{moran_fdr_threshold}_{dataset_name}_pc{pc_label}_nc{lambda_label}_r{res_str_single}.csv",
    )
    moran_scores_raw.to_csv(moran_all_csv, index=False)
    moran_scores_sig.to_csv(moran_sig_csv, index=False)
    print(f"Saved all Moran's I scores to: {moran_all_csv}")
    print(f"Saved significant Moran's I scores to: {moran_sig_csv}")

    top_moran = moran_scores_sig.sort_values(by="pval_sim", ascending=True).head(
        top_moran_n
    )
    for gene in top_moran["Gene"].tolist():
        with rc_context({"figure.figsize": (12, 8)}):
            sq.pl.spatial_scatter(
                adata_moran,
                library_id="spatial",
                color=[gene],
                shape=None,
                size=2,
                img=False,
                show=False,
                save=(
                    f"spatial_scatter_top_moran_I_{gene}_{dataset_name}"
                    f"_pc{pc_label}_nc{lambda_label}_r{res_str_single}.png"
                ),
            )
        plt.close("all")


# -----------------------------------------------------------------------------
# Read config and resolve paths
# -----------------------------------------------------------------------------

args = parse_args()
with open(args.config) as f:
    cfg = json.load(f)

project = cfg["project"]
dataset_name = cfg["dataset_name"]
input_label = cfg.get("input_label", cfg.get("output_label", "filtered_qc_v1"))
recluster_inclusion = cfg.get("recluster_inclusion", "qc_pass_only")
if recluster_inclusion not in ["qc_pass_only", "all_cells"]:
    raise ValueError(
        "recluster_inclusion must be either 'qc_pass_only' or 'all_cells'. "
        f"Received: {recluster_inclusion!r}"
    )
run_label = cfg.get("run_label", f"{input_label}_{recluster_inclusion}")
pc_label = cfg["pc_label"]
lambda_label = cfg["lambda_label"]
res_label = cfg["res_label"]
res_label_list = res_label if isinstance(res_label, list) else [res_label]
res_str = "_".join(res_label_list)
nbr_weight_decay = cfg["nbr_weight_decay"]
coord_keys = tuple(cfg.get("coord_keys", ["x", "y", "xy"]))

base_dir = cfg.get("base_dir", "data/xenium")
processed_path = os.path.join(base_dir, "processed", project, dataset_name)
output_path = cfg.get(
    "output_path",
    os.path.join(
        base_dir,
        "output",
        project,
        dataset_name,
        f"reclustering_{run_label}",
    ),
)
ensure_directory(output_path, "Output reclustering sample")

input_h5ad = cfg.get(
    "input_h5ad",
    os.path.join(
        processed_path,
        f"adata_expression_clean_{dataset_name}_recluster_{run_label}"
        f"_cells_used_for_banksy_with_clusters_{res_str}.h5ad",
    ),
)
if not os.path.isfile(input_h5ad):
    raise FileNotFoundError(f"Script 06 reclustered-cell AnnData does not exist: {input_h5ad}")

output_h5ad = cfg.get(
    "output_h5ad",
    os.path.join(
        processed_path,
        f"adata_expression_clean_{dataset_name}_recluster_{run_label}"
        f"_cells_used_for_banksy_with_clusters_{res_str}_squidpy.h5ad",
    ),
)

# -----------------------------------------------------------------------------
# Load clean passing-cell object and run selected Squidpy resolutions
# -----------------------------------------------------------------------------

log_time(f"Reading script 06 reclustered-cell clean AnnData from: {input_h5ad}")
adata = ad.read_h5ad(input_h5ad)
prepare_spatial_coordinates(adata, coord_keys)

analysis_resolutions = resolve_analysis_resolutions(cfg, res_label_list)
print(f"Running Squidpy analyses for resolutions: {analysis_resolutions}")

for res in analysis_resolutions:
    res_str_single = str(res).replace(".", "p")
    cluster_key = make_recluster_label_col(
        nbr_weight_decay,
        pc_label,
        lambda_label,
        res,
        run_label,
    )
    run_squidpy_for_resolution(
        adata=adata,
        cluster_key=cluster_key,
        output_dir=output_path,
        dataset_name=dataset_name,
        pc_label=pc_label,
        lambda_label=lambda_label,
        res_str_single=res_str_single,
        cfg=cfg,
    )

adata.write_h5ad(output_h5ad)
print(f"Wrote Squidpy-annotated reclustered-cell AnnData to: {output_h5ad}")
log_time(f"Finished Squidpy recluster spatial analysis for {dataset_name}.")
