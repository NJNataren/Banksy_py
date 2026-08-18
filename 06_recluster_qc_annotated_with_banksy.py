#!/usr/bin/env python
# coding: utf-8

"""
Title: QC-Annotated Xenium BANKSY Reclustering
Date: 2026-08-18
Summary: Recluster a script 05 QC-annotated clean expression AnnData object
with BANKSY using either qc_pass_only or all_cells inclusion, copy reclustering
labels back onto the full object, and write clustering QC plots from clean
log-normalized expression.
"""

import argparse
import gzip
import json
import os
import pickle
import random
from datetime import datetime

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

from banksy.cluster_methods import run_Leiden_partition_parallel
from banksy.embed_banksy import generate_banksy_matrix
from banksy.initialize_banksy import initialize_banksy
from banksy.main import concatenate_all, median_dist_to_nearest_neighbour
from banksy.plot_banksy import plot_results
from banksy_utils.cluster_utils import create_spatial_nonspatial_adata, pad_clusters
from banksy_utils.umap_pca import pca_umap


# Match the plotting/logging defaults used by script 00 so the reclustering
# outputs look comparable to the original BANKSY run.
sc.logging.print_header()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 1
plt.rcParams["font.family"] = "Arial"
sns.set_style("white")

# BANKSY itself is deterministic; this seed controls stochastic Scanpy/UMAP and
# Leiden-related steps so repeated runs are easier to compare.
seed = 1234
np.random.seed(seed)
random.seed(seed)


def parse_args():
    """Parse command-line arguments for the reclustering workflow."""
    parser = argparse.ArgumentParser(
        prog="recluster script 05 QC-annotated Xenium AnnData with BANKSY"
    )
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="JSON config with project, dataset_name, input_label, and BANKSY settings.",
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
    """Create enough distinct categorical colours for BANKSY cluster plots."""
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


def make_numeric_cluster_category(labels):
    """Return BANKSY cluster labels as an ordered numeric categorical."""
    label_strings = labels.astype(int).astype(str)
    categories = [str(label) for label in sorted(label_strings.astype(int).unique())]
    return pd.Categorical(label_strings, categories=categories, ordered=True)


def make_banksy_label_col(nbr_weight_decay, pc_label, lambda_label, resolution):
    """Build the native BANKSY cluster-label column name for one resolution."""
    return (
        f"labels_{nbr_weight_decay}"
        f"_pc{pc_label}"
        f"_nc{lambda_label}"
        f"_r{resolution:.2f}"
    )


# Recluster labels get a suffix instead of reusing the script 00 label columns.
# That lets the full AnnData carry original and QC-reclustered labels together.
def make_recluster_label_col(
    nbr_weight_decay,
    pc_label,
    lambda_label,
    resolution,
    input_label,
):
    """Build the provenance-safe reclustered label column name."""
    banksy_col = make_banksy_label_col(
        nbr_weight_decay,
        pc_label,
        lambda_label,
        resolution,
    )
    return f"{banksy_col}_recluster_{input_label}"


def copy_obsm_aligned(source_adata, target_adata, source_key, target_key):
    """Copy an embedding between AnnData objects after aligning by cell ID."""
    # Aligning by cell ID protects against subtle row-order changes when labels
    # or embeddings move between BANKSY and clean expression objects.
    if source_key not in source_adata.obsm:
        print(f"Skipping {target_key}: {source_key!r} not found in source obsm")
        return False

    source_positions = pd.Series(np.arange(source_adata.n_obs), index=source_adata.obs_names)
    missing_cells = target_adata.obs_names.difference(source_positions.index)
    if len(missing_cells) > 0:
        raise ValueError(
            f"{len(missing_cells)} cells in target object are missing from "
            f"source embedding {source_key!r}"
        )

    aligned_idx = source_positions.loc[target_adata.obs_names].to_numpy()
    target_adata.obsm[target_key] = np.asarray(source_adata.obsm[source_key])[aligned_idx, :]
    print(f"Copied {source_key!r} to obsm[{target_key!r}]")
    return True


def resolve_scree_n_pcs(adata, selected_n_pcs, cfg):
    """Choose a valid PCA depth for scree plotting."""
    requested_n_pcs = int(cfg.get("scree_n_pcs", max(selected_n_pcs + 20, 50)))
    max_valid_n_pcs = min(adata.n_obs - 1, adata.n_vars - 1)
    scree_n_pcs = min(requested_n_pcs, max_valid_n_pcs)

    if scree_n_pcs < 2:
        print(f"Skipping PCA scree plot for object shape {adata.shape}.")
        return None

    return scree_n_pcs


def plot_pca_scree(adata, output_dir, dataset_name, selected_n_pcs, scree_n_pcs):
    """Run PCA and save a scree plot plus variance summary."""
    # Scanpy writes PCA results into the AnnData object. This function restores
    # any previous PCA slots so scree QC does not leak into BANKSY metadata.
    had_obsm_pca = "X_pca" in adata.obsm
    existing_obsm_pca = adata.obsm.get("X_pca")
    had_varm_pcs = "PCs" in adata.varm
    existing_varm_pcs = adata.varm.get("PCs")
    had_uns_pca = "pca" in adata.uns
    existing_uns_pca = adata.uns.get("pca")

    try:
        sc.tl.pca(adata, n_comps=scree_n_pcs, svd_solver="arpack", random_state=seed)
        variance_ratio = np.asarray(adata.uns["pca"]["variance_ratio"]).copy()
    finally:
        if had_obsm_pca:
            adata.obsm["X_pca"] = existing_obsm_pca
        else:
            adata.obsm.pop("X_pca", None)

        if had_varm_pcs:
            adata.varm["PCs"] = existing_varm_pcs
        else:
            adata.varm.pop("PCs", None)

        if had_uns_pca:
            adata.uns["pca"] = existing_uns_pca
        else:
            adata.uns.pop("pca", None)

    pcs = np.arange(1, len(variance_ratio) + 1)
    variance_df = pd.DataFrame(
        {
            "pc": pcs,
            "variance_ratio": variance_ratio,
            "variance_percent": variance_ratio * 100,
            "cumulative_variance_percent": np.cumsum(variance_ratio) * 100,
            "selected_for_banksy": pcs <= selected_n_pcs,
        }
    )

    scree_dir = os.path.join(output_dir, "pca_qc")
    ensure_directory(scree_dir, "PCA QC")
    variance_csv = os.path.join(scree_dir, f"{dataset_name}_pca_scree_variance.csv")
    variance_df.to_csv(variance_csv, index=False)

    fig, ax1 = plt.subplots(figsize=(9, 5))
    bars = ax1.bar(pcs, variance_df["variance_percent"], color="#4C78A8", alpha=0.85)
    ax1.set_xlabel("Principal component")
    ax1.set_ylabel("Variance explained (%)")
    ax1.set_title(f"{dataset_name}: PCA scree plot before BANKSY reclustering", fontsize=13)

    ax2 = ax1.twinx()
    line, = ax2.plot(
        pcs,
        variance_df["cumulative_variance_percent"],
        color="#222222",
        marker="o",
        linewidth=1.4,
    )
    ax2.set_ylabel("Cumulative variance explained (%)")

    selected_line = ax1.axvline(selected_n_pcs, color="#D62728", linestyle="--", linewidth=1.2)
    ax1.legend(
        [bars, line, selected_line],
        ["Variance explained per PC", "Cumulative variance explained", f"Selected PC count ({selected_n_pcs})"],
        loc="lower right",
        frameon=False,
    )

    fig.tight_layout()
    scree_png = os.path.join(scree_dir, f"{dataset_name}_pca_scree_plot.png")
    fig.savefig(scree_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved PCA scree plot to: {scree_png}")
    print(f"Saved PCA variance summary to: {variance_csv}")


def determine_max_num_labels(nonspatial_labels, spatial_labels):
    """Return the larger cluster-label count from nonspatial/spatial results."""
    if nonspatial_labels > spatial_labels:
        max_num_labels = nonspatial_labels
        print(
            f"The number of nonspatial labels {nonspatial_labels} is greater "
            f"than spatial labels {spatial_labels}; max_num_labels is "
            f"{max_num_labels}."
        )
    elif nonspatial_labels < spatial_labels:
        max_num_labels = spatial_labels
        print(
            f"The number of spatial labels {spatial_labels} is greater than "
            f"nonspatial labels {nonspatial_labels}; max_num_labels is "
            f"{max_num_labels}."
        )
    else:
        max_num_labels = spatial_labels
        print(
            f"Nonspatial and spatial runs generated the same number of labels; "
            f"max_num_labels is {max_num_labels}."
        )

    return int(max_num_labels)


def resolve_obs_column(adata, candidates, column_label):
    """Return the first available `.obs` column from a list of candidates."""
    for column in candidates:
        if column in adata.obs.columns:
            return column

    raise KeyError(
        f"No {column_label} column found. Checked candidates: {candidates}. "
        f"Available obs columns: {list(adata.obs.columns)}"
    )


def plot_cluster_count_violin(
    clean_recluster_adata,
    output_path,
    dataset_name,
    pc_label,
    lambda_label,
    resolutions,
    nbr_weight_decay,
    input_label,
):
    """Plot transcript/count distributions across reclustered BANKSY clusters."""
    cluster_count_plot_path = os.path.join(output_path, "cluster_count_plot")
    ensure_directory(cluster_count_plot_path, "Cluster count plot")
    sc.settings.figdir = cluster_count_plot_path

    count_obs_key = resolve_obs_column(
        clean_recluster_adata,
        ["total_counts", "nCount_Xenium", "transcript_counts"],
        "total-count/QC",
    )
    print(f"Using {count_obs_key!r} for recluster count violin plots")

    for res in resolutions:
        res_str_single = str(res).replace(".", "p")
        groupby_key = make_recluster_label_col(
            nbr_weight_decay,
            pc_label,
            lambda_label,
            res,
            input_label,
        )
        if groupby_key not in clean_recluster_adata.obs.columns:
            raise KeyError(
                f"{groupby_key!r} was not found in clean_recluster_adata.obs. "
                f"Available obs columns: {list(clean_recluster_adata.obs.columns)}"
            )

        sc.pl.violin(
            clean_recluster_adata,
            count_obs_key,
            groupby=groupby_key,
            show=False,
            save=(
                f"_{dataset_name}_recluster_{input_label}"
                f"_pc{pc_label}_nc{lambda_label}_r{res_str_single}.png"
            ),
        )
        plt.close("all")



# -----------------------------------------------------------------------------
# Read config and resolve sample settings
# -----------------------------------------------------------------------------

args = parse_args()
with open(args.config) as f:
    cfg = json.load(f)

# These config values mirror script 00 so reclustering can use the same PC,
# lambda, resolution, spatial decay, and coordinate settings. The inclusion
# mode lets the same script run the main QC-pass-only analysis and an all-cell
# sensitivity analysis.
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
pc_dims = [int(pc_label)]
lambda_label = cfg["lambda_label"]
lambda_list = [float(lambda_label)]
res_label = cfg["res_label"]
res_label_list = res_label if isinstance(res_label, list) else [res_label]
resolutions = [float(res) for res in res_label_list]
res_str = "_".join(res_label_list)
nbr_weight_decay = cfg["nbr_weight_decay"]
coord_keys = tuple(cfg["coord_keys"])
max_workers = int(
    os.environ.get("SLURM_CPUS_PER_TASK", cfg.get("max_workers", 4))
)

# -----------------------------------------------------------------------------
# Resolve input/output paths
# -----------------------------------------------------------------------------

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
ensure_directory(processed_path, "Processed sample")
ensure_directory(output_path, "Output reclustering sample")

# Input is the full script 05 object. It contains both passing and failed cells,
# plus the boolean `qc_keep_for_reclustering` mask used for qc_pass_only runs.
input_h5ad = cfg.get(
    "input_h5ad",
    os.path.join(
        processed_path,
        f"adata_expression_clean_{dataset_name}_qc_annotated_{input_label}.h5ad",
    ),
)
if not os.path.isfile(input_h5ad):
    raise FileNotFoundError(f"Input QC-annotated AnnData does not exist: {input_h5ad}")

# -----------------------------------------------------------------------------
# Load full provenance object and create the temporary BANKSY subset
# -----------------------------------------------------------------------------

log_time(f"Reading QC-annotated clean AnnData from: {input_h5ad}")
full_adata = ad.read_h5ad(input_h5ad)

if "qc_keep_for_reclustering" not in full_adata.obs.columns:
    raise KeyError("Input AnnData must contain obs['qc_keep_for_reclustering'] from script 05.")

for col in [coord_keys[0], coord_keys[1]]:
    if col not in full_adata.obs.columns:
        raise KeyError(f"Missing coordinate column in full_adata.obs: {col}")

# BANKSY expects spatial coordinates in `.obsm`. Rebuild them
# from the Xenium x/y obs columns so the input file does not need prefilled obsm.
full_coords = full_adata.obs[[coord_keys[0], coord_keys[1]]].to_numpy()
full_adata.obsm["xy"] = full_coords
full_adata.obsm["spatial"] = full_coords

# Choose the cells BANKSY will see. The default qc_pass_only mode excludes QC
# failures from graph construction, PCA/UMAP, and Leiden clustering. The all_cells
# mode is a sensitivity run that asks whether failed cells alter clustering.
keep_mask = full_adata.obs["qc_keep_for_reclustering"].astype(bool)
if recluster_inclusion == "qc_pass_only":
    recluster_mask = keep_mask
else:
    recluster_mask = pd.Series(True, index=full_adata.obs_names)

print(f"Full object contains {full_adata.n_obs:,} cells.")
print(f"Reclustering inclusion mode: {recluster_inclusion}")
print(f"Cells used for BANKSY reclustering: {int(recluster_mask.sum()):,}")
print(f"QC-failed cells in full object: {int((~keep_mask).sum()):,}")

clean_recluster_adata = full_adata[recluster_mask].copy()
recluster_coords = clean_recluster_adata.obs[[coord_keys[0], coord_keys[1]]].to_numpy()
clean_recluster_adata.obsm["xy"] = recluster_coords
clean_recluster_adata.obsm["spatial"] = recluster_coords

# Preserve the clean expression state for marker ranking and downstream
# expression-based analyses. These objects should already be normalized/log1p
# from script 00/01/05.
if clean_recluster_adata.raw is None:
    clean_recluster_adata.raw = clean_recluster_adata.copy()
if full_adata.raw is None:
    full_adata.raw = full_adata.copy()

# -----------------------------------------------------------------------------
# PCA scree QC before BANKSY feature expansion
# -----------------------------------------------------------------------------

scree_n_pcs = resolve_scree_n_pcs(clean_recluster_adata, pc_dims[0], cfg)
if scree_n_pcs is not None:
    plot_pca_scree(clean_recluster_adata, output_path, dataset_name, pc_dims[0], scree_n_pcs)

# -----------------------------------------------------------------------------
# BANKSY graph and matrix construction on reclustered cells only
# -----------------------------------------------------------------------------

k_geom = int(cfg.get("k_geom", 15))
max_m = int(cfg.get("max_m", 1))

log_time(f"Start generating spatial weights graph for {dataset_name}.")
nbrs = median_dist_to_nearest_neighbour(clean_recluster_adata, key=coord_keys[2])
print(nbrs)
log_time(f"Finished generating spatial weights graph for {dataset_name}.")

log_time(f"Start generating spatial weights from distance for {dataset_name}.")
banksy_dict = initialize_banksy(
    clean_recluster_adata,
    coord_keys,
    k_geom,
    nbr_weight_decay=nbr_weight_decay,
    max_m=max_m,
    plt_edge_hist=False,
    plt_nbr_weights=False,
    plt_agf_angles=False,
    plt_theta=False,
)
log_time(f"Finished generating spatial weights from distance for {dataset_name}.")

log_time(f"Start generating BANKSY matrix for {dataset_name}.")
banksy_dict, banksy_matrix = generate_banksy_matrix(
    clean_recluster_adata,
    banksy_dict,
    lambda_list,
    max_m,
)
log_time(f"Finished generating BANKSY matrix for {dataset_name}.")

# Add the nonspatial branch used by script 00 for comparison in BANKSY outputs.
banksy_dict["nonspatial"] = {
    0.0: {"adata": concatenate_all([clean_recluster_adata.X], 0, adata=clean_recluster_adata)}
}

# -----------------------------------------------------------------------------
# BANKSY PCA/UMAP and Leiden clustering across all requested resolutions
# -----------------------------------------------------------------------------

log_time(f"Start PCA and UMAP embedding for {dataset_name}.")
pca_umap(banksy_dict, pca_dims=pc_dims, add_umap=True)
log_time(f"Finished PCA and UMAP embedding for {dataset_name}.")

log_time(f"Start Leiden clustering for {dataset_name}.")
results_df, max_num_labels = run_Leiden_partition_parallel(
    banksy_dict,
    resolutions,
    num_nn=int(cfg.get("num_nn", 50)),
    max_workers=max_workers,
)
log_time(f"Finished Leiden clustering for {dataset_name}.")

# -----------------------------------------------------------------------------
# BANKSY result plots and per-resolution spatial AnnData outputs
# -----------------------------------------------------------------------------

c_map = "nipy_spectral"
cluster_palette = make_cluster_palette(max_num_labels)
weights_graph = banksy_dict[f"{nbr_weight_decay}"]["weights"][0]
plot_results(
    results_df,
    weights_graph,
    c_map,
    match_labels=False,
    coord_keys=coord_keys,
    max_num_labels=max_num_labels,
    save_fig=True,
    save_seperate_fig=True,
    dataset_name=f"{dataset_name}_recluster_{run_label}",
    save_fullfig=True,
    save_path=output_path,
    color_list=cluster_palette,
)
print(results_df)

# Determine the maximum number of labels using the same explicit script 00
# pattern. These values define placeholder cluster annotation dictionaries.
nonspatial_results = results_df[results_df["decay"] == "nonspatial"]
spatial_results = results_df[results_df["decay"] == nbr_weight_decay]
nonspatial_labels = int(nonspatial_results["num_labels"].max())
spatial_labels = int(spatial_results["num_labels"].max())
max_num_labels = determine_max_num_labels(nonspatial_labels, spatial_labels)

# Generate placeholder labels for every observed cluster ID. Manual cell-type
# annotation can replace these later, but numeric placeholders keep all cluster
# outputs valid and easy to inspect immediately after reclustering.
cluster2annotation_spatial = {}
for i in range(spatial_labels):
    cluster2annotation_spatial[str(i)] = str(i)
print(cluster2annotation_spatial)

cluster2annotation_nonspatial = {}
for i in range(nonspatial_labels):
    cluster2annotation_nonspatial[str(i)] = str(i)
print(cluster2annotation_nonspatial)

pad_clusters(cluster2annotation_spatial, list(range(max_num_labels)))
pad_clusters(cluster2annotation_nonspatial, list(range(max_num_labels)))

# `create_spatial_nonspatial_adata` pulls each resolution's clustered BANKSY
# object out of `results_df`, matching the script 00 workflow.
spatial_adatas = {}
for res in resolutions:
    print(f"Creating spatial/nonspatial AnnData objects for resolution {res:.2f}")
    adata_spatial, adata_nonspatial = create_spatial_nonspatial_adata(
        results_df=results_df,
        pca_dims=pc_dims,
        lambda_list=lambda_list,
        resolutions=[res],
        cluster2annotation_spatial=cluster2annotation_spatial,
        cluster2annotation_nonspatial=cluster2annotation_nonspatial,
        weights_scheme=nbr_weight_decay,
    )
    source_col = make_banksy_label_col(nbr_weight_decay, pc_label, lambda_label, res)
    recluster_col = make_recluster_label_col(
        nbr_weight_decay,
        pc_label,
        lambda_label,
        res,
        run_label,
    )
    adata_spatial.obs[recluster_col] = make_numeric_cluster_category(adata_spatial.obs[source_col])
    spatial_adatas[res] = adata_spatial

    res_str_single = str(res).replace(".", "p")
    spatial_output_path = os.path.join(
        processed_path,
        f"adata_spatial_{dataset_name}_recluster_{run_label}_{res_str_single}.h5ad",
    )
    adata_spatial.write_h5ad(spatial_output_path)
    print(f"Wrote spatial AnnData to: {spatial_output_path}")

# Save the heavy BANKSY objects/results so reclustering can be inspected later
# without rerunning the graph/matrix/Leiden steps.
banksy_dict_path = os.path.join(
    processed_path,
    f"{dataset_name}_recluster_{run_label}_pc{pc_label}_nc{lambda_label}_r{res_str}_banksy_dict.pkl.gz",
)
with gzip.open(banksy_dict_path, "wb") as f:
    pickle.dump(banksy_dict, f)
print(f"Wrote BANKSY dictionary to: {banksy_dict_path}")

results_csv = os.path.join(
    processed_path,
    f"results_df_{dataset_name}_recluster_{run_label}_pc{pc_label}_nc{lambda_label}_r{res_str}.csv",
)
results_df.to_csv(results_csv)
print(f"Wrote BANKSY results table to: {results_csv}")

# -----------------------------------------------------------------------------
# Copy recluster labels back to full and clean reclustered-cell objects
# -----------------------------------------------------------------------------

# The cluster table is full-cell for provenance. In qc_pass_only mode, failed
# cells are carried forward as `excluded_by_qc`; in all_cells mode, every cell
# receives a BANKSY cluster label.
cluster_table = pd.DataFrame({"index": full_adata.obs_names})
for res in resolutions:
    source_col = make_banksy_label_col(nbr_weight_decay, pc_label, lambda_label, res)
    recluster_col = make_recluster_label_col(
        nbr_weight_decay,
        pc_label,
        lambda_label,
        res,
        run_label,
    )

    recluster_labels = spatial_adatas[res].obs[source_col].astype(str)
    full_labels = pd.Series(pd.NA, index=full_adata.obs_names, dtype="object")
    full_labels.loc[recluster_labels.index] = recluster_labels
    categories = sorted(recluster_labels.astype(int).astype(str).unique(), key=lambda x: int(x))

    if recluster_inclusion == "qc_pass_only":
        full_labels = full_labels.fillna("excluded_by_qc")
        categories.append("excluded_by_qc")
    elif full_labels.isna().any():
        raise ValueError("Some full_adata cells did not receive labels in all_cells mode.")

    full_adata.obs[recluster_col] = pd.Categorical(full_labels, categories=categories)
    clean_recluster_adata.obs[recluster_col] = make_numeric_cluster_category(
        recluster_labels.reindex(clean_recluster_adata.obs_names)
    )

    cluster_table[recluster_col] = full_adata.obs[recluster_col].astype(str).to_numpy()
    print(f"Copied reclustering labels to full and reclustered clean objects: {recluster_col}")

cluster_table_csv = os.path.join(
    processed_path,
    f"{dataset_name}_recluster_{run_label}_cell_cluster_id_across_clustering_res_{res_str}.csv",
)
cluster_table.to_csv(cluster_table_csv, index=False)
print(f"Wrote full-cell cluster table to: {cluster_table_csv}")

# Copy BANKSY embeddings into the clean objects. In qc_pass_only mode, full-object
# embeddings are NaN for failed cells because those cells were not embedded by
# BANKSY. In all_cells mode, every cell receives an embedding.
embedding_source = spatial_adatas[resolutions[0]]
umap_source_keys = ["X_umap", f"reduced_pc_{pc_dims[0]}_umap", f"reduced_pc_{pc_label}_umap"]
pca_source_keys = ["X_pca", "pca", "X_pca_banksy", f"reduced_pc_{pc_dims[0]}", f"reduced_pc_{pc_label}"]
umap_target_key = f"X_umap_{nbr_weight_decay}_pc{pc_label}_nc{lambda_label}_recluster_{run_label}"
pca_target_key = f"X_pca_{nbr_weight_decay}_pc{pc_label}_nc{lambda_label}_recluster_{run_label}"

for source_key in umap_source_keys:
    if copy_obsm_aligned(embedding_source, clean_recluster_adata, source_key, umap_target_key):
        full_embedding = np.full(
            (full_adata.n_obs, clean_recluster_adata.obsm[umap_target_key].shape[1]),
            np.nan,
        )
        full_pos = pd.Series(np.arange(full_adata.n_obs), index=full_adata.obs_names)
        full_embedding[full_pos.loc[clean_recluster_adata.obs_names].to_numpy(), :] = (
            clean_recluster_adata.obsm[umap_target_key]
        )
        full_adata.obsm[umap_target_key] = full_embedding
        break
else:
    print("No UMAP embedding found to copy")

for source_key in pca_source_keys:
    if copy_obsm_aligned(embedding_source, clean_recluster_adata, source_key, pca_target_key):
        full_embedding = np.full(
            (full_adata.n_obs, clean_recluster_adata.obsm[pca_target_key].shape[1]),
            np.nan,
        )
        full_pos = pd.Series(np.arange(full_adata.n_obs), index=full_adata.obs_names)
        full_embedding[full_pos.loc[clean_recluster_adata.obs_names].to_numpy(), :] = (
            clean_recluster_adata.obsm[pca_target_key]
        )
        full_adata.obsm[pca_target_key] = full_embedding
        break
else:
    print("No PCA embedding found to copy")

# -----------------------------------------------------------------------------
# Save clean expression AnnData objects with reclustering metadata
# -----------------------------------------------------------------------------

clean_pass_h5ad = os.path.join(
    processed_path,
    f"adata_expression_clean_{dataset_name}_recluster_{run_label}_cells_used_for_banksy_with_clusters_{res_str}.h5ad",
)
full_clustered_h5ad = os.path.join(
    processed_path,
    f"adata_expression_clean_{dataset_name}_qc_annotated_{input_label}_with_banksy_reclusters_{run_label}_{res_str}.h5ad",
)

# -----------------------------------------------------------------------------
# Marker ranking on clean expression, not BANKSY-expanded matrices
# -----------------------------------------------------------------------------

marker_method = cfg.get("marker_method", "wilcoxon")
n_marker_genes = int(cfg.get("n_marker_genes", 20))
marker_table_path = os.path.join(output_path, "top_marker_tables")
ensure_directory(marker_table_path, "Top marker table")
ranked_marker_keys = {}

for res in resolutions:
    res_str_single = str(res).replace(".", "p")
    groupby_key = make_recluster_label_col(
        nbr_weight_decay,
        pc_label,
        lambda_label,
        res,
        run_label,
    )
    markers_key = f"{groupby_key}_markers_{marker_method}"

    sc.tl.rank_genes_groups(
        clean_recluster_adata,
        groupby=groupby_key,
        method=marker_method,
        key_added=markers_key,
        use_raw=False,
    )
    marker_df = sc.get.rank_genes_groups_df(
        clean_recluster_adata,
        group=None,
        key=markers_key,
    )
    marker_csv = os.path.join(
        marker_table_path,
        f"rank_genes_groups_{dataset_name}_recluster_{run_label}_pc{pc_label}_nc{lambda_label}_r{res_str_single}_{marker_method}.csv",
    )
    marker_df.to_csv(marker_csv, index=False)
    ranked_marker_keys[res] = markers_key
    print(f"res {res:.2f}: wrote marker table {marker_csv}")

clean_recluster_adata.write_h5ad(clean_pass_h5ad)
print(f"Wrote clean expression AnnData used for BANKSY with reclustering metadata: {clean_pass_h5ad}")

full_adata.write_h5ad(full_clustered_h5ad)
print(f"Wrote full provenance AnnData with reclustering metadata: {full_clustered_h5ad}")

# Plot the marker rankings for visual review, matching the script 00 output type.
top_marker_plot_path = os.path.join(output_path, "top_marker_plot")
ensure_directory(top_marker_plot_path, "Top marker plot")
for res, markers_key in ranked_marker_keys.items():
    res_str_single = str(res).replace(".", "p")
    output_png = os.path.join(
        top_marker_plot_path,
        f"rank_genes_groups_{dataset_name}_recluster_{run_label}_pc{pc_label}_nc{lambda_label}_r{res_str_single}_{marker_method}.png",
    )
    sc.pl.rank_genes_groups(
        clean_recluster_adata,
        key=markers_key,
        n_genes=n_marker_genes,
        fontsize=15,
        show=False,
    )
    fig = plt.gcf()
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"res {res:.2f}: wrote marker plot {output_png}")

# -----------------------------------------------------------------------------
# Cluster-count QC plots on clean reclustered-cell expression object
# -----------------------------------------------------------------------------

# These violin plots mirror script 00 and help assess whether reclustered groups
# are driven by strong count/transcript differences rather than biology.
plot_cluster_count_violin(
    clean_recluster_adata=clean_recluster_adata,
    output_path=output_path,
    dataset_name=dataset_name,
    pc_label=pc_label,
    lambda_label=lambda_label,
    resolutions=resolutions,
    nbr_weight_decay=nbr_weight_decay,
    input_label=run_label,
)


log_time(f"Finished QC-annotated BANKSY reclustering for {dataset_name}.")
