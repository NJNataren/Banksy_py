#!/usr/bin/env python
# coding: utf-8

# In[11]:


#!/usr/bin/env python
# coding: utf-8

"""
Title: Xenium BANKSY Clustering With Clean Expression Outputs
Date: 2026-07-30
Summary: Run config-driven BANKSY clustering for Xenium samples, save clean
log-normalized expression AnnData before BANKSY feature expansion, and copy
BANKSY metadata back onto the clean object for marker analysis.
"""


# In[12]:


# # Xenium spatial clustering with BANKSY
# ### **Author:** Nathalie Nataren
# ### **Date:** 17/04/2026
# 
# **Description:** The purpose of this analysis is to perform clustering of Xenium spatial data as part of the metastatic melanoma ICI therapy response study (VBCT lab).
# 




# In[13]:


import anndata as ad
import os
import numpy as np
import pandas as pd
from IPython.display import display

import scipy.sparse as sparse
from scipy.io import mmread
from scipy.stats import pearsonr, pointbiserialr

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib import rc_context
from datetime import datetime


import seaborn as sns
import scanpy as sc
sc.logging.print_header()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 1 # errors (0), warnings (1), info (2), hints (3)
plt.rcParams["font.family"] = "Arial"
sns.set_style("white")

import random
import warnings
warnings.filterwarnings("ignore") 

# Note that BANKSY itself is deterministic, here the seeds affect the umap clusters and leiden partition
seed = 1234
np.random.seed(seed)
random.seed(seed)


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
        extra_colors = [
            mpl.colors.to_hex(extra_cmap(i / max(extra_count - 1, 1)))
            for i in range(extra_count)
        ]
        colors.extend(extra_colors)

    return colors[:n_colors]


def make_numeric_cluster_category(labels):
    """Return BANKSY cluster labels as an ordered numeric categorical."""
    label_strings = labels.astype(int).astype(str)
    categories = [str(label) for label in sorted(label_strings.astype(int).unique())]
    return pd.Categorical(label_strings, categories=categories, ordered=True)


def make_banksy_label_col(nbr_weight_decay, pc_label, lambda_label, resolution):
    """Build the BANKSY cluster-label column name for one resolution."""
    return (
        f"labels_{nbr_weight_decay}"
        f"_pc{pc_label}"
        f"_nc{lambda_label}"
        f"_r{resolution:.2f}"
    )


def copy_obsm_aligned(source_adata, target_adata, source_key, target_key):
    """Copy an embedding between AnnData objects after aligning by cell ID."""
    if source_key not in source_adata.obsm:
        print(f"Skipping {target_key}: {source_key!r} not found in source obsm")
        return False

    source_positions = pd.Series(
        np.arange(source_adata.n_obs),
        index=source_adata.obs_names,
    )
    missing_cells = target_adata.obs_names.difference(source_positions.index)
    if len(missing_cells) > 0:
        raise ValueError(
            f"{len(missing_cells)} cells in clean_adata are missing from "
            f"source embedding {source_key!r}"
        )

    aligned_idx = source_positions.loc[target_adata.obs_names].to_numpy()
    target_adata.obsm[target_key] = np.asarray(source_adata.obsm[source_key])[aligned_idx, :]
    print(f"Copied {source_key!r} to clean_adata.obsm[{target_key!r}]")
    return True


def ensure_directory(path, label):
    """Create `path` and any missing parent directories, then log its location."""
    os.makedirs(path, exist_ok=True)
    print(f"{label} directory ready: {os.path.abspath(path)}")


def resolve_obs_column(adata, candidates, column_label):
    """Return the first available `.obs` column from a list of candidates."""
    for column in candidates:
        if column in adata.obs.columns:
            return column

    raise KeyError(
        f"No {column_label} column found. Checked candidates: {candidates}. "
        f"Available obs columns: {list(adata.obs.columns)}"
    )


def resolve_scree_n_pcs(adata, selected_n_pcs, cfg):
    """Choose a valid PCA depth for scree plotting beyond the selected PC count."""
    requested_n_pcs = int(cfg.get("scree_n_pcs", max(selected_n_pcs + 20, 50)))
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


def plot_pca_scree(adata, output_dir, dataset_name, selected_n_pcs, scree_n_pcs):
    """Run PCA and save scree plot plus variance summary for clean expression data."""
    # Scanpy stores PCA results on the AnnData object. Keep this scree-only PCA
    # from leaking into later BANKSY embedding transfer by restoring prior slots.
    had_obsm_pca = "X_pca" in adata.obsm
    existing_obsm_pca = adata.obsm.get("X_pca")
    had_varm_pcs = "PCs" in adata.varm
    existing_varm_pcs = adata.varm.get("PCs")
    had_uns_pca = "pca" in adata.uns
    existing_uns_pca = adata.uns.get("pca")

    try:
        sc.tl.pca(
            adata,
            n_comps=scree_n_pcs,
            svd_solver="arpack",
            random_state=seed,
        )
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




# In[14]:


# ##################################
# #   ** LOCAL TESTING BLOCK **    #
# ##################################

# ## These values mirror config/00_clustering/vbct/small/CK_skin_res.json.
# ## Leave this block commented when using --config from Slurm/the shell.
# dataset_name = "CK_skin_res" # sample name
# pc_label = "30" # Label for the number of principal components used for filenames
# pc_dims = [int(pc_label)] # Number of principal components stored as a list
# lambda_label = "0.20" # BANKSY lambda label
# lambda_list = [float(lambda_label)] # BANKSY lambda as a float
# res_label = ["0.70", "0.80", "0.90", "1.00"] # BANKSY clustering resolutions
# resolutions = [float(res) for res in res_label]
# nbr_weight_decay = "scaled_gaussian"
# coord_keys = ("x", "y", "xy") # Coordinate keys in .obs/.obsm
# project = "vbct" # Project/study output folder
# raw_subdir = "vbct" # Subdirectory under data/xenium/raw_data
# max_workers = 8 # Maximum CPUs for Leiden clustering




# In[15]:


########################
#   PARSE ARGUMENTS    #
########################
# This block of code feeds arguments to this python script from a config files found in /config

## Import argparse and json packages to read in variables from the per sample .json config files 
import argparse
import json

parser = argparse.ArgumentParser(prog="used to parse arguments from 00_xenium_clustering_clean_adata.py") # Initialise the parser
parser.add_argument("--config", type=str, help="Optional JSON config file for each Xenium sample", required=False) # This defines the flag and tells the script to look for a JSON config
# parse_known_args keeps notebook/kernel arguments from breaking local runs.
args, _unknown_args = parser.parse_known_args()

if args.config:
    with open(args.config) as f: # Opens the file path provided by the user
        cfg = json.load(f) # Converts the json config into a python dictionary called "cfg"
else:
    print("No --config supplied; using CK_skin_res local testing defaults.")
    cfg = {
        "project": "vbct",
        "raw_subdir": "vbct",
        "dataset_name": "CK_skin_res",
        "pc_label": "30",
        "lambda_label": "0.20",
        "res_label": ["0.70", "0.80", "0.90", "1.00"],
        "nbr_weight_decay": "scaled_gaussian",
        "coord_keys": ["x", "y", "xy"],
        "max_workers": 8,
        "scree_n_pcs": 75,
    }

## Set the dataset_name and related settings to use during this analysis by taking the argument values from the "cfg" dictionary read in from the JSON config
dataset_name = cfg["dataset_name"] # sample name
pc_label = cfg["pc_label"] # Label for the number of principal components used for the purpose of filenames
pc_dims = [int(pc_label)] # The number of principal components stored a list for analyses
lambda_label = cfg["lambda_label"] # File name label for Lambda setting, see comment below. 
lambda_list = [float(lambda_label)] # Lambda setting to tune BANKSY clustering, lambda = 0 is non-spatial, 0.2 is for cell typing, 0.8 if for domain segmentation 
res_label = cfg["res_label"] # BANKSY clustering resolution label for resolution chosen to produce plots
resolutions = [float(res) for res in cfg["res_label"]] # BANSY can take a list of resolutions and perform clustering at each which is saved in the BANKSY dictionary
nbr_weight_decay = cfg["nbr_weight_decay"] # This parameter dictates how much neighbouring cells impact to the neighbourhood expression calculations. Using scaled gaussian, the 
# close neigbours contribute more and this decays as you move out to cells further away in the neighbourhood window. It is scaled for local cell density so that weighting doesn't change
# across regions if cells are pack more closely or loosely in different regions
coord_keys = tuple(cfg["coord_keys"]) # Keys to specify coordinate indexes in the anndata Object
raw_subdir = cfg.get("raw_subdir", "") # Optional subdirectory under data/xenium/raw_data, e.g. ptmt/Run_1
project = cfg.get("project", "")
max_workers = int(os.environ.get("SLURM_CPUS_PER_TASK", cfg.get("max_workers", 4))) # Parameter for the run_Leiden_partition_parallel() clustering function 




# In[16]:


########################
#       SET PATHS      #
########################
## Set file paths and read in xenium data

## Create a base path
base_dir = "data/xenium"

## Create a path to the raw data e.g., unprocessed anndata files.
## `raw_subdir` lets each config point to a run-specific folder while keeping
## older configs valid when the raw files live directly under raw_data.
raw_path = os.path.join(base_dir, "raw_data", raw_subdir)

if not os.path.isdir(raw_path):
    raise FileNotFoundError(
        f"Raw input directory '{raw_path}' does not exist. "
        "Check the config raw_subdir value and mounted data location."
    )
    
else:
    print(f"Directory '{raw_path} exists.")

## Create sample-specific processed/output folders before any writes. This
## protects HPC and notebook runs where project or QC_testing parents are absent.
processed_path = os.path.join(base_dir, "processed", project, dataset_name)
ensure_directory(processed_path, "Processed sample")

## Create a path for output data, if it does not already exist
output_path = os.path.join(base_dir, "output", project, dataset_name)
ensure_directory(output_path, "Output sample")

## Create a path for QC results, if it does not already exist
qc_path = os.path.join(base_dir, "output", project, "QC_testing", dataset_name)
ensure_directory(qc_path, "QC testing sample")




# In[17]:


## Function to log sub task start time
def log_time(step):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {step}")




# In[18]:


###########################
#       LOAD ANNDATA      #
###########################

## Read in the raw AnnData file
raw_adata_path = os.path.join(raw_path, f"{dataset_name}_raw.h5ad")
if not os.path.isfile(raw_adata_path):
    raise FileNotFoundError(
        f"Raw AnnData file '{raw_adata_path}' does not exist. "
        "Check dataset_name and raw_subdir in the config."
    )

adata = ad.read_h5ad(raw_adata_path)
log_time(f"Loading in data for {dataset_name}")

## Create 'xy' spatial coordinates from adata.obs
adata.obsm['xy'] = np.vstack([adata.obs['x'], adata.obs['y']]).T




# In[19]:


res_label




# In[20]:


############################################
#       PLOT ZERO COUNT / EMPTY CELLS      #
############################################

zero_count_mask = adata.obs["nCount_Xenium"] <= 0
adata.obs["zero_count_cell_cat"] = (
    zero_count_mask
    .map({True: "Fail", False: "Pass"})
    .astype("category")
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

# Draw Pass cells first and failed empty cells second so sparse failures remain visible.
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


# In[ ]:


#####################################
#       FILTER ZERO COUNT CELLS     #
#####################################

## Filter out cells with zero counts
adata = adata[adata.obs['nCount_Xenium'] > 0].copy()




# In[ ]:


####################################
#       DOWNCAST ANNDATA FILES     #
####################################

### Down cast the anndata files to 32-bit float type to reduce RAM load
## Function to check the float size of the adata.obsm data
from banksy_utils.object_downcasting_utils import check_float, downcast_float
check_float(adata)

## Function to downcast adata.obsm to 32 bit and check that the resulting array is an N-dimensional array
downcast_float(adata, "float32")

# Write the downcast file to AnnData format 
adata.write_h5ad(
    filename = os.path.join(processed_path, f"{dataset_name}_float_32.h5ad"),
    compression="gzip"
    )

# Define float adata file
float_adata = f"{dataset_name}_float_32.h5ad"




# In[ ]:


#######################################
#       LOAD DATA AND COORDINATES     #
#######################################

from banksy_utils.load_data import load_adata, display_adata

## To either load data from .h5ad directly or convert raw data to .h5ad format
load_adata_directly = True

## Keys to specify coordinate indexes in the anndata Object
coord_keys = coord_keys

raw_y, raw_x, adata = load_adata(filepath=processed_path, adata_filename=float_adata, load_adata_directly=True, coord_keys=coord_keys)




# In[ ]:


# ---------------------------------------------------------------------------- #
#                        STORE RAW COUNTS AND NORMALISE                        #
# ---------------------------------------------------------------------------- #

## Save the raw UMI integer counts
adata.layers["counts"] = adata.X.copy()

## Inspect the counts
print(adata.layers["counts"][:5,:5])




# In[ ]:


from banksy_utils.filter_utils import normalize_total, filter_hvg, print_max_min

## Normalizes the anndata dataset using the BANKSY normalise_total() function
normalize_total(adata)
print(adata.X)




# In[ ]:


## Perform log-transformation and save the normalised and log-transformed data in adata.raw
sc.pp.log1p(adata)
print(adata.X)

selected_n_pcs = int(pc_label)
scree_n_pcs = resolve_scree_n_pcs(adata, selected_n_pcs, cfg)
if scree_n_pcs is not None:
    plot_pca_scree(
        adata=adata,
        output_dir=output_path,
        dataset_name=dataset_name,
        selected_n_pcs=selected_n_pcs,
        scree_n_pcs=scree_n_pcs,
    )

adata.raw = adata.copy() ## This needs to be moved to after the creation of adata_spatial to truly save it in the BANKSY created anndata object




# In[ ]:


###############################################################
#       SAVE CLEAN EXPRESSION ANNDATA BEFORE BANKSY MATRIX     #
###############################################################

## Save a clean cell-by-gene expression object after total-count
## normalization and log1p transformation, before BANKSY feature expansion.
## This keeps expression values suitable for downstream marker summaries.
clean_expression_h5ad = os.path.join(
    processed_path,
    f"adata_expression_clean_{dataset_name}_normalized_log1p.h5ad",
)

adata.write_h5ad(clean_expression_h5ad)
print(f"Wrote clean expression AnnData: {clean_expression_h5ad}")
print(adata)




# In[ ]:


# ## Generate spatial weights graph
# #### In BANKSY, we imagine edges / connections in the graph to represent neighbour relationships between cells
# 
# In doing so, the banksy algorithm requires the following specifications in the main BANKSY algorithm:
# 
# 1. The number of spatial neighbours `num_neighbours`, this is also known as the $k_{geom}$ parameter in the manuscript.
#     
# 2. Assigning weights (dependent on inter-cell spatial distances) to edges in the connected spatial graph. By default, we use the `gaussian decay` option, where weights decay as a function of distnace to the index cell with $\sigma$ = `sigma`. As default, we set $\sigma$ to be the median distance between cells, and do not prescribe any cutoff-radius `p_outside` (no cutoff is conducted).
# 
# 3. The Azumithal Gabor Filter parameter `max_m` which indicates whether to use the AGF (`max_m = 1`) or just the mean expression (`max_m = 0`). By default, we set `max_m = 1`.
#     
# ### Construction of the spatial $k_{geom}$ Nearest-Neighbour graph
# 
# We represent connections between cells and its neighbours in a graph $G = \{N,E,W\}$, comprising of a set of nodes $n \in N$. edges representing connectivity between cells $e \in E$, the edges are be weighted $w \in W$ as a function of the spatial distance between cells.
#     
# Weight of edges can be represented by uniform distance (i.e., a closer neighbour will have a higher weight), or using `reciprocal` ($\frac{1}{r})$. As mentioned above, BANKSY by default applies a gaussian envelope function to map the distance (between nodes) to the `weights` the connection of cell-to-neighbor

from banksy.main import median_dist_to_nearest_neighbour

# set params
# ==========
plot_graph_weights = True
k_geom = 15 # only for fixed type
max_m = 1 # azumithal transform up to kth order
nbr_weight_decay = nbr_weight_decay # can also be "reciprocal", "uniform" or "ranked"

# Find median distance to closest neighbours, the median distance will be `sigma`
log_time(f"Start generating spatial weights graph for {dataset_name}.")
nbrs = median_dist_to_nearest_neighbour(adata, key = coord_keys[2])
log_time(f"Finished generating spatil weights graph for {dataset_name}.")





# In[ ]:


# ### Generate spatial weights from distance
# 
# Here, we generate the spatial weights using the gaussian decay function from the median distance to the k-th nearest neighbours as specified earlier.
# 
# The utility functions `plot_edge_histograms`, `plot_weights`, `plot_theta_graph` can be used to visualize the characteristics of the edges, weights and theta (from the AGF) respectively.
# 
# ### Optional Visualization of Weights and edge graphs
# 
# (1) **Visualize the edge histogram** to show the histogram of distances between cells and the weights between cells by setting `plt_edge_hist = True`
# 
# (2) **Visualize weights** by plotting the connections. Line thickness is proportional to edge weight (normalized by highest weight across all edges) by setting `plt_weights = True`
# 
# (3) **Visualize weights with Azimuthal angles**. Visualizing the azimuthal angles computed by the AGF by colour, the azimuthal connectivities are annotated in red. Warning: this plot many take a few minutes to compute for large datasets. By default, `plt_agf_angles = False`
# 
# (4) **Visualize angles around random cell**. Plot points around a random index cell, annotated with angles from the index cell. `plot_theta = True`

from banksy.initialize_banksy import initialize_banksy
log_time(f"Start generating spatial weights from distance for {dataset_name}.")
banksy_dict = initialize_banksy(
    adata,
    coord_keys,
    k_geom,
    nbr_weight_decay=nbr_weight_decay,
    max_m=max_m,
    plt_edge_hist=False,
    plt_nbr_weights=False,
    plt_agf_angles=False, # takes long time to plot
    plt_theta=False,
)
log_time(f"Finished generating spatial weights from distance for {dataset_name}.")




# In[ ]:


# ## Generate BANKSY matrix
# 
# To generate the BANKSY matrix, we proceed with the following:
# 
# 1. Matrix multiply sparse CSR weights matrix with cell-gene matrix to get **neighbour matrix** and the **AGF matrix** if `max_m` > 1
# 2. Z-score both matrices along **genes**
# 3. Multiply each matrix by a weighting factor $\lambda$ (We refer to this parameter as `lambda` in our manuscript and code)
# 4. Concatenate the matrices along the genes dimension in the form -> `horizontal_concat(cell_mat, nbr_mat, agf_mat)`
# 
# Here, we save all the results in the dictionary (`banksy_dict`), which contains the results from the subsequent operations for BANKSY. 

from banksy.embed_banksy import generate_banksy_matrix
banksy_dict, banksy_matrix = generate_banksy_matrix(adata, banksy_dict, lambda_list, max_m)
log_time(f"Finished generating BANKSY matrix for {dataset_name}.")




# In[ ]:


# ### Append Non-spatial results to the `banksy_dict` for comparsion

from banksy.main import concatenate_all

banksy_dict["nonspatial"] = {
    # Here we simply append the nonspatial matrix (adata.X) to obtain the nonspatial clustering results
    0.0: {"adata": concatenate_all([adata.X], 0, adata=adata), }
}

print(banksy_dict['nonspatial'][0.0]['adata'])




# In[ ]:


## Perform UMAP embedding
from banksy_utils.umap_pca import pca_umap
log_time(f"Start PCA and UMAP embedding for {dataset_name}.")
pca_umap(banksy_dict,
         pca_dims = pc_dims,
         add_umap = True,
         #plt_remaining_var = True,
         )
log_time(f"Finish PCA and UMAP embedding for {dataset_name}.")




# In[ ]:


# ### Cluster cells using a partition algorithm
# 
# For the purpose of this dataset, we cluster cells using the `leiden` algorithm (use `!pip install leidenalg` if package missing) partition methods. Other clustering algorithms include `louvain` (another resolution based clustering algorithm), or `mclust` (a clustering based on gaussian mixture model). 
# 
# Note that by default, we recommend resolution-based clustering (i.e., `leiden` or `louvain`) if no prior information on the number of clusters known. However, if the number of clusters is known *a priori*, the user can use `mclust` (gaussian-mixture model) by specifying the number of clusters beforehand.

from banksy.cluster_methods import run_Leiden_partition_parallel
log_time(f"Start Leiden clustering for {dataset_name}.")
## Parallelised leiden partitioning function for clustering, if a computer has n cores, set max_workers to n-1
results_df, max_num_labels = run_Leiden_partition_parallel(
    
    banksy_dict,
    resolutions,
    num_nn=50,
    max_workers=max_workers # 8-core laptop
)
log_time(f"Finished Leiden clustering for {dataset_name}.")




# In[ ]:


# ## Dynamically extract the number of principal components from the results_df
# pc_dims = results_df[results_df['decay'] == nbr_weight_decay]
# pc_dims = pc_label['num_pcs'].iloc[0]




# In[ ]:


# ## Plot results
# 
# ### Visualize the clustering results from BANKSY, including the clusters from the Umap embbedings

from banksy.plot_banksy import plot_results

## Use enough categorical colours for all cluster IDs. The continuous cmap is
## still used by helper plots such as UMAP/PCA/connectivity, while color_list
## drives the main spatial cluster map and per-cluster panels.
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
    #save_path = os.path.join(file_path, 'tmp_png'),
    save_fig=True, # save the spatial map of all clusters
    save_seperate_fig=True, # save the figure of all clusters plotted seperately
    dataset_name=f"{dataset_name}",
    save_fullfig=True,
    save_path=output_path,
    color_list=cluster_palette,
)

print(results_df)




# In[ ]:


max_num_labels




# In[ ]:


##########################################################
#       DETERMINE MAX LABELS FOR CLUSTER ANNOTATION      #
##########################################################

## Identify the largest label number and set that to max_num_labels
## extract the seperate spatial and non-spatial results
nonspatial_results = results_df[results_df["decay"] == "nonspatial"]
spatial_results = results_df[results_df["decay"] == "scaled_gaussian"]

nonspatial_labels = nonspatial_results["num_labels"].max()
spatial_labels = spatial_results["num_labels"].max()

def determine_max_num_labels(nonspatial_labels, spatial_labels):
## extract the number of labels 

    if nonspatial_labels > spatial_labels:
        max_num_labels = nonspatial_labels
        print(f"The number of nonspatial labels {nonspatial_labels} is greater than spatial labels {spatial_labels}, so {nonspatial_labels} is assigned to max_num_labels which {max_num_labels}.")

    elif nonspatial_labels < spatial_labels:
        max_num_labels = spatial_labels
        print(f"The number of spatial labels {spatial_labels} is greater than nonspatial labels {nonspatial_labels}, so {spatial_labels} is assigned to max_num_labels which is {max_num_labels}.")

    else: 
        max_num_labels = spatial_labels
        print(f"nonspatial and spatial decay have generated the same number of clusters, so {nonspatial_labels} is assigned to max_num_labels which is {max_num_labels}.")

    return max_num_labels


max_num_labels=determine_max_num_labels(nonspatial_labels, spatial_labels)




# In[ ]:


print(max_num_labels)




# In[ ]:


########################################################
#       GENERATE PLACE HOLDER LABELS FOR CLUSTERS      #
########################################################

## Generate place holder labels for clusters
from banksy_utils.cluster_utils import pad_clusters, create_spatial_nonspatial_adata
## Spatial cluster place holders
cluster2annotation_spatial = {}

for i in range(spatial_labels):
    cluster2annotation_spatial[str(i)] = str(i)
print(cluster2annotation_spatial)

## Non-spatial cluster place holders
cluster2annotation_nonspatial = {}

for i in range(nonspatial_labels):
    cluster2annotation_nonspatial[str(i)] = str(i)
print(cluster2annotation_nonspatial)

pad_clusters(cluster2annotation_spatial, list(range(max_num_labels)))




# In[ ]:


print("results_df.index:")
for idx in results_df.index.tolist():
    print(f"  {repr(idx)}")
print(f"\nLooking for: 'scaled_gaussian_pc{pc_dims[0]}_nc{lambda_list[0]:0.2f}_r0.50'")
print(f"pc_dims = {pc_dims}")
print(f"lambda_list = {lambda_list}")
print(f"resolutions = {resolutions}")




# In[ ]:


#########################################################
#       CREATE SPATIAL AND NONSPATIAL ADATA OBJECTS     #
#########################################################
# Loop through the list of desired clustering resolutions and generate spatial and non-spatial object in a dictionary
adata_dict = {}
for resolution in resolutions:
    print(f"\n--- Creating adata_spatial object and processing resolution: {resolution} ---")
    adata_spatial, adata_nonspatial = create_spatial_nonspatial_adata(
        results_df=results_df,
        pca_dims = pc_dims,
        lambda_list= lambda_list,
        resolutions= [resolution],
        cluster2annotation_spatial=cluster2annotation_spatial,
        cluster2annotation_nonspatial=cluster2annotation_nonspatial,
      )
    adata_dict[resolution] = {"spatial": adata_spatial, "nonspatial": adata_nonspatial}




# In[ ]:


res_label


# 



# In[ ]:


## Create a flat dictionary of the spatial and non-spatial 
# spatial
spatial_adatas = {}
for res, adatas in adata_dict.items():
    spatial_adatas[res]  = adatas["spatial"]




# In[ ]:


# Save individaul anndata objects at each resolution
for res in resolutions:
    res_str  = str(res).replace(".", "p")
    spatial_output_path = os.path.join(
        processed_path,
        f"adata_spatial_{dataset_name}_{res_str}.h5ad",
    )
    spatial_output_dir = os.path.dirname(spatial_output_path)
    print(f"Writing spatial AnnData to: {spatial_output_path}")
    print(f"Spatial output parent directory exists: {os.path.isdir(spatial_output_dir)}")
    if not os.path.isdir(spatial_output_dir):
        raise FileNotFoundError(
            f"Spatial output parent directory does not exist: {spatial_output_dir}"
        )
    spatial_adatas[res].write_h5ad(spatial_output_path)




# In[ ]:


import gzip
import pickle

## Export dictionary
res_str = "_".join(res_label)

## Save the banksy dict as a pickle file to load it in later and avoid having to initialize repeatedly to save time
## Use gzip to save it with compression
#for res in resolutions:
with gzip.open(os.path.join(processed_path, f"{dataset_name}_pc{pc_label}_nc{lambda_label}_r{res_str}_banksy_dict.pkl.gz"), "wb") as f:
    pickle.dump(banksy_dict, f)

## Export the results_df data frame

results_df.to_csv(os.path.join(processed_path,f"results_df_{dataset_name}_pc{pc_label}_nc{lambda_label}_r{res_str}.csv"))




# In[ ]:


#############################################
#   Generate a cell cluster identity table  #
#############################################

## Grab the cell ids from the raw anndata object 
raw_obs = adata.obs
raw_cell_ids = raw_obs.reset_index()
raw_cell_ids = raw_cell_ids[['index']]

## Interate over all the spatial adatas and merge the .obs cluster ids to the index of cell ids
merged = raw_cell_ids.copy()

for res in resolutions:
    cluster_label_col = make_banksy_label_col(
        nbr_weight_decay,
        pc_label,
        lambda_label,
        res,
    )
    df = spatial_adatas[res].obs.reset_index()
    df = df[["index", cluster_label_col]]
    merged = pd.merge(merged, df, on="index", how="left")
print(merged)

merged.to_csv(os.path.join(processed_path, f"{dataset_name}_cell_cluster_id_across_clustering_res_{res_str}.csv"))




# In[ ]:


#####################################################################
#       COPY BANKSY METADATA TO CLEAN EXPRESSION ANNDATA              #
#####################################################################

## Read the clean expression object saved before BANKSY feature expansion.
## Cluster labels and embeddings are aligned by cell ID so metadata remains
## correct even if object row order changes. Clean expression values are not
## replaced with BANKSY-expanded features.
try:
    clean_expression_h5ad
except NameError:
    clean_expression_h5ad = os.path.join(
        processed_path,
        f"adata_expression_clean_{dataset_name}_normalized_log1p.h5ad",
    )

clean_adata = ad.read_h5ad(clean_expression_h5ad)
clean_shape_before = clean_adata.shape
clean_var_names_before = clean_adata.var_names.copy()

cluster_label_table = merged.set_index("index")
cluster_label_cols = [
    make_banksy_label_col(nbr_weight_decay, pc_label, lambda_label, res)
    for res in resolutions
]

for cluster_label_col in cluster_label_cols:
    if cluster_label_col not in cluster_label_table.columns:
        raise KeyError(
            f"{cluster_label_col!r} was not found in the cell cluster identity table. "
            f"Available columns: {list(cluster_label_table.columns)}"
        )

    labels = cluster_label_table[cluster_label_col].reindex(clean_adata.obs_names)
    missing_labels = int(labels.isna().sum())
    if missing_labels > 0:
        raise ValueError(
            f"{missing_labels} cells in the clean expression object did not receive "
            f"a BANKSY label from {cluster_label_col!r}"
        )

    # Store cluster IDs as an ordered categorical sorted numerically so Scanpy
    # plots groups as 0, 1, 2, ..., 10 instead of lexical order 0, 1, 10, ... .
    clean_adata.obs[cluster_label_col] = make_numeric_cluster_category(labels)
    print(f"Copied {cluster_label_col} to clean expression AnnData obs")

embedding_source = spatial_adatas[resolutions[0]]
umap_target_key = f"X_umap_{nbr_weight_decay}_pc{pc_label}_nc{lambda_label}"
umap_source_keys = [
    "X_umap",
    f"reduced_pc_{pc_dims[0]}_umap",
    f"reduced_pc_{pc_label}_umap",
]
for source_key in umap_source_keys:
    if copy_obsm_aligned(
        source_adata=embedding_source,
        target_adata=clean_adata,
        source_key=source_key,
        target_key=umap_target_key,
    ):
        break
else:
    print("No UMAP embedding found to copy")

pca_target_key = f"X_pca_{nbr_weight_decay}_pc{pc_label}_nc{lambda_label}"
pca_source_keys = [
    "X_pca",
    "pca",
    "X_pca_banksy",
    f"reduced_pc_{pc_dims[0]}",
    f"reduced_pc_{pc_label}",
]
for source_key in pca_source_keys:
    if copy_obsm_aligned(
        source_adata=embedding_source,
        target_adata=clean_adata,
        source_key=source_key,
        target_key=pca_target_key,
    ):
        break
else:
    print("No PCA embedding found to copy")

required_spatial_cols = ["x", "y"]
missing_spatial_cols = [
    col for col in required_spatial_cols
    if col not in clean_adata.obs.columns
]
if missing_spatial_cols:
    raise KeyError(
        f"Cannot create spatial coordinates. Missing obs columns: {missing_spatial_cols}"
    )

spatial_coords = clean_adata.obs[["x", "y"]].to_numpy()
clean_adata.obsm["xy"] = spatial_coords
clean_adata.obsm["spatial"] = spatial_coords
print("Stored clean spatial coordinates in obsm['xy'] and obsm['spatial']")

recommended_qc_cols = [
    "nCount_Xenium",
    "nFeature_Xenium",
    "transcript_counts",
    "cell_area",
    "nucleus_area",
    "total_counts",
]
control_patterns = ["control", "codeword", "negative", "blank"]
present_qc_cols = [col for col in recommended_qc_cols if col in clean_adata.obs.columns]
missing_qc_cols = [col for col in recommended_qc_cols if col not in clean_adata.obs.columns]
control_qc_cols = [
    col for col in clean_adata.obs.columns
    if any(pattern in col.lower() for pattern in control_patterns)
]
print(f"QC columns present: {present_qc_cols}")
print(f"QC columns missing: {missing_qc_cols}")
print(f"Control/codeword-like columns present: {control_qc_cols}")

if clean_adata.shape != clean_shape_before:
    raise RuntimeError(
        "clean_adata shape changed during BANKSY metadata transfer"
    )
if not clean_adata.var_names.equals(clean_var_names_before):
    raise RuntimeError(
        "clean_adata variables changed during BANKSY metadata transfer"
    )

clean_clustered_h5ad = os.path.join(
    processed_path,
    f"adata_expression_clean_{dataset_name}_with_banksy_clusters_{res_str}.h5ad",
)
clean_adata.write_h5ad(clean_clustered_h5ad)
print(f"Wrote clean expression AnnData with BANKSY metadata: {clean_clustered_h5ad}")
print(clean_adata)




# In[ ]:


#banksy_dict[f"{nbr_weight_decay}"][0.2]['adata'] #lambda = 0.20




# In[ ]:


# ---------------------------------------------------------------------------- #
#              IDENTIFY TOP CLUSTER GENE MARKERS ON CLEAN ADATA                #
# ---------------------------------------------------------------------------- #

## Rank marker genes for each BANKSY clustering resolution using clean_adata.
## The group labels come from BANKSY, but expression values come from the clean
## log-normalized cell-by-gene matrix rather than the BANKSY spatial feature matrix.
marker_method = "wilcoxon"
n_marker_genes = 20

marker_table_path = os.path.join(
    output_path,
    "top_marker_tables",
)
os.makedirs(marker_table_path, exist_ok=True)

ranked_marker_keys = {}

for res in resolutions:
    res_str_single = str(res).replace(".", "p")
    groupby_key = make_banksy_label_col(nbr_weight_decay, pc_label, lambda_label, res)
    markers_key = f"{groupby_key}_markers_{marker_method}"

    if groupby_key not in clean_adata.obs.columns:
        raise KeyError(
            f"{groupby_key!r} was not found in clean_adata.obs. "
            f"Available obs columns: {list(clean_adata.obs.columns)}"
        )

    sc.tl.rank_genes_groups(
        clean_adata,
        groupby=groupby_key,
        method=marker_method,
        key_added=markers_key,
        use_raw=False,
    )

    marker_df = sc.get.rank_genes_groups_df(
        clean_adata,
        group=None,
        key=markers_key,
    )
    marker_csv = os.path.join(
        marker_table_path,
        f"rank_genes_groups_{dataset_name}_pc{pc_label}_nc{lambda_label}_r{res_str_single}_{marker_method}.csv",
    )
    marker_df.to_csv(marker_csv, index=False)

    ranked_marker_keys[res] = markers_key
    print(f"res {res:.2f}: wrote {markers_key}")
    print(f"res {res:.2f}: wrote marker table {marker_csv}")

## Preserve marker-ranking results in the clean AnnData object.
clean_adata.write_h5ad(clean_clustered_h5ad)
print(f"Updated clean expression AnnData with marker rankings: {clean_clustered_h5ad}")




# In[ ]:


# ---------------- Inspect marker keys saved on clean_adata ---------------- #
for res, markers_key in ranked_marker_keys.items():
    print(f"res {res:.2f}: {markers_key}")


# Old BANKSY-spatial marker inspection cell retired in this test notebook. Marker ranking now runs on clean_adata.
# 

# Old BANKSY-spatial marker inspection cell retired in this test notebook. Marker ranking now runs on clean_adata.
# 

# Old BANKSY-spatial marker inspection cell retired in this test notebook. Marker ranking now runs on clean_adata.
# 

# Old BANKSY-spatial marker inspection cell retired in this test notebook. Marker ranking now runs on clean_adata.
# 

# Old BANKSY-spatial marker inspection cell retired in this test notebook. Marker ranking now runs on clean_adata.
# 



# In[ ]:


# ---------------------------------------------------------------------------- #
#                              SET PATHS FOR PLOTS                             #
# ---------------------------------------------------------------------------- #

## Create a path for top marker plots generated from the clean expression object.
top_marker_plot_path = os.path.join(
    output_path,
    "top_marker_plot",
)
os.makedirs(top_marker_plot_path, exist_ok=True)
print(f"Top marker plots will be written to: {top_marker_plot_path}")

## Create a path for cluster counts plots.
cluster_count_plot_path = os.path.join(
    output_path,
    "cluster_count_plot",
)
os.makedirs(cluster_count_plot_path, exist_ok=True)
print(f"Cluster count plots will be written to: {cluster_count_plot_path}")

## Create a path for cluster heatmaps.
heatmap_path = os.path.join(
    output_path,
    "heatmap_path",
)
os.makedirs(heatmap_path, exist_ok=True)
print(f"Heatmaps will be written to: {heatmap_path}")




# In[ ]:


# ------------------ Plot marker genes for each clean cluster set ------------------ #

for res, markers_key in ranked_marker_keys.items():
    res_str_single = str(res).replace(".", "p")
    output_png = os.path.join(
        top_marker_plot_path,
        f"rank_genes_groups_{dataset_name}_pc{pc_label}_nc{lambda_label}_r{res_str_single}_{marker_method}.png",
    )

    sc.pl.rank_genes_groups(
        clean_adata,
        key=markers_key,
        n_genes=n_marker_genes,
        fontsize=15,
        show=False,
    )
    fig = plt.gcf()
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"res {res:.2f}: wrote marker plot {output_png}")




# In[ ]:


clean_adata.obs




# In[ ]:


# ------------------ Total counts across clean BANKSY clusters ------------------ #

sc.settings.figdir = cluster_count_plot_path
count_obs_key = resolve_obs_column(
    clean_adata,
    ["total_counts", "nCount_Xenium", "transcript_counts"],
    "total-count/QC",
)
print(f"Using {count_obs_key!r} for cluster count violin plots")

for res in resolutions:
    res_str_single = str(res).replace(".", "p")
    groupby_key = make_banksy_label_col(nbr_weight_decay, pc_label, lambda_label, res)

    if groupby_key not in clean_adata.obs.columns:
        raise KeyError(
            f"{groupby_key!r} was not found in clean_adata.obs. "
            f"Available obs columns: {list(clean_adata.obs.columns)}"
        )

    sc.pl.violin(
        clean_adata,
        count_obs_key,
        groupby=groupby_key,
        save=f"_{dataset_name}_pc{pc_label}_nc{lambda_label}_r{res_str_single}.png",
    )




# In[ ]:


# ##############################
# #       CLUSTER HEATMAP      #
# ##############################
# # Set the figure directory to your desired location
# sc.settings.figdir = heatmap_path

# ## Generate heat map for top markers per cluster
# for res in resolutions:
#     ad_res = spatial_adatas[res]
#     #markers_key = f"labels_scaled_gaussian_pc{pc_label}_nc{lambda_label}_r{res:.2f}_raw"

#     groupby_key = make_banksy_label_col(nbr_weight_decay, pc_label, lambda_label, res)
#     markers_key = f"{groupby_key}_raw"
#     print(markers_key)
#     sc.tl.dendrogram(
#         ad_res,
#         groupby=groupby_key,
#         key_added=markers_key
#     )

#     ## Determine vmin and vmax dynamically based on gene expression
#     expr_values = adata_spatial.X
#     vmax_dynamic = np.percentile(expr_values[expr_values > 0], 99)
#     vmin_dynamic = np.percentile(expr_values, 1)

#     sc.pl.rank_genes_groups_heatmap(
#         ad_res, 
#         key=markers_key,
#         cmap = "plasma",
#         vmin=vmin_dynamic,
#         vmax=vmax_dynamic,
#         show_gene_labels = True,
#         save = f"_{dataset_name}_pc{pc_label}_nc{lambda_label}_r{res_label}_raw_markers.png",
#         figsize=(20, 15)
#     )




# In[ ]:


# ########################################################
# #       GENERATE DENDROGRAM AND CORRELATION MATRIX     #
# ########################################################
# groupby_key = f'banksy_cluster_pc{pc_label}_nc{lambda_label}_r{res_label}_ann'

# sc.tl.dendrogram(adata_spatial, groupby=groupby_key)

# sc.pl.dendrogram(
#     adata_spatial,
#     groupby=groupby_key,
#     save=f"_{dataset_name}_pc{pc_label}_nc{lambda_label}_r{res_label}.png"
# )

# sc.pl.correlation_matrix(
#     adata_spatial,
#     groupby = groupby_key,
#     save=f"{dataset_name}_pc{pc_label}_nc{lambda_label}_r{res_label}_cluster_correlation_plot.png",
#     figsize=(5, 3.5)
# )




# In[ ]:


# #####################################
# #       EXPORT CLUSTER MARKERS      #
# #####################################
# ### Export the top 20 clusters in wide format to .csv

# from banksy_utils.annotation_utils import export_clusters_wide

# key= f"banksy_cluster_pc{pc_label}_nc{lambda_label}_r{res_label}_markers_raw"

# export_clusters_wide(
#     adata= adata_spatial,
#     key= key,
#     gene_type= "raw",
#     top_n= 20,
#     dataset_name= dataset_name,
#     file_path= processed_path
# )




# In[ ]:


# ### Export the top 20 clusters in long format with scores to .csv

# from banksy_utils.annotation_utils import export_cluster_markers

# export_cluster_markers(
#     adata = adata_spatial,
#     key = key,
#     top_n = 20,
#     dataset_name = dataset_name,
#     file_path = processed_path
# )




# In[ ]:


# ################################################
# #       ADD CELL TYPE INFO TO TOP MARKERS      #
# ################################################

# ## Read in the exported cluster markers generated with "export_clusters_wide()"
# ## and merge with the cell type annotation master document
# n_genes_label =  20 # Set this to the number of genes that were exorted as top markers for each cluster
# top_cluster_file = f"cluster_top_{n_genes_label}_genes_with_scores_{dataset_name}_{key}.csv"
# top_markers =pd.read_csv(os.path.join(processed_path, top_cluster_file), index_col=0)#index_col = 0 prevents insertion of an unnamed col
# top_markers = top_markers.sort_values(by='cluster', ascending=False)
# print(top_markers)

# ## Read in the master annotation file which contains all the curated cell type labels for each marker gene
# master_markers_file = pd.read_csv(raw_path, "xenium_gene_list_annotation_master_v1_March_2026_Tier1_v2.csv")

# ## Create a data frame of the genes and their corresponding primary and secondary cell type annotations 
# cell_annotations = master_markers_file[["Gene", "primary_annotation", "secondary_annotation"]]
# #cell_annotations = cell_annotations.sort_index()


# merged_markers = top_markers.reset_index().merge(
#     cell_annotations,
#     left_on="gene",
#     right_on="Gene",
#     how='left'
# )
# merged_markers['cluster'] = merged_markers['cluster'].astype(int)
# merged_markers = merged_markers.sort_values('cluster', ascending=True).reset_index(drop=True)
# merged_markers = merged_markers.drop("Gene", axis=1)
# print(merged_markers)

# ## Export the cell type annotated top markers to .csv to use in manual annotation and ChatGPT annotation
# merged_markers.to_csv(os.path.join(processed_path, f"cell_type_cluster_top_{n_genes_label}_genes_{dataset_name}_{key}.csv"))


# ### Non-spatial clustering results




# In[ ]:


# # banksy_dict['scaled_gaussian'][0.2]['adata'] #lambda = 0.2
# adata_nonspatial.obs['banksy_cluster_nonspatial']= adata_nonspatial.obs[f"labels_nonspatial_pc{pc_label}_nc0.00_r{res_label}"].astype(str)


