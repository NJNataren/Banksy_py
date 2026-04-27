#!/usr/bin/env python
# coding: utf-8

# In[38]:


# # Xenium spatial clustering with BANKSY
# ### **Author:** Nathalie Nataren
# ### **Date:** 17/04/2026
# 
# **Description:** The purpose of this analysis is to perform clustering of Xenium spatial data as part of the metastatic melanoma ICI therapy response study (VBCT lab).
# 


# In[1]:


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


# In[ ]:


# ##################################
# #   ** LOCAL TESTING BLOCK **    #
# ##################################

# ## Set the dataset_name and related settings to use during this analysis
# dataset_name = "CK_skin_res" # sample name
# #dataset_name = "BE_brain_non_res" # sample name
# pc_label = "20" # Label for the number of principal components used for the purpose of filenames
# pc_dims = [20] # The number of principal components stored a list for analyses
# lambda_label = "0.20" # File name label for Lambda setting, see comment below. 
# lambda_list = [float(lambda_label)] # Lambda setting to tune BANKSY clustering, lambda = 0 is non-spatial, 0.2 is for cell typing, 0.8 if for domain segmentation 
# res_label = ["0.50", "0.60", "0.70"] # BANKSY clustering resolution label for resolution chosen to produce plots

# resolutions = [float(res) for res in res_label] # BANSY can take a list of resolutions and perform clustering at each which is saved in the BANKSY dictionary
# #resolutions = [float(res_label)] # BANSY can take a list of resolutions and perform clustering at each which is saved in the BANKSY dictionary
# nbr_weight_decay = "scaled_gaussian" # This parameter dictates how much neighbouring cells impact to the neighbourhood expression calculations. Using scaled gaussian, the 
# # close neigbours contribute more and this decays as you move out to cells further away in the neighbourhood window. It is scaled for local cell density so that weighting doesn't change
# # across regions if cells are pack more closely or loosely in different regions
# coord_keys = ('x', 'y', 'xy') # Keys to specify coordinate indexes in the anndata Object
# max_workers=8 # maximum CPUs for Leiden clustering


# In[ ]:


########################
#   PARSE ARGUMENTS    #
########################
# This block of code feeds arguments to this python script from a config files found in /config

## Import argparse and json packages to read in variables from the per sample .json config files 
import argparse
import json

parser = argparse.ArgumentParser(prog="used to parse arguments form xenium_clustering.py to run on slurm") # Initialise the parser
parser.add_argument("--config", type=str, help="Provide a JSON config file for each Xenium sample", required=True) # This defines the flag and tells the script to look for a JSON config
# in the form of a string config file path
args = parser.parse_args() # Looks at what is passed throught the terminal (in the slurm script in this case) after --config and stores it 

with open(args.config) as f: # Opens the file path provided by the user
    cfg = json.load(f) # Converts the json config into a python dictionary called "cfg"

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
max_workers = int(os.environ.get("SLURM_CPUS_PER_TASK", 4)) # Parameter for the run_Leiden_partition_parallel() clustering function 


# In[19]:


########################
#       SET PATHS      #
########################
## Set file paths and read in xenium data

## Create a base path
base_dir = "data/xenium"

## Create a path to the raw data e.g., unprocessed anndata files, if it does not already exist
raw_path = os.path.join(base_dir, "raw_data")

if not os.path.isdir(raw_path):
    os.makedirs(raw_path)
    print(f"Directory '{raw_path} successfully.")
    
else:
    print(f"Directory '{raw_path} exists.")

## Create path for processed data e.g., the pre-clustered but unfiltered anndata files, if it does not already exist
processed_path = os.path.join(base_dir, "processed", f"{dataset_name}")

if not os.path.isdir(processed_path):
    os.makedirs(processed_path)
    print(f"Directory '{processed_path}' created successfully.")
    
else:
    print(f"Directory '{processed_path}' already exists.")

## Create a path for output data, if it does not already exist
output_path = os.path.join(base_dir, "output", f"{dataset_name}")

if not os.path.isdir(output_path):
    os.makedirs(output_path)
    print(f"Directory '{output_path}' created successfully.")
    
else:
    print(f"Directory '{output_path}' already exists.")

## Create a path for QC results, if it does not already exist
qc_path = os.path.join(base_dir, "output", "QC_testing", f"{dataset_name}")

if not os.path.isdir(qc_path):
    os.makedirs(qc_path)
    print(f"Directory '{qc_path}' created successfully.")
else:
    print(f"Directory '{qc_path}' already exists.")


# In[20]:


## Function to log sub task start time
def log_time(step):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {step}")


# In[21]:


###########################
#       LOAD ANNDATA      #
###########################

## Read in the raw AnnData file
adata = ad.read_h5ad(os.path.join(raw_path, f"{dataset_name}_raw.h5ad"))
log_time(f"Loading in data for {dataset_name}")

## Create 'xy' spatial coordinates from adata.obs
adata.obsm['xy'] = np.vstack([adata.obs['x'], adata.obs['y']]).T


# In[22]:


res_label


# In[23]:


#####################################
#       FILTER ZERO COUNT CELLS     #
#####################################

## Filter out cells with zero counts
adata = adata[adata.obs['nCount_Xenium'] > 0].copy()


# In[24]:


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


# In[25]:


#######################################
#       LOAD DATA AND COORDINATES     #
#######################################

from banksy_utils.load_data import load_adata, display_adata

## To either load data from .h5ad directly or convert raw data to .h5ad format
load_adata_directly = True

## Keys to specify coordinate indexes in the anndata Object
coord_keys = coord_keys

raw_y, raw_x, adata = load_adata(filepath=processed_path, adata_filename=float_adata, load_adata_directly=True, coord_keys=coord_keys)


# In[26]:


# ## Normalize and log transform data
# 
# Normalize data matrix using `normalize_total`.

## Save the raw counts
adata.layers["counts"] = adata.X.copy()

## Inspect the counts
print(adata.layers["counts"][:5,:5])


# In[27]:


from banksy_utils.filter_utils import normalize_total, filter_hvg, print_max_min

## Normalizes the anndata dataset
normalize_total(adata)
print(adata.X)


# In[28]:


## Perform log-transformation and save the log-normalised in adata.raw
sc.pp.log1p(adata)
print(adata.X)

adata.raw = adata.copy()


# In[29]:


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



# In[30]:


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


# In[31]:


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


# In[32]:


# ### Append Non-spatial results to the `banksy_dict` for comparsion

from banksy.main import concatenate_all

banksy_dict["nonspatial"] = {
    # Here we simply append the nonspatial matrix (adata.X) to obtain the nonspatial clustering results
    0.0: {"adata": concatenate_all([adata.X], 0, adata=adata), }
}

print(banksy_dict['nonspatial'][0.0]['adata'])


# In[33]:


## Perform UMAP embedding
from banksy_utils.umap_pca import pca_umap
log_time(f"Start PCA and UMAP embedding for {dataset_name}.")
pca_umap(banksy_dict,
         pca_dims = pc_dims,
         add_umap = True,
         #plt_remaining_var = True,
         )
log_time(f"Finish PCA and UMAP embedding for {dataset_name}.")


# In[34]:


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


# In[35]:


# ## Dynamically extract the number of principal components from the results_df
# pc_dims = results_df[results_df['decay'] == nbr_weight_decay]
# pc_dims = pc_label['num_pcs'].iloc[0]


# In[36]:


# ## Plot results
# 
# ### Visualize the clustering results from BANKSY, including the clusters from the Umap embbedings

from banksy.plot_banksy import plot_results

c_map =  'tab20' # specify color map
weights_graph =  banksy_dict[f"{nbr_weight_decay}"]['weights'][0]

plot_results(
    results_df,
    weights_graph,
    c_map,
    match_labels = False,
    coord_keys = coord_keys,
    max_num_labels  =  max_num_labels, 
    #save_path = os.path.join(file_path, 'tmp_png'),
    save_fig = True, # save the spatial map of all clusters
    save_seperate_fig = True, # save the figure of all clusters plotted seperately
    dataset_name = f"{dataset_name}",
    save_fullfig=True,
    save_path = output_path
)

print(results_df)


# In[64]:


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


# In[38]:


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


# In[39]:


print("results_df.index:")
for idx in results_df.index.tolist():
    print(f"  {repr(idx)}")
print(f"\nLooking for: 'scaled_gaussian_pc{pc_dims[0]}_nc{lambda_list[0]:0.2f}_r0.50'")
print(f"pc_dims = {pc_dims}")
print(f"lambda_list = {lambda_list}")
print(f"resolutions = {resolutions}")


# In[40]:


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


# In[41]:


res_label


# 

# In[42]:


## Create a flat dictionary of the spatial and non-spatial 
# spatial
spatial_adatas = {}
for res, adatas in adata_dict.items():
    spatial_adatas[res]  = adatas["spatial"]


# In[43]:


# Save individaul anndata objects at each resolution
for res in resolutions:
    res_str  = str(res).replace(".", "p")
    spatial_adatas[res].write_h5ad(os.path.join(processed_path, f"adata_spatial_{dataset_name}_{res_str}.h5ad"))


# In[44]:


## Export dictionary
res_str = "_".join(res_label)
## Save the banksy dict as a pickle file to load it in later and avoid having to initialize repeatedly to save time
## Use gzip to save it with compression
import gzip
import pickle
#for res in resolutions:
with gzip.open(os.path.join(processed_path, f"{dataset_name}_pc{pc_label}_nc{lambda_label}_r{res_str}_banksy_dict.pkl.gz"), "wb") as f:
    pickle.dump(banksy_dict, f)

## Export the results_df data frame

results_df.to_csv(os.path.join(processed_path,f"results_df_{dataset_name}_pc{pc_label}_nc{lambda_label}.csv"))


# In[45]:


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
    df = spatial_adatas[res].obs.reset_index()
    #res_ = float(res)
    #df = df[['index', f'labels_scaled_gaussian_pc{pc_label}_nc{lambda_label}_r{res_label}0']]
    df = df[['index', f'labels_scaled_gaussian_pc{pc_label}_nc{lambda_label}_r{res:.2f}']]
    merged = pd.merge(merged, df, on = 'index', how = 'left')
print(merged)

merged.to_csv(os.path.join(processed_path, f"{dataset_name}_cell_cluster_id_across_clustering_res_{res_str}.csv"))


# In[46]:


resolutions


# In[47]:


res_label


# In[48]:


# #####################################################
# #       STORECLUSTER LABELS AS A STRING IN .OBS     #
# #####################################################

# banksy_dict[f"{nbr_weight_decay}"][0.2]['adata'] #lambda = 0.20
# adata_spatial.obs[f"banksy_cluster_pc{pc_label}_nc{lambda_label}_r{res_label}"]= adata_spatial.obs[f"labels_scaled_gaussian_pc{pc_label}_nc{lambda_label}_r{res_label}"].astype(str)


# In[49]:


# print(adata.shape)          # (n_cells, 304)
# print(adata_spatial.shape)  # (n_cells, 912)

# # Check if the genes overlap
# common_genes = adata.var.index.intersection(adata_spatial.var.index)
# print(f"Common genes: {len(common_genes)}")

# # Check what the extra genes are
# print(adata_spatial.var.index[:20])


# In[50]:


# ## Store the cluster labels as a string in .obs
# banksy_dict[f"{nbr_weight_decay}"][0.2]['adata'] #lambda = 0.20
# adata_spatial.obs[f"banksy_cluster_pc{pc_label}_nc{lambda_label}_r{res_label}"]= adata_spatial.obs[f"labels_scaled_gaussian_pc{pc_label}_nc{lambda_label}_r{res_label}"].astype(str)


# In[51]:


# ## Pull out genes as data frame
# genes = adata.var_names
# gene_df = pd.DataFrame(adata.X.toarray(), columns=genes, index=adata.obs_names)
# #print(gene_df)

# gene_count_across_cells = gene_df.sum(axis=0)
# gene_count_across_cells.to_frame()
# gene_count_across_cells.columns = ["gene", "counts"]

# print(gene_count_across_cells)
# gene_count_across_cells.to_csv(f"{processed_path}{dataset_name}_gene_count_across_cells_not_SCTransformed.csv")
# type(gene_count_across_cells)


# In[52]:


# print(adata_spatial.__repr__())


# # Top expressed genes per cluster

# ### Spatial clustering results


# In[53]:


# ##############################################
# #       TOP EXPRESSED GENES PER CLUSTER      #
# ##############################################
# ### Spatial clustering results

# import scanpy as sc

# ## Rank genes and identify cluster gene markers
# sc.tl.rank_genes_groups(
#     adata_spatial,
#     groupby=f"banksy_cluster_pc{pc_label}_nc{lambda_label}_r{res_label}",
#     method='wilcoxon',
#     key_added=f"banksy_cluster_pc{pc_label}_nc{lambda_label}_r{res_label}_markers"
# )


# In[54]:


# ## Using the filter_ranked_genes_by_type to create a key for 
# from banksy_utils.annotation_utils import filter_ranked_genes_by_type

# key_filtered = filter_ranked_genes_by_type(
#     adata_spatial,
#     key=f"banksy_cluster_pc{pc_label}_nc{lambda_label}_r{res_label}_markers",
#     top_n=100,
#     gene_type = "raw",
#     new_key_suffix="_raw",
# )


# In[55]:


# # Plot the gene markers for each cluster
# sc.pl.rank_genes_groups(
#     adata_spatial, 
#     key=key_filtered,
#     method='wilcoxon',
#     n_genes=20, 
#     fontsize=15,
#     save= f"_{dataset_name}_pc{pc_label}_nc{lambda_label}_r{res_label}_spatial_top_gene_raw.png")


# In[56]:


# ## Total counts across unlabeled clustered
# sc.pl.violin(
#     adata_spatial, 
#     "total_counts", 
#     groupby=f"banksy_cluster_pc{pc_label}_nc{lambda_label}_r{res_label}", 
#     save= f"_{dataset_name}_pc{pc_label}_nc{lambda_label}_r{res_label}.png"
#     )


# In[57]:


# ##############################
# #       CLUSTER HEATMAP      #
# ##############################
# # Set the figure directory to your desired location
# sc.settings.figdir = output_path

# ## Generate heat map for top markers per cluster
# sc.tl.dendrogram(
#     adata_spatial,
#     groupby=f'banksy_cluster_pc{pc_label}_nc{lambda_label}_r{res_label}_raw',
#     key_added=f'banksy_cluster_pc{pc_label}_nc{lambda_label}_r{res_label}_raw'
# )

# ## Determine vmin and vmax dynamically based on gene expression
# expr_values = adata_spatial.X
# vmax_dynamic = np.percentile(expr_values[expr_values > 0], 99)
# vmin_dynamic = np.percentile(expr_values, 1)

# sc.pl.rank_genes_groups_heatmap(
#     adata_spatial, 
#     key=f'banksy_cluster_pc{pc_label}_nc{lambda_label}_r{res_label}_markers_raw',
#     #key= 'raw_test',
#     n_genes=5,
#     cmap = "plasma",
#     vmin=vmin_dynamic,
#     vmax=vmax_dynamic,
#     show_gene_labels = True,
#     save = f"_{dataset_name}_pc{pc_label}_nc{lambda_label}_r{res_label}_raw_markers.png",
#     figsize=(20, 15)
# )


# In[58]:


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


# In[59]:


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


# In[60]:


# ### Export the top 20 clusters in long format with scores to .csv

# from banksy_utils.annotation_utils import export_cluster_markers

# export_cluster_markers(
#     adata = adata_spatial,
#     key = key,
#     top_n = 20,
#     dataset_name = dataset_name,
#     file_path = processed_path
# )


# In[61]:


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


# In[62]:


# # banksy_dict['scaled_gaussian'][0.2]['adata'] #lambda = 0.2
# adata_nonspatial.obs['banksy_cluster_nonspatial']= adata_nonspatial.obs[f"labels_nonspatial_pc{pc_label}_nc0.00_r{res_label}"].astype(str)


# In[63]:


# # Save AnnData objects with clustering results

# adata_spatial.write_h5ad(
#     filename = f"data/xenium/processed/{dataset_name}_clustered_spatial_pc{pc_label}_nc{lambda_label}_r{res_label}.h5ad",
#     compression="gzip"
#     )


# # In[ ]:


# ## Export dictionary

# ## Save the banksy dict as a pickle file to load it in later and avoid having to initialize repeatedly to save time
# ## Use gzip to save it with compression
# import gzip
# import pickle
# with gzip.open(f"data/xenium/processed/{dataset_name}_pc{pc_label}_nc{lambda_label}_r{res_label}_banksy_dict.pkl.gz", "rb") as f:
#     pickle.dump(banksy_dict, f)

# ## Export the results_df data frame

# results_df.to_csv(f"data/xenium/processed/results_df_{dataset_name}_pc{pc_label}_nc{lambda_label}_r{res_label}.csv")

