#!/usr/bin/env python
# coding: utf-8

#  # Squipy Xenium spatial melanoma neighbourhood and gene of interest analysis

# ## Import packages

# In[75]:


import anndata as ad
import os
import numpy as np
import pandas as pd

import spatialdata as sd
from spatialdata_io import xenium

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import squidpy as sq

from pathlib import Path
import shutil

import squidpy as sq
from IPython.display import display
import warnings
warnings.filterwarnings("ignore") 

import scipy.sparse as sparse
from scipy.io import mmread

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib import rc_context

import random
# Note that BANKSY itself is deterministic, here the seeds affect the umap clusters and leiden partition
seed = 1234
np.random.seed(seed)
random.seed(seed)

#import dask
#dask.config.set({"dataframe.query-planning": True})


# ## Set file paths and read in xenium data

# In[ ]:


# #############
# #   TEST    #
# #############

# ## Set the dataset_name and related settings to use during this analysis
# dataset_name = "GR_lung_non_res_roi"

# pca_label = "35"
# pca_dims = int(pca_label)
# pca_dims = [pca_dims]

# lambda_label = "0.20"
# lambda_list = [float(lambda_label)]

# res_label = "0.50"
# resolutions = [float(res_label)]

# nbr_weight_decay = "scaled_gaussian"

# # Keys to specify coordinate indexes in the anndata Object
# coord_keys = ('x', 'y', 'xy')


# ### Create directories for processed and output data and set the dataset_name

# In[77]:


# #############
# #   TEST    #
# #############

# #  Create a path for raw data
# raw_data = f"data/xenium/raw_data/"

# # Create path for processed data
# processed_path = f"data/xenium/processed/anndata_conversion_test/{dataset_name}/"

# if not os.path.isdir(processed_path):
#     os.makedirs(processed_path)
#     print(f"Directory '{processed_path}' created successfully.")
# else:
#     print(f"Directory '{processed_path}' already exists.")


# # Create path for output data
# output_path = f"data/xenium/output/anndata_conversion_test/{dataset_name}/"

# if not os.path.isdir(output_path):
#     os.makedirs(output_path)
#     print(f"Directory '{output_path}' created successfully.")
# else:
#     print(f"Directory '{output_path}' already exists.")

# # Create a path for QC results
# qc_path = f"data/xenium/output/anndata_conversion_test/{dataset_name}/QC"

# if not os.path.isdir(qc_path):
#     os.makedirs(qc_path)
#     print(f"Directory '{qc_path}' created successfully.")
# else:
#     print(f"Directory '{qc_path}' already exists.")


# In[ ]:


# ## Set the dataset_name and related settings to use during this analysis
import json
import argparse

parser = argparse.ArgumentParser(prog="used to parse arguments form 00_QC_xenium_spatial.py to run on slurm")
parser.add_argument("--config", type=str, help="Provide a JSON config file for each Xenium sample", required=True)
args = parser.parse_args()

with open(args.config) as f:
    cfg = json.load(f)

# Derive the values from the config JSON

dataset_name = cfg["dataset_name"]
pca_label = cfg["pca_label"]
pca_dims = [int(pca_label)]
lambda_label = cfg["lambda_label"]
lambda_list = [float(lambda_label)]
res_label = cfg["res_label"]
resolutions = [float(res_label)]
nbr_weight_decay = cfg["nbr_weight_decay"]

# Keys to specify coordinate indexes in the anndata Object
coord_keys = tuple(cfg["coord_keys"])


# In[78]:


## Create a base path
base_dir = "data/xenium"

## Create a path to the raw data e.g., unprocessed anndata files
raw_path = os.path.join(base_dir, "raw_data")

if not os.path.isdir(raw_path):
    os.makedirs(raw_path)
    print(f"Directory '{raw_path} successfully.")
    
else:
    print(f"Directory '{raw_path} exists.")

## Create path for processed data e.g., the pre-clustered but unfiltered anndata files
processed_path = os.path.join(base_dir,f"processed/{dataset_name}")

if not os.path.isdir(processed_path):
    os.makedirs(processed_path)
    print(f"Directory '{processed_path}' created successfully.")
    
else:
    print(f"Directory '{processed_path}' already exists.")

## Create a path for output data
output_path = os.path.join(base_dir,f"output/{dataset_name}")

if not os.path.isdir(output_path):
    os.makedirs(output_path)
    print(f"Directory '{output_path}' created successfully.")
    
else:
    print(f"Directory '{output_path}' already exists.")

## Create a path for QC results
qc_path = f"data/xenium/output/QC_testing/{dataset_name}"

if not os.path.isdir(qc_path):
    os.makedirs(qc_path)
    print(f"Directory '{qc_path}' created successfully.")
else:
    print(f"Directory '{qc_path}' already exists.")


# ## Load Xenium data
# Allocate path to folder and the file name of the designated AnnaData Object (in `.h5ad` format) <br>

# In[38]:


## Read in the raw AnnData file
adata = ad.read_h5ad(f"data/xenium/raw_data/{dataset_name}_raw.h5ad")
#output_path = f"data/xenium/processed/{dataset_name}_xy.h5ad"
adata.obs_keys()

# Create 'xy' spatial coordinates from adata.obs
import numpy as np
adata.obsm['xy'] = np.vstack([adata.obs['x'], adata.obs['y']]).T


# In[79]:


# Filter out cells with zero counts
adata = adata[adata.obs['nCount_Xenium'] > 0].copy()


# In[80]:


# Function to check the float size of the adata.obsm data
from banksy_utils.object_downcasting_utils import check_float, downcast_float
check_float(adata)

# Function to downcast adata.obsm to 32 bit and check that the resulting array is an N-dimensional array
downcast_float(adata, "float32")


# In[81]:


# Check the float size of the adata.obsm data after downcasting to make sure the float type is correct
check_float(adata)


# In[86]:


# Write the downcast file to AnnData format 
import anndata as ad
adata.write_h5ad(
    #filename = f"data/xenium/processed/{dataset_name}_float_32.h5ad",
    filename = f"{processed_path}/{dataset_name}_raw_sopa_float_32.h5ad",
    compression="gzip"
    )

# Define float adata file
float_adata = f"{dataset_name}_raw_sopa_float_32.h5ad"


# In[87]:


f"{processed_path}{dataset_name}_raw_sopa_float_32.h5ad"


# In[88]:


from banksy_utils.load_data import load_adata, display_adata

# To either load data from .h5ad directly or convert raw data to .h5ad format
load_adata_directly = True

# Keys to specify coordinate indexes in the anndata Object
coord_keys = ('x', 'y', 'xy')

raw_y, raw_x, adata = load_adata(filepath=processed_path, adata_filename=float_adata, load_adata_directly=True, coord_keys=('x','y','xy'))


# In[89]:


adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("MT-")

# Calulates QC metrics and put them in place to the adata object
sc.pp.calculate_qc_metrics(adata, 
                           qc_vars=["mt"], 
                           log1p=True, 
                           inplace=True,
                           percent_top=[10,20]
                           )


# In[90]:


from banksy_utils.plot_utils import plot_qc_hist, plot_cell_positions

# bin options for fomratting histograms
# Here, we set 'auto' for 1st figure, 80 bins for 2nd figure. and so on
hist_bin_options = ['auto', 80, 80, 100]

plot_qc_hist(adata, 
         total_counts_cutoff = 200, # for visualization
         n_genes_high_cutoff = 1000, # for visualization
         n_genes_low_cutoff = 0, # for visualization
         bin_options = hist_bin_options)

plt.savefig(os.path.join(qc_path, f"banksy_counts_and_genes_plot_{dataset_name}.png"), dpi =300, bbox_inches='tight')
plt.show()


# ## Filter 1 - minimum number of transcripts filter

# ### Plot knee plot (log10 transcripts per cell)

# In[91]:


## from https://pachterlab.github.io/kallistobustools/tutorials/kb_getting_started/python/kb_intro_2_python/

threshold = 10
knee = np.sort((np.array(adata.X.sum(axis=1))).flatten())[::-1]
#knee = np.sort(adata.obs['nCount_Xenium'].values)[::-1]
#knee = np.sort(np.asarray(adata.X.sum(axis=1)).flatten())[::-1]

cell_set = np.arange(len(knee))
#num_cells = cell_set[knee > threshold] [::-1][0]
num_cells = (knee > threshold).sum()

fig, ax = plt.subplots(figsize=(10, 7))

#ax.loglog(cell_set, knee, linewidth=5, color="g")
ax.semilogy(knee, linewidth =5, color="g")
#ax.loglog(cell_set, knee, linewidth=5, color="g")
ax.axhline(y=threshold, linewidth=3, color="b")
ax.axvline(x=num_cells, linewidth=3, color="r")

ax.set_xlabel("Cells (ranked)")
ax.set_ylabel("Total counts per cell")

ax.set_yticks([10, 20, 40, 60, 80, 100, 120, 140, 160, 200, 250, 300, 400, 500, 600, 800, 1000, 2000])
#ax.set_yticks([10, 20, 40, 60, 80, 100, 120, 140, 160,200])
ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())

plt.grid(True, which="both")
plt.title(f"Knee Plot | Threshold = {threshold} counts | {num_cells:,} cells pass")
plt.savefig(os.path.join(qc_path, f"transcript_knee_plot_{dataset_name}.png"), dpi =300, bbox_inches='tight')
plt.show()

print(f"Cells passing threshold of {threshold} counts: {num_cells:,} / {len(knee):,} ({num_cells/len(knee)*100:.1f}%)")


# In[92]:


# ### Read in the raw matrix here!!
# raw_matrix = pd.read_csv(os.path.join(raw_data, "GR_lung_non_res_cell_feature_matrix.csv"))
# raw_matrix.rename(columns={raw_matrix.columns[0]: "cell_ids"}, inplace=True)

# # Defince cells IDs to check
# cell_ids = [
#     'akljhmnp-1',
#     'akljiabi-1', 
#     'aklkomhk-1',
#     'akllenhp-1',
#     'aklobaoh-1',
#     'aklocplk-1'
# ]

# # Subset the raw matric by the cells above to compare to the adata.X data for these same cells

# raw_matrix_subset = raw_matrix[raw_matrix['cell_ids'].isin(cell_ids)]

# raw_matrix_subset.to_csv(os.path.join(qc_path, 'raw_matrix_subset.csv'))



# In[93]:


# # Define your cell IDs
# cell_ids = [
#     'akljhmnp-1',
#     'akljiabi-1', 
#     'aklkomhk-1',
#     'akllenhp-1',
#     'aklobaoh-1',
#     'aklocplk-1'
# ]

# # Subset adata.X for these cells
# adata_subset = adata[cell_ids, :]

# # Get the X matrix for these cells
# X_subset = adata_subset.X.toarray()  # convert sparse to dense

# # Sum counts across cells (per gene)
# gene_sums = X_subset.sum(axis=0)  # sum down rows = per gene total
# print(pd.DataFrame(gene_sums.reshape(1, -1), columns=adata.var_names))

# # Sum counts per cell (across genes)
# cell_sums = X_subset.sum(axis=1)  # sum across columns = per cell total
# print(pd.DataFrame(cell_sums, index=cell_ids, columns=['total_counts']))

# # Subset obsm for these cells
# obsm_subset = adata_subset.obsm['xy']  # or whatever key you need
# print(obsm_subset)

# # Output the full X matrix as a readable DataFrame
# X_df = pd.DataFrame(
#     X_subset,
#     index=cell_ids,
#     columns=adata.var_names
# )
# print(X_df)

# # Optionally save to CSV
# X_df.to_csv(os.path.join(qc_path, 'subset_X.csv'))


# #### Apply the minimum transcripts per cell mask 

# In[94]:


## Apply a mask to retrieve cells that don't pass the filter
mask = adata.obs['nCount_Xenium'] <= threshold

# Turn the masked cells into a data frame
mask.to_frame()
mask_cells = mask.reset_index()
mask_cells = mask_cells.rename(columns={"index":"cell_id", "nCount_Xenium" : "threshold_passed"})

# Filter for the poor cells and keep these ids for future plotting
poor_cells = mask_cells[mask_cells['threshold_passed']]

# Put results for the cells passing the threshold into the object directly
adata.obs['threshold_passed'] = adata.obs['nCount_Xenium'] > threshold #True = passed


# #### Plot location of poor cells on tissue

# In[95]:


adata.uns["spatial"] = {dataset_name: {}}

with rc_context({"figure.figsize": (12,8)}):
    sq.pl.spatial_scatter(
        adata,
        library_id="dataset_name",
        spatial_key="xy",
        color="threshold_passed",
        shape=None,
        size=2,
        img=False
    )


plt.savefig(
    os.path.join(qc_path, f"tissue_spatial_scatter_threshold_passed_cells_qc_{dataset_name}.png"),
    dpi=300,
    bbox_inches='tight'
    )


# In[96]:


adata.obs["threshold_passed_cat"] = adata.obs["threshold_passed"].map({True: "Pass", False: "Fail"}).astype("category")

with rc_context({"figure.figsize": (12, 8)}):
    sq.pl.spatial_scatter(
        adata,
        library_id="dataset_name",
        spatial_key="xy",
        color="threshold_passed_cat",
        shape=None,
        size=2,
        img=False
    )
    plt.legend(fontsize=20)


# #### Plot of total transripts per cell across sample (rainbow in viridis)

# In[97]:


with rc_context({"figure.figsize": (12, 8)}):
    sq.pl.spatial_scatter(
        adata,
        library_id="dataset_name",
        spatial_key="xy",
        color="nCount_Xenium",
        shape=None,
        size=2,
        img=False,
        vmax = adata.obs['nCount_Xenium'].max(),
        #vmax= 300,
        palette="gist_stern"
    )

#
plt.savefig(
    os.path.join(qc_path, f"tissue_spatial_scatter_threshold_passed_cells_qc_{dataset_name}.png"),
    dpi=300,
    bbox_inches='tight'
    )


# In[98]:


# obs_data= adata.obs
# obs_data.to_csv(f"{processed_path}/obs_data_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}.csv")

# var_data= adata.var
# var_data.to_csv(f"{processed_path}/var_data_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}.csv")

# uns_data= adata.uns
# uns_data.to_csv(f"{processed_path}/uns_data_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}.csv")

# obsm_data= adata.obsm
# obsm_data.to_csv(f"{processed_path}/obsm_data_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}.csv")



# ### Check a slice of cells to inspect counts

# In[99]:


# # Define cells to look at
# cell_ids = [
#     'akljhmnp-1',
#     'akljiabi-1', 
#     'aklkomhk-1',
#     'akllenhp-1',
#     'aklobaoh-1',
#     'aklocplk-1'
# ]

# # Subset obs for these cells
# obs_subset = adata.obs.loc[cell_ids]
# print(obs_subset)

# # Check specific columns
# print(obs_subset[['nCount_Xenium', 'nCount_SCT', 'total_counts']])


# obs_subset.to_csv(f"{qc_path}/obs_subset_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}.csv")


# In[100]:


# # Define your cell IDs
# cell_ids = [
#     'akljhmnp-1',
#     'akljiabi-1', 
#     'aklkomhk-1',
#     'akllenhp-1',
#     'aklobaoh-1',
#     'aklocplk-1'
# ]

# # Subset adata.X for these cells
# adata_subset = adata[cell_ids, :]

# # Get the X matrix for these cells
# X_subset = adata_subset.X.toarray()  # convert sparse to dense

# # Sum counts across cells (per gene)
# gene_sums = X_subset.sum(axis=0)  # sum down rows = per gene total
# print(pd.DataFrame(gene_sums.reshape(1, -1), columns=adata.var_names))

# # Sum counts per cell (across genes)
# cell_sums = X_subset.sum(axis=1)  # sum across columns = per cell total
# print(pd.DataFrame(cell_sums, index=cell_ids, columns=['total_counts']))

# # Subset obsm for these cells
# obsm_subset = adata_subset.obsm['xy']  # or whatever key you need
# print(obsm_subset)


# In[101]:


# # Define your cell IDs
# cell_ids = [
#     'akljhmnp-1',
#     'akljiabi-1', 
#     'aklkomhk-1',
#     'akllenhp-1',
#     'aklobaoh-1',
#     'aklocplk-1'
# ]

# # Subset adata.X for these cells
# adata_subset = adata[cell_ids, :]

# # Get the X matrix for these cells
# X_subset = adata_subset.X.toarray()

# # Sum counts across cells (per gene) and save
# gene_sums_df = pd.DataFrame(gene_sums.reshape(1, -1), columns=adata.var_names)
# gene_sums_df.to_csv(os.path.join(qc_path, "gene_sums_subset.csv"), index=False)

# # Sum counts per cell (across genes) and save
# cell_sums_df = pd.DataFrame(cell_sums, index=cell_ids, columns=['total_counts'])
# cell_sums_df.to_csv(os.path.join(qc_path, "cell_sums_subset.csv"))

# # Full X matrix for these cells and save
# X_subset_df = pd.DataFrame(X_subset, index=cell_ids, columns=adata.var_names)
# X_subset_df.to_csv(os.path.join(qc_path, "X_subset.csv"))

# # Subset obsm and save
# obsm_subset_df = pd.DataFrame(obsm_subset, index=cell_ids, columns=['x', 'y'])
# obsm_subset_df.to_csv(os.path.join(qc_path, "obsm_subset.csv"))

# print("All files saved to:", qc_path)


# In[102]:


## Convert the the whole adata.X into a matrix and export to .csv
pd.DataFrame(
    adata.X.toarray(),  # convert sparse matrix to dense
    index=adata.obs_names, 
    columns=adata.var_names
).to_csv(os.path.join(qc_path, f"adata_X_inspection_{dataset_name}.csv"))


# In[103]:


obs= adata.obs
obs.to_csv(os.path.join(qc_path, f"obs_{dataset_name}.csv"))


# ### Inspect gene counts for the cells with the highest counts

# In[104]:


# cell_ids = [
# "blmnefce-1",
# "ienlfejo-1",
# "mcmnjhco-1",
# "hpdkdmgb-1",
# "hbjkojbd-1",
# "aomhfiaa-1",
# "hdgdbfcn-1",
# "goojoceo-1"
# ]

# # Subset adata.X for these cells
# adata_subset = adata[cell_ids, :]

# adata_subset_X = pd.DataFrame(
#     adata_subset.X.toarray(),
#     index=cell_ids,
#     columns=adata.var_names
# )

# adata_subset_X.to_csv(os.path.join(qc_path, f"adata_subset_X_inspection_{dataset_name}.csv"))

# # # Get the X matrix for these cells
# # X_subset = adata_subset.X.toarray()  # convert sparse to dense

# # Sum counts across cells (per gene)
# gene_sums = X_subset.sum(axis=0)  # sum down rows = per gene total
# print(pd.DataFrame(gene_sums.reshape(1, -1), columns=adata.var_names))

# # # # Sum counts per cell (across genes)
# # # cell_sums = X_subset.sum(axis=1)  # sum across columns = per cell total
# # # print(pd.DataFrame(cell_sums, index=cell_ids, columns=['total_counts']))

# # # Subset obsm for these cells
# # obsm_subset = adata_subset.obsm['xy']  # or whatever key you need
# # print(obsm_subset)




# ### Sort the adata.obs by total_counts to see which genes have exceedingly high counts

# In[105]:


obs = adata.obs

obs.sort_values(by=["transcript_counts"], ascending=False)

obs.to_csv(os.path.join(qc_path, f"obs_{dataset_name}.csv"))

print("All files saved to:", qc_path)


# In[106]:


# adata_X = pd.DataFrame(
#     adata.X.toarray(),
#     index=adata.obs_names,
#     columns=adata.var_names
# )

# ## Sum counts per gene across your subset cells
# gene_totals = adata_subset_X.sum(axis=0)  # sum down rows

# # Sort descending and display
# top_genes = gene_totals.sort_values(ascending=False)

# # print(top_genes.head(20))  # top 20 genes

# # # Optional: see all genes with counts > 0
# # print(top_genes[top_genes > 0])

# # # Optional: save to CSV
# # top_genes.to_csv("top_genes_subset.csv", header=["total_counts"])


# In[108]:


adata_X = pd.DataFrame(
    adata.X.toarray(),
    index=adata.obs_names,
    columns=adata.var_names
)

gene_total_object = adata_X.sum(axis=0)


top_genes_object = gene_total_object.sort_values(ascending=False).reset_index()
top_genes_object.columns = ['gene', 'total_counts']
top_genes_object.to_csv(os.path.join(qc_path, f"gene_counts_across_all_cells_{dataset_name}_test.csv"))


# In[110]:


## Plot the top 50 genes
top_genes_object = top_genes_object.iloc[0:101]

plt.figure(figsize=(25, 5))
plt.scatter(top_genes_object['gene'], top_genes_object['total_counts'])
plt.xticks(rotation=90)
plt.tight_layout()

plt.savefig(
    os.path.join(qc_path, f"top_50_genes_by_counts_{dataset_name}.png"),
    dpi=300,
    bbox_inches='tight'
    )


# ## Read in pre-labelled object

# In[111]:


#file_path =f"{output_path}/{dataset_name}_clustered_spatial_pc{pca_label}_nc{lambda_label}_r{res_label}.h5ad"
file_path =f"{raw_path}/{dataset_name}/{dataset_name}_clustered_spatial_pc{pca_label}_nc{lambda_label}_r{res_label}.h5ad"

print("Exists?", os.path.exists(file_path))


# In[115]:


## Read in the processed AnnData file
#adata_lab = ad.read_h5ad(f"{output_path}/{dataset_name}_clustered_spatial_pc{pca_label}_nc{lambda_label}_r{res_label}.h5ad")
adata_lab = ad.read_h5ad(f"{raw_path}/{dataset_name}/{dataset_name}_clustered_spatial_pc{pca_label}_nc{lambda_label}_r{res_label}.h5ad")
adata_lab.obs_keys()


## Read in the banksy_dict dictionary as a .pkl file
## Use gzip to save it with compression
import gzip
import pickle
#with gzip.open(f"{output_path}/{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_banksy_dict.pkl.gz", "rb") as f:
with gzip.open(f"{raw_path}/{dataset_name}/{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_banksy_dict.pkl.gz", "rb") as f:
    banksy_dict = pickle.load(f)


## Read in the results_df data frame as a .pkl file
#with gzip.open(f"{output_path}/results_df_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}.pkl.gz", "rb") as f:
with gzip.open(f"{raw_path}/{dataset_name}/results_df_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}.pkl.gz", "rb") as f:
    results_df = pickle.load(f)


# In[116]:


import scipy.sparse as sp
import numpy as np

# Check if values are integers
if sp.issparse(adata.X):
    values = adata.X.data  # non-zero values only
else:
    values = adata.X.flatten()

print(values[:20])  # look at actual values
print(f"Are all integers: {np.all(values == values.astype(int))}")
print(f"Max value: {values.max()}")
print(f"Min value: {values.min()}")


# In[117]:


# Check the ratio
print(adata.obs['nCount_Xenium'].max())        # ~1800 - ALL transcripts
print(np.asarray(adata.X.sum(axis=1)).max())   # ~160  - gene panel only


# #### PCA dimensions

# In[119]:


import scipy.sparse as sp

if sp.issparse(adata_lab.X):
    X = adata_lab.X.toarray()
else:
    X = np.array(adata_lab.X) 

## Import the noise_equiv_singular_value() function 
from banksy_utils.pca import plot_singular_values
import numpy as np

# Transpose the matirx so that the columns are the genes (features) not the cell IDs
X = X.T
X.shape

# Calculate the higher singular value 
from banksy_utils.pca import noise_equiv_singular_value

noise_sv, all_permuted_svs = noise_equiv_singular_value(
    data = X,
    num_permutations = 50,
    #average_type = "mean",
    verbose = True
    )
print(f"Noise threshold (mean top singular value of permuted data): {noise_sv:.4f}")


# In[120]:


## Perform PCA 
from sklearn.decomposition import PCA

pca = PCA(n_components=50).fit(X)


## Plot singular values vs noise to determine the optimal number of PC above the noise threshold
from banksy_utils.pca import plot_singular_values

plot_singular_values(
    pca, 
    noise_highest_sv=noise_sv, 
    title=f"Singular values vs noise_{dataset_name}",
    figsize=(20,20),
    #save_path = f"{output_path}singular_values_vs_noise_{dataset_name}"
    )

plt.xlabel("Number of principal components", fontsize=20)
plt.ylabel("Number of singular values", fontsize=20)

plt.savefig(
    os.path.join(qc_path, f"singular_values_vs_noise_{dataset_name}.png"),
    dpi=300,
    bbox_inches='tight'
    )


# ## Display previously generated UMAP plots

# ##

# In[121]:


# Identify the largest label number and set that to max_num_labels

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

    if nonspatial_labels < spatial_labels:
        max_num_labels = spatial_labels
        print(f"The number of spatial labels {nonspatial_labels} is greater than nonspatial labels {spatial_labels}, so {spatial_labels} is assigned to max_num_labels which {max_num_labels}.")

    elif nonspatial_labels == spatial_labels:
        print(f"nonspatial and spatial decay have generated the same number of clusters, so {nonspatial_labels} is assigned to max_num_labels which is {max_num_labels}.")

    return(max_num_labels)


max_num_labels=determine_max_num_labels(nonspatial_labels, spatial_labels)


# In[122]:


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


# In[123]:


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


# In[124]:


## Save annotations in two different anndata objects (adata_spatial and adata_nonspatial)

adata_spatial, adata_nonspatial = create_spatial_nonspatial_adata(results_df,
                                    pca_dims,
                                    lambda_list, 
                                    resolutions,
                                    cluster2annotation_spatial,
                                    cluster2annotation_nonspatial)


# In[125]:


## Store the cluster labels as a string in .obs
banksy_dict[f"{nbr_weight_decay}"][0.2]['adata'] #lambda = 0.20
adata_spatial.obs[f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}"]= adata_spatial.obs[f"labels_scaled_gaussian_pc{pca_label}_nc{lambda_label}_r{res_label}"].astype(str)


# In[126]:


## The top 20 gene markers for each cluster are assessed for cell cluster annotation

new_labels = {
    '0' : 'Invasive_Neural_crest-like_melanoma_0',
    '1' : 'M2-like_Tumour_Associated_Macrophages_1',
    '2' : 'Inflamed_Bronchial_Epithelial_Cells_(low_confidence)_2',
    '3' : 'Cycling_melanoma_cells_3',
    '4' : 'Cancer_Associated_Fibroblasts_4',
    '5' : 'Vascular_Endothelial_Cells_5', 
    '6' : 'Cytotoxic_CD8+_and_CD4+_T_Cells_6',
    '7' : 'Immature_Myeloid_Cells_or_Tumour_Myeloid_Cells_7',
    '8' : 'Mesenchymal-like_Melanoma_Cells_8',
    '9' : 'Vascular_smooth_muscle_or_pericytes_9'
}

# cluster_col = f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}"
adata_spatial.obs[f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann"] = adata_spatial.obs["cell type"].map(new_labels)

new_labels


# #### Cluster level violin plot of cells that pass the minimum trancript count threshhold "threshold_passed_str"

# In[128]:


## Transfer the minimum threshold values from the adata.obs to the adata_spatial.obs
adata_spatial.obs['threshold_passed'] = adata.obs['threshold_passed']
## Transfer the minimum threshold category labels from the adata.obs to the adata_spatial.obs
adata_spatial.obs['threshold_passed_cat'] = adata.obs['threshold_passed_cat']


# In[129]:


obs = adata_spatial.obs
cell_label = "banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann"

plt.figure(figsize=(8, 6))
sns.violinplot(data=obs, y=obs["nCount_Xenium"],  x=obs[f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann"], log_scale=10, color="lightblue")

sns.stripplot(data=obs, y=obs["nCount_Xenium"], x=obs[f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann"], hue =obs["threshold_passed"], size=4)

plt.xlabel("Cluster", fontsize=15)
plt.ylabel("Total transcript counts per cell (log10)", fontsize=15)
plt.xticks(rotation=45, fontsize=12, ha='right', va='top')


plt.savefig(
    os.path.join(qc_path, f"threshold_passed_cells_per_cluster_violin_plot_{dataset_name}.png"),
    dpi=300,
    bbox_inches='tight'
    )


# # Top expressed genes per cluster

# ### Spatial clustering results

# In[130]:


import scanpy as sc

## Rank genes and identify cluster gene markers
sc.tl.rank_genes_groups(
    adata_spatial,
    groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}",
    method='wilcoxon',
    key_added=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_markers"
)


# In[131]:


## Check which cluster keys exist
adata_spatial.uns[f'banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_markers']['names'].dtype


# In[132]:


## Using the filter_ranked_genes_by_type to create a key for 
from banksy_utils.annotation_utils import filter_ranked_genes_by_type

key_filtered = filter_ranked_genes_by_type(
    adata_spatial,
    key=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_markers",
    top_n=100,
    gene_type = "raw",
    new_key_suffix="_raw",
)


# In[133]:


# Plot the gene markers for each cluster
sc.pl.rank_genes_groups(
    adata_spatial, 
    key=key_filtered,
    method='wilcoxon',
    #key=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_markers", 
    #key=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_markers_raw", 
    #key=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_markers_nbr_0", 
    #key=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_markers_nbr_1", 
    #key_added = 'banksy_cluster_pc35_nc0.20_r0.50_markers_raw',
    #key_added= "raw_test",
    n_genes=20, 
    fontsize=15,
    save= f"_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_spatial_top_gene_raw.png")


# In[134]:


# # Plot the gene markers for each cluster as a violin plot that conpares to every other cluster
# with rc_context({"figure.figsize": (9, 1.5)}):
#     sc.pl.rank_genes_groups_violin(
#         adata_spatial, 
#         n_genes=20, 
#         jitter=False, 
#         key = key_filtered
#         )


# In[135]:


sc.pl.violin(
    adata_spatial, 
    "total_counts", 
    groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}", 
    save= f"_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}.png"
    )


# In[136]:


# Set the figure directory to your desired location
sc.settings.figdir = output_path

## Generate heat map for top markers per cluster
sc.tl.dendrogram(
    adata_spatial,
    groupby=f'banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_raw',
    key_added=f'banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_raw'
)

## Determine vmin and vmax dynamically based on gene expression
expr_values = adata_spatial.X
vmax_dynamic = np.percentile(expr_values[expr_values > 0], 99)
vmin_dynamic = np.percentile(expr_values, 1)

sc.pl.rank_genes_groups_heatmap(
    adata_spatial, 
    key=f'banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_markers_raw',
    #key= 'raw_test',
    n_genes=5,
    cmap = "plasma",
    vmin=vmin_dynamic,
    vmax=vmax_dynamic,
    show_gene_labels = True,
    save = f"_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_raw_markers.png",
    figsize=(20, 15)
)


# In[137]:

groupby_key = f'banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann'

sc.tl.dendrogram(adata_spatial, groupby=groupby_key)

sc.pl.dendrogram(
    adata_spatial,
    groupby=groupby_key,
    save=f"_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}.png"
)

sc.pl.correlation_matrix(
    adata_spatial,
    groupby_key = f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann",
    save=f"{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_cluster_correlation_plot.png",
    figsize=(5, 3.5)
)


# ### Export the top 20 clusters in wide format to .csv

# In[138]:


from banksy_utils.annotation_utils import export_clusters_wide

key= f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_markers_raw"

export_clusters_wide(
    adata= adata_spatial,
    key= key,
    gene_type= "raw",
    top_n= 20,
    dataset_name= dataset_name,
    file_path= processed_path
)


# ### Export the top 20 clusters in long format with scores to .csv

# In[139]:


from banksy_utils.annotation_utils import export_cluster_markers

export_cluster_markers(
    adata = adata_spatial,
    key = key,
    top_n = 30,
    dataset_name = dataset_name,
    file_path = processed_path
)


# ### Add cell_type information to top clusters for manual annotation

# In[144]:


## Read in the exported cluster markers generated with "export_clusters_wide()"
## and merge with the cell type annotation master document
n_genes_label =  20 # Set this to the number of genes that were exorted as top markers for each cluster
top_markers = pd.read_csv(f"{processed_path}/cluster_top_{n_genes_label}_genes_with_scores_{dataset_name}_{key}.csv", index_col=0) #index_col = 0 prevents insertion of an unnamed col
top_markers = top_markers.sort_index()
#print(top_markers)

## Read in the master annotation file which contains all the curated cell type labels for each marker gene
#master_markers_file = pd.read_csv(f"{raw_path}/xenium_gene_list_annotation_master.csv")
master_markers_file = pd.read_csv(f"{raw_path}/xenium_gene_list_annotation_master_v1_March_2026_Tier1.csv")

## Create a data frame of the genes and their corresponding primary and secondary cell type annotations 
cell_annotations = master_markers_file[["Gene", "primary_annotation", "secondary_annotation"]]
#cell_annotations = cell_annotations.sort_index()

## Merge the top marker results and the cell type annotations
merged_markers = (
    pd.merge(
        top_markers.reset_index(), 
        cell_annotations, 
        left_on="gene", 
        right_on="Gene", 
        how="left", 
        sort="False",
        validate="m:1"
    )
    .set_index("index")
    .loc[top_markers.index]
)

merged_markers = merged_markers.drop("Gene", axis=1)
print(merged_markers)

## Export the cell type annotated top markers to .csv to use in manual annotation and ChatGPT annotation
merged_markers.to_csv(f"{processed_path}cell_type_cluster_top_{n_genes_label}_genes_{dataset_name}_{key}.csv")


# In[145]:


raw_path


# In[146]:


# # Save AnnData objects with cluster placeholder resultsn for downstream annotation

# adata_spatial.write_h5ad(
#     filename = f"{output_path}{dataset_name}_clustered_spatial_pc{pca_label}_nc{lambda_label}_r{res_label}.h5ad",
#     compression="gzip"
#     )

# adata_nonspatial.write_h5ad(
#     filename = f"{output_path}{dataset_name}_clustered_nonspatial_pc{pca_label}_r{res_label}.h5ad"
#)


# In[147]:


adata_spatial.obs


# In[148]:


## The top 20 gene markers for each cluster are assessed for cell cluster annotation

new_labels = {
    '0' : 'Invasive_Neural_crest-like_melanoma_0',
    '1' : 'M2-like_Tumour_Associated_Macrophages_1',
    '2' : 'Inflamed_Bronchial_Epithelial_Cells_(low_confidence)_2',
    '3' : 'Cycling_melanoma_cells_3',
    '4' : 'Cancer_Associated_Fibroblasts_4',
    '5' : 'Vascular_Endothelial_Cells_5', 
    '6' : 'Cytotoxic_CD8+_and_CD4+_T_Cells_6',
    '7' : 'Immature_Myeloid_Cells_or_Tumour_Myeloid_Cells_7',
    '8' : 'Mesenchymal-like_Melanoma_Cells_8',
    '9' : 'Vascular_smooth_muscle_or_pericytes_9'
}

cluster_col = f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}"
adata_spatial.obs[f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann"] = adata_spatial.obs["cell type"].map(new_labels)

new_labels


# In[149]:


## Store the UMAP results in the .obsm slot
adata_spatial.obsm['X_umap'] = banksy_dict['scaled_gaussian'][0.2]['adata'].obsm[f"reduced_pc_{pca_label}_umap"]


# ## Visualisations

# ### Cluster Annotated UMAP

# In[150]:


# Plot the umap with bulk labels
import scanpy as sc
from matplotlib.pyplot import rc_context

color_vars = [
f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann"
]

with rc_context({"figure.figsize": (3,3)}):
    sc.pl.umap(
        adata_spatial, 
        color=color_vars, 
        s=50, 
        frameon=True, 
        vmax="p99",
        add_outline=True,
        legend_fontsize=10,
        title=f"{dataset_name} clusters",
        save = f"_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}.png"
        )


# In[151]:


# Plot the umap with bulk labels

color_vars = [
f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann"
]

fig, axes = plt.subplots(1, 2, figsize=(20,10))

#with rc_context({"figure.figsize": (3,3)}):
sc.pl.umap(
    adata_spatial, 
    color=color_vars, 
    s=5, 
    frameon=True, 
    vmax="p99",
    add_outline=True,
    legend_fontsize=10,
    title=f"{dataset_name} clusters",
    ax=axes[0],
    show=False
    #save = f"_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}.png",
    )
plt.legend(fontsize=20, title_fontsize=50)

color_vars = [
"threshold_passed_cat"
]

#with rc_context({"figure.figsize": (3,3)}):
sc.pl.umap(
    adata_spatial, 
    color=color_vars, 
    s=5, 
    frameon=True, 
    vmax="p99",
    add_outline=True,
    legend_fontsize=10,
    title=f"{dataset_name} clusters",
    ax=axes[1],
    show=False
    #save = f"_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}.png",
    )

plt.tight_layout()
plt.legend(fontsize=20, title_fontsize=50)
plt.savefig(os.path.join(qc_path, f"umap_{dataset_name}_threshold.png"), dpi=300, bbox_inches='tight')
plt.show()


# ### Gene of interest UMAP and violin plotting

# In[152]:


## Create a path for gene of interest UMAP plots
umap_path = os.path.join(base_dir,f"output/{dataset_name}/umap/")

if not os.path.isdir(umap_path):
    os.makedirs(umap_path)
    print(f"Directory '{umap_path}' created successfully.")
    
else:
    print(f"Directory '{umap_path}' already exists.")


## Create a path for gene of interest violin plots
violin_path = os.path.join(base_dir,f"output/{dataset_name}/violin/")

if not os.path.isdir(violin_path):
    os.makedirs(violin_path)
    print(f"Directory '{violin_path}' created successfully.")
    
else:
    print(f"Directory '{violin_path}' already exists.")

## Create a path for gene of interest violin plots
dotplot_path = os.path.join(base_dir,f"output/{dataset_name}/dotplot/")

if not os.path.isdir(dotplot_path):
    os.makedirs(dotplot_path)
    print(f"Directory '{dotplot_path}' created successfully.")
    
else:
    print(f"Directory '{dotplot_path}' already exists.")

## Create a path for gene of interest UMAP plots
neighbour_path = os.path.join(base_dir,f"output/{dataset_name}/neighbourhood_centrality_co-occurence/")

if not os.path.isdir(neighbour_path):
    os.makedirs(neighbour_path)
    print(f"Directory '{neighbour_path}' created successfully.")
    
else:
    print(f"Directory '{neighbour_path}' already exists.")


## Create a path for gene of interest spatial scatter plots
spatial_scatter_path = os.path.join(base_dir,f"output/{dataset_name}/spatial_scatter/")

if not os.path.isdir(spatial_scatter_path):
    os.makedirs(spatial_scatter_path)
    print(f"Directory '{spatial_scatter_path}' created successfully.")
    
else:
    print(f"Directory '{spatial_scatter_path}' already exists.")


# In[153]:


output_path


# ## UMAPs

# In[154]:


# Set the figure directory to your desired location
sc.settings.figdir = umap_path


# #### DSG2

# In[155]:


color_vars = [
"DSG2"
]

with rc_context({"figure.figsize": (8,5)}):
    sc.pl.umap(
        adata_spatial, 
        color=color_vars, 
        s=20, 
        frameon=True, 
        vmax="p99",
        add_outline=True,
        legend_fontsize=9,
        title=[f"DSG2 expression in {dataset_name}"],
        save = f"_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_DSG2_response_genes.png"
        )


# ### SERPINE1

# In[156]:


color_vars = [
"SERPINE1"
]

with rc_context({"figure.figsize": (8,5)}):
    sc.pl.umap(
        adata_spatial, 
        color=color_vars, 
        s=20, 
        frameon=True, 
        vmax="p99",
        add_outline=True,
        legend_fontsize=9,
        title=[f"SERPINE1 expression in {dataset_name}"],
        save = f"_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_SERPINE1_response_genes.png"
        )


# ### DKK1

# In[157]:


color_vars = [
"DKK1"
]

with rc_context({"figure.figsize": (8,5)}):
    sc.pl.umap(
        adata_spatial, 
        color=color_vars, 
        s=20, 
        frameon=True, 
        vmax="p99",
        add_outline=True,
        legend_fontsize=9,
        title=[f"DKK1 expression in {dataset_name}"],
        save = f"_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_DKK1_response_genes.png"
        )


# ### T cell markers

# In[158]:


color_vars = [
"FOXP3", "GZMB"
]

with rc_context({"figure.figsize": (8,5)}):
    sc.pl.umap(
        adata_spatial, 
        color=color_vars, 
        s=20, 
        frameon=True, 
        vmax="p99",
        add_outline=True,
        legend_fontsize=9,
        title=[f"FOXP3 expression in {dataset_name}", f"GZMB expression in {dataset_name}"],
        save = f"_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_T_cell_markers_response_genes.png"
        )


# ### Integrin genes

# In[159]:


color_vars = [
"ITGAE", 
"ITGA1", 
"ITGA4",
"ITGA5"
]

with rc_context({"figure.figsize": (8,5)}):
    sc.pl.umap(
        adata_spatial, 
        color=color_vars, 
        s=20, 
        frameon=True, 
        vmax="p99",
        add_outline=True,
        legend_fontsize=9,
        ncols=2,
        title=[
            f"ITGAE expression in {dataset_name}", 
            f"ITGA1 expression in {dataset_name}", 
            f"ITGA4 expression in {dataset_name}", 
            f"ITGA5 expression in {dataset_name}"],
        save = f"_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_integrin_genes.png"
        )


# ### TPD52 family

# In[160]:


color_vars = [
"TPD52", "TPD52L1", "TPD52L2"
]

with rc_context({"figure.figsize": (8,5)}):
    sc.pl.umap(
        adata_spatial, 
        color=color_vars, 
        s=20, 
        frameon=True, 
        vmax="p99",
        add_outline=True,
        legend_fontsize=9,
        ncols=2,
        title=[f"TPD52 expression in {dataset_name}", f"TPD52L1 expression in {dataset_name}", f"TPD52L2 expression in {dataset_name}"],
        save = f"_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_TPD52-like_genes.png"
        )


# #### Selectins

# In[161]:


color_vars = [
"SELE", "SELL", "SELP"
]

with rc_context({"figure.figsize": (8,5)}):
    sc.pl.umap(
        adata_spatial, 
        color=color_vars, 
        s=20, 
        frameon=True, 
        vmax="p99",
        add_outline=True,
        legend_fontsize=9,
        ncols=2,
        title=[f"SELE expression in {dataset_name}", f"SELL expression in {dataset_name}", f"SELP expression in {dataset_name}"],
        save = f"_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_selectins.png"
        )


# ### Fucosyltransferase genes

# In[162]:


color_vars = [
"FUT6", "FUT7"
]

with rc_context({"figure.figsize": (8,5)}):
    sc.pl.umap(
        adata_spatial, 
        color=color_vars, 
        s=20, 
        frameon=True, 
        vmax="p99",
        add_outline=True,
        legend_fontsize=9,
        title=[f"FUT6 expression in {dataset_name}", f"FUT7 expression in {dataset_name}"],
        save = f"_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_FUT6_FUT7.png"
        )


# ### Adhesion molecules

# In[163]:


color_vars = [
"VCAM1", "ICAM1"
]

with rc_context({"figure.figsize": (8,5)}):
    sc.pl.umap(
        adata_spatial, 
        color=color_vars, 
        s=20, 
        frameon=True, 
        vmax="p99",
        add_outline=True,
        legend_fontsize=9,
        title=[f"ICAM1 expression in {dataset_name}", f"VCAM1 expression in {dataset_name}"],
        save = f"_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_adhesion_genes.png"
        )


# ### ICI response genes

# In[164]:


color_vars = [
"PDCD1", "CD274", "CTLA4"
]

with rc_context({"figure.figsize": (8,5)}):
    sc.pl.umap(
        adata_spatial, 
        color=color_vars, 
        s=20, 
        frameon=True, 
        vmax="p99",
        add_outline=True,
        legend_fontsize=9,
        ncols=2,
        title=[f"PDCD1 expression in {dataset_name}", f"CD274 expression in {dataset_name}", f"CTLA4 expression in {dataset_name}"],
        save = f"_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_ICI_response_genes.png"
        )


# ## Gene expression violin plots

# In[165]:


# Set the figure directory to your desired location
sc.settings.figdir = violin_path


# In[166]:


with rc_context({"figure.figsize": (10, 8)}):
    ax= sc.pl.violin(
        adata_spatial,
        ["nFeature_Xenium"],
        groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann",
        rotation=270,
        stripplot=False,  # remove the internal dots
        inner="box",  # adds a boxplot inside violins
        save=f"{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_gene_count_per_cluster_violin_plot.png"
    )


# ### DSG2

# In[167]:


# with rc_context({"figure.figsize": (10,8)}):
#     sc.pl.violin(
#         adata_spatial, 
#         [
#             "DSG2"
#             ], 
#         #groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}",
#         groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann",
#         rotation=270,
#         save=f"{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_DSG2_family_violin_plot.png"
#         )

with rc_context({"figure.figsize": (10,8)}):
   ax = sc.pl.violin(
        adata_spatial, 
        keys = [
            "DSG2"
            ], 
        #groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}",
        groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann",
        rotation=270,
        #save=f"{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_DSG2_family_violin_plot.png"
        show=False
        )
ax.set_yscale('log')
ax.set_xlabel('Cell type', fontsize=12)
ax.set_ylabel("DSG2 expression (log)")

plt.tight_layout()
plt.savefig(
    os.path.join(qc_path, "{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_DSG2_family_violin_plot.png"),
    dpi=300,
    bbox_inches='tight'
)


# ### SERPINE1

# In[168]:


# with rc_context({"figure.figsize": (10, 8)}):
#     sc.pl.violin(
#         adata_spatial, 
#         [
#             "SERPINE1"
#             ], 
#         #groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}",
#         groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann",
#         rotation=270,
#         save=f"{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_SERPINE1_family_violin_plot.png"
#         )

with rc_context({"figure.figsize": (10,8)}):
   ax = sc.pl.violin(
        adata_spatial, 
        keys = [
            "SERPINE1"
            ], 
        #groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}",
        groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann",
        rotation=270,
        #save=f"{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_DSG2_family_violin_plot.png"
        show=False
        )
ax.set_yscale('log')
ax.set_xlabel('Cell type', fontsize=12)
ax.set_ylabel("SERPINE1 expression (log)")

plt.tight_layout()
plt.savefig(
    os.path.join(qc_path, "{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_SERPINE1_family_violin_plot.png"),
    dpi=300,
    bbox_inches='tight'
)


# ### DKK1

# In[169]:


# with rc_context({"figure.figsize": (10,8)}):
#     sc.pl.violin(
#         adata_spatial, 
#         [
#             "DKK1"
#             ], 
#         #groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}",
#         groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann",
#         rotation=270,
#         save=f"{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_DKK1_family_violin_plot.png"
#         )

with rc_context({"figure.figsize": (10,8)}):
   ax = sc.pl.violin(
        adata_spatial, 
        keys = [
            "DKK1"
            ], 
        #groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}",
        groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann",
        rotation=270,
        #save=f"{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_DSG2_family_violin_plot.png"
        show=False
        )
ax.set_yscale('log')
ax.set_xlabel('Cell type', fontsize=12)
ax.set_ylabel("DKK1 expression (log)")

plt.tight_layout()
plt.savefig(
    os.path.join(qc_path, "{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_DKK1_family_violin_plot.png"),
    dpi=300,
    bbox_inches='tight'
)


# ### T cell markers

# In[170]:


# with rc_context({"figure.figsize": (10, 8)}):
#     sc.pl.violin(
#         adata_spatial, 
#         [
#             "FOXP3", 
#             "GZMB"
#             ], 
#         #groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}",
#         groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann",
#         rotation=270,
#         save=f"{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_T_cell_markers_violin_plot.png"
#         )

with rc_context({"figure.figsize": (10,8)}):
   axes = sc.pl.violin(
        adata_spatial, 
        keys = [
             "FOXP3", 
             "GZMB"
            ], 
        #groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}",
        groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann",
        rotation=270,
        #save=f"{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_DSG2_family_violin_plot.png"
        show=False
        )

# Apply customisation to each graph axis
for ax in axes:
    ax.set_yscale('log')
    ax.set_xlabel('Cell type', fontsize=12)

axes[0].set_ylabel("T cell markers expression (log)")

plt.tight_layout()
plt.savefig(
    os.path.join(qc_path, "{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_T_cell_markers_violin_plot.png"),
    dpi=300,
    bbox_inches='tight'
)


# ### Integrins

# In[171]:


with rc_context({"figure.figsize": (6, 5)}):
    sc.pl.violin(
        adata_spatial, 
        [
            "ITGAE", 
            "ITGA1", 
            "ITGA4",
            "ITGA5"
            ], 
        #groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}",
        groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann",
        rotation=270,
        save=f"{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_integrin_genes_violin_plot.png"
        )


# ### TPD52 genes

# In[172]:


with rc_context({"figure.figsize": (6, 5)}):
    sc.pl.violin(
        adata_spatial, 
        [
            "TPD52",
            "TPD52L1",
            "TPD52L2"
            ], 
        #groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}",
        groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann",
        rotation=270,
        save=f"{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_TPD52_genes_violin_plot.png"
        )


# #### Selectins

# In[173]:


with rc_context({"figure.figsize": (6,5)}):
    sc.pl.violin(
        adata_spatial, 
        [
        "SELE", 
        "SELL", 
        "SELP"
        ], 
        #groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}",
        groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann",
        rotation=270,
        save=f"{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_selectin_genes_violin_plot.png"
        )


# ### Fucosyltransferase genes

# In[174]:


with rc_context({"figure.figsize": (6, 5)}):
    sc.pl.violin(
        adata_spatial, 
        [
            "FUT6", 
            "FUT7"
        ], 
        #groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}",
        groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann",
        rotation=270,
        save=f"{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_fucosyltransferase_genes_violin_plot.png"
        )


# ### Adhesion genes

# In[175]:


with rc_context({"figure.figsize": (6,5 )}):
    sc.pl.violin(
        adata_spatial, 
        [
            "VCAM1",
            "ICAM1"
        ], 
        #groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}",
        groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann",
        rotation=270,
        save=f"{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_adhesion_genes_violin_plot.png"
        )


# ### ICI response genes

# In[176]:


with rc_context({"figure.figsize": (6, 5)}):
    sc.pl.violin(
        adata_spatial, 
        [
            "PDCD1", 
            "CD274", 
            "CTLA4"
        ], 
        #groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}",
        groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann",
        rotation=270,
        save=f"{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_adhesion_genes_violin_plot.png"
        )


# ## Gene expression dot plots

# In[177]:


# Set the figure directory to your desired location
sc.settings.figdir = dotplot_path


# In[178]:


adata_spatial.obs[f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann"]


# In[179]:


# {'0': 'Differentiated_melanoma_cells_0',
#  '1': 'Invasive_melanoma_cells_1',
#  '2': 'Cytotoxic_CD8+_T_cells_2',
#  '3': 'Cancer_associated_fibroblasts_3',
#  '4': 'IFN-responsive_dedifferentiated_melanoma_cells_4',
#  '5': 'M2-like_tumour_associated_macrophages_5',
#  '6': 'Cycling_melanoma_cells_6',
#  '7': 'CD4+_naive_or_memory_T_cells_(TCF7+ SELL+)_7',
#  '8': 'Vascular_endothelial_cells_8',
#  '9': 'Plasma_B_cells_9',
#  '10': 'Tumor–Epidermal Interface Melanoma Cells_(low_confidence)_10',
#  '11': 'Invasive_melanoma_cells_(low_confidence)_11',
#  '12': 'Melanoma_cells_lymphatic_interface_(MITF+)_12',
#  '13': 'Proliferative–Invasive_melanoma_cells_(MITF+ VIM+)_13',
#  '14': 'Immune-Interacting_melanoma_cells_MITF+_(low_confidence)_14',
#  '15': 'Proliferative–Invasive_melanoma_cells_MITF+ VIM+ MYH11+_(low_confidence)_15',
#  '16': 'Interface_melanoma_cells_keratinocyte-immune-Interacting_(low_confidence)_16',
#  '17': 'Immune-Interacting_stromal_cells_(low_confidence)_17'}


# In[180]:


## cycle through the cluster labels and create dotplot for each element
from banksy_utils.annotation_utils import extract_marker_genes_dict

## Ensure that the patch line edge colopur is set to black in global matplotlib rcparams so that the bracket to group genes per cluster is shown
plt.rcParams['patch.edgecolor'] = 'black'

cluster_keys = list(new_labels.keys())
values = list(new_labels.values())

## Compute the dendrogram for clusters
sc.tl.dendrogram(
    adata_spatial,
    groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann"
    )

for i in range(len(cluster_keys)):
    key = cluster_keys[i]
    subset_cluster= new_labels[key]
    cluster_label = values[i]


    marker_genes_dict=extract_marker_genes_dict(
    adata=adata_spatial, 
    filtered_key=key_filtered, 
    gene_type='raw', 
    groupby = f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann",
    subset_cluster= subset_cluster,
    top_n=20
    )
    marker_genes_dict

    # Skip if no marker genes found
    if not marker_genes_dict or all(len(v) == 0 for v in marker_genes_dict.values()):
        print(f"Skipping {cluster_label} — no marker genes found")
        continue

## Generate dot plot
    with rc_context({"figure.figsize": (10,10)}):
        sc.pl.dotplot(
            adata_spatial, 
            #marker_genes_list,
            marker_genes_dict,
            groupby=f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann",
            #dendrogram=True,
            save=f"{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}_{cluster_label}_dot_plot.png"
            )


# In[181]:


cluster_keys


# ## Neighbourhood enrichment analysis

# In[182]:


# Set the figure directory to your desired location
sc.settings.figdir = neighbour_path


# In[183]:


# Adding spatial coordinates to .obsm
adata_spatial.obsm["spatial"] = adata.obs[["x", "y"]].to_numpy()


# In[184]:


sq.gr.spatial_neighbors(adata_spatial)


# In[185]:


cluster_key = f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_ann"
sq.gr.nhood_enrichment(adata_spatial, cluster_key=cluster_key)


# In[186]:


with rc_context({"figure.figsize": (10,10)}):
    sq.pl.nhood_enrichment(
        adata_spatial, 
        cluster_key=cluster_key,
        save=f"neighbourhood_enrichment_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}.png")


# ### Moran I's score

# In[187]:


# ## Can't create sdata if there are strings in the .uns, the code below will check the object for strings

# bad_uns    = [k for k in adata_spatial.uns.keys()    if not isinstance(k, str)]
# bad_obsm   = [k for k in adata_spatial.obsm.keys()   if not isinstance(k, str)]
# bad_varm   = [k for k in adata_spatial.varm.keys()   if not isinstance(k, str)]
# bad_layers = [k for k in adata_spatial.layers.keys() if not isinstance(k, str)]

# print("uns:", bad_uns)
# print("obsm:", bad_obsm)
# print("varm:", bad_varm)
# print("layers:", bad_layers)


# In[188]:


# ## Convert non strings to strings
# if getattr(adata_spatial, "is_view", False):
#     adata_spatial = adata_spatial.copy()

# # 1) uns: remap non-string keys to strings
# if bad_uns:
#     new_uns = {}
#     for k, v in adata_spatial.uns.items():
#         new_uns[str(k)] = v
#     adata_spatial.uns = new_uns  # safe reassignment

# # 2) obsm / varm / layers: rename non-string keys in-place
# def rename_nonstring_mapping(m):
#     # collect renames first to avoid modifying while iterating
#     to_rename = [(k, str(k)) for k in list(m.keys()) if not isinstance(k, str)]
#     for old, new in to_rename:
#         m[new] = m[old]
#         del m[old]

# rename_nonstring_mapping(adata_spatial.obsm)
# rename_nonstring_mapping(adata_spatial.varm)
# rename_nonstring_mapping(adata_spatial.layers)


# In[189]:


from spatialdata import SpatialData
adata_spatial.obs.rename(columns={"cell type": "cell_type"}, inplace=True)

# Create SpatialData with a single table (the AnnData)
sdata = SpatialData(tables={"cells": adata_spatial})

print(sdata)


# In[190]:


sdata.tables["subsample"] = sc.pp.subsample(adata_spatial, fraction=0.5, copy=True)


# In[191]:


adata_subsample = sdata.tables["subsample"]


# In[192]:


adata_subsample


# In[193]:


sq.gr.spatial_neighbors(adata_subsample, coord_type="generic", delaunay=True)
sq.gr.spatial_autocorr(
    adata_subsample,
    mode="moran",
    n_perms=100,
    n_jobs=1,
)
adata_subsample.uns["moranI"].head(10)


# In[194]:


## Store moran I scores as a data frame
moran_scores = adata_subsample.uns["moranI"]
moran_scores = moran_scores.reset_index(names=["Gene"])

## Exclude _nbr_0 and _nbr_1
moran_scores_raw = moran_scores[~moran_scores["Gene"].str.contains("_nbr_0", case=False, na=False)]
moran_scores_raw = moran_scores_raw[~moran_scores_raw["Gene"].str.contains("_nbr_1", case=False, na=False)]

## Filter for genes that have an FDR if less than 10%
moran_scores_raw = moran_scores_raw[moran_scores_raw["pval_sim_fdr_bh"] <= 0.10]

## Write Moran I scores to .csv
moran_scores_raw.to_csv(f"{processed_path}/moran_scores_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}.csv")


# In[195]:


## Plot spatial scatter plots for genes with the top 10 highest Moran I scores

## Create a list of the top scoring genes
top_moran = moran_scores_raw.sort_values(by="pval_sim", ascending=True)
top_moran = top_moran.iloc[:10]
top_moran = list(top_moran['Gene'])
top_moran

for i in top_moran:
    with rc_context({"figure.figsize": (12,8)}):
        sq.pl.spatial_scatter(
            adata_subsample,
            library_id="spatial",
            color=[i],
            shape=None,
            size=2,
            img=False,
            save= f"spatial_scatter_top_moran_I_{i}_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}.png"
)


# ## Compute centrality scores

# In[196]:


sq.gr.centrality_scores(adata_subsample, cluster_key=cluster_key)


# In[197]:


sq.pl.centrality_scores(adata_subsample, cluster_key=cluster_key, figsize=(22, 4))


# ## Compute co-occurence probability

# In[198]:


from spatialdata import SpatialData
sdata = SpatialData()


# In[199]:


sdata.tables["subsample"] = sc.pp.subsample(adata_subsample, fraction=0.5, copy=True)


# In[200]:


#adata_subsample = sdata.tables["subsample"]
adata_subsample


# In[201]:


adata.obs.keys()


# In[202]:


cluster_keys = list(new_labels.keys())


# In[203]:


print(new_labels)
print(cluster_keys)


# In[204]:


cluster_key_label = f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_raw"

## Calculate and plot the co-occurence of clusters
sq.gr.co_occurrence(
    adata_subsample,
    cluster_key=cluster_key_label
    )


for i in cluster_keys:     
    cluster_number = str(i)
     
    with rc_context({"figure.figsize": (12,12)}):
        sq.pl.co_occurrence(
            adata_subsample,
            cluster_key=cluster_key_label,
            clusters=i,
            figsize=(10, 10),
            save=f"co-occurence_cluster_{cluster_number}_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}.png"
            )

# ['banksy_cluster_pc35_nc0.20_r0.50_colors',
#  'banksy_cluster_pc35_nc0.20_r0.50_ann_colors',
#  'banksy_cluster_pc35_nc0.20_r0.50_raw_colors']

with rc_context({"figure.figsize": (12,8)}):
    sq.pl.spatial_scatter(
        adata_subsample,
        library_id="dataset_name",
        color=cluster_key_label,
        shape=None,
        size=2,
        save=f"spatial_cluster_scatter_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}.png"
        )


# In[205]:


# See all color related keys
[key for key in adata_subsample.uns.keys() if 'color' in key.lower()]

adata_subsample.uns["banksy_cluster_pc35_nc0.20_r0.50_raw_colors"]


# In[206]:


[col for col in adata_subsample.obs.columns if 'banksy' in col.lower()]


# ## Spatial expression for genes of interest

# In[208]:


# Set the figure directory to your desired location
sc.settings.figdir = spatial_scatter_path

goi_list = [
    "DSG2",
    "SERPINE1",
    "DKK1",
    "FOXP3", 
    "GZMB",
    "ITGAE",
    "ITGA1",
    "ITGA4",
    "ITGA5",
    "TPD52",
    "TPD52L1",
    "TPD52L2",
    "SELE", 
    "SELL", 
    "SELP",
    "FUT6", 
    "FUT7",
    "VCAM1",
    "ICAM1",
    "PDCD1", 
    "CD274", 
    "CTLA4",
]

# Set the figure directory to your desired location
sc.settings.figdir = spatial_scatter_path

for gene in goi_list:
    sq.pl.spatial_scatter(
        adata_subsample,
        library_id="spatial",
        color=[gene],
        shape=None,
        size=2,
        img=False,
        save= f"spatial_scatter_{gene}_{dataset_name}_pc{pca_label}_nc{lambda_label}_r{res_label}.png"
)

