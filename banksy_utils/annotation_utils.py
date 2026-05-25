# Author: Nathalie Nataren
# Date: 30/09/2025

# These helper functions can be used to make cluster gene marker extraction eaiser to aid manual cluster annotation

import copy
import pandas as pd
import os
import numpy as np

def filter_ranked_genes_by_type(
    adata, 
    key, # i.e., f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_markers_raw"
    gene_type : str, # "raw", "_nbr_0" or "_nbr_1"
    top_n : int, 
    new_key_suffix : str):
    """

    This function filters the ranked gene results (after running sc.tl.rank_genes_groups()) to include
    only genes of specific annotation and saves the results under a new key in 
    adata_spatial.uns, where:
    "raw" -> the gene's won expression in each cell (this captures cell identity)
    "_nbr_0" -> mean neighbour expression of a given gene, averaged over cell's spatial neighbour
    "_nbr_1" -> Azimuthal Gabor filter output for a given gene, this measures the local expression gradient around each cell (is a cell in a tissue boundary)

    Parameters
    ----------

    adata : AnnData object
        BANKSY clustered anndata object
    key : str
        The adata.obs column name created by `sc.tl.rank_genes_groups()` that indicates top markers for each cluster.
    gene_type : str
        select the gene type from those described above, "raw", "_nbr_0", and "_nbr_1".
    top_n : int
        The number of top genes to extract for a specified gene type.
    new_key_suffix: str
        The name for the new key that will be entered into adata.uns which can then be called by scanpy functions such as "sc.pl.rank_genes_groups()"

    """

    # Copy the full rank_genes_groups results
    original= adata.uns[key]
    filtered= copy.deepcopy(original) 

    # Get cluster names
    clusters = original['names'].dtype.names

    # Define all the fields to capture from .uns
    fields = ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']

    # Build empty array for each field
    new_arrays = {}

    for f in fields:
        # Gene names are strings,all other values are floats
        if f == "names":
            dtype = object
        else:
            dtype = float
        
        # Generate one column per cluster all with the same dtype
        column_spec = []
        for cluster_name in clusters:
            column_spec.append((cluster_name, dtype))
            #print(f"This is the cluster name '{cluster_name}'.")

        # Empty table: top_n rows, one column per cluster
        new_arrays[f] = np.recarray(shape=(top_n,), dtype=column_spec) 


    # # Create new structured arrays for all the .uns data
    # new_names = np.recarray(shape=(top_n,), dtype = [(c, object) for c in clusters])
    # new_scores = np.recarray(shape=(top_n,), dtype = [(c, float) for c in clusters])

    for c in clusters:
        # Get the full list of genes and scores
        genes = original['names'][c]
        #scores = original['scores'][c]

        # Filter the genes based on their gene_type annotation i.e., 'raw' for no suffix
        if gene_type == 'raw':
            mask=[not(g.endswith("_nbr_0") or g.endswith("_nbr_1")) for g in genes]
        elif gene_type == "nbr_0":
            mask=[g.endswith("_nbr_0") for g in genes]
        elif gene_type == "nbr_1":
            mask=[g.endswith("_nbr_1") for g in genes]
        else:
            raise ValueError("Invalid gene_type. Use 'raw', 'nbr_0' or 'nbr_1'.")

        # Indices of kept gene, top_n only, array of original positions of our top_n genes
        idx = np.where(mask)[0][:top_n]

        # Loop through and slice every field with the same indices to that we take the top_n across all, then pad to top_n if needed

        for f in fields:
            # Get the original values for field f and cluster c
            full_column = original[f][c]
            # Select only the positions for top_n that we want to keep
            kept_values = full_column[idx]
            #Convert to a list so we can append padding
            values = list(kept_values)

            # Decide on padding depnding on the type of missing data

            if f == "names":
                pad_value =""
            else:
                pad_value = np.nan
            
            # Calculate bow many entries short by
            n_missing = top_n -len(values)

            # Pad by n_missing
            values = values + [pad_value] * n_missing

            # Write into the correct column of the new array
            new_arrays[f][c] = values

    # Replace each filed in the filtered copy
    for f in fields:
        filtered[f] = new_arrays[f]

    new_key = key + new_key_suffix
    adata.uns[new_key] = filtered


    # Copy group labels to .obs for downstream plotting e.g., sc.pl.rank_genes_groups_heatmap(
    original_groupby = key.replace("_markers", "")
    new_groupby = original_groupby + new_key_suffix

    if original_groupby in adata.obs.columns:
        adata.obs[new_groupby] = adata.obs[original_groupby].copy()
        print(f"Copied cluster labels to adata.obs['{new_groupby}']")
    else:
         print(f"Could not find original groupby column '{original_groupby}' in adata.obs.")

    print(f"Stored filtered {gene_type} markers in adata.uns['{new_key}']")
    return new_key

        # Patch params so sc.pl.rank_genes_groups_* functions can find the right groupby column
    if 'params' in original and isinstance(original['params'], dict):
        filtered['params'] = dict(original['params'])
    else:
        filtered['params'] = {}
    filtered['params']['groupby'] = new_groupby

    print(f"Stored filtered {gene_type} markers in adata.uns['{new_key}']")
    return new_key

        # # Apply the filters and store the top_n genes and scores
        # # initialise empty lists
        # filtered_genes = []
        # filtered_scores = []

        # # Iterate through each gene and its corresponding mask value
        # for g, m, score in zip(genes, mask, scores):
        #     if m:
        #         filtered_genes.append(g)
        #         filtered_scores.append(score)
    
        # # Keep only the top N genes and scores
        # filtered_genes = filtered_genes[:top_n]
        # filtered_scores = filtered_scores[:top_n]

        # # pad the top genes if necessary
        # filtered_genes += [''] * (top_n - len(filtered_genes))
        # filtered_scores += [np.nan] * (top_n - len(filtered_scores))

        # # Assign to the structured arrays
        # new_names[c] = filtered_genes
        # new_scores[c] = filtered_scores

    # # store the filtered results under a new key

    # new_key = key + new_key_suffix
    # filtered['names'] = new_names
    # filtered['scores'] = new_scores
    # adata.uns[new_key] = filtered

    # # Copy group labels to .obs for downstream plotting e.g., sc.pl.rank_genes_groups_heatmap(
    # original_groupby = key.replace("_markers", "")
    # new_groupby = original_groupby + new_key_suffix

    # if original_groupby in adata.obs.columns:
    #     adata.obs[new_groupby] = adata.obs[original_groupby].copy()
    #     print(f"Copied cluster labels to adata.obs['{new_groupby}']")
    # else:
    #      print(f"Could not find original groupby column '{original_groupby}' in adata.obs.")

    # print(f"Stored filtered {gene_type} markers in adata.uns['{new_key}']")
    # return new_key

#######################################################################################

def export_clusters_wide(
    adata,
    key, # i.e., f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_markers_raw"
    gene_type,
    top_n,
    dataset_name,
    file_path
    ):

    """
    This function exports the top ranked genes per cluster after running sc.pl.rank_genes_groups(). 
    The genes are exported to a .csv in wide format, with the top marker genes with clusters being listed column-wise.
    """

    # Extract cluster genes
    cluster_genes = adata.uns[key]
    
    #Extract cluster gene data frame
    clusters = cluster_genes['names'].dtype.names

    # Create cluster gene data frame and pad genes
    top_genes = pd.DataFrame({
        cluster_group: (
            list(cluster_genes['names'][cluster_group][:top_n]) + 
            [None] * (top_n - len(cluster_genes['names'][cluster_group])) # Add padding if there are fewer genes than top_n
        )
        for cluster_group in clusters
    })

    print(top_genes)

    # Write the top genes to .csv
    top_genes.to_csv(f"{file_path}{dataset_name}_cluster_top_{top_n}_genes_{gene_type}.csv")

    ######################################################################################

def export_cluster_markers(
    adata,
    key, # i.e., f"banksy_cluster_pc{pca_label}_nc{lambda_label}_r{res_label}_markers_raw"
    top_n,
    dataset_name, 
    file_path
    ):
    
    """
    This function exports the top ranked genes per cluster with Wilcoxon z-score after running sc.pl.rank_genes_groups(),
    which can be used for manual annotation.
    The top marker genes for each cluster are exported to a .csv with the following columns:
        "cluster" => cluster number
        "gene" => HGNC gene symbol
        "score"  => the Wilcoxon z-score generated by sc.pl.rank_genes_groups()
    """
    os.makedirs(file_path, exist_ok=True)
    # Extract cluster genes
    cluster_genes = adata.uns[key]
    cluster_ids = cluster_genes['names'].dtype.names
    marker_data = adata.uns[key]
    
    all_cluster_genes = []

    n_genes_label = top_n-1
    for cluster in cluster_ids:
        genes = list(marker_data['names'][cluster][:top_n]) + \
            [None] * (top_n - len(marker_data['names'][cluster]))
        scores = list(marker_data['scores'][cluster][:top_n]) + \
            [None] * (top_n - len(marker_data['scores'][cluster]))

    ## Create a data frame of the genes and their z-scores
        cluster_df = pd.DataFrame({
            "cluster": int(cluster),
            "gene": genes, 
            "score": scores
            })

        all_cluster_genes.append(cluster_df)

    ## Combine into one long data frame
    all_cluster_genes_df = pd.concat(all_cluster_genes, ignore_index=True)
    all_cluster_genes_df.head(40)

    # Write to .csv
    all_cluster_genes_df.to_csv(f"{file_path}/cluster_top_{top_n}_genes_with_scores_{dataset_name}_{key}.csv")
    try:
        all_cluster_genes_df.to_csv(f"{file_path}/cluster_top_{top_n}_genes_with_scores_{dataset_name}_{key}.csv", index=False)
        print(f"Saved to {file_path}")
    except Exception as e:
        print(f"Failed to save: {e}")

    return(all_cluster_genes_df)

    #####################################################################################################

def extract_marker_genes_dict(
    adata, 
    filtered_key, 
    groupby, 
    subset_cluster='None',
    gene_type='raw', 
    top_n=20
    ):
    """
    Build a marker_genes_dict of raw top markers from the filtered rank_genes_groups in .uns to use for dotplot.
    Returns a dictionary that uses the manually annotated cluster labels stored in "banksy_cluster_pc{pca_label}_nc{res_label}_r{res_label}_ann"
    """
    results = adata.uns[filtered_key] # stores the results stored in the .uns slot
    clusters = results['names'].dtype.names # stores the clusters
    categories = adata.obs[groupby].cat.categories

    ## Check that the luster labels provided by subset_cluster
    assert len(clusters) == len(categories), (
        f"Mismatch between number of clusters in '{filtered_key}' ({len(clusters)})"
        f"and number of categories in '{groupby}' ({len(categories)}).\n"
        f"clusters found: {list(clusters)}\n"
        f"Categories found: {categories.tolist()}\n"
        f"Check that your new_labels dict matches the number of clusters you have provided."
    )

    marker_genes_dict = {} # initialise the marker_genes_dict

    for c in clusters:
        cluster_label = categories[int(c)] # this line maps the cluster number (first converting it to an int) to the human readable label stored in adata.obs[groupby]
        
        # skip subsetting if subsetting is not requested "None", but if a subset list is provided, 
        # this checks if the current cluster being assessed in the for loop is in the subset_cluster list, and if it is not, then the results for this cluster are not added to the dictionary
        if subset_cluster is not None and cluster_label not in subset_cluster:
            continue

        all_genes = results['names'][c]

        if gene_type == "raw":
            genes = [g for g in all_genes if not(g.endswith("_nbr_0") or g.endswith("_nbr_1")) and g != '']
        elif gene_type == "_nbr_0":
            genes = [g for g in all_genes if g.endswith("_nbr_0")]
        elif gene_type == "_nbr_1":
            genes = [g for f in all_genes if g.endswith("_nbr_1")]
        else:
            raise ValueError("gene_type must be 'raw', 'nbr_0' or 'nbr_1'.")


        # map numeric cluster (e.g., "0") to descriptive label (e.g., "CD4+ T cells")
        marker_genes_dict[cluster_label] = genes[:top_n]

    return marker_genes_dict