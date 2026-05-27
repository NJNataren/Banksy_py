#!/usr/bin/env python
# coding: utf-8

"""Export dotplot-ready gene expression summaries from configured AnnData objects."""

import argparse
import json
import os



parser = argparse.ArgumentParser(
    prog="export dotplot data from configured Xenium AnnData objects"
)
parser.add_argument(
    "--config",
    type=str,
    help="JSON config defining AnnData objects, markers, and output CSV.",
    required=True,
)
args = parser.parse_args()

with open(args.config) as f:
    cfg = json.load(f)

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

marker_file = cfg["marker_file"]
gene_column = cfg.get("gene_column", "Gene")
marker_group_column = cfg.get("marker_group_column")
expression_source = cfg.get("expression_source", "raw")
objects = cfg["objects"]
output_csv = cfg["output_csv"]

markers = pd.read_csv(marker_file)
marker_genes = markers[gene_column].dropna().astype(str).drop_duplicates().tolist()

if marker_group_column and marker_group_column in markers.columns:
    marker_groups = (
        markers[[gene_column, marker_group_column]]
        .dropna(subset=[gene_column])
        .drop_duplicates(subset=[gene_column])
        .set_index(gene_column)[marker_group_column]
        .astype(str)
        .to_dict()
    )
else:
    marker_groups = {gene: "marker" for gene in marker_genes}


def to_dense(x):
    if sparse.issparse(x):
        return x.toarray()
    return np.asarray(x)


rows = []

for obj in objects:
    sample = obj["sample"]
    resolution = str(obj.get("resolution", ""))
    adata_path = obj["adata_path"]
    groupby = obj["groupby"]

    print(f"Reading {adata_path}")
    adata = ad.read_h5ad(adata_path)

    if expression_source == "raw" and adata.raw is not None:
        expr_adata = adata.raw
        var_names = pd.Index(adata.raw.var_names)
        source_used = "raw"
    else:
        expr_adata = adata
        var_names = pd.Index(adata.var_names)
        source_used = "X"

    present_genes = [gene for gene in marker_genes if gene in var_names]
    missing_genes = sorted(set(marker_genes) - set(present_genes))

    if missing_genes:
        print(f"{sample} r{resolution}: missing {len(missing_genes)} marker genes")

    if not present_genes:
        print(f"{sample} r{resolution}: no marker genes found, skipping")
        continue

    expr = to_dense(expr_adata[:, present_genes].X)
    clusters = adata.obs[groupby].astype(str)

    for cluster_id in sorted(clusters.unique()):
        mask = (clusters == cluster_id).to_numpy()
        cluster_expr = expr[mask, :]
        mean_expression = np.asarray(cluster_expr.mean(axis=0)).ravel()
        percent_expressing = np.asarray((cluster_expr > 0).mean(axis=0)).ravel() * 100

        for i, gene in enumerate(present_genes):
            rows.append(
                {
                    "sample": sample,
                    "resolution": resolution,
                    "cluster_id": cluster_id,
                    "sample_cluster": f"{sample}__r{resolution}__cluster_{cluster_id}",
                    "groupby": groupby,
                    "gene": gene,
                    "marker_group": marker_groups.get(gene, "marker"),
                    "mean_expression": mean_expression[i],
                    "percent_expressing": percent_expressing[i],
                    "n_cells": int(mask.sum()),
                    "expression_source": source_used,
                    "adata_path": adata_path,
                }
            )

out = pd.DataFrame(rows)
os.makedirs(os.path.dirname(output_csv), exist_ok=True)
out.to_csv(output_csv, index=False)
print(f"Wrote {output_csv}")
