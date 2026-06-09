#!/usr/bin/env python
# coding: utf-8

"""
Title: Create Expression AnnData With BANKSY Cluster Labels
Date: 2026-06-09
Summary: Build a clean gene-expression AnnData object for downstream plotting by
normalizing the original Xenium expression matrix and copying existing BANKSY
cluster labels into `.obs`. The output preserves raw counts in a layer,
log-normalized expression in `.X` and `.raw`, can attach archived cell-type
labels from a cell-level CSV, and avoids using BANKSY-expanded features as
marker-expression values.
"""

import argparse
import json
import os


# Import the scientific stack after argument parsing inside main() so `--help`
# still works even if the active shell is not in the project conda environment.


def parse_args():
    """Parse command-line arguments for the expression-object builder."""
    parser = argparse.ArgumentParser(
        prog="create expression AnnData with BANKSY cluster labels"
    )
    parser.add_argument(
        "--config",
        type=str,
        help="JSON config defining expression AnnData, BANKSY objects, and output path.",
        required=True,
    )
    return parser.parse_args()


def load_config(config_path):
    """Load the JSON configuration file.

    Args:
        config_path: Path to the JSON config file.

    Returns:
        Dictionary containing input paths, output path, and cluster-label settings.
    """
    with open(config_path) as f:
        return json.load(f)


def normalize_expression_adata(expr_adata, filter_obs_key, filter_min, counts_layer):
    """Filter, preserve counts, and log-normalize an expression AnnData object.

    Args:
        expr_adata: Original Xenium expression AnnData object.
        filter_obs_key: Observation column used to remove empty cells.
        filter_min: Minimum value required in `filter_obs_key`.
        counts_layer: Name of the layer used to store raw count-like values.

    Returns:
        AnnData object with raw counts in `layers[counts_layer]`, log-normalized
        expression in `.X`, and log-normalized expression copied to `.raw`.
    """
    import scanpy as sc

    # Reapply the same empty-cell filter used before clustering so cell IDs line
    # up with the saved BANKSY spatial objects.
    if filter_obs_key in expr_adata.obs.columns:
        print(f"Filtering cells where {filter_obs_key} > {filter_min}")
        expr_adata = expr_adata[expr_adata.obs[filter_obs_key] > filter_min].copy()
    else:
        print(f"Filter column {filter_obs_key!r} not found; skipping cell filter")

    # Preserve raw count-like values before normalization. Layers are valid here
    # because this object is still the original cell-by-gene expression matrix.
    print("Saving raw counts to layer:", counts_layer)
    expr_adata.layers[counts_layer] = expr_adata.X.copy()

    print("Normalizing total counts")
    sc.pp.normalize_total(expr_adata, inplace=True)

    print("Applying log1p transform")
    sc.pp.log1p(expr_adata)

    print("Saving log-normalized expression to adata.raw")
    expr_adata.raw = expr_adata.copy()

    return expr_adata


def copy_banksy_labels(expr_adata, cluster_adata, groupby, cluster_adata_path):
    """Copy one BANKSY cluster-label column into expression AnnData `.obs`.

    Args:
        expr_adata: Normalized expression AnnData receiving cluster labels.
        cluster_adata: Existing BANKSY spatial AnnData containing cluster labels.
        groupby: Cluster-label column to copy from `cluster_adata.obs`.
        cluster_adata_path: Path used only for informative error messages.

    Raises:
        KeyError: If `groupby` is absent from `cluster_adata.obs`.
        ValueError: If any expression cells do not receive a cluster label.
    """
    if groupby not in cluster_adata.obs.columns:
        raise KeyError(
            f"{groupby!r} was not found in {cluster_adata_path}. "
            f"Available obs columns: {list(cluster_adata.obs.columns)}"
        )

    # Align by cell ID rather than row order so labels remain correct even if an
    # object was filtered or read back in a different order.
    labels = cluster_adata.obs[groupby].reindex(expr_adata.obs_names)
    missing_labels = int(labels.isna().sum())
    if missing_labels > 0:
        raise ValueError(
            f"{missing_labels} cells in the expression object did not receive "
            f"a label from {groupby!r} in {cluster_adata_path}"
        )

    expr_adata.obs[groupby] = labels.astype("category")
    print(f"Copied {groupby} to expression AnnData obs")


def copy_optional_annotation(
    expr_adata,
    cluster_adata,
    annotation_column,
    annotation_output,
):
    """Copy a readable annotation column when it is complete and available.

    Args:
        expr_adata: Normalized expression AnnData receiving annotations.
        cluster_adata: Existing BANKSY spatial AnnData containing annotations.
        annotation_column: Source column in `cluster_adata.obs`.
        annotation_output: Destination column in `expr_adata.obs`.
    """
    if annotation_column not in cluster_adata.obs.columns:
        return

    annotations = cluster_adata.obs[annotation_column].reindex(expr_adata.obs_names)
    missing_annotations = int(annotations.isna().sum())
    if missing_annotations == 0:
        expr_adata.obs[annotation_output] = annotations.astype("category")
        print(f"Copied {annotation_column!r} to obs[{annotation_output!r}]")
    else:
        print(f"Skipped {annotation_column!r}: {missing_annotations} missing annotations")


def attach_cluster_metadata(expr_adata, objects):
    """Attach cluster labels from all configured BANKSY objects.

    Args:
        expr_adata: Normalized expression AnnData to annotate.
        objects: List of config entries with BANKSY object paths and groupby keys.

    Returns:
        The same expression AnnData with cluster metadata added to `.obs`.
    """
    import anndata as ad

    # Copy cluster labels from each configured BANKSY spatial object into obs.
    # This keeps the clustering result but leaves behind the BANKSY-expanded
    # matrix, which should not be used as marker-expression signal.
    for obj in objects:
        cluster_adata_path = obj["adata_path"]
        groupby = obj["groupby"]
        annotation_column = obj.get("annotation_column", "cell type")
        annotation_output = obj.get("annotation_output", f"{groupby}_annotation")

        print(f"Reading BANKSY AnnData: {cluster_adata_path}")
        cluster_adata = ad.read_h5ad(cluster_adata_path)

        copy_banksy_labels(expr_adata, cluster_adata, groupby, cluster_adata_path)
        copy_optional_annotation(
            expr_adata,
            cluster_adata,
            annotation_column,
            annotation_output,
        )

    return expr_adata


def resolve_archive_label_config(cfg):
    """Return archive-label settings from either new or legacy config keys.

    Args:
        cfg: Script 02 config dictionary.

    Returns:
        Archive-label config dictionary, or `None` when no archive labels should
        be attached.
    """
    if cfg.get("archive_labels"):
        return cfg["archive_labels"]

    if cfg.get("archive_label_csv"):
        return {
            "label_csv": cfg["archive_label_csv"],
            "sample": cfg.get("archive_label_sample", cfg.get("sample")),
            "columns": cfg.get("archive_label_columns"),
        }

    return None


def attach_archive_labels(expr_adata, archive_cfg):
    """Attach archived cell-level labels to expression AnnData `.obs`.

    Args:
        expr_adata: Normalized expression AnnData receiving archived metadata.
        archive_cfg: Config block describing the archived label CSV and columns.

    Returns:
        The same expression AnnData with archived metadata columns added.

    Raises:
        ValueError: If required config or CSV columns are missing, or if cell IDs
            are duplicated within the selected archived sample.
    """
    import pandas as pd

    if not archive_cfg:
        return expr_adata

    label_csv = archive_cfg["label_csv"]
    sample = archive_cfg.get("sample")
    cell_id_column = archive_cfg.get("cell_id_column", "cell_id")
    sample_column = archive_cfg.get("sample_column", "sample")
    columns = archive_cfg.get(
        "columns",
        {
            "cluster_id": "archive_cluster_id",
            "cell_type_label": "archive_cell_type_label",
        },
    )
    missing_values = archive_cfg.get(
        "missing_values",
        {
            "archive_cluster_id": "missing_archive_cluster",
            "archive_cell_type_label": "missing_archive_label",
        },
    )
    add_has_label_column = archive_cfg.get("add_has_label_column", True)
    has_label_column = archive_cfg.get("has_label_column", "has_archive_label")

    print(f"Reading archived label table: {label_csv}")
    archive_labels = pd.read_csv(label_csv, low_memory=False)

    required_columns = {cell_id_column, *columns.keys()}
    if sample:
        required_columns.add(sample_column)
    missing_columns = required_columns - set(archive_labels.columns)
    if missing_columns:
        raise ValueError(
            f"{label_csv} is missing required archived label columns: "
            f"{sorted(missing_columns)}"
        )

    if sample:
        archive_labels = archive_labels[
            archive_labels[sample_column].astype(str) == str(sample)
        ].copy()
        if archive_labels.empty:
            raise ValueError(f"No archived labels for sample {sample!r} in {label_csv}")

    if archive_labels[cell_id_column].duplicated().any():
        duplicate_count = int(archive_labels[cell_id_column].duplicated().sum())
        raise ValueError(
            f"{label_csv} has {duplicate_count} duplicated cell IDs after sample "
            "filtering; cannot safely align archived labels"
        )

    archive_labels = archive_labels.set_index(cell_id_column)
    matched_cells = expr_adata.obs_names.isin(archive_labels.index)
    print(
        "Archived label coverage: "
        f"{int(matched_cells.sum())}/{expr_adata.n_obs} expression cells matched"
    )

    for source_column, output_column in columns.items():
        values = archive_labels[source_column].reindex(expr_adata.obs_names)
        fill_value = missing_values.get(output_column, "missing_archive_label")
        values = values.fillna(fill_value).astype(str)
        expr_adata.obs[output_column] = pd.Categorical(values)
        print(f"Copied archived {source_column!r} to obs[{output_column!r}]")

    if add_has_label_column:
        expr_adata.obs[has_label_column] = matched_cells
        print(f"Added obs[{has_label_column!r}]")

    return expr_adata


def main():
    """Create and save a clean expression AnnData with BANKSY labels in `.obs`."""
    args = parse_args()
    cfg = load_config(args.config)

    import anndata as ad

    expression_adata_path = cfg["expression_adata_path"]
    output_h5ad = cfg["output_h5ad"]
    objects = cfg["objects"]

    # These optional config values mirror choices made in the clustering workflow.
    filter_obs_key = cfg.get("filter_obs_key", "nCount_Xenium")
    filter_min = cfg.get("filter_min", 0)
    counts_layer = cfg.get("counts_layer", "counts")

    # Start from the original expression object, not the BANKSY feature object.
    # This keeps the output object gene-focused and suitable for marker dotplots.
    print(f"Reading expression AnnData: {expression_adata_path}")
    expr_adata = ad.read_h5ad(expression_adata_path)

    expr_adata = normalize_expression_adata(
        expr_adata,
        filter_obs_key,
        filter_min,
        counts_layer,
    )
    expr_adata = attach_cluster_metadata(expr_adata, objects)
    expr_adata = attach_archive_labels(
        expr_adata,
        resolve_archive_label_config(cfg),
    )

    # Save a clean downstream object: log-normalized expression in .X/.raw, raw
    # counts in a layer, and BANKSY cluster labels in obs.
    os.makedirs(os.path.dirname(output_h5ad), exist_ok=True)
    expr_adata.write_h5ad(output_h5ad)
    print(f"Wrote {output_h5ad}")
    print(expr_adata)


if __name__ == "__main__":
    main()
