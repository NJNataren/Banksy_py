#!/usr/bin/env python
# coding: utf-8

"""
Title: Export Cluster Labels From AnnData Spatial Objects
Date: 2026-06-03
Summary: Read archived or current BANKSY spatial AnnData objects and export
cell-level cluster assignments to CSV. The config defines input objects,
optional cluster/annotation columns, output paths, and optional `.obs` metadata
columns to preserve.
"""

import argparse
import json
import os
import re

import pandas as pd


DEFAULT_OPTIONAL_OBS_COLUMNS = [
    "x",
    "y",
    "cell type",
    "nCount_Xenium",
    "nFeature_Xenium",
]


SPATIAL_FILENAME_RE = re.compile(
    r"(?P<sample>.+)_clustered_spatial_pc(?P<pc_label>[^_]+)_"
    r"nc(?P<lambda_label>[^_]+)_r(?P<resolution>[^.]+(?:\.\d+)?)\.h5ad$"
)


def parse_args():
    """Parse command-line arguments for the cluster-label exporter."""
    parser = argparse.ArgumentParser(
        prog="export cell-level cluster labels from spatial AnnData objects"
    )
    parser.add_argument(
        "--config",
        required=True,
        help="JSON config defining archived AnnData inputs and output CSV paths.",
    )
    return parser.parse_args()


def load_config(config_path):
    """Load the JSON exporter config.

    Args:
        config_path: Path to a JSON config file.

    Returns:
        Config dictionary.
    """
    with open(config_path) as f:
        return json.load(f)


def parse_spatial_filename(adata_path):
    """Infer sample and clustering parameters from an archived spatial filename.

    Args:
        adata_path: Path to a `*_clustered_spatial_pc*_nc*_r*.h5ad` file.

    Returns:
        Dictionary with inferred metadata. Missing values are returned as empty
        strings when the filename does not match the expected archive pattern.
    """
    match = SPATIAL_FILENAME_RE.search(os.path.basename(adata_path))
    if not match:
        return {
            "sample": "",
            "pc_label": "",
            "lambda_label": "",
            "resolution": "",
        }
    return match.groupdict()


def infer_cluster_column(obs_columns, pc_label, lambda_label, resolution):
    """Infer the numeric BANKSY cluster column from `.obs` columns.

    Args:
        obs_columns: Available observation columns.
        pc_label: Principal-component label from config or filename.
        lambda_label: BANKSY lambda label from config or filename.
        resolution: BANKSY resolution label from config or filename.

    Returns:
        Best matching cluster column name, or `None` if no match is found.
    """
    preferred = f"banksy_cluster_pc{pc_label}_nc{lambda_label}_r{resolution}"
    if preferred in obs_columns:
        return preferred

    candidates = [
        column
        for column in obs_columns
        if column.startswith("banksy_cluster_pc")
        and not column.endswith("_raw")
        and not column.endswith("_ann")
    ]
    if len(candidates) == 1:
        return candidates[0]

    label_candidates = [
        column for column in obs_columns if column.startswith("labels_scaled_gaussian_pc")
    ]
    if len(label_candidates) == 1:
        return label_candidates[0]

    return None


def infer_annotation_column(obs_columns, cluster_column):
    """Infer a readable annotation column for a numeric cluster column.

    Args:
        obs_columns: Available observation columns.
        cluster_column: Numeric cluster column.

    Returns:
        Annotation column name, or `None` if unavailable.
    """
    if not cluster_column:
        return None

    annotation_column = f"{cluster_column}_ann"
    if annotation_column in obs_columns:
        return annotation_column

    if "cell type" in obs_columns:
        return "cell type"

    return None


def resolve_object_metadata(obj):
    """Resolve per-object metadata from config and filename values.

    Args:
        obj: Config entry for one AnnData object.

    Returns:
        Dictionary containing sample, resolution, pc, lambda, and path values.
    """
    adata_path = obj["adata_path"]
    inferred = parse_spatial_filename(adata_path)
    return {
        "adata_path": adata_path,
        "sample": str(obj.get("sample") or inferred["sample"]),
        "resolution": str(obj.get("resolution") or inferred["resolution"]),
        "pc_label": str(obj.get("pc_label") or inferred["pc_label"]),
        "lambda_label": str(obj.get("lambda_label") or inferred["lambda_label"]),
    }


def export_one_object(obj, optional_obs_columns):
    """Export cell-level labels from one configured spatial AnnData object.

    Args:
        obj: Config entry for one AnnData object.
        optional_obs_columns: Metadata columns to copy from `.obs` when present.

    Returns:
        DataFrame with one row per cell.

    Raises:
        ValueError: If the cluster column cannot be found or inferred.
    """
    import anndata as ad

    metadata = resolve_object_metadata(obj)
    adata_path = metadata["adata_path"]
    print(f"Reading {adata_path}")
    adata = ad.read_h5ad(adata_path, backed="r")
    obs_columns = list(adata.obs.columns)

    cluster_column = obj.get("cluster_column") or infer_cluster_column(
        obs_columns,
        metadata["pc_label"],
        metadata["lambda_label"],
        metadata["resolution"],
    )
    if not cluster_column or cluster_column not in obs_columns:
        adata.file.close()
        raise ValueError(
            f"Could not resolve cluster column for {adata_path}. Available obs "
            f"columns: {obs_columns}"
        )

    annotation_column = obj.get("annotation_column") or infer_annotation_column(
        obs_columns,
        cluster_column,
    )

    obs = adata.obs
    out = pd.DataFrame(
        {
            "cell_id": obs.index.astype(str),
            "sample": metadata["sample"],
            "resolution": metadata["resolution"],
            "pc_label": metadata["pc_label"],
            "lambda_label": metadata["lambda_label"],
            "cluster_id": obs[cluster_column].astype(str).to_numpy(),
            "cluster_column": cluster_column,
            "annotation_column": annotation_column or "",
            "adata_path": adata_path,
        }
    )

    if annotation_column and annotation_column in obs_columns:
        out["cluster_annotation"] = obs[annotation_column].astype(str).to_numpy()
    else:
        out["cluster_annotation"] = ""

    # Copy optional metadata only when present so one config can handle archived
    # objects with slightly different `.obs` schemas.
    for column in optional_obs_columns:
        if column in obs_columns and column not in out.columns:
            out[column] = obs[column].to_numpy()

    adata.file.close()
    print(f"Exported {len(out)} cells from {metadata['sample']} r{metadata['resolution']}")
    return out


def write_outputs(combined, cfg):
    """Write combined labels and an optional cluster-size summary.

    Args:
        combined: Combined cell-level label table.
        cfg: Exporter config dictionary.
    """
    output_csv = cfg["output_csv"]
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    combined.to_csv(output_csv, index=False)
    print(f"Wrote {output_csv}")

    summary_csv = cfg.get("summary_csv")
    if summary_csv:
        summary = (
            combined.groupby(
                [
                    "sample",
                    "resolution",
                    "pc_label",
                    "lambda_label",
                    "cluster_id",
                    "cluster_annotation",
                ],
                dropna=False,
            )
            .size()
            .reset_index(name="n_cells")
            .sort_values(["sample", "resolution", "cluster_id"])
        )
        os.makedirs(os.path.dirname(summary_csv), exist_ok=True)
        summary.to_csv(summary_csv, index=False)
        print(f"Wrote {summary_csv}")


def main():
    """Run the cluster-label export workflow from a JSON config."""
    args = parse_args()
    cfg = load_config(args.config)
    optional_obs_columns = cfg.get("optional_obs_columns", DEFAULT_OPTIONAL_OBS_COLUMNS)

    tables = [export_one_object(obj, optional_obs_columns) for obj in cfg["objects"]]
    combined = pd.concat(tables, ignore_index=True)
    write_outputs(combined, cfg)


if __name__ == "__main__":
    main()
