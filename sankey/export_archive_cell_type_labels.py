#!/usr/bin/env python
# coding: utf-8

"""
Title: Export Archive Cell Type Labels
Date: 2026-06-04
Summary: Read archived BANKSY spatial AnnData objects from a JSON config and
export a slim cell-level table containing cell IDs, sample metadata, cluster
IDs, and readable cell-type labels.
"""

import argparse
import json
import os
import re

import pandas as pd


SPATIAL_FILENAME_RE = re.compile(
    r"(?P<sample>.+)_clustered_spatial_pc(?P<pc_label>[^_]+)_"
    r"nc(?P<lambda_label>[^_]+)_r(?P<resolution>[^.]+(?:\.\d+)?)\.h5ad$"
)


def parse_args():
    """Parse command-line arguments for the archive cell-type label exporter."""
    parser = argparse.ArgumentParser(
        prog="export archived cell IDs and cell-type labels from AnnData objects"
    )
    parser.add_argument(
        "--config",
        required=True,
        help=(
            "JSON config with an 'objects' list. Each object must include "
            "'adata_path' and may include sample, resolution, pc_label, "
            "lambda_label, cluster_column, and annotation_column."
        ),
    )
    parser.add_argument(
        "--output-csv",
        default=None,
        help=(
            "Output CSV path. Defaults to config key 'cell_type_label_output_csv' "
            "or the config output_csv path with '_cell_type_labels' appended."
        ),
    )
    return parser.parse_args()


def load_config(config_path):
    """Load a JSON label-export config.

    Args:
        config_path: Path to the JSON config.

    Returns:
        Config dictionary.
    """
    with open(config_path) as f:
        return json.load(f)


def parse_spatial_filename(adata_path):
    """Infer sample and clustering metadata from a spatial AnnData filename.

    Args:
        adata_path: Path to a `*_clustered_spatial_pc*_nc*_r*.h5ad` file.

    Returns:
        Dictionary with sample, pc, lambda, and resolution values. Values are
        empty strings when the filename does not match the expected pattern.
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


def resolve_object_metadata(obj):
    """Resolve per-object metadata from explicit config and filename values.

    Args:
        obj: One object entry from the config.

    Returns:
        Dictionary with AnnData path, sample, resolution, pc, and lambda labels.
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


def infer_cluster_column(obs_columns, pc_label, lambda_label, resolution):
    """Infer the archived BANKSY cluster column from `.obs` column names.

    Args:
        obs_columns: Available observation columns.
        pc_label: Principal-component label from config or filename.
        lambda_label: BANKSY lambda label from config or filename.
        resolution: BANKSY resolution label from config or filename.

    Returns:
        Best matching cluster column name, or `None` when no clear match exists.
    """
    preferred = f"banksy_cluster_pc{pc_label}_nc{lambda_label}_r{resolution}"
    if preferred in obs_columns:
        return preferred

    banksy_candidates = [
        column
        for column in obs_columns
        if column.startswith("banksy_cluster_pc")
        and not column.endswith("_raw")
        and not column.endswith("_ann")
    ]
    if len(banksy_candidates) == 1:
        return banksy_candidates[0]

    label_candidates = [
        column for column in obs_columns if column.startswith("labels_scaled_gaussian_pc")
    ]
    if len(label_candidates) == 1:
        return label_candidates[0]

    return None


def infer_annotation_column(obs_columns, cluster_column):
    """Infer the readable annotation column paired with a cluster column.

    Args:
        obs_columns: Available observation columns.
        cluster_column: Numeric cluster column.

    Returns:
        Annotation column name, or `None` if no label column can be inferred.
    """
    if cluster_column:
        annotation_column = f"{cluster_column}_ann"
        if annotation_column in obs_columns:
            return annotation_column

    if "cell type" in obs_columns:
        return "cell type"

    return None


def export_one_object(obj):
    """Export cell IDs and readable labels from one archived AnnData object.

    Args:
        obj: One object entry from the config.

    Returns:
        DataFrame with one row per cell.

    Raises:
        ValueError: If no cluster or annotation column can be resolved.
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
            f"Could not resolve cluster column for {adata_path}. "
            f"Available obs columns: {obs_columns}"
        )

    annotation_column = (
        obj.get("annotation_column")
        or obj.get("cell_type_column")
        or infer_annotation_column(obs_columns, cluster_column)
    )
    if not annotation_column or annotation_column not in obs_columns:
        adata.file.close()
        raise ValueError(
            f"Could not resolve cell-type label column for {adata_path}. "
            f"Available obs columns: {obs_columns}"
        )

    obs = adata.obs
    # Align labels by the AnnData observation index; this is the stable cell ID
    # used later when joining old annotations to current clean expression objects.
    out = pd.DataFrame(
        {
            "cell_id": obs.index.astype(str),
            "sample": metadata["sample"],
            "resolution": metadata["resolution"],
            "pc_label": metadata["pc_label"],
            "lambda_label": metadata["lambda_label"],
            "cluster_id": obs[cluster_column].astype(str).to_numpy(),
            "cell_type_label": obs[annotation_column].astype(str).to_numpy(),
            "cluster_column": cluster_column,
            "cell_type_column": annotation_column,
            "adata_path": adata_path,
        }
    )
    adata.file.close()
    print(f"Exported {len(out)} cells from {metadata['sample']}")
    return out


def resolve_output_csv(cfg, output_csv_arg):
    """Resolve the output path from CLI, config, or the main export path.

    Args:
        cfg: Label-export config dictionary.
        output_csv_arg: Optional output path from the command line.

    Returns:
        Output CSV path.
    """
    if output_csv_arg:
        return output_csv_arg

    if cfg.get("cell_type_label_output_csv"):
        return cfg["cell_type_label_output_csv"]

    if cfg.get("output_csv"):
        stem, ext = os.path.splitext(cfg["output_csv"])
        return f"{stem}_cell_type_labels{ext or '.csv'}"

    return "archive_spatial_cell_type_labels.csv"


def write_output(combined, output_csv):
    """Write the combined cell-type label table.

    Args:
        combined: Combined cell-level label table.
        output_csv: Destination CSV path.
    """
    output_dir = os.path.dirname(output_csv)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    combined.to_csv(output_csv, index=False)
    print(f"Wrote {output_csv}")


def main():
    """Run the archive cell-type label export workflow."""
    args = parse_args()
    cfg = load_config(args.config)
    tables = [export_one_object(obj) for obj in cfg["objects"]]
    combined = pd.concat(tables, ignore_index=True)
    output_csv = resolve_output_csv(cfg, args.output_csv)
    write_output(combined, output_csv)


if __name__ == "__main__":
    main()
