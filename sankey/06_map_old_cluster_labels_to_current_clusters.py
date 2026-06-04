#!/usr/bin/env python
# coding: utf-8

"""
Title: Map Old Cluster Labels To Current Clusters
Date: 2026-06-03
Summary: Compare current BANKSY cluster labels in clean expression AnnData
objects with archived cell-level cluster labels exported from old spatial
AnnData objects. Write a per-current-cluster CSV showing the majority old
cluster annotation by cell overlap.
"""

import argparse
import json
import os

import pandas as pd


def parse_args():
    """Parse command-line arguments for the old-to-current cluster mapper."""
    parser = argparse.ArgumentParser(
        prog="map archived cluster labels to current cluster labels by cell overlap"
    )
    parser.add_argument(
        "--config",
        required=True,
        help="JSON config defining old label CSV, current AnnData objects, and outputs.",
    )
    return parser.parse_args()


def load_config(config_path):
    """Load a JSON mapping config.

    Args:
        config_path: Path to the JSON config.

    Returns:
        Config dictionary.
    """
    with open(config_path) as f:
        return json.load(f)


def load_old_labels(old_labels_csv):
    """Load old cell-level labels exported by script 05.

    Args:
        old_labels_csv: CSV with at least sample, cell_id, cluster_id, and
            cluster_annotation columns.

    Returns:
        DataFrame with string-normalized key columns.
    """
    old = pd.read_csv(old_labels_csv, low_memory=False)
    required = {"sample", "cell_id", "cluster_id", "cluster_annotation"}
    missing = required - set(old.columns)
    if missing:
        raise ValueError(f"{old_labels_csv} is missing required columns: {missing}")

    for column in ["sample", "cell_id", "cluster_id", "cluster_annotation"]:
        old[column] = old[column].astype(str)
    if "resolution" in old.columns:
        old["resolution"] = old["resolution"].astype(str)
    else:
        old["resolution"] = ""

    return old


def load_current_labels(obj):
    """Load current cell-level cluster labels from one clean expression AnnData.

    Args:
        obj: Config entry with sample, current AnnData path, cluster column, and
            resolution metadata.

    Returns:
        DataFrame with one row per current cell.

    Raises:
        KeyError: If the configured current cluster column is absent.
    """
    import anndata as ad

    sample = str(obj["sample"])
    current_resolution = str(obj.get("current_resolution", ""))
    current_adata_path = obj["current_adata_path"]
    current_cluster_column = obj["current_cluster_column"]

    print(f"Reading current labels: {current_adata_path}")
    adata = ad.read_h5ad(current_adata_path, backed="r")
    if current_cluster_column not in adata.obs.columns:
        obs_columns = list(adata.obs.columns)
        adata.file.close()
        raise KeyError(
            f"{current_cluster_column!r} was not found in {current_adata_path}. "
            f"Available obs columns: {obs_columns}"
        )

    current = pd.DataFrame(
        {
            "sample": sample,
            "cell_id": adata.obs.index.astype(str),
            "current_resolution": current_resolution,
            "current_cluster_id": adata.obs[current_cluster_column].astype(str).to_numpy(),
            "current_cluster_column": current_cluster_column,
            "current_adata_path": current_adata_path,
        }
    )
    adata.file.close()
    return current


def summarize_current_cluster(
    joined,
    current_group,
    old_sample,
    min_overlap_fraction,
    mapping_denominator,
    old_labels_are_roi_subset,
):
    """Summarize old-label overlap for one current cluster.

    Args:
        joined: Current-label table joined to old labels for one sample.
        current_group: Tuple identifying the current cluster group.
        old_sample: Sample name used to select archived labels.
        min_overlap_fraction: Fraction below which the majority old label is
            flagged as mixed.
        mapping_denominator: Denominator for majority/mixed decisions. Use
            `"current_cluster"` for whole-sample old labels or `"matched_cells"`
            when old labels cover only an ROI subset.
        old_labels_are_roi_subset: Whether old labels are known to represent an
            ROI/subset of the current object.

    Returns:
        Dictionary describing the majority old label and overlap quality.
    """
    (
        sample,
        current_resolution,
        current_cluster_id,
        current_cluster_column,
        current_adata_path,
    ) = current_group

    new_cluster_n_cells = len(joined)
    matched = joined[joined["old_cluster_id"].notna()].copy()
    matched_n_cells = len(matched)

    base = {
        "sample": sample,
        "old_sample": old_sample,
        "current_resolution": current_resolution,
        "current_cluster_id": current_cluster_id,
        "current_cluster_column": current_cluster_column,
        "current_adata_path": current_adata_path,
        "new_cluster_n_cells": new_cluster_n_cells,
        "matched_old_n_cells": matched_n_cells,
        "matched_old_fraction": matched_n_cells / new_cluster_n_cells
        if new_cluster_n_cells
        else 0.0,
    }

    if matched.empty:
        base.update(
            {
                "old_resolution": "",
                "old_cluster_id_majority": "",
                "old_cluster_annotation_majority": "no old match",
                "majority_old_n_cells": 0,
                "majority_old_fraction_of_new_cluster": 0.0,
                "majority_old_fraction_of_matched_cells": 0.0,
                "mapping_denominator": mapping_denominator,
                "old_labels_are_roi_subset": old_labels_are_roi_subset,
                "n_old_clusters_overlapping": 0,
                "mapping_status": "no_old_match",
                "old_cluster_mix": "",
            }
        )
        return base

    counts = (
        matched.groupby(
            ["old_resolution", "old_cluster_id", "old_cluster_annotation"],
            dropna=False,
        )
        .size()
        .reset_index(name="n_cells")
        .sort_values("n_cells", ascending=False)
    )
    top = counts.iloc[0]
    majority_fraction_new_cluster = top["n_cells"] / new_cluster_n_cells
    majority_fraction_matched_cells = top["n_cells"] / matched_n_cells
    if mapping_denominator == "matched_cells":
        decision_fraction = majority_fraction_matched_cells
    elif mapping_denominator == "current_cluster":
        decision_fraction = majority_fraction_new_cluster
    else:
        raise ValueError(
            "mapping_denominator must be 'current_cluster' or 'matched_cells'; "
            f"got {mapping_denominator!r}"
        )

    status_prefix = "roi_subset_" if old_labels_are_roi_subset else ""
    status = f"{status_prefix}majority_old_label"
    if decision_fraction < min_overlap_fraction:
        status = f"{status_prefix}mixed_old_labels"

    old_cluster_mix = "; ".join(
        f"{row.old_cluster_annotation} ({row.n_cells})"
        for row in counts.head(5).itertuples(index=False)
    )

    base.update(
        {
            "old_resolution": str(top["old_resolution"]),
            "old_cluster_id_majority": str(top["old_cluster_id"]),
            "old_cluster_annotation_majority": str(top["old_cluster_annotation"]),
            "majority_old_n_cells": int(top["n_cells"]),
            "majority_old_fraction_of_new_cluster": majority_fraction_new_cluster,
            "majority_old_fraction_of_matched_cells": majority_fraction_matched_cells,
            "mapping_denominator": mapping_denominator,
            "old_labels_are_roi_subset": old_labels_are_roi_subset,
            "n_old_clusters_overlapping": len(counts),
            "mapping_status": status,
            "old_cluster_mix": old_cluster_mix,
        }
    )
    return base


def map_one_object(obj, old_labels, min_overlap_fraction, default_mapping_denominator):
    """Map old labels onto current clusters for one sample.

    Args:
        obj: Current-object config entry.
        old_labels: Combined old labels from script 05.
        min_overlap_fraction: Threshold for majority versus mixed labels.
        default_mapping_denominator: Global denominator used when an object does
            not define `mapping_denominator`.

    Returns:
        Tuple of `(cluster_mapping, cell_mapping)`, where the first table has
        one row per current cluster and the second has one row per current cell.
    """
    current = load_current_labels(obj)
    sample = str(obj["sample"])
    old_sample = str(obj.get("old_sample", sample))
    old_labels_are_roi_subset = bool(obj.get("old_labels_are_roi_subset", False))
    mapping_denominator = obj.get(
        "mapping_denominator",
        "matched_cells" if old_labels_are_roi_subset else default_mapping_denominator,
    )
    old_sample_df = old_labels[old_labels["sample"] == old_sample].copy()
    old_sample_df["old_sample"] = old_sample

    old_sample_df = old_sample_df.rename(
        columns={
            "resolution": "old_resolution",
            "cluster_id": "old_cluster_id",
            "cluster_annotation": "old_cluster_annotation",
        }
    )
    old_columns = [
        "cell_id",
        "old_sample",
        "old_resolution",
        "old_cluster_id",
        "old_cluster_annotation",
    ]
    joined = current.merge(old_sample_df[old_columns], on="cell_id", how="left")

    group_columns = [
        "sample",
        "current_resolution",
        "current_cluster_id",
        "current_cluster_column",
        "current_adata_path",
    ]
    rows = []
    for current_group, group_df in joined.groupby(group_columns, sort=True):
        rows.append(
            summarize_current_cluster(
                group_df,
                current_group,
                old_sample,
                min_overlap_fraction,
                mapping_denominator,
                old_labels_are_roi_subset,
            )
        )

    mapped = pd.DataFrame(rows)
    joined["has_old_label"] = joined["old_cluster_id"].notna()
    print(
        f"Mapped {len(mapped)} current clusters for {sample} using old sample {old_sample}; "
        f"{int(mapped['matched_old_n_cells'].sum())} cells overlapped old labels"
    )
    return mapped, joined


def add_best_cluster_labels_to_cells(mapping, cell_mapping):
    """Add majority old-cluster labels for each current cluster to cell rows.

    Args:
        mapping: Per-current-cluster mapping table.
        cell_mapping: Per-current-cell mapping table.

    Returns:
        Cell-level DataFrame with both each cell's old label and the best-fit
        old label assigned to its current cluster.
    """
    join_columns = [
        "sample",
        "current_resolution",
        "current_cluster_id",
        "current_cluster_column",
        "current_adata_path",
    ]
    best_label_columns = join_columns + [
        "old_cluster_id_majority",
        "old_cluster_annotation_majority",
        "majority_old_n_cells",
        "majority_old_fraction_of_new_cluster",
        "majority_old_fraction_of_matched_cells",
        "matched_old_fraction",
        "mapping_denominator",
        "old_labels_are_roi_subset",
        "mapping_status",
        "old_cluster_mix",
    ]
    best_labels = mapping[best_label_columns].rename(
        columns={
            "old_cluster_id_majority": "best_fit_old_cluster_id",
            "old_cluster_annotation_majority": "best_fit_old_cluster_annotation",
            "majority_old_n_cells": "best_fit_old_n_cells",
            "majority_old_fraction_of_new_cluster": (
                "best_fit_old_fraction_of_new_cluster"
            ),
            "majority_old_fraction_of_matched_cells": (
                "best_fit_old_fraction_of_matched_cells"
            ),
            "matched_old_fraction": "current_cluster_old_label_coverage",
            "mapping_status": "best_fit_mapping_status",
            "old_cluster_mix": "best_fit_old_cluster_mix",
        }
    )
    return cell_mapping.merge(best_labels, on=join_columns, how="left")


def build_dotplot_cluster_label_guide(mapping):
    """Build a compact cluster-level label guide for dotplot y-axis annotation.

    Args:
        mapping: Per-current-cluster mapping table.

    Returns:
        DataFrame with one row per current dotplot cluster and reviewer-friendly
        best-fit archived label columns.
    """
    guide = mapping.copy()
    guide["sample_cluster"] = (
        guide["sample"].astype(str)
        + "__r"
        + guide["current_resolution"].astype(str)
        + "__cluster_"
        + guide["current_cluster_id"].astype(str)
    )
    guide["dotplot_old_label"] = guide["old_cluster_annotation_majority"].astype(str)
    guide.loc[
        guide["old_cluster_annotation_majority"].astype(str).eq(""),
        "dotplot_old_label",
    ] = "no old match"

    # Keep both fractions: whole-cluster coverage matters for full old objects,
    # while matched-cell fraction is the fairer confidence measure for ROI labels.
    output_columns = [
        "sample_cluster",
        "sample",
        "current_resolution",
        "current_cluster_id",
        "old_sample",
        "old_resolution",
        "old_cluster_id_majority",
        "old_cluster_annotation_majority",
        "dotplot_old_label",
        "new_cluster_n_cells",
        "matched_old_n_cells",
        "matched_old_fraction",
        "majority_old_n_cells",
        "majority_old_fraction_of_new_cluster",
        "majority_old_fraction_of_matched_cells",
        "mapping_denominator",
        "old_labels_are_roi_subset",
        "n_old_clusters_overlapping",
        "mapping_status",
        "old_cluster_mix",
        "current_cluster_column",
        "current_adata_path",
    ]
    guide = guide[output_columns].rename(
        columns={
            "current_resolution": "resolution",
            "current_cluster_id": "cluster_id",
            "old_cluster_id_majority": "best_fit_old_cluster_id",
            "old_cluster_annotation_majority": "best_fit_old_cluster_annotation",
            "matched_old_fraction": "current_cluster_old_label_coverage",
            "majority_old_fraction_of_new_cluster": (
                "best_fit_old_fraction_of_new_cluster"
            ),
            "majority_old_fraction_of_matched_cells": (
                "best_fit_old_fraction_of_matched_cells"
            ),
            "majority_old_n_cells": "best_fit_old_n_cells",
            "old_cluster_mix": "best_fit_old_cluster_mix",
        }
    )
    return guide.sort_values(["sample", "resolution", "cluster_id"])


def write_outputs(mapping, cell_mapping, cfg):
    """Write cluster-level and optional cell-level mapping outputs.

    Args:
        mapping: Combined per-current-cluster mapping table.
        cell_mapping: Combined per-current-cell mapping table.
        cfg: Mapping config dictionary.
    """
    output_csv = cfg["output_csv"]
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    mapping.to_csv(output_csv, index=False)
    print(f"Wrote {output_csv}")

    summary_csv = cfg.get("summary_csv")
    if summary_csv:
        summary = (
            mapping.groupby(["sample", "mapping_status"], dropna=False)
            .size()
            .reset_index(name="n_current_clusters")
            .sort_values(["sample", "mapping_status"])
        )
        os.makedirs(os.path.dirname(summary_csv), exist_ok=True)
        summary.to_csv(summary_csv, index=False)
        print(f"Wrote {summary_csv}")

    dotplot_label_output_csv = cfg.get("dotplot_label_output_csv")
    if dotplot_label_output_csv:
        guide = build_dotplot_cluster_label_guide(mapping)
        os.makedirs(os.path.dirname(dotplot_label_output_csv), exist_ok=True)
        guide.to_csv(dotplot_label_output_csv, index=False)
        print(f"Wrote {dotplot_label_output_csv}")

    cell_output_csv = cfg.get("cell_output_csv")
    if cell_output_csv:
        os.makedirs(os.path.dirname(cell_output_csv), exist_ok=True)
        cell_mapping.to_csv(cell_output_csv, index=False)
        print(f"Wrote {cell_output_csv}")

    cell_best_label_output_csv = cfg.get("cell_best_label_output_csv")
    if cell_best_label_output_csv:
        cell_best_labels = add_best_cluster_labels_to_cells(mapping, cell_mapping)
        os.makedirs(os.path.dirname(cell_best_label_output_csv), exist_ok=True)
        cell_best_labels.to_csv(cell_best_label_output_csv, index=False)
        print(f"Wrote {cell_best_label_output_csv}")


def main():
    """Run old-to-current cluster label mapping from a JSON config."""
    args = parse_args()
    cfg = load_config(args.config)
    old_labels = load_old_labels(cfg["old_cluster_labels_csv"])
    min_overlap_fraction = float(cfg.get("min_overlap_fraction", 0.5))
    default_mapping_denominator = cfg.get("mapping_denominator", "current_cluster")

    mapped_outputs = [
        map_one_object(
            obj,
            old_labels,
            min_overlap_fraction,
            default_mapping_denominator,
        )
        for obj in cfg["objects"]
    ]
    cluster_mappings = [item[0] for item in mapped_outputs]
    cell_mappings = [item[1] for item in mapped_outputs]
    combined = pd.concat(cluster_mappings, ignore_index=True)
    combined_cells = pd.concat(cell_mappings, ignore_index=True)
    write_outputs(combined, combined_cells, cfg)


if __name__ == "__main__":
    main()
