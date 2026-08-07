#!/usr/bin/env python
# coding: utf-8

"""
Title: Export Dotplot Data From Script 00 Clean Objects
Date: 2026-08-06
Summary: Export long-format dotplot-ready gene expression summaries directly
from clean expression AnnData objects created by script 00. The config defines
one clean object per sample plus multiple `.obs` grouping columns, validates
that the requested groupings are present, and writes combined or split CSV
outputs for downstream multi-sample dotplots. When QC configs or external
cell-by-cluster CSVs are available, the exporter can attach readable labels or
cluster columns for the requested resolutions while preserving clean expression
values.
"""

import argparse
from datetime import datetime
import json
import os
import shutil


# Import the scientific stack inside main() so `--help` works even when the
# active shell is not currently in the project conda environment.


def parse_args():
    """Parse command-line arguments for the dotplot summary exporter."""
    parser = argparse.ArgumentParser(
        prog="export dotplot data from script 00 clean Xenium AnnData objects"
    )
    parser.add_argument(
        "--config",
        type=str,
        help="JSON config defining clean AnnData objects, groupings, markers, and output CSV.",
        required=True,
    )
    return parser.parse_args()


def load_config(config_path):
    """Load a JSON dotplot export config.

    Args:
        config_path: Path to the JSON config file.

    Returns:
        Dictionary containing marker settings, clean AnnData objects, grouping
        columns, and output path.
    """
    with open(config_path) as f:
        return json.load(f)


def resolve_qc_config_path(cfg, sample_cfg):
    """Return the QC config path for a sample when annotations are enabled.

    Args:
        cfg: Dotplot export config dictionary.
        sample_cfg: One sample entry from `cfg["samples"]`.

    Returns:
        Path to a QC config, or `None` when no candidate exists.
    """
    if sample_cfg.get("qc_config") is not None:
        return sample_cfg["qc_config"]

    if cfg.get("use_qc_annotations", True) is False:
        return None

    qc_config_dir = sample_cfg.get("qc_config_dir", cfg.get("qc_config_dir"))
    project = sample_cfg.get("project", cfg.get("project"))
    if not qc_config_dir and project:
        qc_config_dir = os.path.join("config", "01_QC", project)

    if not qc_config_dir:
        return None

    sample = sample_cfg["sample"]
    candidate = os.path.join(qc_config_dir, f"{sample}.json")
    if os.path.exists(candidate):
        return candidate

    return None


def load_qc_annotation(qc_config_path):
    """Load `cluster_col` and `new_labels` from a script 01 QC config.

    Args:
        qc_config_path: Path to a sample QC JSON config.

    Returns:
        Dictionary with `cluster_col`, `new_labels`, and `qc_config_path`, or
        `None` when the config does not define usable annotations.
    """
    if not qc_config_path:
        return None

    with open(qc_config_path) as f:
        qc_cfg = json.load(f)

    cluster_col = qc_cfg.get("cluster_col")
    new_labels = qc_cfg.get("new_labels") or {}
    if not cluster_col or not new_labels:
        return None

    return {
        "cluster_col": cluster_col,
        "new_labels": {str(key): str(value) for key, value in new_labels.items()},
        "qc_config_path": qc_config_path,
    }


def to_dense(x):
    """Return `x` as a dense NumPy array."""
    import numpy as np
    from scipy import sparse

    if sparse.issparse(x):
        return x.toarray()
    return np.asarray(x)


def load_marker_settings(cfg):
    """Resolve marker genes and marker groups from the export config.

    Args:
        cfg: Dotplot export config dictionary.

    Returns:
        Tuple of `(use_all_genes, marker_genes, marker_groups)`.

    Raises:
        ValueError: If marker mode is requested without `marker_file`.
    """
    import pandas as pd

    use_all_genes = cfg.get("use_all_genes", False)
    marker_file = cfg.get("marker_file")
    gene_column = cfg.get("gene_column", "Gene")
    marker_group_column = cfg.get("marker_group_column")

    if use_all_genes:
        return use_all_genes, None, {}

    if marker_file is None:
        raise ValueError("marker_file is required unless use_all_genes is true")

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

    return use_all_genes, marker_genes, marker_groups


def select_expression_source(adata, expression_source):
    """Select the configured expression matrix and matching gene index.

    Args:
        adata: Clean expression AnnData object.
        expression_source: Requested source, usually `"raw"` or `"X"`.

    Returns:
        Tuple of `(expr_adata, var_names, source_used)`.
    """
    import pandas as pd

    if expression_source == "raw" and adata.raw is not None:
        return adata.raw, pd.Index(adata.raw.var_names), "raw"

    return adata, pd.Index(adata.var_names), "X"


def resolve_present_genes(sample, resolution, var_names, marker_genes, use_all_genes):
    """Choose genes that are present in the selected expression source.

    Args:
        sample: Sample label used for status messages.
        resolution: Resolution label used for status messages.
        var_names: Gene names available in the selected expression source.
        marker_genes: Requested marker genes, or `None` in all-gene mode.
        use_all_genes: Whether to export all genes instead of marker genes.

    Returns:
        List of genes to summarize.
    """
    if use_all_genes:
        present_genes = var_names.astype(str).tolist()
        print(f"{sample} r{resolution}: exporting all {len(present_genes)} genes")
        return present_genes

    present_genes = [gene for gene in marker_genes if gene in var_names]
    missing_genes = sorted(set(marker_genes) - set(present_genes))

    if missing_genes:
        print(f"{sample} r{resolution}: missing {len(missing_genes)} marker genes")

    if not present_genes:
        print(f"{sample} r{resolution}: no marker genes found, skipping")

    return present_genes


def expand_sample_groupings(cfg):
    """Expand compact sample configs into one summary task per grouping column.

    Args:
        cfg: Dotplot export config dictionary containing `samples` entries.

    Returns:
        List of task dictionaries with sample, AnnData path, resolution, and
        groupby metadata.

    Raises:
        ValueError: If neither `samples` nor legacy `objects` are configured.
    """
    if "samples" not in cfg:
        if "objects" in cfg:
            return cfg["objects"]
        raise ValueError("Config must contain `samples` or legacy `objects`")

    tasks = []
    for sample_cfg in cfg["samples"]:
        sample = sample_cfg["sample"]
        adata_path = sample_cfg["adata_path"]
        default_groupby_label = sample_cfg.get("groupby_label")
        default_missing_group_label = sample_cfg.get("missing_group_label")
        qc_annotation = load_qc_annotation(resolve_qc_config_path(cfg, sample_cfg))

        for grouping in sample_cfg["groupbys"]:
            groupby = grouping["groupby"]
            grouping_qc_annotation = qc_annotation
            if grouping.get("qc_config") is not None:
                grouping_qc_annotation = load_qc_annotation(grouping["qc_config"])

            task = {
                "sample": sample,
                "adata_path": adata_path,
                "resolution": str(grouping.get("resolution", "")),
                "groupby": groupby,
                "groupby_label": grouping.get(
                    "groupby_label",
                    default_groupby_label or groupby,
                ),
                "label_csvs": grouping.get(
                    "label_csvs",
                    sample_cfg.get("label_csvs", []),
                ),
                "label_csv_cell_id_col": grouping.get(
                    "label_csv_cell_id_col",
                    sample_cfg.get("label_csv_cell_id_col", "index"),
                ),
            }
            if (
                grouping_qc_annotation
                and grouping_qc_annotation["cluster_col"] == groupby
            ):
                task["annotation_labels"] = grouping_qc_annotation["new_labels"]
                task["annotation_source"] = grouping_qc_annotation["qc_config_path"]
            if default_missing_group_label is not None:
                task["missing_group_label"] = default_missing_group_label
            if "missing_group_label" in grouping:
                task["missing_group_label"] = grouping["missing_group_label"]
            tasks.append(task)

    return tasks


def normalize_label_csv_configs(label_csvs, default_cell_id_col):
    """Return normalized label CSV config dictionaries.

    Args:
        label_csvs: List of paths or dictionaries from a sample/grouping config.
        default_cell_id_col: Default column containing cell IDs.

    Returns:
        List of dictionaries with `path` and `cell_id_col`.
    """
    if isinstance(label_csvs, (str, os.PathLike)):
        label_csvs = [label_csvs]

    normalized = []
    for item in label_csvs or []:
        if isinstance(item, dict):
            normalized.append({
                "path": item["path"],
                "cell_id_col": item.get("cell_id_col", default_cell_id_col),
            })
        else:
            normalized.append({
                "path": str(item),
                "cell_id_col": default_cell_id_col,
            })

    return normalized


def drop_csv_index_columns(data):
    """Drop accidental row-number columns from a CSV table."""
    index_like_cols = [
        col for col in data.columns
        if str(col) in {"", "...1", "X", "X1", "Unnamed: 0"}
    ]
    if index_like_cols:
        data = data.drop(columns=index_like_cols)
    return data


def attach_external_label_csvs(adata, adata_path, tasks):
    """Attach configured cell-level cluster-label CSV columns to `.obs`.

    Args:
        adata: Clean expression AnnData receiving extra obs columns.
        adata_path: Path to the clean AnnData object, for error messages.
        tasks: Summary tasks that use this object.

    Returns:
        The same AnnData object with external label columns added to `.obs`.
    """
    import pandas as pd

    csv_configs = []
    for task in tasks:
        csv_configs.extend(
            normalize_label_csv_configs(
                task.get("label_csvs", []),
                task.get("label_csv_cell_id_col", "index"),
            )
        )

    seen = set()
    for csv_cfg in csv_configs:
        key = (csv_cfg["path"], csv_cfg["cell_id_col"])
        if key in seen:
            continue
        seen.add(key)

        label_csv = csv_cfg["path"]
        cell_id_col = csv_cfg["cell_id_col"]
        print(f"Reading external cluster labels: {label_csv}")
        labels = pd.read_csv(label_csv, low_memory=False)
        labels = drop_csv_index_columns(labels)

        if cell_id_col not in labels.columns:
            raise KeyError(
                f"Cell ID column {cell_id_col!r} was not found in {label_csv}. "
                f"Available columns: {list(labels.columns)}"
            )
        if labels[cell_id_col].duplicated().any():
            raise ValueError(f"Duplicate cell IDs found in {label_csv}")

        labels = labels.set_index(cell_id_col)
        label_columns = [col for col in labels.columns if col not in adata.obs.columns]
        if not label_columns:
            continue

        aligned = labels[label_columns].reindex(adata.obs_names)
        missing_counts = aligned.isna().sum()
        missing_counts = missing_counts[missing_counts > 0]
        if not missing_counts.empty:
            raise ValueError(
                f"External label CSV {label_csv} does not cover all cells in "
                f"{adata_path}: {missing_counts.to_dict()}"
            )

        for column in label_columns:
            adata.obs[column] = aligned[column].astype("category")
            print(f"Attached external label column: {column}")

    return adata


def validate_clean_adata(adata, adata_path, tasks):
    """Validate one script 00 clean object before summarizing it.

    Args:
        adata: Clean expression AnnData read from `adata_path`.
        adata_path: Path to the clean AnnData object.
        tasks: Summary tasks that use this object.

    Raises:
        ValueError: If the object is empty or a grouping column has no labels.
        KeyError: If required `.obs` grouping columns are missing.
    """
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError(
            f"{adata_path} is empty or has no genes: shape={adata.shape}"
        )

    required_groupbys = sorted({task["groupby"] for task in tasks})
    missing_groupbys = [
        groupby for groupby in required_groupbys if groupby not in adata.obs.columns
    ]
    if missing_groupbys:
        raise KeyError(
            f"{adata_path} is missing required obs columns: {missing_groupbys}. "
            f"Available obs columns: {list(adata.obs.columns)}"
        )

    empty_groupbys = [
        groupby for groupby in required_groupbys if adata.obs[groupby].notna().sum() == 0
    ]
    if empty_groupbys:
        raise ValueError(
            f"{adata_path} has grouping columns with no non-null labels: "
            f"{empty_groupbys}"
        )


def load_clean_adatas(tasks):
    """Read and validate each unique clean AnnData object once.

    Args:
        tasks: Expanded grouping tasks.

    Returns:
        Dictionary mapping AnnData paths to loaded AnnData objects.
    """
    import anndata as ad

    tasks_by_path = {}
    for task in tasks:
        tasks_by_path.setdefault(task["adata_path"], []).append(task)

    adatas = {}
    for adata_path, path_tasks in tasks_by_path.items():
        print(f"Reading clean script 00 AnnData: {adata_path}")
        adata = ad.read_h5ad(adata_path)
        attach_external_label_csvs(adata, adata_path, path_tasks)
        validate_clean_adata(adata, adata_path, path_tasks)
        adatas[adata_path] = adata

    return adatas


def summarize_grouping(
    obj,
    adata,
    expression_source,
    marker_genes,
    marker_groups,
    use_all_genes,
    default_missing_group_label,
):
    """Summarize one grouping column from a clean script 00 AnnData object.

    Args:
        obj: Expanded grouping task with sample, resolution, path, and groupby.
        adata: Clean expression AnnData containing expression and cluster labels.
        expression_source: Requested expression source, usually `"raw"`.
        marker_genes: Marker genes to summarize, or `None` in all-gene mode.
        marker_groups: Mapping from gene to marker group label.
        use_all_genes: Whether to summarize every available gene.
        default_missing_group_label: Label used when cells have missing group
            metadata.

    Returns:
        List of output rows, where each row represents one dotplot dot.
    """
    import numpy as np

    sample = obj["sample"]
    resolution = str(obj.get("resolution", ""))
    adata_path = obj["adata_path"]
    groupby = obj["groupby"]
    groupby_label = obj.get("groupby_label", groupby)
    annotation_labels = obj.get("annotation_labels", {})
    annotation_source = obj.get("annotation_source", "")
    missing_group_label = obj.get(
        "missing_group_label",
        default_missing_group_label,
    )

    expr_adata, var_names, source_used = select_expression_source(
        adata,
        expression_source,
    )
    present_genes = resolve_present_genes(
        sample,
        resolution,
        var_names,
        marker_genes,
        use_all_genes,
    )
    if not present_genes:
        return []

    # Script 00 clean objects retain biological expression values in .X/.raw and
    # carry all requested BANKSY cluster labels in .obs. Here we summarize only
    # those clean values, never the BANKSY-expanded spatial matrices.
    expr = to_dense(expr_adata[:, present_genes].X)
    group_values = adata.obs[groupby].astype("object")
    missing_mask = group_values.isna() | group_values.astype(str).str.strip().isin(
        ["", "nan", "None", "NA", "<NA>"]
    )
    groups = group_values.where(~missing_mask, missing_group_label).astype(str)

    rows = []
    for group_id in sorted(groups.unique()):
        mask = (groups == group_id).to_numpy()
        group_expr = expr[mask, :]
        mean_expression = np.asarray(group_expr.mean(axis=0)).ravel()
        percent_expressing = np.asarray((group_expr > 0).mean(axis=0)).ravel() * 100
        sample_group = f"{sample}__{groupby_label}__{group_id}"
        group_label = annotation_labels.get(group_id, group_id)
        has_annotation = group_id in annotation_labels

        for i, gene in enumerate(present_genes):
            rows.append(
                {
                    "sample": sample,
                    "resolution": resolution,
                    "cluster_id": group_id,
                    "group_id": group_id,
                    "group_label": group_label,
                    "sample_cluster": sample_group,
                    "sample_group": sample_group,
                    "groupby": groupby,
                    "groupby_label": groupby_label,
                    "has_group_annotation": has_annotation,
                    "annotation_source": annotation_source,
                    "gene": gene,
                    "marker_group": marker_groups.get(
                        gene,
                        "all_genes" if use_all_genes else "marker",
                    ),
                    "mean_expression": mean_expression[i],
                    "percent_expressing": percent_expressing[i],
                    "n_cells": int(mask.sum()),
                    "expression_source": source_used,
                    "adata_path": adata_path,
                }
            )

    return rows


def sanitize_output_label(value):
    """Return a filename-safe label for a split output value."""
    return str(value).replace(".", "p").replace("/", "-").replace(" ", "_")


def default_split_output_path(output_csv, split_by, split_value):
    """Build a split-output path when no explicit template is configured.

    Args:
        output_csv: Base CSV path from the config.
        split_by: Column used to split the output.
        split_value: Value of `split_by` for this output table.

    Returns:
        CSV path with a split label inserted before `_dotplot_summary` when
        possible, otherwise before the `.csv` suffix.
    """
    label = sanitize_output_label(split_value)
    suffix = f"_{split_by}_{label}"

    if output_csv.endswith("_dotplot_summary.csv"):
        return output_csv.replace("_dotplot_summary.csv", f"{suffix}_dotplot_summary.csv")

    root, ext = os.path.splitext(output_csv)
    return f"{root}{suffix}{ext or '.csv'}"


def resolve_output_path(cfg, split_by=None, split_value=None):
    """Resolve the CSV output path for combined or split export modes.

    Args:
        cfg: Dotplot export config dictionary.
        split_by: Optional column used to split outputs.
        split_value: Optional value for the current split.

    Returns:
        Path to the CSV that should be written.
    """
    output_csv = cfg["output_csv"]
    if split_by is None:
        return output_csv

    template = cfg.get("output_csv_template")
    if template:
        return template.format(
            split_by=split_by,
            split_value=str(split_value),
            split_label=sanitize_output_label(split_value),
            resolution=str(split_value),
            resolution_label=sanitize_output_label(split_value),
        )

    return default_split_output_path(output_csv, split_by, split_value)


def archive_existing_output(output_path, cfg):
    """Archive an existing output file before overwriting it.

    Args:
        output_path: Destination path that is about to be written.
        cfg: Dotplot export config dictionary.

    Returns:
        Path to the archived file, or `None` if no archive was created.
    """
    if not cfg.get("archive_previous_outputs", True):
        return None

    if not os.path.exists(output_path):
        return None

    archive_dir = cfg.get(
        "archive_output_dir",
        os.path.join(os.path.dirname(output_path), "archive"),
    )
    os.makedirs(archive_dir, exist_ok=True)

    root, ext = os.path.splitext(os.path.basename(output_path))
    timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    archived_path = os.path.join(archive_dir, f"{root}__{timestamp}{ext}")
    shutil.move(output_path, archived_path)
    print(f"Archived previous output to {archived_path}")
    return archived_path


def write_dotplot_table(out, output_csv, cfg):
    """Write one dotplot summary table to CSV and return the path."""
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    archive_existing_output(output_csv, cfg)
    out.to_csv(output_csv, index=False)
    print(f"Wrote {output_csv}")
    return output_csv


def write_dotplot_outputs(out, cfg):
    """Write combined or split dotplot output tables.

    Args:
        out: Long-format dotplot summary DataFrame.
        cfg: Dotplot export config dictionary.

    Returns:
        List of written CSV paths.
    """
    split_by = cfg.get("split_by")
    if not split_by:
        return [write_dotplot_table(out, resolve_output_path(cfg), cfg)]

    if split_by not in out.columns:
        raise ValueError(
            f"split_by column {split_by!r} is not available. "
            f"Available columns: {list(out.columns)}"
        )

    written_paths = []
    for split_value, split_df in out.groupby(split_by, sort=True):
        output_csv = resolve_output_path(cfg, split_by=split_by, split_value=split_value)
        written_paths.append(write_dotplot_table(split_df, output_csv, cfg))

    return written_paths


def export_dotplot_summary(cfg):
    """Export long-format dotplot summary CSVs from script 00 clean objects.

    Args:
        cfg: Dotplot export config dictionary.

    Returns:
        List of written CSV paths.
    """
    import pandas as pd

    expression_source = cfg.get("expression_source", "raw")
    default_missing_group_label = cfg.get(
        "missing_group_label",
        "missing_group",
    )
    tasks = expand_sample_groupings(cfg)
    adatas = load_clean_adatas(tasks)
    use_all_genes, marker_genes, marker_groups = load_marker_settings(cfg)

    rows = []
    for task in tasks:
        rows.extend(
            summarize_grouping(
                task,
                adatas[task["adata_path"]],
                expression_source,
                marker_genes,
                marker_groups,
                use_all_genes,
                default_missing_group_label,
            )
        )

    out = pd.DataFrame(rows)
    return write_dotplot_outputs(out, cfg)


def main():
    """Run the script 00 clean-object dotplot export workflow from JSON."""
    args = parse_args()
    cfg = load_config(args.config)
    export_dotplot_summary(cfg)


if __name__ == "__main__":
    main()
