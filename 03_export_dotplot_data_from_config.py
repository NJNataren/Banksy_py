#!/usr/bin/env python
# coding: utf-8

"""
Title: Export Dotplot Data From Config
Date: 2026-06-09
Summary: Export long-format dotplot-ready gene expression summaries from clean
Xenium expression AnnData objects. The config defines marker selection, the
expression source, grouping columns, and either one combined CSV output or
split CSV outputs by fields such as resolution. Missing group labels can be
kept as an explicit group for archive-label workflows.
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
        prog="export dotplot data from configured Xenium AnnData objects"
    )
    parser.add_argument(
        "--config",
        type=str,
        help="JSON config defining AnnData objects, markers, and output CSV.",
        required=True,
    )
    return parser.parse_args()


def load_config(config_path):
    """Load a JSON dotplot export config.

    Args:
        config_path: Path to the JSON config file.

    Returns:
        Dictionary containing marker settings, AnnData objects, and output path.
    """
    with open(config_path) as f:
        return json.load(f)


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

    # Marker mode keeps exports focused and can preserve broad marker classes
    # for later plot ordering or faceting.
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


def summarize_object(
    obj,
    expression_source,
    marker_genes,
    marker_groups,
    use_all_genes,
    default_missing_group_label,
):
    """Summarize one configured AnnData object for dotplot rendering.

    Args:
        obj: Config entry with sample, resolution, AnnData path, and groupby key.
        expression_source: Requested expression source, usually `"raw"`.
        marker_genes: Marker genes to summarize, or `None` in all-gene mode.
        marker_groups: Mapping from gene to marker group label.
        use_all_genes: Whether to summarize every available gene.
        default_missing_group_label: Label used when cells have missing group
            metadata.

    Returns:
        List of output rows, where each row represents one dotplot dot.

    Raises:
        KeyError: If the configured `groupby` column is missing.
    """
    import anndata as ad
    import numpy as np

    sample = obj["sample"]
    resolution = str(obj.get("resolution", ""))
    adata_path = obj["adata_path"]
    groupby = obj["groupby"]
    groupby_label = obj.get("groupby_label", groupby)
    missing_group_label = obj.get(
        "missing_group_label",
        default_missing_group_label,
    )

    print(f"Reading {adata_path}")
    adata = ad.read_h5ad(adata_path)

    if groupby not in adata.obs.columns:
        raise KeyError(
            f"{groupby!r} was not found in {adata_path}. "
            f"Available obs columns: {list(adata.obs.columns)}"
        )

    # Clean expression objects created by script 02 store log-normalized
    # expression in both .X and .raw. Prefer .raw when requested because it is
    # stable against later subsetting or transformations.
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

    # Subset once per object, then summarize each group. Each output row is one
    # dot in a dotplot: color = mean expression, size = percent expressing.
    # Missing archive labels are retained as their own group so no cells vanish
    # from old-label summaries.
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

        for i, gene in enumerate(present_genes):
            rows.append(
                {
                    "sample": sample,
                    "resolution": resolution,
                    "cluster_id": group_id,
                    "group_id": group_id,
                    "group_label": group_id,
                    "sample_cluster": sample_group,
                    "sample_group": sample_group,
                    "groupby": groupby,
                    "groupby_label": groupby_label,
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
    """Export long-format dotplot summary CSVs.

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
    objects = cfg["objects"]
    use_all_genes, marker_genes, marker_groups = load_marker_settings(cfg)

    rows = []
    for obj in objects:
        rows.extend(
            summarize_object(
                obj,
                expression_source,
                marker_genes,
                marker_groups,
                use_all_genes,
                default_missing_group_label,
            )
        )

    # Write one combined table by default, or split tables when the config asks
    # for one CSV per resolution/sample/etc.
    out = pd.DataFrame(rows)
    return write_dotplot_outputs(out, cfg)


def main():
    """Run the dotplot export workflow from a JSON config."""
    args = parse_args()
    cfg = load_config(args.config)
    export_dotplot_summary(cfg)


if __name__ == "__main__":
    main()
