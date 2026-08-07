#!/usr/bin/env python
# coding: utf-8

"""
Title: Plot Multi-Sample Dotplot From Config Local
Date: 2026-08-07
Summary: Local-working copy of the multi-sample dotplot renderer. Read one or
more long-format dotplot summary CSVs produced by
03_export_dotplot_data_from_config.py and render a combined multi-sample
dotplot from a JSON config using a compact local layout. Rows can represent
BANKSY clusters or any other exported grouping, including archived cell-type
labels. Per-sample filters can keep one resolution, several resolutions, or all
available resolutions for cross-resolution marker review.
"""

import argparse
from datetime import datetime
import json
import os
import shutil


import matplotlib.pyplot as plt
import pandas as pd


REQUIRED_COLUMNS = {
    "sample",
    "resolution",
    "cluster_id",
    "sample_cluster",
    "gene",
    "marker_group",
    "mean_expression",
    "percent_expressing",
}


DEFAULT_FIGURE = {
    "min_width": 12,
    "min_height": 8,
    "width_per_gene": 0.35,
    "height_per_cluster": 0.42,
    "max_dot_size": 220,
    "dpi": 100,
    "cmap": "viridis",
    "color_vmin": None,
    "color_vmax": None,
    "x_tick_fontsize": 16,
    "y_tick_fontsize": 14,
    "axis_label_fontsize": 16,
    "title_fontsize": 18,
    "title_y": 0.88,
    "colorbar_fontsize": 16,
    "colorbar_pad": 0.01,
    "size_legend_bbox_to_anchor": [1.16, 1.0],
    "size_legend_loc": "upper left",
    "gene_group_fontsize": 16,
    "show_gene_group_labels_top": True,
    "show_gene_group_labels_bottom": True,
    "gene_group_label_top_y": 1.03,
    "gene_group_label_bottom_y": -0.03,
    "x_tick_labeltop": True,
    "x_tick_labelbottom": True,
    "highlight_label_color": "#d62728",
    "highlight_label_weight": "bold",
    "tight_layout_rect": [0.03, 0.04, 0.98, 0.88],
}


# Local plotting is intentionally dense so large panels such as PTMT can show
# many genes without requiring config-only micromanagement of dot spacing.
DOTPLOT_WIDTH_COMPRESSION = 0.8
DOTPLOT_HEIGHT_COMPRESSION = 0.82
DOTPLOT_SIZE_BOOST = 1.12
GENE_DENDROGRAM_HEIGHT = 1.9
GENE_DENDROGRAM_HSPACE = 0.24


def parse_args():
    """Parse command-line arguments for the multi-sample dotplot renderer."""
    parser = argparse.ArgumentParser(
        prog="plot multi-sample dotplot from exported expression summaries"
    )
    parser.add_argument(
        "--config",
        required=True,
        help="JSON config defining input CSVs, plot options, and output path.",
    )
    return parser.parse_args()


def load_config(config_path):
    """Load the JSON plotting config.

    Args:
        config_path: Path to a JSON plot config.

    Returns:
        Dictionary containing input CSV paths and plotting options.
    """
    with open(config_path) as f:
        return json.load(f)


def validate_export_columns(df, csv_path):
    """Validate that a script 03 export CSV has the columns needed for plotting.

    Args:
        df: DataFrame read from a dotplot summary CSV.
        csv_path: Path used in the error message.

    Raises:
        ValueError: If required plotting columns are missing.
    """
    missing_columns = REQUIRED_COLUMNS - set(df.columns)
    if missing_columns:
        raise ValueError(f"{csv_path} is missing required columns: {missing_columns}")


def load_export_tables(input_csvs):
    """Read and concatenate one or more exported dotplot summary CSVs.

    Args:
        input_csvs: List of `_dotplot_summary.csv` paths produced by script 03.

    Returns:
        Combined DataFrame containing all requested samples/resolutions/clusters.
    """
    if not input_csvs:
        raise ValueError(
            "Config must define at least one CSV in dotplot_summary_csvs or input_csvs"
        )

    tables = []
    for csv_path in input_csvs:
        print(f"Reading {csv_path}")
        df = pd.read_csv(csv_path)
        validate_export_columns(df, csv_path)
        df["source_csv"] = csv_path
        tables.append(df)

    combined = pd.concat(tables, ignore_index=True)
    if "group_id" not in combined.columns:
        combined["group_id"] = combined["cluster_id"]
    if "group_label" not in combined.columns:
        combined["group_label"] = combined["group_id"]
    if "sample_group" not in combined.columns:
        combined["sample_group"] = combined["sample_cluster"]
    if "groupby_label" not in combined.columns:
        combined["groupby_label"] = combined.get("groupby", "cluster")

    for column in [
        "sample",
        "resolution",
        "cluster_id",
        "sample_cluster",
        "group_id",
        "group_label",
        "sample_group",
        "groupby_label",
        "gene",
    ]:
        combined[column] = combined[column].astype(str)

    return combined


def normalize_resolution_value(value):
    """Return a comparable resolution label while preserving non-numeric values."""
    value = str(value).strip()
    if not value:
        return value

    try:
        return f"{float(value):.6g}"
    except ValueError:
        return value


def filter_values(df, column, allowed_values):
    """Filter a DataFrame column when an allowed-value list is configured."""
    if not allowed_values:
        return df

    allowed_values = [str(value) for value in allowed_values]
    if column == "resolution":
        allowed_normalized = {normalize_resolution_value(value) for value in allowed_values}
        value_series = df[column].astype(str).map(normalize_resolution_value)
        filtered = df[value_series.isin(allowed_normalized)].copy()
    else:
        filtered = df[df[column].astype(str).isin(allowed_values)].copy()

    if filtered.empty:
        raise ValueError(f"Filtering {column!r} to {allowed_values} removed all rows")
    return filtered


def resolve_sample_filter_resolutions(item):
    """Return the configured resolutions for one sample filter.

    Args:
        item: One `sample_filters` entry.

    Returns:
        `None` when all resolutions should be kept, otherwise a list of
        requested resolution labels.
    """
    if item.get("include_all_resolutions", False):
        return None

    if "resolutions" in item:
        resolutions = item["resolutions"]
        if isinstance(resolutions, str):
            resolutions = [resolutions]
        return [str(value).strip() for value in resolutions if str(value).strip()]

    resolution = str(item.get("resolution", "")).strip()
    if not resolution or resolution.lower() == "all":
        return None

    return [resolution]


def apply_sample_filters(df, sample_filters):
    """Apply sample-specific filters such as chosen resolutions per sample.

    Args:
        df: Combined dotplot summary table.
        sample_filters: List of config entries with `sample` plus optional
            `resolution`, `resolutions`, or `include_all_resolutions`.

    Returns:
        Filtered DataFrame containing only requested sample/resolution rows.

    Raises:
        ValueError: If a requested sample is absent, if a sample is listed more
            than once, or if requested resolutions are unavailable.
    """
    if not sample_filters:
        return df

    filtered_parts = []
    seen_samples = set()

    for item in sample_filters:
        sample = str(item["sample"])
        requested_resolutions = resolve_sample_filter_resolutions(item)

        if sample in seen_samples:
            raise ValueError(f"Sample {sample!r} appears more than once in sample_filters")
        seen_samples.add(sample)

        sample_df = df[df["sample"].astype(str) == sample].copy()
        if sample_df.empty:
            raise ValueError(f"Sample {sample!r} was not found in the input CSVs")

        if requested_resolutions is not None:
            requested_normalized = {
                normalize_resolution_value(value) for value in requested_resolutions
            }
            available_normalized = sample_df["resolution"].astype(str).map(
                normalize_resolution_value
            )
            sample_df = sample_df[available_normalized.isin(requested_normalized)].copy()
            if sample_df.empty:
                raise ValueError(
                    f"Resolutions {requested_resolutions!r} were not found for "
                    f"sample {sample!r}"
                )

        filtered_parts.append(sample_df)

    if not filtered_parts:
        raise ValueError("sample_filters removed all rows")

    return pd.concat(filtered_parts, ignore_index=True)


def load_gene_file(gene_file, gene_column, gene_group_column):
    """Load an optional gene filter and grouping table.

    Args:
        gene_file: Optional CSV containing genes to plot.
        gene_column: Column containing gene symbols.
        gene_group_column: Optional column containing gene group labels.

    Returns:
        Tuple of `(requested_genes, gene_groups)`.
    """
    if not gene_file:
        return None, None

    gene_df = pd.read_csv(gene_file)
    if gene_column not in gene_df.columns:
        raise ValueError(
            f"Gene column {gene_column!r} not found in {gene_file}. "
            f"Available columns: {list(gene_df.columns)}"
        )
    if gene_group_column and gene_group_column not in gene_df.columns:
        raise ValueError(
            f"Gene group column {gene_group_column!r} not found in {gene_file}. "
            f"Available columns: {list(gene_df.columns)}"
        )

    gene_df = gene_df.dropna(subset=[gene_column]).copy()
    gene_df[gene_column] = gene_df[gene_column].astype(str)

    gene_groups = None
    if gene_group_column:
        gene_df[gene_group_column] = (
            gene_df[gene_group_column].fillna("unannotated").astype(str)
        )
        gene_df = gene_df.sort_values([gene_group_column, gene_column])
        gene_groups = (
            gene_df.drop_duplicates(gene_column)
            .set_index(gene_column)[gene_group_column]
            .to_dict()
        )

    requested_genes = gene_df[gene_column].drop_duplicates().tolist()
    return requested_genes, gene_groups


def is_true_flag_value(value):
    """Return whether a table value is an explicit TRUE-style flag."""
    if pd.isna(value):
        return False

    return str(value).strip().lower() in {"true", "t", "yes", "y", "1"}


def load_highlight_genes(highlight_gene_file, highlight_gene_column, gene_column="Gene"):
    """Load genes whose x-axis labels should be visually highlighted.

    Args:
        highlight_gene_file: Optional CSV containing genes to highlight.
        highlight_gene_column: Column containing either gene symbols or a
            TRUE/FALSE-style flag marking rows to highlight.
        gene_column: Column containing gene symbols when `highlight_gene_column`
            is a boolean-style flag column.

    Returns:
        Set of gene symbols to highlight.
    """
    if not highlight_gene_file:
        return set()

    highlight_df = pd.read_csv(highlight_gene_file)
    if highlight_gene_column not in highlight_df.columns:
        raise ValueError(
            f"Highlight gene column {highlight_gene_column!r} not found in "
            f"{highlight_gene_file}. Available columns: {list(highlight_df.columns)}"
        )

    highlight_flags = highlight_df[highlight_gene_column].map(is_true_flag_value)
    if highlight_flags.any():
        if gene_column not in highlight_df.columns:
            raise ValueError(
                f"Gene column {gene_column!r} is needed when "
                f"{highlight_gene_column!r} contains TRUE/FALSE highlight flags. "
                f"Available columns: {list(highlight_df.columns)}"
            )

        highlighted_genes = set(
            highlight_df.loc[highlight_flags, gene_column]
            .dropna()
            .astype(str)
            .drop_duplicates()
            .tolist()
        )
        print(
            f"Highlighted {len(highlighted_genes)} genes using "
            f"{highlight_gene_column!r} flags from {highlight_gene_file}"
        )
        return highlighted_genes

    return set(
        highlight_df[highlight_gene_column]
        .dropna()
        .astype(str)
        .drop_duplicates()
        .tolist()
    )


def is_keep_gene_value(value):
    """Return whether a review-table value means the gene should be plotted."""
    if pd.isna(value):
        return True

    normalized = str(value).strip().lower()
    if not normalized:
        return True

    return normalized not in {"false", "f", "no", "n", "0", "drop", "exclude"}


def apply_gene_review_filter(df, gene_order, gene_groups, cfg):
    """Remove genes marked for exclusion in an optional review CSV.

    Args:
        df: Filtered long-format dotplot summary table.
        gene_order: Ordered genes currently planned for the x-axis.
        gene_groups: Optional mapping from gene name to group label.
        cfg: Plot config dictionary.

    Returns:
        Tuple of filtered DataFrame, filtered gene order, and filtered gene groups.
    """
    review_file = cfg.get("gene_review_file")
    if not review_file:
        return df, gene_order, gene_groups

    gene_column = cfg.get("gene_review_gene_column", cfg.get("gene_column", "Gene"))
    keep_column = cfg.get("gene_review_keep_column", "keep_for_dotplot")
    review_df = pd.read_csv(review_file)

    for column in [gene_column, keep_column]:
        if column not in review_df.columns:
            raise ValueError(
                f"Gene review column {column!r} not found in {review_file}. "
                f"Available columns: {list(review_df.columns)}"
            )

    review_df = review_df.dropna(subset=[gene_column]).copy()
    review_df[gene_column] = review_df[gene_column].astype(str)
    review_df["keep_for_dotplot_resolved"] = review_df[keep_column].map(
        is_keep_gene_value
    )

    dropped_genes = set(
        review_df.loc[~review_df["keep_for_dotplot_resolved"], gene_column]
        .drop_duplicates()
        .tolist()
    )
    if not dropped_genes:
        return df, gene_order, gene_groups

    df = df[~df["gene"].astype(str).isin(dropped_genes)].copy()
    gene_order = [gene for gene in gene_order if gene not in dropped_genes]
    if gene_groups:
        gene_groups = {
            gene: group for gene, group in gene_groups.items() if gene not in dropped_genes
        }

    if df.empty or not gene_order:
        raise ValueError("gene_review_file removed all genes from the plot")

    print(f"Dropped {len(dropped_genes)} genes using {review_file}")
    return df, gene_order, gene_groups


def apply_expression_zscore(df, cfg):
    """Optionally z-score mean expression values within each gene.

    Args:
        df: Filtered long-format dotplot summary table after unwanted genes have
            been removed.
        cfg: Plot config dictionary.

    Returns:
        DataFrame with `mean_expression_zscore` added when configured.
    """
    if not cfg.get("z_score_expression", False):
        return df

    expression_column = cfg.get("z_score_expression_column", "mean_expression")
    if expression_column not in df.columns:
        raise ValueError(
            f"Z-score expression column {expression_column!r} was not found. "
            f"Available columns: {list(df.columns)}"
        )

    df = df.copy()
    expression = pd.to_numeric(df[expression_column], errors="coerce")
    gene_means = expression.groupby(df["gene"]).transform("mean")
    gene_stds = expression.groupby(df["gene"]).transform(lambda values: values.std(ddof=0))
    z_scores = (expression - gene_means) / gene_stds.replace(0, pd.NA)
    z_scores = z_scores.fillna(0.0)

    z_score_clip = cfg.get("z_score_clip")
    if z_score_clip is not None:
        z_scores = z_scores.clip(lower=-float(z_score_clip), upper=float(z_score_clip))

    df["mean_expression_zscore"] = z_scores
    cfg.setdefault("color_value_column", "mean_expression_zscore")
    cfg.setdefault("colorbar_label", "Mean expression z-score")
    cfg.setdefault("figure", {})
    cfg["figure"].setdefault("cmap", "RdBu_r")
    cfg["figure"].setdefault("color_vmin", -float(z_score_clip or 2.5))
    cfg["figure"].setdefault("color_vmax", float(z_score_clip or 2.5))
    print("Z-scored mean_expression within each gene for dotplot colors")
    return df


def cluster_gene_order(df, gene_order, cfg):
    """Optionally reorder genes by hierarchical clustering of expression summaries.

    Args:
        df: Filtered long-format dotplot summary table.
        gene_order: Current ordered gene list from the marker/review files.
        cfg: Plot config dictionary.

    Returns:
        Tuple of `(gene_order, linkage_matrix, clustered_genes)`. The linkage
        output is `None` when clustering is disabled or cannot be computed.
    """
    if not cfg.get("cluster_genes", False):
        return gene_order, None, []

    import numpy as np
    from scipy.cluster.hierarchy import leaves_list, linkage
    from scipy.spatial.distance import pdist

    value_column = cfg.get(
        "gene_clustering_value_column",
        cfg.get("color_value_column", "mean_expression"),
    )
    if value_column not in df.columns:
        raise ValueError(
            f"Gene clustering value column {value_column!r} was not found. "
            f"Available columns: {list(df.columns)}"
        )

    matrix = df.pivot_table(
        index="sample_group",
        columns="gene",
        values=value_column,
        aggfunc="mean",
    )
    matrix = matrix.reindex(columns=gene_order)
    matrix = matrix.apply(pd.to_numeric, errors="coerce")

    fill_value = cfg.get("gene_clustering_fill_value", 0.0)
    matrix = matrix.fillna(float(fill_value))

    gene_matrix = matrix.T
    missing_genes = [gene for gene in gene_order if gene not in gene_matrix.index]
    if missing_genes:
        print(f"Gene clustering skipped {len(missing_genes)} genes missing from matrix")

    variance_threshold = float(cfg.get("gene_clustering_variance_threshold", 0.0))
    gene_variance = gene_matrix.var(axis=1)
    variable_genes = gene_variance[gene_variance > variance_threshold].index.tolist()
    constant_genes = [
        gene for gene in gene_order
        if gene in gene_matrix.index and gene not in variable_genes
    ]

    if len(variable_genes) < 2:
        print(
            "Gene clustering skipped: fewer than two variable genes after "
            "removing constant expression profiles"
        )
        return gene_order, None, []

    variable_matrix = gene_matrix.loc[variable_genes]
    metric = cfg.get("gene_clustering_metric", "correlation")
    method = cfg.get("gene_clustering_method", "average")
    fallback_metric = cfg.get("gene_clustering_fallback_metric")

    try:
        distances = pdist(variable_matrix.to_numpy(), metric=metric)
        if not np.isfinite(distances).all():
            raise ValueError(f"{metric!r} produced non-finite distances")
    except Exception as error:
        if not fallback_metric:
            raise ValueError(
                "Gene clustering failed. Consider setting "
                "`gene_clustering_fallback_metric`, for example 'euclidean'."
            ) from error
        print(
            f"Gene clustering metric {metric!r} failed ({error}); "
            f"falling back to {fallback_metric!r}"
        )
        distances = pdist(variable_matrix.to_numpy(), metric=fallback_metric)
        if not np.isfinite(distances).all():
            raise ValueError(
                f"Fallback gene clustering metric {fallback_metric!r} produced "
                "non-finite distances"
            )

    linkage_matrix = linkage(distances, method=method)
    leaf_order = leaves_list(linkage_matrix)
    ordered_variable_genes = variable_matrix.index[leaf_order].tolist()

    constant_gene_position = cfg.get("constant_gene_position", "end")
    if constant_gene_position == "start":
        ordered_genes = constant_genes + ordered_variable_genes
    elif constant_gene_position == "drop":
        ordered_genes = ordered_variable_genes
    else:
        ordered_genes = ordered_variable_genes + constant_genes

    remaining_genes = [gene for gene in gene_order if gene not in ordered_genes]
    ordered_genes.extend(remaining_genes)

    print(
        f"Clustered {len(ordered_variable_genes)} variable genes using "
        f"{method} linkage and {metric} distance"
    )
    if constant_genes:
        print(
            f"Kept {len(constant_genes)} constant genes at "
            f"{constant_gene_position!r} of gene order"
        )
    return ordered_genes, linkage_matrix, ordered_variable_genes


def resolve_gene_order(df, cfg):
    """Filter genes and determine x-axis order.

    Args:
        df: Combined dotplot summary table.
        cfg: Plot config dictionary.

    Returns:
        Tuple of filtered DataFrame, ordered gene list, and gene-group mapping.
    """
    requested_genes, gene_groups = load_gene_file(
        cfg.get("gene_file"),
        cfg.get("gene_column", "Gene"),
        cfg.get("gene_group_column"),
    )

    if requested_genes:
        exported_genes = set(df["gene"].astype(str))
        missing_genes = [gene for gene in requested_genes if gene not in exported_genes]
        if missing_genes:
            print(f"Warning: {len(missing_genes)} requested genes were not found")

        df = df[df["gene"].astype(str).isin(requested_genes)].copy()
        if df.empty:
            raise ValueError("No requested genes were found in the export CSVs")

        gene_order = [gene for gene in requested_genes if gene in exported_genes]
        return df, gene_order, gene_groups

    gene_order = cfg.get("gene_order")
    if gene_order:
        gene_order = [str(gene) for gene in gene_order]
        df = df[df["gene"].astype(str).isin(gene_order)].copy()
        if df.empty:
            raise ValueError("Configured gene_order removed all rows")
        return df, [gene for gene in gene_order if gene in set(df["gene"])], gene_groups

    # Without an external gene list, keep marker groups together using labels
    # already present in the exported summary table.
    gene_order = (
        df[["marker_group", "gene"]]
        .drop_duplicates()
        .sort_values(["marker_group", "gene"])["gene"]
        .tolist()
    )
    return df, gene_order, gene_groups


def numeric_aware_sort(values):
    """Sort labels numerically when possible, then lexicographically."""
    def sort_key(value):
        text_value = str(value)
        try:
            return (0, float(text_value), text_value)
        except ValueError:
            return (1, text_value)

    return sorted([str(value) for value in values], key=sort_key)


def sort_with_optional_order(values, configured_order, normalize_values=False):
    """Sort values using a configured prefix order followed by remaining values."""
    values = [str(value) for value in values]
    if not configured_order:
        return numeric_aware_sort(values)

    configured_order = [str(value) for value in configured_order]
    if normalize_values:
        value_lookup = {normalize_resolution_value(value): value for value in values}
        configured_present = [
            value_lookup[normalize_resolution_value(value)]
            for value in configured_order
            if normalize_resolution_value(value) in value_lookup
        ]
        configured_present_set = set(configured_present)
    else:
        configured_present = [value for value in configured_order if value in values]
        configured_present_set = set(configured_present)

    remaining = numeric_aware_sort(
        value for value in values if value not in configured_present_set
    )
    return configured_present + remaining


def build_cluster_order(df, cfg):
    """Build y-axis order from sample, resolution, and cluster identifiers.

    Args:
        df: Combined dotplot summary table.
        cfg: Plot config dictionary.

    Returns:
        Ordered list of `sample_cluster` identifiers for the y-axis.
    """
    sample_order = sort_with_optional_order(
        df["sample"].drop_duplicates(), cfg.get("sample_order")
    )
    resolution_order = sort_with_optional_order(
        df["resolution"].drop_duplicates(),
        cfg.get("resolution_order"),
        normalize_values=True,
    )

    ordered_rows = []
    unique_clusters = df[
        ["sample", "resolution", "group_id", "group_label", "sample_group"]
    ].drop_duplicates()
    for sample in sample_order:
        sample_df = unique_clusters[unique_clusters["sample"] == sample]
        for resolution in resolution_order:
            resolution_mask = sample_df["resolution"].astype(str).map(
                normalize_resolution_value
            ) == normalize_resolution_value(resolution)
            res_df = sample_df[resolution_mask].copy()
            if res_df.empty:
                continue
            res_df["cluster_sort"] = pd.to_numeric(
                res_df["group_id"], errors="coerce"
            )
            res_df = res_df.sort_values(
                ["cluster_sort", "group_label"], na_position="last"
            )
            ordered_rows.extend(res_df["sample_group"].tolist())

    return ordered_rows


def make_cluster_labels(df, cluster_order, cfg):
    """Create readable y-axis labels for sample/resolution/group rows."""
    label_df = (
        df[[
            "sample_group",
            "sample",
            "resolution",
            "group_id",
            "group_label",
            "groupby_label",
        ]]
        .drop_duplicates("sample_group")
        .set_index("sample_group")
    )
    label_template = cfg.get(
        "y_label_template",
        "{sample} | r{resolution} | {group_label}",
    )
    labels = []
    for sample_group in cluster_order:
        row = label_df.loc[sample_group]
        labels.append(
            label_template.format(
                sample=row["sample"],
                resolution=row["resolution"],
                group_id=row["group_id"],
                cluster_id=row["group_id"],
                group_label=row["group_label"],
                groupby_label=row["groupby_label"],
            )
        )
    return labels


def add_gene_group_labels(ax, gene_order, gene_groups, figure_cfg):
    """Draw vertical separators and labels for grouped genes.

    Args:
        ax: Matplotlib axes containing the dotplot.
        gene_order: Ordered list of genes along the x-axis.
        gene_groups: Optional mapping from gene name to group label.
        figure_cfg: Plot sizing and style options from the config.
    """
    label_artists = []
    if not gene_groups or not gene_order:
        return label_artists

    current_group = gene_groups.get(gene_order[0], "unannotated")
    start = 0
    group_runs = []

    for idx, gene_name in enumerate(gene_order[1:], start=1):
        group_name = gene_groups.get(gene_name, "unannotated")
        if group_name != current_group:
            group_runs.append((current_group, start, idx - 1))
            ax.axvline(idx - 0.5, color="#9a9a9a", linewidth=0.8)
            current_group = group_name
            start = idx

    group_runs.append((current_group, start, len(gene_order) - 1))
    show_top = figure_cfg.get("show_gene_group_labels_top", True)
    show_bottom = figure_cfg.get("show_gene_group_labels_bottom", True)
    top_y = figure_cfg.get("gene_group_label_top_y", 1.03)
    bottom_y = figure_cfg.get("gene_group_label_bottom_y", -0.03)

    for group_name, start, end in group_runs:
        midpoint = (start + end) / 2
        if show_top:
            top_label = ax.text(
                midpoint,
                top_y,
                group_name,
                ha="center",
                va="bottom",
                rotation=90,
                fontsize=figure_cfg["gene_group_fontsize"],
                color="#333333",
                transform=ax.get_xaxis_transform(),
                clip_on=False,
            )
            # Keep outside-axis group labels from forcing tight_layout to shrink
            # the dot panel when their configured offset changes.
            top_label.set_in_layout(False)
            label_artists.append(top_label)
        if show_bottom:
            bottom_label = ax.text(
                midpoint,
                bottom_y,
                group_name,
                ha="center",
                va="top",
                rotation=90,
                fontsize=figure_cfg["gene_group_fontsize"],
                color="#333333",
                transform=ax.get_xaxis_transform(),
                clip_on=False,
            )
            bottom_label.set_in_layout(False)
            label_artists.append(bottom_label)

    return label_artists


def style_highlighted_gene_labels(ax, highlight_genes, figure_cfg):
    """Highlight selected gene tick labels on the x-axis.

    Args:
        ax: Matplotlib axes containing the dotplot.
        highlight_genes: Set of gene symbols to highlight.
        figure_cfg: Plot style options from the config.
    """
    if not highlight_genes:
        return

    for tick_label in ax.xaxis.get_ticklabels():
        if tick_label.get_text() in highlight_genes:
            tick_label.set_color(figure_cfg["highlight_label_color"])
            tick_label.set_fontweight(figure_cfg["highlight_label_weight"])


def add_group_separators(ax, df, cluster_order, cfg):
    """Draw horizontal separators between sample and resolution blocks."""
    label_df = (
        df[["sample_group", "sample", "resolution"]]
        .drop_duplicates("sample_group")
        .set_index("sample_group")
    )
    previous_sample = None
    previous_resolution = None
    for idx, sample_group in enumerate(cluster_order):
        row = label_df.loc[sample_group]
        sample = row["sample"]
        resolution = row["resolution"]
        if previous_sample is not None and sample != previous_sample:
            ax.axhline(idx - 0.5, color="#555555", linewidth=1.0)
        elif (
            cfg.get("show_resolution_separators", True)
            and previous_resolution is not None
            and resolution != previous_resolution
        ):
            ax.axhline(idx - 0.5, color="#9a9a9a", linewidth=0.7, linestyle="--")
        previous_sample = sample
        previous_resolution = resolution


def plot_dotplot(
    df,
    gene_order,
    cluster_order,
    gene_groups,
    highlight_genes,
    cfg,
    gene_linkage=None,
    clustered_genes=None,
):
    """Render and save the multi-sample dotplot.

    Args:
        df: Filtered long-format dotplot summary table.
        gene_order: Ordered genes for the x-axis.
        cluster_order: Ordered cluster identifiers for the y-axis.
        gene_groups: Optional mapping from gene to group label.
        highlight_genes: Set of genes whose x-axis labels should be highlighted.
        cfg: Plot config dictionary.
        gene_linkage: Optional hierarchical clustering linkage matrix for genes.
        clustered_genes: Genes represented in `gene_linkage`, in dendrogram leaf order.
    """
    figure_cfg = {**DEFAULT_FIGURE, **cfg.get("figure", {})}
    output_png = cfg["output_png"]

    gene_to_x = {gene: i for i, gene in enumerate(gene_order)}
    cluster_to_y = {cluster: i for i, cluster in enumerate(cluster_order)}

    color_value_column = cfg.get("color_value_column", "mean_expression")
    if color_value_column not in df.columns:
        raise ValueError(
            f"Color value column {color_value_column!r} was not found. "
            f"Available columns: {list(df.columns)}"
        )

    plot_df = df.copy()
    plot_df["x"] = plot_df["gene"].map(gene_to_x)
    plot_df["y"] = plot_df["sample_group"].map(cluster_to_y)
    plot_df = plot_df.dropna(subset=["x", "y"])
    plot_df["dot_size"] = (
        plot_df["percent_expressing"] / 100
    ) * figure_cfg["max_dot_size"] * DOTPLOT_SIZE_BOOST

    compact_min_width = figure_cfg["min_width"] * DOTPLOT_WIDTH_COMPRESSION
    compact_width_per_gene = figure_cfg["width_per_gene"] * DOTPLOT_WIDTH_COMPRESSION
    compact_min_height = figure_cfg["min_height"] * DOTPLOT_HEIGHT_COMPRESSION
    compact_height_per_cluster = (
        figure_cfg["height_per_cluster"] * DOTPLOT_HEIGHT_COMPRESSION
    )

    fig_width = max(compact_min_width, len(gene_order) * compact_width_per_gene)
    fig_height = max(
        compact_min_height,
        len(cluster_order) * compact_height_per_cluster,
    )
    show_gene_dendrogram = bool(cfg.get("show_gene_dendrogram", False))
    dendrogram_enabled = (
        show_gene_dendrogram and gene_linkage is not None and bool(clustered_genes)
    )
    if dendrogram_enabled:
        from scipy.cluster.hierarchy import dendrogram

        dendrogram_height = GENE_DENDROGRAM_HEIGHT
        dendrogram_hspace = GENE_DENDROGRAM_HSPACE
        fig, (dendro_ax, ax) = plt.subplots(
            2,
            1,
            figsize=(fig_width, fig_height + dendrogram_height),
            gridspec_kw={
                "height_ratios": [dendrogram_height, fig_height],
                "hspace": dendrogram_hspace,
            },
            sharex=True,
        )
        dendrogram_data = dendrogram(
            gene_linkage,
            orientation="top",
            no_plot=True,
        )
        clustered_x_offset = gene_to_x[clustered_genes[0]]
        for xs, ys in zip(dendrogram_data["icoord"], dendrogram_data["dcoord"]):
            # SciPy places leaves at 5, 15, 25, ...; transform those coordinates
            # onto the dotplot's integer gene positions so branches align with ticks.
            aligned_xs = [((x - 5.0) / 10.0) + clustered_x_offset for x in xs]
            dendro_ax.plot(aligned_xs, ys, color="#555555", linewidth=0.8)
        dendro_ax.set_axis_off()
    else:
        dendro_ax = None
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    scatter = ax.scatter(
        plot_df["x"],
        plot_df["y"],
        s=plot_df["dot_size"],
        c=plot_df[color_value_column],
        cmap=figure_cfg["cmap"],
        vmin=figure_cfg.get("color_vmin"),
        vmax=figure_cfg.get("color_vmax"),
        edgecolors="black",
        linewidths=0.25,
    )

    ax.set_xticks(range(len(gene_order)))
    ax.set_xticklabels(
        gene_order, rotation=90, fontsize=figure_cfg["x_tick_fontsize"]
    )
    ax.tick_params(
        axis="x",
        top=figure_cfg.get("x_tick_labeltop", True),
        labeltop=figure_cfg.get("x_tick_labeltop", True),
        bottom=figure_cfg.get("x_tick_labelbottom", True),
        labelbottom=figure_cfg.get("x_tick_labelbottom", True),
    )
    ax.set_yticks(range(len(cluster_order)))
    ax.set_yticklabels(
        make_cluster_labels(df, cluster_order, cfg),
        fontsize=figure_cfg["y_tick_fontsize"],
    )
    ax.set_ylim(-0.75, len(cluster_order) - 0.25)
    ax.set_xlim(-0.75, len(gene_order) - 0.35)
    if dendro_ax is not None:
        dendro_ax.set_xlim(ax.get_xlim())
    ax.set_xlabel("Gene", fontsize=figure_cfg["axis_label_fontsize"])
    ax.set_ylabel(
        cfg.get("y_axis_label", "Sample / resolution / group"),
        fontsize=figure_cfg["axis_label_fontsize"],
    )
    ax.tick_params(axis="x", labelsize=figure_cfg["x_tick_fontsize"])
    ax.tick_params(axis="y", labelsize=figure_cfg["y_tick_fontsize"])
    style_highlighted_gene_labels(ax, highlight_genes, figure_cfg)
    title_artist = fig.suptitle(
        cfg.get("title", "Multi-sample dotplot summary"),
        fontsize=figure_cfg["title_fontsize"],
        y=figure_cfg["title_y"],
    )

    gene_group_label_artists = add_gene_group_labels(
        ax,
        gene_order,
        gene_groups,
        figure_cfg,
    )
    add_group_separators(ax, df, cluster_order, cfg)

    ax.grid(axis="both", color="#e6e6e6", linewidth=0.5)
    ax.set_axisbelow(True)

    cbar = fig.colorbar(
        scatter,
        ax=ax,
        pad=figure_cfg.get("colorbar_pad", 0.01),
    )
    cbar.set_label(
        cfg.get("colorbar_label", "Mean expression"),
        fontsize=figure_cfg["colorbar_fontsize"],
    )
    cbar.ax.tick_params(labelsize=figure_cfg["colorbar_fontsize"])

    for pct in cfg.get("size_legend_percentages", [25, 50, 75, 100]):
        ax.scatter(
            [],
            [],
            s=(pct / 100) * figure_cfg["max_dot_size"] * DOTPLOT_SIZE_BOOST,
            c="white",
            edgecolors="black",
            linewidths=0.25,
            label=f"{pct}%",
        )
    size_legend_artist = ax.legend(
        title="Percent expressing",
        bbox_to_anchor=figure_cfg.get("size_legend_bbox_to_anchor", [1.16, 1.0]),
        loc=figure_cfg.get("size_legend_loc", "upper left"),
        borderaxespad=0,
        frameon=False,
    )

    if dendrogram_enabled:
        rect = figure_cfg["tight_layout_rect"]
        fig.subplots_adjust(
            left=rect[0],
            bottom=rect[1],
            right=rect[2],
            top=rect[3],
            hspace=dendrogram_hspace,
        )
        # The colorbar narrows the dotplot axes but not the dendrogram axes.
        # Match the final axes geometry so dendrogram leaves sit over the grid.
        dotplot_position = ax.get_position()
        dendrogram_position = dendro_ax.get_position()
        dendro_ax.set_position([
            dotplot_position.x0,
            dendrogram_position.y0,
            dotplot_position.width,
            dendrogram_position.height,
        ])
        dendro_ax.set_xlim(ax.get_xlim())
    else:
        fig.tight_layout(rect=figure_cfg["tight_layout_rect"])
    os.makedirs(os.path.dirname(output_png), exist_ok=True)
    archive_existing_output(output_png, cfg)
    fig.savefig(
        output_png,
        dpi=figure_cfg["dpi"],
        bbox_inches="tight",
        bbox_extra_artists=[title_artist, size_legend_artist, *gene_group_label_artists],
    )
    print(f"Wrote {output_png}")


def archive_existing_output(output_path, cfg):
    """Archive an existing plot file before overwriting it.

    Args:
        output_path: Destination image path that is about to be written.
        cfg: Plot config dictionary.

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


def get_dotplot_summary_csvs(cfg):
    """Return the configured script 03 summary CSV paths.

    Args:
        cfg: Plot config dictionary.

    Returns:
        List of `_dotplot_summary.csv` paths.
    """
    return cfg.get("dotplot_summary_csvs", cfg.get("input_csvs", []))


def main():
    """Read a plot config, combine export CSVs, and write one dotplot image."""
    args = parse_args()
    cfg = load_config(args.config)

    df = load_export_tables(get_dotplot_summary_csvs(cfg))
    df = apply_sample_filters(df, cfg.get("sample_filters"))
    df = filter_values(df, "sample", cfg.get("samples"))
    df = filter_values(df, "resolution", cfg.get("resolutions"))

    df, gene_order, gene_groups = resolve_gene_order(df, cfg)
    df, gene_order, gene_groups = apply_gene_review_filter(
        df=df,
        gene_order=gene_order,
        gene_groups=gene_groups,
        cfg=cfg,
    )
    df = apply_expression_zscore(df, cfg)
    gene_order, gene_linkage, clustered_genes = cluster_gene_order(df, gene_order, cfg)
    highlight_genes = load_highlight_genes(
        cfg.get("highlight_gene_file"),
        cfg.get("highlight_gene_column", cfg.get("gene_column", "Gene")),
        cfg.get("gene_column", "Gene"),
    )
    cluster_order = build_cluster_order(df, cfg)
    if not cluster_order:
        raise ValueError("No clusters remained after filtering")

    plot_dotplot(
        df=df,
        gene_order=gene_order,
        cluster_order=cluster_order,
        gene_groups=gene_groups,
        highlight_genes=highlight_genes,
        cfg=cfg,
        gene_linkage=gene_linkage,
        clustered_genes=clustered_genes,
    )


if __name__ == "__main__":
    main()
