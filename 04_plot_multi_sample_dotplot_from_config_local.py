#!/usr/bin/env python
# coding: utf-8

"""
Title: Plot Multi-Sample Dotplot From Config Local
Date: 2026-06-02
Summary: Local-working copy of the multi-sample dotplot renderer. Read one or
more long-format dotplot summary CSVs produced by
03_export_dotplot_data_from_config.py and render a combined multi-sample
dotplot from a JSON config.
"""

import argparse
import json
import os

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
    "dpi": 300,
    "cmap": "viridis",
}


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
    for column in ["sample", "resolution", "cluster_id", "sample_cluster", "gene"]:
        combined[column] = combined[column].astype(str)

    return combined


def filter_values(df, column, allowed_values):
    """Filter a DataFrame column when an allowed-value list is configured."""
    if not allowed_values:
        return df

    allowed_values = [str(value) for value in allowed_values]
    filtered = df[df[column].astype(str).isin(allowed_values)].copy()
    if filtered.empty:
        raise ValueError(f"Filtering {column!r} to {allowed_values} removed all rows")
    return filtered


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


def sort_with_optional_order(values, configured_order):
    """Sort values using a configured prefix order followed by remaining values."""
    values = [str(value) for value in values]
    if not configured_order:
        return sorted(values)

    configured_order = [str(value) for value in configured_order]
    configured_present = [value for value in configured_order if value in values]
    remaining = sorted(value for value in values if value not in configured_present)
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
    )

    ordered_rows = []
    unique_clusters = df[
        ["sample", "resolution", "cluster_id", "sample_cluster"]
    ].drop_duplicates()
    for sample in sample_order:
        sample_df = unique_clusters[unique_clusters["sample"] == sample]
        for resolution in resolution_order:
            res_df = sample_df[sample_df["resolution"] == resolution].copy()
            if res_df.empty:
                continue
            res_df["cluster_sort"] = pd.to_numeric(
                res_df["cluster_id"], errors="coerce"
            )
            res_df = res_df.sort_values(
                ["cluster_sort", "cluster_id"], na_position="last"
            )
            ordered_rows.extend(res_df["sample_cluster"].tolist())

    return ordered_rows


def make_cluster_labels(df, cluster_order):
    """Create readable y-axis labels for sample/resolution/cluster rows."""
    label_df = (
        df[["sample_cluster", "sample", "resolution", "cluster_id"]]
        .drop_duplicates("sample_cluster")
        .set_index("sample_cluster")
    )
    labels = []
    for sample_cluster in cluster_order:
        row = label_df.loc[sample_cluster]
        labels.append(f"{row['sample']} | r{row['resolution']} | c{row['cluster_id']}")
    return labels


def add_gene_group_labels(ax, gene_order, cluster_order, gene_groups):
    """Draw vertical separators and labels for grouped genes.

    Args:
        ax: Matplotlib axes containing the dotplot.
        gene_order: Ordered list of genes along the x-axis.
        cluster_order: Ordered list of y-axis cluster identifiers.
        gene_groups: Optional mapping from gene name to group label.
    """
    if not gene_groups or not gene_order:
        return

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

    top_y = len(cluster_order) - 0.05
    for group_name, start, end in group_runs:
        midpoint = (start + end) / 2
        ax.text(
            midpoint,
            top_y,
            group_name,
            ha="center",
            va="bottom",
            rotation=90,
            fontsize=8,
            color="#333333",
        )


def add_sample_separators(ax, df, cluster_order):
    """Draw horizontal separators between sample blocks on the y-axis."""
    label_df = (
        df[["sample_cluster", "sample"]]
        .drop_duplicates("sample_cluster")
        .set_index("sample_cluster")
    )
    previous_sample = None
    for idx, sample_cluster in enumerate(cluster_order):
        sample = label_df.loc[sample_cluster, "sample"]
        if previous_sample is not None and sample != previous_sample:
            ax.axhline(idx - 0.5, color="#555555", linewidth=1.0)
        previous_sample = sample


def plot_dotplot(df, gene_order, cluster_order, gene_groups, cfg):
    """Render and save the multi-sample dotplot.

    Args:
        df: Filtered long-format dotplot summary table.
        gene_order: Ordered genes for the x-axis.
        cluster_order: Ordered cluster identifiers for the y-axis.
        gene_groups: Optional mapping from gene to group label.
        cfg: Plot config dictionary.
    """
    figure_cfg = {**DEFAULT_FIGURE, **cfg.get("figure", {})}
    output_png = cfg["output_png"]

    gene_to_x = {gene: i for i, gene in enumerate(gene_order)}
    cluster_to_y = {cluster: i for i, cluster in enumerate(cluster_order)}

    plot_df = df.copy()
    plot_df["x"] = plot_df["gene"].map(gene_to_x)
    plot_df["y"] = plot_df["sample_cluster"].map(cluster_to_y)
    plot_df = plot_df.dropna(subset=["x", "y"])
    plot_df["dot_size"] = (
        plot_df["percent_expressing"] / 100
    ) * figure_cfg["max_dot_size"]

    fig_width = max(
        figure_cfg["min_width"], len(gene_order) * figure_cfg["width_per_gene"]
    )
    fig_height = max(
        figure_cfg["min_height"],
        len(cluster_order) * figure_cfg["height_per_cluster"],
    )
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    scatter = ax.scatter(
        plot_df["x"],
        plot_df["y"],
        s=plot_df["dot_size"],
        c=plot_df["mean_expression"],
        cmap=figure_cfg["cmap"],
        edgecolors="black",
        linewidths=0.25,
    )

    ax.set_xticks(range(len(gene_order)))
    ax.set_xticklabels(gene_order, rotation=90)
    ax.set_yticks(range(len(cluster_order)))
    ax.set_yticklabels(make_cluster_labels(df, cluster_order))
    ax.set_ylim(-0.75, len(cluster_order) - 0.25)
    ax.set_xlabel("Gene")
    ax.set_ylabel("Sample / resolution / cluster")
    ax.set_title(cfg.get("title", "Multi-sample dotplot summary"))

    add_gene_group_labels(ax, gene_order, cluster_order, gene_groups)
    add_sample_separators(ax, df, cluster_order)

    ax.grid(axis="both", color="#e6e6e6", linewidth=0.5)
    ax.set_axisbelow(True)

    cbar = fig.colorbar(scatter, ax=ax, pad=0.01)
    cbar.set_label(cfg.get("colorbar_label", "Mean expression"))

    for pct in cfg.get("size_legend_percentages", [25, 50, 75, 100]):
        ax.scatter(
            [],
            [],
            s=(pct / 100) * figure_cfg["max_dot_size"],
            c="white",
            edgecolors="black",
            linewidths=0.25,
            label=f"{pct}%",
        )
    ax.legend(
        title="Percent expressing",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        borderaxespad=0,
        frameon=False,
    )

    fig.tight_layout()
    os.makedirs(os.path.dirname(output_png), exist_ok=True)
    fig.savefig(output_png, dpi=figure_cfg["dpi"], bbox_inches="tight")
    print(f"Wrote {output_png}")


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
    df = filter_values(df, "sample", cfg.get("samples"))
    df = filter_values(df, "resolution", cfg.get("resolutions"))

    df, gene_order, gene_groups = resolve_gene_order(df, cfg)
    cluster_order = build_cluster_order(df, cfg)
    if not cluster_order:
        raise ValueError("No clusters remained after filtering")

    plot_dotplot(
        df=df,
        gene_order=gene_order,
        cluster_order=cluster_order,
        gene_groups=gene_groups,
        cfg=cfg,
    )


if __name__ == "__main__":
    main()
