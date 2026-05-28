#!/usr/bin/env python
# coding: utf-8

"""
Title: Plot Dotplot From Exported Expression Summary
Date: 2026-05-28
Summary: Read the long-format cluster-level expression CSV produced by
03_export_dotplot_data_from_config.py and render a dotplot where dot color
represents mean expression and dot size represents percent expressing. Optionally
filter and group genes using a marker annotation CSV.
"""

import argparse
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


def parse_args():
    """Parse command-line arguments for the dotplot renderer."""
    parser = argparse.ArgumentParser(
        prog="plot dotplot from exported cluster-level expression summary"
    )
    parser.add_argument(
        "--input-csv",
        default="data/xenium/processed/cross_sample_dotplot_exports/local_test_CK_skin_res_0p5_dotplot_summary.csv",
        help="CSV created by 03_export_dotplot_data_from_config.py.",
    )
    parser.add_argument(
        "--output-png",
        default="testing/local_test_CK_skin_res_0p5_dotplot.png",
        help="Path for the output dotplot image.",
    )
    parser.add_argument(
        "--max-dot-size",
        type=float,
        default=260,
        help="Maximum matplotlib marker size for 100 percent expressing.",
    )
    parser.add_argument(
        "--gene-file",
        default=None,
        help="Optional CSV containing genes to plot.",
    )
    parser.add_argument(
        "--gene-column",
        default="Gene",
        help="Column in --gene-file containing gene symbols.",
    )
    parser.add_argument(
        "--gene-group-column",
        default=None,
        help="Optional column in --gene-file used to group/order genes on the x-axis.",
    )
    return parser.parse_args()


def validate_export_columns(df):
    """Raise an error if the exported dotplot summary is missing required columns.

    Args:
        df: DataFrame read from the long-format dotplot summary CSV.

    Raises:
        ValueError: If required columns for plotting are absent.
    """
    missing_columns = REQUIRED_COLUMNS - set(df.columns)
    if missing_columns:
        raise ValueError(f"Input CSV is missing required columns: {missing_columns}")


def filter_and_order_genes(df, gene_file, gene_column, gene_group_column):
    """Filter genes and determine x-axis order for the dotplot.

    Args:
        df: Exported dotplot summary table.
        gene_file: Optional CSV containing a gene list for filtering.
        gene_column: Column in `gene_file` containing gene symbols.
        gene_group_column: Optional column in `gene_file` used to group genes.

    Returns:
        Tuple of filtered DataFrame, ordered gene list, and optional gene-to-group map.
    """
    gene_groups = None

    if gene_file:
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

        # When a grouping column is supplied, sort by group first so the final
        # plot visually separates major annotation classes.
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
        exported_genes = set(df["gene"].astype(str))
        missing_genes = [gene for gene in requested_genes if gene not in exported_genes]
        if missing_genes:
            print(
                f"Warning: {len(missing_genes)} requested genes were not found in the export CSV"
            )

        df = df[df["gene"].astype(str).isin(requested_genes)].copy()
        if df.empty:
            raise ValueError("No requested genes were found in the export CSV")

        gene_order = [gene for gene in requested_genes if gene in exported_genes]
        return df, gene_order, gene_groups

    # Without an external gene list, keep marker groups together using the group
    # labels already present in the exported summary table.
    gene_order = (
        df[["marker_group", "gene"]]
        .drop_duplicates()
        .sort_values(["marker_group", "gene"])
        ["gene"]
        .tolist()
    )
    return df, gene_order, gene_groups


def get_cluster_order(df):
    """Return y-axis cluster order from sample, resolution, and cluster IDs."""
    return (
        df[["sample", "resolution", "cluster_id", "sample_cluster"]]
        .drop_duplicates()
        .sort_values(["sample", "resolution", "cluster_id"])
        ["sample_cluster"]
        .tolist()
    )


def add_gene_group_labels(ax, gene_order, cluster_order, gene_groups):
    """Draw vertical separators and labels for grouped genes.

    Args:
        ax: Matplotlib axes containing the dotplot.
        gene_order: Ordered list of genes along the x-axis.
        cluster_order: Ordered list of y-axis cluster labels.
        gene_groups: Mapping from gene name to group label.
    """
    if not gene_groups or not gene_order:
        return

    group_runs = []
    current_group = gene_groups.get(gene_order[0], "unannotated")
    start = 0

    for idx, gene_name in enumerate(gene_order[1:], start=1):
        group_name = gene_groups.get(gene_name, "unannotated")
        if group_name != current_group:
            group_runs.append((current_group, start, idx - 1))
            ax.axvline(idx - 0.5, color="#9a9a9a", linewidth=0.8)
            current_group = group_name
            start = idx

    group_runs.append((current_group, start, len(gene_order) - 1))

    # Place group labels above the uppermost cluster row. The y-axis padding set
    # in plot_dotplot keeps these labels from being clipped.
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


def plot_dotplot(df, gene_order, cluster_order, gene_groups, output_png, max_dot_size):
    """Render and save the dotplot from a dot-level summary table.

    Args:
        df: Filtered long-format dotplot summary table.
        gene_order: Ordered genes for the x-axis.
        cluster_order: Ordered sample/resolution/cluster labels for the y-axis.
        gene_groups: Optional mapping from gene to group label.
        output_png: Path to save the rendered PNG.
        max_dot_size: Marker size corresponding to 100 percent expressing.
    """
    gene_to_x = {gene: i for i, gene in enumerate(gene_order)}
    cluster_to_y = {cluster: i for i, cluster in enumerate(cluster_order)}

    plot_df = df.copy()
    plot_df["x"] = plot_df["gene"].map(gene_to_x)
    plot_df["y"] = plot_df["sample_cluster"].map(cluster_to_y)
    plot_df["dot_size"] = (plot_df["percent_expressing"] / 100) * max_dot_size

    # Scale the canvas with the number of genes and clusters so dense plots do
    # not collapse rows or clip large dots.
    fig_width = max(10, len(gene_order) * 0.35)
    fig_height = max(6, len(cluster_order) * 0.85)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    scatter = ax.scatter(
        plot_df["x"],
        plot_df["y"],
        s=plot_df["dot_size"],
        c=plot_df["mean_expression"],
        cmap="viridis",
        edgecolors="black",
        linewidths=0.25,
    )

    ax.set_xticks(range(len(gene_order)))
    ax.set_xticklabels(gene_order, rotation=90)
    ax.set_yticks(range(len(cluster_order)))
    ax.set_yticklabels(cluster_order)
    ax.set_ylim(-0.75, len(cluster_order) - 0.25)
    ax.set_xlabel("Gene")
    ax.set_ylabel("Sample / resolution / cluster")
    ax.set_title("Dotplot summary from exported expression values")

    add_gene_group_labels(ax, gene_order, cluster_order, gene_groups)

    ax.grid(axis="both", color="#e6e6e6", linewidth=0.5)
    ax.set_axisbelow(True)

    cbar = fig.colorbar(scatter, ax=ax, pad=0.01)
    cbar.set_label("Mean expression")

    # Dummy scatter points create a dot-size legend without adding extra data to
    # the plotting area.
    for pct in [25, 50, 75, 100]:
        ax.scatter(
            [],
            [],
            s=(pct / 100) * max_dot_size,
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
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    print(f"Wrote {output_png}")


def main():
    """Read the export CSV, optionally filter genes, and write a dotplot PNG."""
    args = parse_args()
    df = pd.read_csv(args.input_csv)
    validate_export_columns(df)

    df, gene_order, gene_groups = filter_and_order_genes(
        df,
        args.gene_file,
        args.gene_column,
        args.gene_group_column,
    )
    cluster_order = get_cluster_order(df)

    plot_dotplot(
        df=df,
        gene_order=gene_order,
        cluster_order=cluster_order,
        gene_groups=gene_groups,
        output_png=args.output_png,
        max_dot_size=args.max_dot_size,
    )


if __name__ == "__main__":
    main()
