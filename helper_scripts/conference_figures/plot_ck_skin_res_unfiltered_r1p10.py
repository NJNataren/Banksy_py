#!/usr/bin/env python
# coding: utf-8

"""
Title: Plot CK Skin Unfiltered R1.10 Conference Figures
Date: 2026-08-27
Summary: Generate presentation-ready spatial, UMAP, marker, dotplot, and cluster
abundance figures for the unfiltered CK_skin_res script 00 clean expression
AnnData object at BANKSY resolution 1.10.
"""

import argparse
from pathlib import Path


DEFAULT_ADATA = (
    "data/xenium/processed/vbct/CK_skin_res/"
    "adata_expression_clean_CK_skin_res_with_banksy_clusters_"
    "0.50_0.60_0.70_0.80_0.90_1.00_1.10.h5ad"
)
DEFAULT_CLUSTER_COL = "labels_scaled_gaussian_pc30_nc0.20_r1.10"
DEFAULT_OUTPUT_DIR = "figures/conference/CK_skin_res_unfiltered_r1p10"
DEFAULT_ENDOTHELIAL_MARKERS = [
    "AQP1",
    "CALCRL",
    "CDH5",
    "ECSCR",
    "PLVAP",
    "SELP",
    "VWF",
    "TFPI",
]
DEFAULT_UMAP_KEYS = [
    "X_umap_scaled_gaussian_pc30_nc0.20",
    "X_umap",
]


pd = None
plt = None
sc = None


def load_plotting_stack():
    """Import plotting dependencies after CLI parsing."""
    global pd, plt, sc

    if sc is not None:
        return

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as pyplot
    import pandas as pandas
    import scanpy as scanpy

    pd = pandas
    plt = pyplot
    sc = scanpy


def parse_args():
    """Parse command-line arguments for the CK skin conference figure helper."""
    parser = argparse.ArgumentParser(
        description=(
            "Create conference-oriented CK_skin_res r1.10 figures from the "
            "unfiltered/script00 clean expression AnnData object."
        )
    )
    parser.add_argument("--adata", default=DEFAULT_ADATA, help="Input clean h5ad path.")
    parser.add_argument(
        "--cluster-col",
        default=DEFAULT_CLUSTER_COL,
        help="Observation column containing the BANKSY cluster labels to plot.",
    )
    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help="Directory where PNG/PDF figures will be written.",
    )
    parser.add_argument(
        "--markers",
        nargs="*",
        default=DEFAULT_ENDOTHELIAL_MARKERS,
        help="Marker genes to use for spatial expression panels and dotplot.",
    )
    parser.add_argument(
        "--point-size",
        type=float,
        default=3.0,
        help="Point size for spatial and UMAP scatter plots.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Figure resolution for saved raster images.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Show plots interactively instead of only saving them.",
    )
    return parser.parse_args()


def save_current_figure(output_dir, stem, dpi, save_pdf=True):
    """Save the current Matplotlib figure as PNG and optionally PDF."""
    png_path = output_dir / f"{stem}.png"
    plt.savefig(png_path, dpi=dpi, bbox_inches="tight")
    print(f"Wrote {png_path}")

    if save_pdf:
        pdf_path = output_dir / f"{stem}.pdf"
        plt.savefig(pdf_path, bbox_inches="tight")
        print(f"Wrote {pdf_path}")

    plt.close()


def prepare_adata(adata_path, cluster_col):
    """Load AnnData and add plotting-friendly spatial and cluster columns."""
    print(f"Reading {adata_path}")
    adata = sc.read_h5ad(adata_path)

    if cluster_col not in adata.obs.columns:
        raise KeyError(
            f"Cluster column {cluster_col!r} was not found. "
            f"Available obs columns include: {list(adata.obs.columns[:25])}"
        )

    # Scanpy's generic embedding plot expects X_<basis>, while Xenium objects
    # commonly store tissue coordinates under obsm['spatial'] or obsm['xy'].
    if "X_spatial" not in adata.obsm:
        if "spatial" in adata.obsm:
            adata.obsm["X_spatial"] = adata.obsm["spatial"]
        elif "xy" in adata.obsm:
            adata.obsm["X_spatial"] = adata.obsm["xy"]
        else:
            raise KeyError("No spatial coordinates found in obsm['spatial'] or obsm['xy'].")

    adata.obs["clusters_r1p10"] = adata.obs[cluster_col].astype(str).astype("category")
    return adata


def available_genes(adata, requested_genes):
    """Return requested genes present in `adata`, preserving requested order."""
    present = [gene for gene in requested_genes if gene in adata.var_names]
    missing = [gene for gene in requested_genes if gene not in adata.var_names]

    if missing:
        print(f"Skipping missing marker genes: {', '.join(missing)}")
    if not present:
        raise ValueError("None of the requested marker genes were found in adata.var_names.")

    return present


def plot_spatial_clusters(adata, output_dir, point_size, dpi, show):
    """Plot tissue coordinates colored by r1.10 cluster labels."""
    sc.pl.embedding(
        adata,
        basis="spatial",
        color="clusters_r1p10",
        size=point_size,
        frameon=False,
        legend_loc="right margin",
        title="CK_skin_res BANKSY clusters, r1.10",
        show=show,
    )
    save_current_figure(output_dir, "CK_skin_res_spatial_clusters_r1p10", dpi)

    sc.pl.embedding(
        adata,
        basis="spatial",
        color="clusters_r1p10",
        size=point_size,
        frameon=False,
        legend_loc="on data",
        title="CK_skin_res BANKSY clusters, r1.10",
        show=show,
    )
    save_current_figure(output_dir, "CK_skin_res_spatial_clusters_r1p10_on_data", dpi)


def plot_umap_clusters(adata, output_dir, point_size, dpi, show):
    """Plot BANKSY UMAP coordinates when a copied UMAP embedding is available."""
    for obsm_key in DEFAULT_UMAP_KEYS:
        if obsm_key in adata.obsm:
            adata.obsm["X_conference_umap"] = adata.obsm[obsm_key]
            print(f"Using {obsm_key} for conference UMAP plot.")
            break
    else:
        print("Skipping UMAP cluster plot: no known UMAP embedding found.")
        return

    sc.pl.embedding(
        adata,
        basis="conference_umap",
        color="clusters_r1p10",
        size=point_size,
        frameon=False,
        legend_loc="right margin",
        title="CK_skin_res BANKSY UMAP, r1.10 clusters",
        show=show,
    )
    save_current_figure(output_dir, "CK_skin_res_umap_clusters_r1p10", dpi)


def plot_marker_spatial_panels(adata, markers, output_dir, point_size, dpi, show):
    """Plot spatial expression panels for selected marker genes."""
    sc.pl.embedding(
        adata,
        basis="spatial",
        color=markers,
        size=point_size,
        frameon=False,
        cmap="viridis",
        ncols=4,
        title=markers,
        show=show,
    )
    save_current_figure(output_dir, "CK_skin_res_spatial_endothelial_markers", dpi)


def plot_marker_dotplot(adata, markers, output_dir, dpi, show):
    """Plot marker expression summarized across r1.10 BANKSY clusters."""
    sc.pl.dotplot(
        adata,
        var_names=markers,
        groupby="clusters_r1p10",
        standard_scale="var",
        cmap="RdBu_r",
        dendrogram=False,
        show=show,
    )
    save_current_figure(output_dir, "CK_skin_res_endothelial_marker_dotplot_r1p10", dpi)


def plot_cluster_abundance(adata, output_dir, dpi):
    """Plot the number of cells assigned to each r1.10 BANKSY cluster."""
    counts = (
        adata.obs["clusters_r1p10"]
        .value_counts()
        .rename_axis("cluster")
        .reset_index(name="n_cells")
    )
    counts["cluster_sort"] = pd.to_numeric(counts["cluster"], errors="coerce")
    counts = counts.sort_values(["cluster_sort", "cluster"], na_position="last")

    fig_width = max(8, 0.35 * len(counts))
    plt.figure(figsize=(fig_width, 4.5))
    plt.bar(counts["cluster"], counts["n_cells"], color="#4c78a8")
    plt.xlabel("BANKSY cluster, r1.10")
    plt.ylabel("Cells")
    plt.title("CK_skin_res r1.10 cluster abundance")
    plt.xticks(rotation=90)
    save_current_figure(output_dir, "CK_skin_res_cluster_abundance_r1p10", dpi)


def main():
    """Run all CK skin r1.10 conference figure exports."""
    args = parse_args()
    load_plotting_stack()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = prepare_adata(Path(args.adata), args.cluster_col)
    markers = available_genes(adata, args.markers)

    print(f"Using {adata.n_obs:,} cells and {adata.n_vars:,} genes.")
    print(f"Plotting markers: {', '.join(markers)}")

    plot_spatial_clusters(adata, output_dir, args.point_size, args.dpi, args.show)
    plot_umap_clusters(adata, output_dir, args.point_size, args.dpi, args.show)
    plot_marker_spatial_panels(adata, markers, output_dir, args.point_size, args.dpi, args.show)
    plot_marker_dotplot(adata, markers, output_dir, args.dpi, args.show)
    plot_cluster_abundance(adata, output_dir, args.dpi)


if __name__ == "__main__":
    main()
