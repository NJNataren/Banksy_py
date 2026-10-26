"""
Title: Annotate selected BANKSY clusters and plot shared-colour UMAPs
Date: 2026-10-26
Summary: Validate a reviewed cluster-to-cell-type mapping, attach annotations to
clean script 06 AnnData objects, and render matching Matplotlib and Scanpy UMAPs.
The original expression matrix, cell identifiers, cluster columns, and saved
embeddings are preserved. Figure formats, typography, titles, legends, and point
rasterization are config-driven.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scanpy as sc


# Keep text editable in Inkscape while allowing dense point clouds to be rasterized.
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Apply reviewed cell-type mappings to clean BANKSY recluster objects "
            "and render Scanpy cell-type UMAPs."
        )
    )
    parser.add_argument("--config", required=True, help="Path to the JSON run config.")
    parser.add_argument(
        "--sample",
        action="append",
        help="Process only this sample; repeat to select multiple samples.",
    )
    parser.add_argument(
        "--input-root",
        help="Override the config input root, for example an HPC repository root.",
    )
    parser.add_argument(
        "--output-root",
        help="Override the config output root. Relative output paths use this root.",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate mappings, cluster coverage, and embeddings without writing outputs.",
    )
    parser.add_argument(
        "--no-write-annotated",
        action="store_true",
        help="Create plots and summaries without writing annotated H5AD objects.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing annotated H5AD and figure outputs.",
    )
    return parser.parse_args()


def load_json(path: Path) -> dict:
    """Read a JSON object from disk."""
    with path.open() as handle:
        return json.load(handle)


def resolve_path(root: Path, value: str) -> Path:
    """Resolve an absolute path or a path relative to `root`."""
    path = Path(value)
    return path if path.is_absolute() else root / path


def load_palette(path: Path) -> tuple[list[str], dict[str, str]]:
    """Load and validate an ordered shared cell-type palette."""
    payload = load_json(path)
    order = payload["cell_type_order"]
    colors = payload["colors"]
    if len(order) != len(set(order)):
        raise ValueError("cell_type_order contains duplicate labels")
    if set(order) != set(colors):
        raise ValueError("Palette labels and cell_type_order do not match exactly")
    for label, color in colors.items():
        if not matplotlib.colors.is_color_like(color):
            raise ValueError(f"Invalid colour for {label!r}: {color!r}")
    return order, colors


def load_mapping(path: Path) -> pd.DataFrame:
    """Load the reviewed cluster annotation table with cluster IDs as strings."""
    mapping = pd.read_csv(path, dtype={"sample": str, "resolution": str, "cluster_id": str})
    required = {
        "sample",
        "branch",
        "resolution",
        "groupby",
        "cluster_id",
        "cell_type",
        "annotation_status",
    }
    missing = required - set(mapping.columns)
    if missing:
        raise KeyError(f"Mapping table is missing columns: {sorted(missing)}")
    duplicate = mapping.duplicated(["sample", "branch", "resolution", "groupby", "cluster_id"])
    if duplicate.any():
        rows = mapping.loc[duplicate, sorted(required)].to_dict("records")
        raise ValueError(f"Duplicate cluster mappings found: {rows}")
    if mapping[list(required)].isna().any().any():
        raise ValueError("Mapping table contains missing values in required columns")
    return mapping


def select_samples(cfg: dict, requested: list[str] | None) -> list[dict]:
    """Return configured samples, optionally restricted by command-line selection."""
    samples = cfg["samples"]
    configured = [sample["sample"] for sample in samples]
    if len(configured) != len(set(configured)):
        raise ValueError("The run config contains duplicate sample names")
    if not requested:
        return samples
    unknown = sorted(set(requested) - set(configured))
    if unknown:
        raise ValueError(f"Requested samples are not configured: {unknown}")
    requested_set = set(requested)
    return [sample for sample in samples if sample["sample"] in requested_set]


def validate_sample_mapping(
    adata: ad.AnnData,
    sample_cfg: dict,
    mapping: pd.DataFrame,
    cell_type_order: list[str],
) -> pd.DataFrame:
    """Validate one mapping against the selected object and cluster grouping."""
    sample = sample_cfg["sample"]
    groupby = sample_cfg["groupby"]
    resolution = str(sample_cfg["resolution"])
    branch = sample_cfg["branch"]

    if groupby not in adata.obs:
        raise KeyError(f"{sample}: missing cluster column {groupby!r}")
    embedding_key = sample_cfg["umap_key"]
    if embedding_key not in adata.obsm:
        raise KeyError(f"{sample}: missing UMAP embedding adata.obsm[{embedding_key!r}]")
    embedding = np.asarray(adata.obsm[embedding_key])
    if embedding.shape != (adata.n_obs, 2):
        raise ValueError(f"{sample}: expected a two-dimensional UMAP; found {embedding.shape}")
    if not np.isfinite(embedding).all():
        raise ValueError(f"{sample}: UMAP coordinates contain non-finite values")

    selected = mapping[
        (mapping["sample"] == sample)
        & (mapping["branch"] == branch)
        & (mapping["resolution"].astype(float) == float(resolution))
        & (mapping["groupby"] == groupby)
    ].copy()
    if selected.empty:
        raise ValueError(f"{sample}: no mapping rows match branch/resolution/groupby")

    object_clusters = set(adata.obs[groupby].dropna().astype(str).unique())
    mapped_clusters = set(selected["cluster_id"].astype(str))
    missing = sorted(object_clusters - mapped_clusters)
    extra = sorted(mapped_clusters - object_clusters)
    if missing or extra:
        raise ValueError(
            f"{sample}: mapping coverage mismatch; unmapped={missing}, absent_from_object={extra}"
        )
    unknown_types = sorted(set(selected["cell_type"]) - set(cell_type_order))
    if unknown_types:
        raise ValueError(f"{sample}: cell types missing from shared palette: {unknown_types}")
    if set(selected["annotation_status"]) != {"reviewed"}:
        raise ValueError(f"{sample}: all applied annotations must have status 'reviewed'")
    return selected


def apply_annotations(
    adata: ad.AnnData,
    sample_cfg: dict,
    sample_mapping: pd.DataFrame,
    cell_type_order: list[str],
    palette: dict[str, str],
    mapping_path: Path,
) -> None:
    """Attach reviewed labels and reproducible annotation metadata in place."""
    groupby = sample_cfg["groupby"]
    cluster_to_type = sample_mapping.set_index("cluster_id")["cell_type"].to_dict()
    source_clusters = adata.obs[groupby].astype(str)
    labels = source_clusters.map(cluster_to_type)
    if labels.isna().any():
        missing = sorted(source_clusters[labels.isna()].unique())
        raise ValueError(f"Unmapped clusters remained after validation: {missing}")

    # Keep the project-wide category order even when a type is absent from a sample.
    adata.obs["cell_type"] = pd.Categorical(labels, categories=cell_type_order, ordered=True)
    adata.obs["cell_type_annotation_status"] = pd.Categorical(
        ["reviewed"] * adata.n_obs,
        categories=["reviewed", "unresolved"],
    )
    adata.obs["cell_type_source_cluster"] = source_clusters.astype("category")
    adata.uns["cell_type_colors"] = [palette[label] for label in cell_type_order]
    adata.uns["cell_type_annotation"] = {
        "branch": sample_cfg["branch"],
        "resolution": str(sample_cfg["resolution"]),
        "groupby": groupby,
        "umap_key": sample_cfg["umap_key"],
        "mapping_file": str(mapping_path),
        "cell_type_order": cell_type_order,
    }


def shared_legend(cell_type_order: list[str], palette: dict[str, str]) -> list[Line2D]:
    """Build legend handles in project-wide cell-type order."""
    return [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=palette[label],
            markeredgecolor="none",
            markersize=7,
            label=label,
        )
        for label in cell_type_order
    ]


def add_centroid_labels(
    ax: plt.Axes,
    coords: np.ndarray,
    labels: pd.Series,
    present_types: list[str],
    fontsize: float,
) -> None:
    """Place one readable cell-type label at each population's median location."""
    label_array = labels.astype(str).to_numpy()
    for cell_type in present_types:
        location = np.nanmedian(coords[label_array == cell_type], axis=0)
        ax.text(
            location[0],
            location[1],
            cell_type,
            fontsize=fontsize,
            fontweight="bold",
            ha="center",
            va="center",
            color="#111111",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 1.5},
            zorder=3,
        )


def style_umap_axis(ax: plt.Axes, title: str, plot_cfg: dict) -> None:
    """Apply config-driven title and axis styling to a UMAP panel."""
    ax.set_title(title if plot_cfg.get("show_sample_titles", True) else "")
    if plot_cfg.get("show_axis_labels", True):
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
    else:
        ax.set_xlabel("")
        ax.set_ylabel("")
    ax.set_aspect("equal", adjustable="datalim")
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)


def plot_matplotlib_umap(
    ax: plt.Axes,
    adata: ad.AnnData,
    sample_cfg: dict,
    cell_type_order: list[str],
    palette: dict[str, str],
    plot_cfg: dict,
) -> None:
    """Plot one annotated UMAP directly with Matplotlib."""
    coords = np.asarray(adata.obsm[sample_cfg["umap_key"]])
    labels = adata.obs["cell_type"]
    colors = labels.astype(str).map(palette).to_numpy()
    ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=colors,
        s=float(plot_cfg["point_size"]),
        alpha=float(plot_cfg["alpha"]),
        linewidths=0,
        rasterized=bool(plot_cfg.get("rasterized_points", True)),
    )
    present = [label for label in cell_type_order if (labels == label).any()]
    if plot_cfg.get("label_on_data", True):
        add_centroid_labels(ax, coords, labels, present, float(plot_cfg["label_fontsize"]))
    style_umap_axis(ax, sample_cfg["title"], plot_cfg)


def plot_scanpy_umap(
    ax: plt.Axes,
    adata: ad.AnnData,
    sample_cfg: dict,
    cell_type_order: list[str],
    palette: dict[str, str],
    plot_cfg: dict,
) -> None:
    """Plot the same saved embedding through `scanpy.pl.umap`."""
    had_generic = "X_umap" in adata.obsm
    generic = np.asarray(adata.obsm["X_umap"]).copy() if had_generic else None
    had_colors = "cell_type_colors" in adata.uns
    stored_colors = np.asarray(adata.uns["cell_type_colors"]).copy() if had_colors else None
    adata.obsm["X_umap"] = np.asarray(adata.obsm[sample_cfg["umap_key"]])
    try:
        sc.pl.umap(
            adata,
            color="cell_type",
            palette=palette,
            size=float(plot_cfg["point_size"]),
            alpha=float(plot_cfg["alpha"]),
            legend_loc="on data" if plot_cfg.get("label_on_data", True) else None,
            legend_fontsize=float(plot_cfg["label_fontsize"]),
            legend_fontweight="bold",
            frameon=False,
            title=sample_cfg["title"] if plot_cfg.get("show_sample_titles", True) else "",
            ax=ax,
            show=False,
        )
        for collection in ax.collections:
            collection.set_rasterized(bool(plot_cfg.get("rasterized_points", True)))
        style_umap_axis(ax, sample_cfg["title"], plot_cfg)
    finally:
        if had_generic:
            adata.obsm["X_umap"] = generic
        else:
            del adata.obsm["X_umap"]
        # Scanpy normalizes palette strings to RGBA values in `.uns`. Restore the
        # canonical project hex codes so plotting cannot mutate saved metadata.
        if had_colors:
            adata.uns["cell_type_colors"] = stored_colors
        else:
            adata.uns.pop("cell_type_colors", None)


def save_figure(fig: plt.Figure, path: Path, dpi: int, overwrite: bool, formats: list[str]) -> None:
    """Save a figure in every configured format with explicit overwrite checks."""
    supported = {"png", "svg", "pdf"}
    normalized = [str(value).lower().lstrip(".") for value in formats]
    unknown = sorted(set(normalized) - supported)
    if unknown:
        raise ValueError(f"Unsupported figure formats: {unknown}")
    if len(normalized) != len(set(normalized)):
        raise ValueError(f"Duplicate figure formats requested: {normalized}")
    paths = [path.with_suffix(f".{extension}") for extension in normalized]
    existing = [output_path for output_path in paths if output_path.exists()]
    if existing and not overwrite:
        raise FileExistsError(f"Outputs already exist; use --overwrite: {existing}")
    path.parent.mkdir(parents=True, exist_ok=True)
    for output_path in paths:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight", facecolor="white")
        print(f"Wrote {output_path}")
    plt.close(fig)


def plot_individual_scanpy(
    adata: ad.AnnData,
    sample_cfg: dict,
    cell_type_order: list[str],
    palette: dict[str, str],
    plot_cfg: dict,
    figure_root: Path,
    overwrite: bool,
) -> None:
    """Write a Scanpy UMAP for one sample in every configured format."""
    fig, ax = plt.subplots(1, 1, figsize=tuple(plot_cfg["individual_figsize"]))
    plot_scanpy_umap(ax, adata, sample_cfg, cell_type_order, palette, plot_cfg)
    if plot_cfg.get("show_legend", True):
        fig.legend(
            handles=shared_legend(cell_type_order, palette),
            loc="center left",
            bbox_to_anchor=(0.985, 0.5),
            frameon=False,
            title="Cell type",
        )
    if plot_cfg.get("show_figure_title", True):
        fig.suptitle(
            f"Shared cell-type colours — {sample_cfg['branch_label']}, r{sample_cfg['resolution']}",
            fontsize=float(plot_cfg["suptitle_fontsize"]),
        )
    path = figure_root / "per_sample_scanpy" / f"{sample_cfg['sample']}_cell_type_umap_scanpy.png"
    save_figure(fig, path, int(plot_cfg["dpi"]), overwrite, plot_cfg.get("output_formats", ["png"]))


def make_panel_figure(
    records: list[tuple[dict, np.ndarray, pd.Categorical]],
    method: str,
    cell_type_order: list[str],
    palette: dict[str, str],
    plot_cfg: dict,
    figure_root: Path,
    overwrite: bool,
) -> None:
    """Write one multi-sample panel using Matplotlib or Scanpy rendering."""
    ncols = int(plot_cfg["panel_ncols"])
    nrows = math.ceil(len(records) / ncols)
    width, height = plot_cfg["panel_size_per_sample"]
    fig, axes = plt.subplots(nrows, ncols, figsize=(width * ncols, height * nrows), squeeze=False)
    for ax, (sample_cfg, coords, labels) in zip(axes.flat, records):
        plot_adata = ad.AnnData(obs=pd.DataFrame(index=np.arange(len(labels))))
        plot_adata.obs["cell_type"] = labels.copy()
        plot_adata.uns["cell_type_colors"] = [palette[label] for label in cell_type_order]
        plot_adata.obsm[sample_cfg["umap_key"]] = coords
        if method == "matplotlib":
            plot_matplotlib_umap(ax, plot_adata, sample_cfg, cell_type_order, palette, plot_cfg)
        elif method == "scanpy":
            plot_scanpy_umap(ax, plot_adata, sample_cfg, cell_type_order, palette, plot_cfg)
        else:
            raise ValueError(f"Unknown plotting method: {method}")
    for ax in axes.flat[len(records):]:
        ax.set_visible(False)
    if plot_cfg.get("show_legend", True):
        fig.legend(
            handles=shared_legend(cell_type_order, palette),
            loc="center left",
            bbox_to_anchor=(0.985, 0.5),
            frameon=False,
            title="Cell type",
        )
    if plot_cfg.get("show_figure_title", True):
        method_title = "Matplotlib" if method == "matplotlib" else "Scanpy sc.pl.umap"
        fig.suptitle(
            f"VBCT min15/area95 reviewed cell types — {method_title}",
            fontsize=float(plot_cfg["suptitle_fontsize"]),
        )
    path = figure_root / f"combined_{method}_cell_type_umaps.png"
    save_figure(fig, path, int(plot_cfg["dpi"]), overwrite, plot_cfg.get("output_formats", ["png"]))


def build_count_tables(
    adata: ad.AnnData,
    sample_cfg: dict,
    sample_mapping: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Create per-source-cluster and merged-cell-type count tables."""
    groupby = sample_cfg["groupby"]
    counts = (
        adata.obs.groupby([groupby, "cell_type"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .rename(columns={groupby: "cluster_id"})
    )
    counts["cluster_id"] = counts["cluster_id"].astype(str)
    cluster_table = sample_mapping.merge(counts, on=["cluster_id", "cell_type"], how="left")
    if cluster_table["n_cells"].isna().any():
        raise ValueError(f"{sample_cfg['sample']}: count table did not reconcile with mapping")
    type_table = (
        adata.obs.groupby("cell_type", observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    type_table.insert(0, "sample", sample_cfg["sample"])
    type_table.insert(1, "branch", sample_cfg["branch"])
    type_table.insert(2, "resolution", str(sample_cfg["resolution"]))
    if int(type_table["n_cells"].sum()) != adata.n_obs:
        raise ValueError(f"{sample_cfg['sample']}: cell-type counts do not sum to n_obs")
    return cluster_table, type_table


def write_annotated_object(
    adata: ad.AnnData,
    path: Path,
    compression: str | None,
    overwrite: bool,
) -> None:
    """Write an annotated object separately from its source object."""
    if path.exists() and not overwrite:
        raise FileExistsError(f"Annotated object already exists; use --overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path, compression=compression)
    print(f"Wrote {path}")


def main() -> None:
    """Validate mappings, annotate objects, and create comparison UMAP figures."""
    args = parse_args()
    config_path = Path(args.config).resolve()
    cfg = load_json(config_path)
    project_root = config_path.parents[int(cfg.get("project_root_parent_depth", 3))]
    input_root = Path(args.input_root or cfg.get("input_root", project_root)).resolve()
    output_root = Path(args.output_root or cfg.get("output_root", project_root)).resolve()
    mapping_path = resolve_path(project_root, cfg["mapping_csv"])
    palette_path = resolve_path(project_root, cfg["palette_json"])
    mapping = load_mapping(mapping_path)
    cell_type_order, palette = load_palette(palette_path)
    samples = select_samples(cfg, args.sample)
    figure_root = resolve_path(output_root, cfg["figure_dir"])
    table_root = resolve_path(output_root, cfg["summary_dir"])
    plot_cfg = cfg["plot"]

    cluster_tables = []
    type_tables = []
    panel_records = []
    for sample_cfg in samples:
        sample = sample_cfg["sample"]
        input_path = resolve_path(input_root, sample_cfg["input_h5ad"])
        if not input_path.exists():
            raise FileNotFoundError(f"{sample}: input object does not exist: {input_path}")
        print(f"Reading {sample}: {input_path}")
        adata = ad.read_h5ad(input_path)
        original_shape = adata.shape
        original_obs_names = adata.obs_names.copy()
        sample_mapping = validate_sample_mapping(
            adata, sample_cfg, mapping, cell_type_order
        )
        print(
            f"Validated {sample}: {adata.n_obs:,} cells, "
            f"{len(sample_mapping)} clusters, r{sample_cfg['resolution']}"
        )
        if args.validate_only:
            continue

        apply_annotations(
            adata,
            sample_cfg,
            sample_mapping,
            cell_type_order,
            palette,
            mapping_path,
        )
        if adata.shape != original_shape or not adata.obs_names.equals(original_obs_names):
            raise RuntimeError(f"{sample}: annotation changed shape or cell identifier order")
        cluster_table, type_table = build_count_tables(adata, sample_cfg, sample_mapping)
        cluster_tables.append(cluster_table)
        type_tables.append(type_table)
        plot_individual_scanpy(
            adata,
            sample_cfg,
            cell_type_order,
            palette,
            plot_cfg,
            figure_root,
            args.overwrite,
        )
        if plot_cfg.get("make_combined_figure", True):
            panel_records.append(
                (
                    sample_cfg,
                    np.asarray(adata.obsm[sample_cfg["umap_key"]]).copy(),
                    adata.obs["cell_type"].values.copy(),
                )
            )
        if cfg.get("write_annotated_objects", True) and not args.no_write_annotated:
            output_h5ad = resolve_path(output_root, sample_cfg["output_h5ad"])
            write_annotated_object(
                adata,
                output_h5ad,
                cfg.get("h5ad_compression", "gzip"),
                args.overwrite,
            )
        del adata

    if args.validate_only:
        print(f"Validation passed for {len(samples)} sample(s); no outputs written.")
        return

    if cfg.get("write_summary_tables", True):
        table_root.mkdir(parents=True, exist_ok=True)
        cluster_counts_path = table_root / "reviewed_cluster_to_cell_type_counts.csv"
        type_counts_path = table_root / "reviewed_cell_type_counts.csv"
        for summary_path in (cluster_counts_path, type_counts_path):
            if summary_path.exists() and not args.overwrite:
                raise FileExistsError(f"Summary already exists; use --overwrite: {summary_path}")
        pd.concat(cluster_tables, ignore_index=True).to_csv(cluster_counts_path, index=False)
        pd.concat(type_tables, ignore_index=True).to_csv(type_counts_path, index=False)
        print(f"Wrote {cluster_counts_path}")
        print(f"Wrote {type_counts_path}")
    if plot_cfg.get("make_combined_figure", True):
        make_panel_figure(
            panel_records,
            "scanpy",
            cell_type_order,
            palette,
            plot_cfg,
            figure_root,
            args.overwrite,
        )


if __name__ == "__main__":
    main()
