"""
Title: Annotate selected BANKSY clusters and plot shared-colour UMAPs
Date: 2026-10-26
Updated: 2026-11-19 - samples come from a cohort config (xenium_workflow/cohorts.py);
--input-root/--output-root default to the repository root.
Updated: 2026-11-22 - labels come from a versioned annotation set
(annotation_sets.py, config/annotations/): one obs column per tier
(`cell_type_tier_1`, `cell_type_tier_2`, ...), provenance stamp in
`uns['cell_type_annotation']`, derived output paths, a change table against the
parent set, and a run manifest. The old `cell_type` column is retired.
Summary: Validate an annotation set's cluster-to-label mapping against clean
script 06 AnnData objects, attach one label column per tier, and render Scanpy
UMAPs per tier with the set's fixed colours. The expression matrix, cell
identifiers, cluster columns, and saved embeddings are preserved. Figure
formats, typography, titles, legends, and point rasterization are config-driven.

Inputs: step config (`annotation`, plot settings); script 06 clean objects
(`recluster_h5ad` in the cohort) under --input-root, or, for replot configs
(`input_key: annotated_h5ad`), this set's annotated objects.
Outputs (under --output-root, paths from xenium_workflow/annotation_sets.py):
  processed/.../cell_type_annotation/<cohort>/<version>/<sample>/...annotated_<version>_r<res>.h5ad
  processed/.../cell_type_annotation/<cohort>/<version>/cluster_annotation_counts.csv,
      cell_type_counts_<tier>.csv, changes_vs_<parent>.csv, run_manifest.json
  output/.../cell_type_annotation/<cohort>/<version>/<tier>/[<output_subdir>/] UMAP figures
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scanpy as sc

from xenium_workflow.annotation_sets import (
    add_annotation_arguments, annotated_h5ad_path, check_stamp, cluster_changes, inherited_column,
    label_column, level_title, level_vocabulary, load_annotation, load_mapping, mapping_for_sample,
    resolve_step, sibling_dir, stamp, step_dir,
)
from xenium_workflow.cohorts import REPO_ROOT, add_root_arguments, load_json, resolve_path, select_samples


# Keep text editable in Inkscape while allowing dense point clouds to be rasterized.
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

STEP = "cell_type_annotation"


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Apply an annotation set's cluster mapping to clean BANKSY recluster objects "
            "and render Scanpy UMAPs per label tier."
        )
    )
    parser.add_argument("--config", required=True, help="Path to the JSON run config.")
    parser.add_argument(
        "--sample",
        action="append",
        help="Process only this sample; repeat to select multiple samples.",
    )
    # Script 08 writes every tier; --label-level is not used here (see plot_levels).
    add_annotation_arguments(parser, label_level=False)
    # Both default to the repository root. Locally, script 06 inputs are on the
    # Phoenix mount, so pass --input-root ~/mnt/phoenix_scratch/Banksy_py.
    add_root_arguments(parser)
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


def validate_sample_mapping(adata: ad.AnnData, sample_cfg: dict, sample_mapping: pd.DataFrame) -> None:
    """Check that the mapping covers exactly the object's clusters and the UMAP is usable."""
    sample = sample_cfg["sample"]
    groupby = sample_cfg["groupby"]
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
    if adata.obs[groupby].isna().any():
        raise ValueError(f"{sample}: cluster column {groupby!r} has missing values")

    object_clusters = set(adata.obs[groupby].astype(str).unique())
    mapped_clusters = set(sample_mapping["cluster_id"])
    missing = sorted(object_clusters - mapped_clusters)
    extra = sorted(mapped_clusters - object_clusters)
    if missing or extra:
        raise ValueError(
            f"{sample}: mapping coverage mismatch; unmapped={missing}, absent_from_object={extra}"
        )


def apply_annotations(adata: ad.AnnData, sample_cfg: dict, sample_mapping: pd.DataFrame,
                      annotation: dict) -> None:
    """Attach one label column per tier plus provenance, in place."""
    groupby = sample_cfg["groupby"]
    source_clusters = adata.obs[groupby].astype(str)
    by_cluster = sample_mapping.set_index("cluster_id")
    for level in annotation["levels"]:
        column = label_column(level)
        order, colors = level_vocabulary(annotation, level)
        labels = source_clusters.map(by_cluster[column])
        if labels.isna().any():
            missing = sorted(source_clusters[labels.isna()].unique())
            raise ValueError(f"Unmapped clusters remained after validation: {missing}")
        # Keep the set-wide category order even when a label is absent from a sample.
        adata.obs[column] = pd.Categorical(labels, categories=order, ordered=True)
        adata.uns[f"{column}_colors"] = [colors[label] for label in order]
        if level != annotation["levels"][0]:
            adata.obs[inherited_column(level)] = source_clusters.map(by_cluster[inherited_column(level)]).astype(bool)
    adata.obs["cell_type_annotation_status"] = pd.Categorical(
        source_clusters.map(by_cluster["annotation_status"]), categories=["reviewed", "proposed"])
    adata.obs["cell_type_source_cluster"] = source_clusters.astype("category")
    adata.uns["cell_type_annotation"] = {
        **stamp(annotation),
        "branch": sample_cfg["branch"],
        "resolution": str(sample_cfg["resolution"]),
        "groupby": groupby,
        "umap_key": sample_cfg["umap_key"],
        "mapping_file": str(annotation["mapping_path"].relative_to(REPO_ROOT)),
    }
    # Retired 2026-11-22; never leave a stale single-level column behind.
    if "cell_type" in adata.obs:
        del adata.obs["cell_type"]
    adata.uns.pop("cell_type_colors", None)


def shared_legend(order: list[str], palette: dict[str, str]) -> list[Line2D]:
    """Build legend handles in set-wide label order."""
    return [
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor=palette[label],
               markeredgecolor="none", markersize=7, label=label)
        for label in order
    ]


def add_centroid_labels(ax: plt.Axes, coords: np.ndarray, labels: pd.Series,
                        present_types: list[str], fontsize: float) -> None:
    """Place one readable label at each population's median location."""
    label_array = labels.astype(str).to_numpy()
    for cell_type in present_types:
        location = np.nanmedian(coords[label_array == cell_type], axis=0)
        ax.text(
            location[0], location[1], cell_type, fontsize=fontsize, fontweight="bold",
            ha="center", va="center", color="#111111",
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


def plot_matplotlib_umap(ax, adata, sample_cfg, column, order, palette, plot_cfg) -> None:
    """Plot one annotated UMAP directly with Matplotlib."""
    coords = np.asarray(adata.obsm[sample_cfg["umap_key"]])
    labels = adata.obs[column]
    colors = labels.astype(str).map(palette).to_numpy()
    ax.scatter(
        coords[:, 0], coords[:, 1], c=colors, s=float(plot_cfg["point_size"]),
        alpha=float(plot_cfg["alpha"]), linewidths=0,
        rasterized=bool(plot_cfg.get("rasterized_points", True)),
    )
    present = [label for label in order if (labels == label).any()]
    if plot_cfg.get("label_on_data", True):
        add_centroid_labels(ax, coords, labels, present, float(plot_cfg["label_fontsize"]))
    style_umap_axis(ax, sample_cfg["title"], plot_cfg)


def plot_scanpy_umap(ax, adata, sample_cfg, column, order, palette, plot_cfg) -> None:
    """Plot the saved embedding through `scanpy.pl.umap`."""
    had_generic = "X_umap" in adata.obsm
    generic = np.asarray(adata.obsm["X_umap"]).copy() if had_generic else None
    color_key = f"{column}_colors"
    had_colors = color_key in adata.uns
    stored_colors = np.asarray(adata.uns[color_key]).copy() if had_colors else None
    adata.obsm["X_umap"] = np.asarray(adata.obsm[sample_cfg["umap_key"]])
    try:
        sc.pl.umap(
            adata,
            color=column,
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
        # canonical hex codes so plotting cannot mutate saved metadata.
        if had_colors:
            adata.uns[color_key] = stored_colors
        else:
            adata.uns.pop(color_key, None)


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
        plt.close(fig)
        raise FileExistsError(f"Outputs already exist; use --overwrite: {existing}")
    path.parent.mkdir(parents=True, exist_ok=True)
    for output_path in paths:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight", facecolor="white")
        print(f"Wrote {output_path}")
    plt.close(fig)


def legend_title(level: str) -> str:
    return f"Cell type ({level_title(level)})"


def plot_individual_scanpy(adata, sample_cfg, annotation, level, plot_cfg, figure_root, overwrite) -> None:
    """Write a Scanpy UMAP for one sample and tier in every configured format."""
    column = label_column(level)
    order, palette = level_vocabulary(annotation, level)
    fig, ax = plt.subplots(1, 1, figsize=tuple(plot_cfg["individual_figsize"]))
    plot_scanpy_umap(ax, adata, sample_cfg, column, order, palette, plot_cfg)
    if plot_cfg.get("show_legend", True):
        present = [label for label in order if (adata.obs[column] == label).any()]
        fig.legend(handles=shared_legend(present, palette), loc="center left",
                   bbox_to_anchor=(0.985, 0.5), frameon=False, title=legend_title(level))
    if plot_cfg.get("show_figure_title", True):
        fig.suptitle(
            f"{sample_cfg['branch_label']}, r{sample_cfg['resolution']} — "
            f"annotation {annotation['version']}, {level_title(level)}",
            fontsize=float(plot_cfg["suptitle_fontsize"]),
        )
    path = figure_root / "per_sample_scanpy" / f"{sample_cfg['sample']}_{column}_umap_scanpy.png"
    save_figure(fig, path, int(plot_cfg["dpi"]), overwrite, plot_cfg.get("output_formats", ["png"]))


def make_panel_figure(records, annotation, level, plot_cfg, figure_root, overwrite) -> None:
    """Write one multi-sample Scanpy panel for one tier."""
    column = label_column(level)
    order, palette = level_vocabulary(annotation, level)
    ncols = int(plot_cfg["panel_ncols"])
    nrows = math.ceil(len(records) / ncols)
    width, height = plot_cfg["panel_size_per_sample"]
    fig, axes = plt.subplots(nrows, ncols, figsize=(width * ncols, height * nrows), squeeze=False)
    present = set()
    for ax, (sample_cfg, coords, labels) in zip(axes.flat, records):
        plot_adata = ad.AnnData(obs=pd.DataFrame(index=np.arange(len(labels)).astype(str)))
        plot_adata.obs[column] = labels.copy()
        plot_adata.uns[f"{column}_colors"] = [palette[label] for label in order]
        plot_adata.obsm[sample_cfg["umap_key"]] = coords
        plot_scanpy_umap(ax, plot_adata, sample_cfg, column, order, palette, plot_cfg)
        present |= set(pd.Series(labels).astype(str))
    for ax in axes.flat[len(records):]:
        ax.set_visible(False)
    if plot_cfg.get("show_legend", True):
        fig.legend(handles=shared_legend([label for label in order if label in present], palette),
                   loc="center left", bbox_to_anchor=(0.985, 0.5), frameon=False, title=legend_title(level))
    if plot_cfg.get("show_figure_title", True):
        fig.suptitle(
            f"{annotation['project'].upper()} {annotation['cohort']['branch']} — "
            f"annotation {annotation['version']}, {level_title(level)}",
            fontsize=float(plot_cfg["suptitle_fontsize"]),
        )
    path = figure_root / f"combined_scanpy_{column}_umaps.png"
    save_figure(fig, path, int(plot_cfg["dpi"]), overwrite, plot_cfg.get("output_formats", ["png"]))


def build_count_tables(adata, sample_cfg, sample_mapping, annotation) -> tuple[pd.DataFrame, dict]:
    """Per-source-cluster counts (all tiers) and per-tier label counts for one sample."""
    groupby = sample_cfg["groupby"]
    n_by_cluster = adata.obs[groupby].astype(str).value_counts()
    cluster_table = sample_mapping.copy()
    cluster_table["n_cells"] = cluster_table["cluster_id"].map(n_by_cluster)
    if cluster_table["n_cells"].isna().any() or int(cluster_table["n_cells"].sum()) != adata.n_obs:
        raise ValueError(f"{sample_cfg['sample']}: cluster counts did not reconcile with the object")
    cluster_table["n_cells"] = cluster_table["n_cells"].astype(int)
    level_tables = {}
    for level in annotation["levels"]:
        column = label_column(level)
        table = adata.obs.groupby(column, observed=True).size().rename("n_cells").reset_index()
        table = table.rename(columns={column: "cell_type"})
        table.insert(0, "label_level", level)
        table.insert(0, "resolution", str(sample_cfg["resolution"]))
        table.insert(0, "branch", sample_cfg["branch"])
        table.insert(0, "sample", sample_cfg["sample"])
        if int(table["n_cells"].sum()) != adata.n_obs:
            raise ValueError(f"{sample_cfg['sample']}: {column} counts do not sum to n_obs")
        level_tables[level] = table
    return cluster_table, level_tables


def write_annotated_object(adata: ad.AnnData, path: Path, compression: str | None, overwrite: bool) -> None:
    """Write an annotated object separately from its source object."""
    if path.exists() and not overwrite:
        raise FileExistsError(f"Annotated object already exists; use --overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path, compression=compression)
    print(f"Wrote {path}")


def write_table(frame: pd.DataFrame, path: Path, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"Summary already exists; use --overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)
    print(f"Wrote {path}")


def change_table(annotation: dict, cluster_counts: pd.DataFrame) -> pd.DataFrame | None:
    """Clusters whose labels changed since the parent set, with their cell counts."""
    if not annotation["parent"]:
        return None
    parent = load_annotation(sibling_dir(annotation, annotation["parent"]))
    changes = cluster_changes(annotation, parent)
    counts = cluster_counts[["sample", "cluster_id", "n_cells"]]
    return changes.merge(counts, on=["sample", "cluster_id"], how="left")


def main() -> None:
    """Validate the set, annotate objects, and create per-tier UMAP figures."""
    args = parse_args()
    cfg = load_json(Path(args.config).resolve())
    input_root = args.input_root.resolve()
    output_root = args.output_root.resolve()
    annotation, _ = resolve_step(cfg, args, use_level=False)
    mapping = load_mapping(annotation)
    samples = select_samples(annotation["cohort"], cfg.get("samples"), args.sample)
    plot_levels = cfg.get("plot_levels", annotation["levels"])
    unknown_levels = sorted(set(plot_levels) - set(annotation["levels"]))
    if unknown_levels:
        raise ValueError(f"plot_levels {unknown_levels} are not in set {annotation['version']}")
    # Annotation reads the script 06 object; replot-only configs read this set's annotated one.
    input_key = cfg.get("input_key", "recluster_h5ad")
    if input_key not in ("recluster_h5ad", "annotated_h5ad"):
        raise ValueError("input_key must be recluster_h5ad or annotated_h5ad")
    replot = input_key == "annotated_h5ad"
    write_objects = cfg.get("write_annotated_objects", True) and not args.no_write_annotated
    if replot and cfg.get("write_annotated_objects", True):
        raise ValueError("Replotting annotated objects requires write_annotated_objects: false")
    table_root = resolve_path(output_root, step_dir(annotation, "processed", STEP))
    plot_cfg = cfg["plot"]
    print(f"Annotation set {annotation['version']} ({annotation['status']}), cohort "
          f"{annotation['cohort_name']}, levels {annotation['levels']}, plotting {plot_levels}")

    cluster_tables, level_tables, panel_records, written = [], {lvl: [] for lvl in annotation["levels"]}, \
        {lvl: [] for lvl in plot_levels}, []
    for sample_cfg in samples:
        sample = sample_cfg["sample"]
        if replot:
            input_path = resolve_path(input_root, annotated_h5ad_path(annotation, sample_cfg))
        else:
            input_path = resolve_path(input_root, sample_cfg["recluster_h5ad"])
        if not input_path.exists():
            raise FileNotFoundError(f"{sample}: input object does not exist: {input_path}")
        print(f"Reading {sample}: {input_path}")
        adata = ad.read_h5ad(input_path)
        original_shape = adata.shape
        original_obs_names = adata.obs_names.copy()
        sample_mapping = mapping_for_sample(mapping, sample_cfg)
        validate_sample_mapping(adata, sample_cfg, sample_mapping)
        if replot:
            check_stamp(dict(adata.uns.get("cell_type_annotation", {})), annotation, str(input_path))
        print(f"Validated {sample}: {adata.n_obs:,} cells, {len(sample_mapping)} clusters, "
              f"r{sample_cfg['resolution']}")
        if args.validate_only:
            continue

        if not replot:
            apply_annotations(adata, sample_cfg, sample_mapping, annotation)
        if adata.shape != original_shape or not adata.obs_names.equals(original_obs_names):
            raise RuntimeError(f"{sample}: annotation changed shape or cell identifier order")
        cluster_table, per_level = build_count_tables(adata, sample_cfg, sample_mapping, annotation)
        cluster_tables.append(cluster_table)
        for level, table in per_level.items():
            level_tables[level].append(table)
        for level in plot_levels:
            figure_root = resolve_path(output_root, step_dir(annotation, "output", STEP, level,
                                                             cfg.get("output_subdir")))
            plot_individual_scanpy(adata, sample_cfg, annotation, level, plot_cfg, figure_root, args.overwrite)
            if plot_cfg.get("make_combined_figure", True):
                panel_records[level].append((
                    sample_cfg,
                    np.asarray(adata.obsm[sample_cfg["umap_key"]]).copy(),
                    adata.obs[label_column(level)].values.copy(),
                ))
        if write_objects and not replot:
            output_h5ad = resolve_path(output_root, annotated_h5ad_path(annotation, sample_cfg))
            write_annotated_object(adata, output_h5ad, cfg.get("h5ad_compression", "gzip"), args.overwrite)
            written.append({"sample": sample, "n_obs": int(adata.n_obs), "path": str(output_h5ad)})
        del adata

    if args.validate_only:
        print(f"Validation passed for {len(samples)} sample(s); no outputs written.")
        return

    if cfg.get("write_summary_tables", True):
        cluster_counts = pd.concat(cluster_tables, ignore_index=True)
        write_table(cluster_counts, table_root / "cluster_annotation_counts.csv", args.overwrite)
        for level, tables in level_tables.items():
            write_table(pd.concat(tables, ignore_index=True),
                        table_root / f"cell_type_counts_{level}.csv", args.overwrite)
        changes = change_table(annotation, cluster_counts)
        if changes is not None:
            write_table(changes, table_root / f"changes_vs_{annotation['parent']}.csv", args.overwrite)
            print(f"{len(changes)} clusters changed label since {annotation['parent']} "
                  f"({int(changes['n_cells'].sum()):,} cells in processed samples)")
        manifest = {
            **stamp(annotation),
            "run_at": datetime.now().isoformat(timespec="seconds"),
            "config": str(Path(args.config)),
            "input_root": str(input_root),
            "samples": [s["sample"] for s in samples],
            "total_cells": int(cluster_counts["n_cells"].sum()),
            "annotated_h5ads": written,
        }
        manifest_path = table_root / "run_manifest.json"
        if manifest_path.exists() and not args.overwrite:
            raise FileExistsError(f"Manifest already exists; use --overwrite: {manifest_path}")
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        print(f"Wrote {manifest_path}")
    if plot_cfg.get("make_combined_figure", True):
        for level in plot_levels:
            figure_root = resolve_path(output_root, step_dir(annotation, "output", STEP, level,
                                                             cfg.get("output_subdir")))
            make_panel_figure(panel_records[level], annotation, level, plot_cfg, figure_root, args.overwrite)


if __name__ == "__main__":
    main()
