#!/usr/bin/env python
# coding: utf-8

"""
Title: Squidpy Spatial Analysis Of Reviewed Cell Types
Date: 2026-10-26 (renamed from 10_squidpy_recluster_spatial_analysis.py 2026-11-20)
Updated: 2026-11-19 - samples, palette, and mapping come from a cohort config
(xenium_workflow/cohorts.py); --input-root/--output-root default to the repository root.
Updated: 2026-11-20 - also writes <stem>_nhood_enrichment.csv beside each result
H5AD (xenium_workflow/spatial_tables.py) for the cross-sample comparison helper. Cell-type mode
writes centrality and Moran CSVs to <result_dir>/<sample>/<branch>_r<res>/,
not the figure folder. Moran's I permutations are seeded. The legacy script 06
cluster mode was archived (archive/2026-11-20_script10_legacy_cluster_mode/).
Updated: 2026-11-22 - reads a released annotation set (xenium_workflow/annotation_sets.py) and
one label level (`label_level`, e.g. tier_1 -> obs `cell_type_tier_1`); output
paths are derived from the set and level; the H5AD stamp is checked.
Summary: Reuse script 08 annotated clean expression objects and reviewed cell
types for Squidpy neighbourhood enrichment, centrality, co-occurrence, and
Moran's I, one cohort sample at a time. Export PNG, SVG, and PDF figures with editable vector text and 300-DPI raster
layers. The workflow is compatible with Squidpy 1.6 spatial-scatter plotting.
"""

import argparse
import json
import os
import random
from datetime import datetime
from matplotlib import rc_context
from matplotlib.collections import PathCollection

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import squidpy as sq

from xenium_workflow.spatial_tables import enrichment_from_labels, enrichment_table_path, write_enrichment_table
from xenium_workflow.annotation_sets import (
    add_annotation_arguments, annotated_h5ad_path, check_stamp, label_column, level_vocabulary,
    load_mapping, mapping_for_sample, require_released, resolve_step, step_dir,
)
from xenium_workflow.cohorts import add_root_arguments, resolve_path, select_samples

STEP = "spatial_analysis"


sc.logging.print_header()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 1
plt.rcParams["font.family"] = "Arial"
sns.set_style("white")

seed = 1234
np.random.seed(seed)
random.seed(seed)


def parse_args():
    """Parse command-line arguments for the cell-type spatial analysis."""
    parser = argparse.ArgumentParser(
        prog="run Squidpy spatial analyses on script 08 reviewed cell types"
    )
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Script 10 JSON config naming an annotation set, label level, and analysis settings.",
    )
    add_annotation_arguments(parser)
    parser.add_argument("--allow-draft", action="store_true",
                        help="Run on a draft annotation set (local tests only)")
    parser.add_argument("--sample", action="append", help="Select a script 08 sample; repeat as needed")
    parser.add_argument("--validate-only", action="store_true", help="Check annotated inputs without running analyses")
    # Annotated inputs resolve under --input-root; figures and results under --output-root.
    add_root_arguments(parser)
    return parser.parse_args()


def log_time(step):
    """Print a timestamped workflow message."""
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {step}")


def ensure_directory(path, label):
    """Create `path` and log the resolved location."""
    os.makedirs(path, exist_ok=True)
    print(f"{label} directory ready: {os.path.abspath(path)}")


def prepare_spatial_coordinates(adata, coord_keys):
    """Store Xenium x/y coordinates in `.obsm['xy']` and `.obsm['spatial']`."""
    for col in [coord_keys[0], coord_keys[1]]:
        if col not in adata.obs.columns:
            raise KeyError(f"Missing coordinate column in adata.obs: {col}")

    coords = adata.obs[[coord_keys[0], coord_keys[1]]].to_numpy()
    adata.obsm["xy"] = coords
    adata.obsm["spatial"] = coords


def figure_export_settings(cfg):
    """Validate export settings before running expensive spatial calculations."""
    formats = cfg.get("figure_formats", ["png", "svg", "pdf"])
    if (not isinstance(formats, list) or not formats
            or any(fmt not in ("png", "svg", "pdf") for fmt in formats)
            or len(set(formats)) != len(formats)):
        raise ValueError("figure_formats must be a unique nonempty list of png, svg, pdf")
    dpi = cfg.get("figure_dpi", 300)
    if isinstance(dpi, bool) or not isinstance(dpi, (int, float)) or not np.isfinite(dpi) or dpi <= 0:
        raise ValueError("figure_dpi must be a finite positive number")
    rasterize = cfg.get("rasterize_spatial_points", True)
    if not isinstance(rasterize, bool):
        raise ValueError("rasterize_spatial_points must be a boolean")
    return formats, dpi, rasterize


def save_spatial_figure(fig, output_stem, cfg, spatial_points=False, transparent=True):
    """Save all requested formats from one figure and close it even on failure.

    Only spatial cell-point collections are rasterized; text, axes, centrality
    markers, and co-occurrence curves remain vector elements in SVG/PDF.
    """
    try:
        formats, dpi, rasterize = figure_export_settings(cfg)
        if spatial_points and rasterize:
            for ax in fig.axes:
                for collection in ax.collections:
                    if isinstance(collection, PathCollection):
                        collection.set_rasterized(True)
        # SVG keeps text nodes; PDF embeds TrueType fonts instead of Type 3 glyphs.
        with rc_context({"svg.fonttype": "none", "pdf.fonttype": 42}):
            for fmt in formats:
                path = f"{output_stem}.{fmt}"
                fig.savefig(path, dpi=dpi, bbox_inches="tight", transparent=transparent)
                print(f"Saved spatial figure: {path}")
    finally:
        plt.close(fig)


def run_squidpy_for_resolution(
    adata,
    cluster_key,
    output_dir,
    dataset_name,
    pc_label,
    lambda_label,
    res_str_single,
    cfg,
    table_dir,
):
    """Run spatial analyses and export each plot in the configured formats.

    Figures go to `<output_dir>/squidpy`; CSV tables go to `table_dir`
    (beside the result H5AD).
    """
    figure_export_settings(cfg)
    if cluster_key not in adata.obs.columns:
        raise KeyError(
            f"{cluster_key!r} was not found in adata.obs. "
            f"Available obs columns: {list(adata.obs.columns)}"
        )

    if not isinstance(adata.obs[cluster_key].dtype, pd.CategoricalDtype):
        adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")

    # Squidpy plotting reads this colour vector; the caller sets it from the
    # shared palette so every sample uses the same cell-type colours.
    color_key = f"{cluster_key}_colors"
    if color_key not in adata.uns:
        raise KeyError(f"Missing {color_key!r}; set it from the shared palette")

    squidpy_path = os.path.join(output_dir, "squidpy")
    ensure_directory(squidpy_path, "Squidpy")
    sc.settings.figdir = squidpy_path
    table_path = str(table_dir)
    ensure_directory(table_path, "Squidpy tables")

    # All Squidpy analyses use the script 08 annotated clean expression object,
    # which holds only qc_pass_only cells, so QC-failed cells do not contribute
    # to any spatial graph, neighbourhood summary, or expression autocorrelation.
    log_time(f"Running Squidpy spatial graph for {cluster_key}")
    sq.gr.spatial_neighbors(adata, coord_type="generic")

    log_time(f"Running neighbourhood enrichment for {cluster_key}")
    # Permutation z-scores vary between runs unless the permutation seed is fixed.
    sq.gr.nhood_enrichment(adata, cluster_key=cluster_key, seed=seed)
    with rc_context({"figure.figsize": (10, 10)}):
        sq.pl.nhood_enrichment(
            adata,
            cluster_key=cluster_key,
        )
        save_spatial_figure(
            plt.gcf(), os.path.join(squidpy_path,
                f"neighbourhood_enrichment_{dataset_name}"
                f"_pc{pc_label}_nc{lambda_label}_r{res_str_single}"
            ), cfg,
        )

    log_time(f"Running centrality scores for {cluster_key}")
    sq.gr.centrality_scores(adata, cluster_key=cluster_key)
    centrality_key = f"{cluster_key}_centrality_scores"
    if centrality_key in adata.uns:
        centrality_csv = os.path.join(
            table_path,
            f"centrality_scores_{dataset_name}_pc{pc_label}_nc{lambda_label}_r{res_str_single}.csv",
        )
        adata.uns[centrality_key].to_csv(centrality_csv)
        print(f"Saved centrality scores to: {centrality_csv}")

    with rc_context({"figure.figsize": (22, 4)}):
        sq.pl.centrality_scores(adata, cluster_key=cluster_key)
        fig = plt.gcf()
        centrality_stem = os.path.join(
            squidpy_path,
            f"centrality_scores_{dataset_name}_pc{pc_label}_nc{lambda_label}_r{res_str_single}",
        )
        save_spatial_figure(fig, centrality_stem, cfg, transparent=False)

    log_time(f"Running co-occurrence for {cluster_key}")
    sq.gr.co_occurrence(adata, cluster_key=cluster_key)
    cluster_categories = list(adata.obs[cluster_key].cat.categories)
    max_cooccurrence_clusters = int(
        cfg.get("max_cooccurrence_clusters", len(cluster_categories))
    )
    for cluster_number in cluster_categories[:max_cooccurrence_clusters]:
        with rc_context({"figure.figsize": (10, 10)}):
            sq.pl.co_occurrence(
                adata,
                cluster_key=cluster_key,
                clusters=cluster_number,
            )
            save_spatial_figure(
                plt.gcf(), os.path.join(squidpy_path,
                    f"co_occurrence_cluster_{cluster_number}_{dataset_name}"
                    f"_pc{pc_label}_nc{lambda_label}_r{res_str_single}"
                ), cfg,
            )

    # Moran's I uses `.X`, so this must be the clean normalized/log1p expression
    # matrix (carried from script 06 into script 08), never a BANKSY-expanded object.
    log_time("Running Moran's I on clean log-normalized expression")
    moran_fraction = float(cfg.get("moran_subsample_fraction", 0.5))
    moran_n_perms = int(cfg.get("moran_n_perms", 100))
    moran_n_jobs = int(cfg.get("moran_n_jobs", 1))
    moran_fdr_threshold = float(cfg.get("moran_fdr_threshold", 0.10))
    top_moran_n = int(cfg.get("top_moran_n", 10))

    if moran_fraction < 1:
        adata_moran = sc.pp.subsample(
            adata,
            fraction=moran_fraction,
            copy=True,
            random_state=seed,
        )
    else:
        adata_moran = adata.copy()

    sq.gr.spatial_neighbors(adata_moran, coord_type="generic", delaunay=True)
    # Seed the permutations so simulated p-values (and the FDR gene list) are
    # reproducible; Moran's I itself is deterministic. Added 2026-11-20.
    sq.gr.spatial_autocorr(
        adata_moran,
        mode="moran",
        n_perms=moran_n_perms,
        n_jobs=moran_n_jobs,
        seed=seed,
    )

    moran_scores = adata_moran.uns["moranI"].reset_index(names=["Gene"])
    moran_scores_raw = moran_scores[
        ~moran_scores["Gene"].str.contains("_nbr_0", case=False, na=False)
    ]
    moran_scores_raw = moran_scores_raw[
        ~moran_scores_raw["Gene"].str.contains("_nbr_1", case=False, na=False)
    ]
    moran_scores_sig = moran_scores_raw[
        moran_scores_raw["pval_sim_fdr_bh"] <= moran_fdr_threshold
    ]

    moran_all_csv = os.path.join(
        table_path,
        f"moran_scores_all_{dataset_name}_pc{pc_label}_nc{lambda_label}_r{res_str_single}.csv",
    )
    moran_sig_csv = os.path.join(
        table_path,
        f"moran_scores_fdr{moran_fdr_threshold}_{dataset_name}_pc{pc_label}_nc{lambda_label}_r{res_str_single}.csv",
    )
    moran_scores_raw.to_csv(moran_all_csv, index=False)
    moran_scores_sig.to_csv(moran_sig_csv, index=False)
    print(f"Saved all Moran's I scores to: {moran_all_csv}")
    print(f"Saved significant Moran's I scores to: {moran_sig_csv}")

    top_moran = moran_scores_sig.sort_values(by="pval_sim", ascending=True).head(
        top_moran_n
    )
    for gene in top_moran["Gene"].tolist():
        with rc_context({"figure.figsize": (12, 8)}):
            sq.pl.spatial_scatter(
                adata_moran,
                library_id="spatial",
                color=[gene],
                shape=None,
                size=2,
                img=False,
                # Squidpy 1.6 forwards unknown kwargs to Matplotlib; `show`
                # therefore raises on PathCollection instead of controlling display.
                return_ax=True,
            )
            save_spatial_figure(
                plt.gcf(), os.path.join(squidpy_path,
                    f"spatial_scatter_top_moran_I_{gene}_{dataset_name}"
                    f"_pc{pc_label}_nc{lambda_label}_r{res_str_single}"
                ), cfg, spatial_points=True,
            )


def run_annotation_config(cfg, args):
    """Run label-level analyses on the script 08 annotated objects of one annotation set.

    Validation reads only backed metadata. Full runs load one sample at a time;
    annotation inputs remain unchanged and spatial results use separate paths.
    """
    annotation, level = resolve_step(cfg, args)
    require_released(annotation, args.allow_draft, "Script 10")
    column = label_column(level)
    order, colors = level_vocabulary(annotation, level)
    mapping = load_mapping(annotation)
    cohort = annotation["cohort"]
    input_root, output_root = args.input_root.resolve(), args.output_root.resolve()
    figure_base = resolve_path(output_root, step_dir(annotation, "output", STEP, level))
    result_base = resolve_path(output_root, step_dir(annotation, "processed", STEP, level))
    samples = select_samples(cohort, cfg.get("samples"), args.sample)
    print(f"Annotation set {annotation['version']} ({annotation['status']}), level {level} ({column})")
    for sample in samples:
        name = sample["sample"]
        source = resolve_path(input_root, annotated_h5ad_path(annotation, sample))
        if not source.is_file():
            raise FileNotFoundError(f"{name}: no script 08 object for set {annotation['version']} at {source}")
        adata = ad.read_h5ad(source, backed="r" if args.validate_only else None)
        try:
            check_stamp(dict(adata.uns.get("cell_type_annotation", {})), annotation, name)
            groupby = sample["groupby"]
            for key in (column, groupby):
                if key not in adata.obs or adata.obs[key].isna().any():
                    raise ValueError(f"{name}: missing labels in {key}")
            if not adata.obs_names.is_unique:
                raise ValueError(f"{name}: duplicate cell identifiers")
            labels = adata.obs[column].astype(str)
            if set(labels) - set(order):
                raise ValueError(f"{name}: labels outside the {level} vocabulary")
            # Verify labels against the set's mapping for this exact grouping.
            rows = mapping_for_sample(mapping, sample)
            expected = adata.obs[groupby].astype(str).map(rows.set_index("cluster_id")[column])
            if expected.isna().any() or not (expected.to_numpy() == labels.to_numpy()).all():
                raise ValueError(f"{name}: {column} disagrees with the annotation set mapping")
            prepare_spatial_coordinates(adata, cohort["coord_keys"])
            if not np.isfinite(adata.obsm["spatial"]).all():
                raise ValueError(f"{name}: nonfinite spatial coordinates")
            # Absent labels must not create empty groups in Squidpy's calculations.
            present = [label for label in order if label in set(labels)]
            print(f"Validated {name}: {adata.n_obs:,} cells; {len(present)} {level} labels; "
                  f"{sample['branch']} r{sample['resolution']}")
            if args.validate_only:
                continue
            adata.obs[column] = pd.Categorical(labels, categories=present)
            adata.uns[f"{column}_colors"] = [colors[label] for label in present]
            tag = f"{sample['branch']}_r{sample['resolution']}"
            output = figure_base / name / tag
            destination = result_base / name / f"{source.stem}_squidpy.h5ad"
            if destination.resolve() == source.resolve():
                raise ValueError("Spatial output must not overwrite the annotated input")
            os.makedirs(destination.parent, exist_ok=True)
            run_cfg = dict(cfg)
            run_cfg.setdefault("max_cooccurrence_clusters", len(present))
            run_squidpy_for_resolution(
                adata, column, str(output), name, "annotated", level,
                sample["resolution"].replace(".", "p"), run_cfg,
                table_dir=destination.parent / tag,
            )
            adata.uns["spatial_analysis_provenance"] = {
                "annotation_version": annotation["version"], "annotation_sha256": annotation["sha256"],
                "label_level": level, "label_column": column,
                "cohort": annotation["cohort_name"], "input_h5ad": str(source),
                "sample": name, "branch": sample["branch"],
                "resolution": sample["resolution"], "source_groupby": groupby,
            }
            adata.write_h5ad(destination)
            print(f"Saved spatial results: {destination}")
            # Small copy of the enrichment so cross-sample helpers need not read this H5AD.
            table = enrichment_table_path(result_base, name, source)
            write_enrichment_table(table, *enrichment_from_labels(
                adata.uns[f"{column}_nhood_enrichment"], adata.obs[column]))
            print(f"Saved neighbourhood enrichment table: {table}")
        finally:
            if adata.isbacked:
                adata.file.close()
        del adata


def main():
    """Run cell-type spatial analyses for the annotation set named in the config."""
    args = parse_args()
    with open(args.config) as handle:
        cfg = json.load(handle)
    figure_export_settings(cfg)
    run_annotation_config(cfg, args)


if __name__ == "__main__":
    main()
