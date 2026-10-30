#!/usr/bin/env python
# coding: utf-8

"""
Title: Squidpy Spatial Analysis For BANKSY Reclusters
Date: 2026-10-26
Summary: Reuse script 08 annotated clean expression objects and reviewed cell
types for Squidpy neighbourhood enrichment, centrality, co-occurrence, and
Moran's I. Legacy script 06 cluster-based configurations remain supported.
Export PNG, SVG, and PDF figures with editable vector text and 300-DPI raster
layers. The workflow is compatible with Squidpy 1.6 spatial-scatter plotting.
"""

import argparse
import json
import os
import random
from pathlib import Path
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


sc.logging.print_header()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 1
plt.rcParams["font.family"] = "Arial"
sns.set_style("white")

seed = 1234
np.random.seed(seed)
random.seed(seed)


def parse_args():
    """Parse command-line arguments for the Squidpy recluster analysis."""
    parser = argparse.ArgumentParser(
        prog="run Squidpy spatial analyses on script 06 BANKSY reclusters"
    )
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="JSON config with project, dataset_name, input_label, and reclustering settings.",
    )
    parser.add_argument("--sample", action="append", help="Select a script 08 sample; repeat as needed")
    parser.add_argument("--validate-only", action="store_true", help="Check annotated inputs without running analyses")
    parser.add_argument("--annotation-root", help="Override the root containing script 08 annotated outputs")
    return parser.parse_args()


def log_time(step):
    """Print a timestamped workflow message."""
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {step}")


def ensure_directory(path, label):
    """Create `path` and log the resolved location."""
    os.makedirs(path, exist_ok=True)
    print(f"{label} directory ready: {os.path.abspath(path)}")




def make_cluster_palette(n_colors):
    """Create enough distinct categorical colours for cluster plots."""
    base_cmaps = ["tab20", "tab20b", "tab20c"]
    colors = []

    for cmap_name in base_cmaps:
        cmap = plt.get_cmap(cmap_name)
        colors.extend([mpl.colors.to_hex(cmap(i)) for i in range(cmap.N)])

    if n_colors > len(colors):
        extra_cmap = plt.get_cmap("nipy_spectral")
        extra_count = n_colors - len(colors)
        colors.extend(
            [
                mpl.colors.to_hex(extra_cmap(i / max(extra_count - 1, 1)))
                for i in range(extra_count)
            ]
        )

    return colors[:n_colors]


def make_banksy_label_col(nbr_weight_decay, pc_label, lambda_label, resolution):
    """Build the native BANKSY cluster-label column name for one resolution."""
    return (
        f"labels_{nbr_weight_decay}"
        f"_pc{pc_label}"
        f"_nc{lambda_label}"
        f"_r{resolution:.2f}"
    )


def make_recluster_label_col(
    nbr_weight_decay,
    pc_label,
    lambda_label,
    resolution,
    input_label,
):
    """Build the provenance-safe reclustered label column name from script 06."""
    banksy_col = make_banksy_label_col(
        nbr_weight_decay,
        pc_label,
        lambda_label,
        resolution,
    )
    return f"{banksy_col}_recluster_{input_label}"


def resolve_analysis_resolutions(cfg, res_label_list):
    """Return the reclustering resolutions that should receive Squidpy analyses."""
    analysis_res_label = cfg.get("analysis_res_label", res_label_list[-1])
    if analysis_res_label == "all":
        return [float(res) for res in res_label_list]
    if isinstance(analysis_res_label, list):
        return [float(res) for res in analysis_res_label]
    return [float(analysis_res_label)]


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
):
    """Run spatial analyses and export each plot in the configured formats."""
    figure_export_settings(cfg)
    if cluster_key not in adata.obs.columns:
        raise KeyError(
            f"{cluster_key!r} was not found in adata.obs. "
            f"Available obs columns: {list(adata.obs.columns)}"
        )

    if not isinstance(adata.obs[cluster_key].dtype, pd.CategoricalDtype):
        adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")

    # Squidpy plotting looks for this colour vector in `.uns`. Supplying it
    # explicitly avoids noisy palette warnings and keeps plots consistent.
    color_key = f"{cluster_key}_colors"
    if color_key not in adata.uns:
        adata.uns[color_key] = make_cluster_palette(
            len(adata.obs[cluster_key].cat.categories)
        )

    squidpy_path = os.path.join(output_dir, "squidpy")
    ensure_directory(squidpy_path, "Squidpy")
    sc.settings.figdir = squidpy_path

    # All Squidpy analyses use the clean expression object used for BANKSY from
    # script 06. In qc_pass_only mode, failed cells are absent here, so they do not contribute to any
    # spatial graph, neighbourhood summary, or expression autocorrelation.
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
            squidpy_path,
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
    # matrix from script 06 rather than any BANKSY-expanded spatial object.
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
    sq.gr.spatial_autocorr(
        adata_moran,
        mode="moran",
        n_perms=moran_n_perms,
        n_jobs=moran_n_jobs,
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
        squidpy_path,
        f"moran_scores_all_{dataset_name}_pc{pc_label}_nc{lambda_label}_r{res_str_single}.csv",
    )
    moran_sig_csv = os.path.join(
        squidpy_path,
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


def resolve_path(root, value):
    """Resolve a configured absolute path or a path relative to its project root."""
    path = Path(value)
    return path if path.is_absolute() else root / path


def run_annotation_config(cfg, args):
    """Reuse script 08 sample definitions and annotated outputs for cell-type analyses.

    Validation reads only backed metadata. Full runs load one sample at a time;
    annotation inputs remain unchanged and spatial results use separate paths.
    """
    config_path = Path(args.config).resolve()
    root = config_path.parents[int(cfg.get("project_root_parent_depth", 4))]
    annotation_path = resolve_path(root, cfg["annotation_config"])
    annotation = json.loads(annotation_path.read_text())
    annotation_project = annotation_path.parents[int(annotation.get("project_root_parent_depth", 3))]
    annotation_root = Path(args.annotation_root or annotation.get("output_root", annotation_project)).resolve()
    palette = json.loads(resolve_path(annotation_project, annotation["palette_json"]).read_text())
    order, colors = palette["cell_type_order"], palette["colors"]
    if len(order) != len(set(order)) or set(order) != set(colors):
        raise ValueError("Shared palette has inconsistent labels")
    if not all(mpl.colors.is_color_like(color) for color in colors.values()):
        raise ValueError("Invalid shared palette colour")
    samples = annotation["samples"]
    names = [sample["sample"] for sample in samples]
    if len(names) != len(set(names)):
        raise ValueError("Duplicate samples in annotation config")
    if args.sample:
        if set(args.sample) - set(names):
            raise ValueError("Requested samples are absent from script 08 config")
        samples = [sample for sample in samples if sample["sample"] in args.sample]
    mapping = pd.read_csv(resolve_path(annotation_project, annotation["mapping_csv"]), dtype={"cluster_id": str, "resolution": str})
    for sample in samples:
        name = sample["sample"]
        source = resolve_path(annotation_root, sample["output_h5ad"])
        adata = ad.read_h5ad(source, backed="r" if args.validate_only else None)
        try:
            groupby = sample["groupby"]
            for key in ("cell_type", groupby):
                if key not in adata.obs or adata.obs[key].isna().any():
                    raise ValueError(f"{name}: missing labels in {key}")
            if not adata.obs_names.is_unique:
                raise ValueError(f"{name}: duplicate cell identifiers")
            labels = adata.obs["cell_type"].astype(str)
            if set(labels) - set(order):
                raise ValueError(f"{name}: unknown cell types")
            # Verify annotations against the reviewed mapping for this exact grouping.
            rows = mapping.loc[(mapping["sample"] == name) & (mapping["branch"] == sample["branch"])
                               & (mapping["resolution"] == sample["resolution"]) & (mapping["groupby"] == groupby)]
            if rows.empty or rows["cluster_id"].duplicated().any():
                raise ValueError(f"{name}: missing or duplicate reviewed mapping")
            expected = adata.obs[groupby].astype(str).map(rows.set_index("cluster_id")["cell_type"])
            if expected.isna().any() or not (expected.to_numpy() == labels.to_numpy()).all():
                raise ValueError(f"{name}: annotations disagree with reviewed cluster mapping")
            prepare_spatial_coordinates(adata, cfg.get("coord_keys", ["x", "y"]))
            if not np.isfinite(adata.obsm["spatial"]).all():
                raise ValueError(f"{name}: nonfinite spatial coordinates")
            # Absent types must not create empty groups in Squidpy's calculations.
            present = [label for label in order if label in set(labels)]
            adata.obs["cell_type"] = pd.Categorical(labels, categories=present)
            adata.uns["cell_type_colors"] = [colors[label] for label in present]
            print(f"Validated {name}: {adata.n_obs:,} cells; {len(present)} types; {sample['branch']} r{sample['resolution']}")
            if args.validate_only:
                continue
            tag = f"{sample['branch']}_r{sample['resolution']}"
            output = resolve_path(root, cfg["figure_dir"]) / name / tag
            destination = resolve_path(root, cfg["result_dir"]) / name / f"{source.stem}_squidpy.h5ad"
            if destination.resolve() == source.resolve():
                raise ValueError("Spatial output must not overwrite the annotated input")
            os.makedirs(destination.parent, exist_ok=True)
            run_cfg = dict(cfg)
            run_cfg.setdefault("max_cooccurrence_clusters", len(present))
            run_squidpy_for_resolution(
                adata, "cell_type", str(output), name, "annotated", "cell_type",
                sample["resolution"].replace(".", "p"), run_cfg,
            )
            adata.uns["spatial_analysis_provenance"] = {
                "annotation_config": str(annotation_path), "input_h5ad": str(source),
                "sample": name, "branch": sample["branch"],
                "resolution": sample["resolution"], "source_groupby": groupby,
            }
            adata.write_h5ad(destination)
            print(f"Saved spatial results: {destination}")
        finally:
            if adata.isbacked:
                adata.file.close()
        del adata


def run_legacy_config(cfg, args):
    """Preserve the existing single-sample BANKSY-cluster workflow."""
    project = cfg["project"]
    dataset_name = cfg["dataset_name"]
    input_label = cfg.get("input_label", cfg.get("output_label", "filtered_qc_v1"))
    recluster_inclusion = cfg.get("recluster_inclusion", "qc_pass_only")
    if recluster_inclusion not in ["qc_pass_only", "all_cells"]:
        raise ValueError(
            "recluster_inclusion must be either 'qc_pass_only' or 'all_cells'. "
            f"Received: {recluster_inclusion!r}"
        )
    run_label = cfg.get("run_label", f"{input_label}_{recluster_inclusion}")
    pc_label = cfg["pc_label"]
    lambda_label = cfg["lambda_label"]
    res_label = cfg["res_label"]
    res_label_list = res_label if isinstance(res_label, list) else [res_label]
    res_str = "_".join(res_label_list)
    nbr_weight_decay = cfg["nbr_weight_decay"]
    coord_keys = tuple(cfg.get("coord_keys", ["x", "y", "xy"]))

    base_dir = cfg.get("base_dir", "data/xenium")
    processed_path = os.path.join(base_dir, "processed", project, dataset_name)
    output_path = cfg.get(
        "output_path",
        os.path.join(
            base_dir,
            "output",
            project,
            dataset_name,
            f"reclustering_{run_label}",
        ),
    )
    if not args.validate_only:
        ensure_directory(output_path, "Output reclustering sample")

    input_h5ad = cfg.get(
        "input_h5ad",
        os.path.join(
            processed_path,
            f"adata_expression_clean_{dataset_name}_recluster_{run_label}"
            f"_cells_used_for_banksy_with_clusters_{res_str}.h5ad",
        ),
    )
    if not os.path.isfile(input_h5ad):
        raise FileNotFoundError(f"Script 06 reclustered-cell AnnData does not exist: {input_h5ad}")

    output_h5ad = cfg.get(
        "output_h5ad",
        os.path.join(
            processed_path,
            f"adata_expression_clean_{dataset_name}_recluster_{run_label}"
            f"_cells_used_for_banksy_with_clusters_{res_str}_squidpy.h5ad",
        ),
    )

    # -----------------------------------------------------------------------------
    # Load clean passing-cell object and run selected Squidpy resolutions
    # -----------------------------------------------------------------------------

    log_time(f"Reading script 06 reclustered-cell clean AnnData from: {input_h5ad}")
    adata = ad.read_h5ad(input_h5ad, backed="r" if args.validate_only else None)
    if args.validate_only:
        try:
            prepare_spatial_coordinates(adata, coord_keys)
            for res in resolve_analysis_resolutions(cfg, res_label_list):
                key = make_recluster_label_col(nbr_weight_decay, pc_label, lambda_label, res, run_label)
                if key not in adata.obs:
                    raise KeyError(key)
            print(f"Validated {dataset_name}: {adata.n_obs:,} cells")
        finally:
            adata.file.close()
        return
    prepare_spatial_coordinates(adata, coord_keys)

    analysis_resolutions = resolve_analysis_resolutions(cfg, res_label_list)
    print(f"Running Squidpy analyses for resolutions: {analysis_resolutions}")

    for res in analysis_resolutions:
        res_str_single = str(res).replace(".", "p")
        cluster_key = make_recluster_label_col(
            nbr_weight_decay,
            pc_label,
            lambda_label,
            res,
            run_label,
        )
        run_squidpy_for_resolution(
            adata=adata,
            cluster_key=cluster_key,
            output_dir=output_path,
            dataset_name=dataset_name,
            pc_label=pc_label,
            lambda_label=lambda_label,
            res_str_single=res_str_single,
            cfg=cfg,
        )

    adata.write_h5ad(output_h5ad)
    print(f"Wrote Squidpy-annotated reclustered-cell AnnData to: {output_h5ad}")
    log_time(f"Finished Squidpy recluster spatial analysis for {dataset_name}.")


def main():
    """Dispatch shared script 08 inputs or legacy cluster-based configurations."""
    args = parse_args()
    with open(args.config) as handle:
        cfg = json.load(handle)
    figure_export_settings(cfg)
    if "annotation_config" in cfg:
        run_annotation_config(cfg, args)
    else:
        if args.sample or args.annotation_root:
            raise ValueError("--sample and --annotation-root require an annotation_config")
        run_legacy_config(cfg, args)


if __name__ == "__main__":
    main()
