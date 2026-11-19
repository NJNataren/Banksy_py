#!/usr/bin/env python
# coding: utf-8

"""
Title: PTMT Marker Expression Spatial QC
Date: 2026-11-10
Summary: Audit selected PTMT marker genes in clean normalized/log1p AnnData
objects before script 03 dotplot export. The script writes sample- and
cluster-level expression summaries, a marker-review scaffold, and spatial and
cluster QC figures without modifying the source AnnData objects or making
automatic marker-inclusion decisions.

The JSON config must define an explicit gene list and sample list. Each sample
entry supplies a clean clustered AnnData path and BANKSY cluster column. Paths
are resolved relative to the repository root unless they are absolute.

Minimal config shape:

{
  "analysis_name": "ptmt_pc55_run1_r1p5_hypoxia_cell_death_qc",
  "project": "ptmt_pc55",
  "genes": [
    {"gene": "HIF1A", "category": "Hypoxia-associated"},
    {"gene": "IL1A", "category": "Damage/inflammation-associated"}
  ],
  "expression_source": "raw",
  "allow_x_fallback": false,
  "output_root": "data/xenium/output/ptmt/QC/ptmt_pc55_run1_r1p5_hypoxia_cell_death_qc",
  "samples": [
    {
      "sample": "10329_run_1_DMSO",
      "adata_path": "data/xenium/processed/ptmt_pc55/10329_run_1_DMSO/<clean>.h5ad",
      "resolution": "1.50",
      "cluster_column": "labels_scaled_gaussian_pc55_nc0.20_r1.50"
    }
  ]
}

In script 00 clean objects, `.raw.X` is a normalized/log1p snapshot rather
than raw counts. Original counts, when retained, are in `layers["counts"]`.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import re
from typing import Any


DEFAULT_GENES = ["HIF1A", "IL1A", "IL1B", "PMAIP1"]
VALID_MODES = ("all", "cluster", "spatial")
VALID_FORMATS = {"png", "pdf", "svg"}


@dataclass
class SampleInspection:
    """Hold validated numerical summaries for one sample."""

    sample_cfg: dict[str, Any]
    audit_rows: list[dict[str, Any]]
    sample_rows: list[dict[str, Any]]
    cluster_rows: list[dict[str, Any]]
    scale_rows: list[dict[str, Any]]
    warnings: list[str]


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for marker-expression QC."""
    parser = argparse.ArgumentParser(
        description=(
            "Audit normalized PTMT marker expression and create pre-dotplot "
            "spatial and cluster QC outputs."
        )
    )
    parser.add_argument(
        "--config",
        required=True,
        help="JSON config defining genes, clean AnnData objects, and outputs.",
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Repository root used to resolve relative paths.",
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        help=(
            "Optional root used only to resolve relative AnnData paths. "
            "Useful for reading inputs from a mounted project checkout while "
            "writing outputs beneath the local --repo-root."
        ),
    )
    parser.add_argument(
        "--sample",
        action="append",
        dest="samples",
        help="Run one named sample; repeat to select multiple samples.",
    )
    parser.add_argument(
        "--mode",
        choices=VALID_MODES,
        default="all",
        help="Create all outputs, cluster QC only, or spatial QC only.",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Audit inputs and print status without writing files or figures.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing QC outputs.",
    )
    return parser.parse_args()


def load_config(path: Path) -> dict[str, Any]:
    """Read one JSON configuration file."""
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def resolve_path(value: str | Path, repo_root: Path) -> Path:
    """Resolve a configured path relative to the repository root."""
    path = Path(value).expanduser()
    return path if path.is_absolute() else repo_root / path


def safe_label(value: Any) -> str:
    """Return a filename-safe representation of a label."""
    value = str(value).replace(".", "p")
    return re.sub(r"[^A-Za-z0-9_-]+", "_", value).strip("_")


def parse_gene_settings(cfg: dict[str, Any]) -> tuple[list[str], dict[str, str]]:
    """Return unique configured genes and their display categories."""
    entries = cfg.get("genes", DEFAULT_GENES)
    if not isinstance(entries, list) or not entries:
        raise ValueError("Config field 'genes' must be a non-empty list")

    genes: list[str] = []
    categories: dict[str, str] = {}
    for entry in entries:
        if isinstance(entry, str):
            gene = entry.strip()
            category = "Marker"
        elif isinstance(entry, dict):
            gene = str(entry.get("gene", "")).strip()
            category = str(entry.get("category", "Marker")).strip() or "Marker"
        else:
            raise TypeError("Each genes entry must be a string or an object")

        if not gene:
            raise ValueError("Configured genes cannot contain blank gene names")
        if gene in categories:
            raise ValueError(f"Gene is configured more than once: {gene}")
        genes.append(gene)
        categories[gene] = category

    return genes, categories


def select_sample_configs(
    cfg: dict[str, Any], selected_samples: list[str] | None
) -> list[dict[str, Any]]:
    """Validate and optionally subset the explicit sample configuration list."""
    sample_cfgs = cfg.get("samples")
    if not isinstance(sample_cfgs, list) or not sample_cfgs:
        raise ValueError("Config field 'samples' must be a non-empty list")

    names = [str(item.get("sample", "")) for item in sample_cfgs]
    if any(not name for name in names):
        raise ValueError("Every sample entry must define a non-empty 'sample'")
    if len(names) != len(set(names)):
        duplicates = sorted({name for name in names if names.count(name) > 1})
        raise ValueError(f"Duplicate sample entries: {duplicates}")

    if not selected_samples:
        return sample_cfgs

    requested = set(selected_samples)
    unknown = sorted(requested - set(names))
    if unknown:
        raise ValueError(f"Requested samples are not in the config: {unknown}")
    return [item for item in sample_cfgs if item["sample"] in requested]


def validate_config(cfg: dict[str, Any]) -> None:
    """Validate global settings that could otherwise produce misleading QC."""
    expression_source = cfg.get("expression_source", "raw")
    if expression_source not in {"raw", "X"}:
        raise ValueError("expression_source must be either 'raw' or 'X'")

    plot_cfg = cfg.get("plot", {})
    quantile = float(plot_cfg.get("within_sample_quantile", 0.99))
    if not 0 < quantile <= 1:
        raise ValueError("plot.within_sample_quantile must be in (0, 1]")

    formats = plot_cfg.get("formats", ["png"])
    invalid_formats = sorted(set(formats) - VALID_FORMATS)
    if invalid_formats:
        raise ValueError(f"Unsupported figure formats: {invalid_formats}")


def to_dense(matrix: Any) -> Any:
    """Return a matrix as a NumPy array without changing its values."""
    import numpy as np
    from scipy import sparse

    if sparse.issparse(matrix):
        return matrix.toarray()
    return np.asarray(matrix)


def expression_matrix(
    adata: Any,
    genes: list[str],
    requested_source: str,
    allow_x_fallback: bool,
) -> tuple[Any, list[str], str]:
    """Extract selected genes from the configured normalized expression source.

    Args:
        adata: Clean clustered AnnData object.
        genes: Requested genes in display order.
        requested_source: Either `raw` or `X`.
        allow_x_fallback: Whether `.X` may replace a missing `.raw` snapshot.

    Returns:
        Dense cell-by-gene matrix, present genes, and source actually used.
    """
    if requested_source == "raw":
        if adata.raw is not None:
            var_names = set(map(str, adata.raw.var_names))
            present = [gene for gene in genes if gene in var_names]
            matrix = to_dense(adata.raw[:, present].X)
            return matrix, present, "raw_normalized_log1p_snapshot"
        if not allow_x_fallback:
            raise ValueError(
                "Configured expression_source is 'raw', but adata.raw is missing "
                "and allow_x_fallback is false"
            )

    var_names = set(map(str, adata.var_names))
    present = [gene for gene in genes if gene in var_names]
    matrix = to_dense(adata[:, present].X)
    source_used = "X_normalized_log1p"
    if requested_source == "raw":
        source_used += "_fallback"
    return matrix, present, source_used


def selected_x_matrix(adata: Any, genes: list[str]) -> tuple[Any, list[str]]:
    """Return selected genes from `.X` for source-parity auditing."""
    var_names = set(map(str, adata.var_names))
    present = [gene for gene in genes if gene in var_names]
    return to_dense(adata[:, present].X), present


def selected_count_matrix(adata: Any, genes: list[str]) -> tuple[Any | None, list[str]]:
    """Return selected genes from `layers['counts']` when available."""
    import numpy as np

    if "counts" not in adata.layers:
        return None, []
    var_index = {str(gene): index for index, gene in enumerate(adata.var_names)}
    present = [gene for gene in genes if gene in var_index]
    indices = [var_index[gene] for gene in present]
    matrix = adata.layers["counts"][:, indices]
    return np.asarray(to_dense(matrix)), present


def coordinate_matrix(adata: Any, sample_cfg: dict[str, Any], cfg: dict[str, Any]) -> tuple[Any, str]:
    """Return finite two-dimensional tissue coordinates and their source."""
    import numpy as np

    requested = sample_cfg.get("coordinate_key", cfg.get("coordinate_key"))
    candidates = [requested] if requested else ["spatial", "xy"]
    for key in candidates:
        if key and key in adata.obsm:
            coords = np.asarray(adata.obsm[key])
            if coords.ndim != 2 or coords.shape[0] != adata.n_obs or coords.shape[1] < 2:
                raise ValueError(f"Invalid obsm['{key}'] coordinate shape: {coords.shape}")
            coords = coords[:, :2]
            if not np.isfinite(coords).all():
                raise ValueError(f"Non-finite values found in obsm['{key}']")
            return coords, f"obsm['{key}']"

    if {"x", "y"}.issubset(adata.obs.columns):
        coords = adata.obs[["x", "y"]].to_numpy()
        if not np.isfinite(coords).all():
            raise ValueError("Non-finite values found in obs x/y coordinates")
        return coords, "obs['x','y']"

    raise KeyError("No coordinates found in obsm['spatial'], obsm['xy'], or obs x/y")


def normalized_cluster_labels(adata: Any, cluster_column: str, missing_label: str) -> Any:
    """Return string cluster labels with missing values made explicit."""
    values = adata.obs[cluster_column].astype("object")
    missing = values.isna() | values.astype(str).str.strip().isin(
        ["", "nan", "None", "NA", "<NA>"]
    )
    return values.where(~missing, missing_label).astype(str).to_numpy()


def cluster_sort_key(value: Any) -> tuple[int, float | str]:
    """Sort numeric cluster labels numerically and other labels lexically."""
    try:
        return 0, float(value)
    except (TypeError, ValueError):
        return 1, str(value)


def robust_positive_vmax(values: Any, quantile: float) -> float:
    """Return a robust positive-cell colour maximum with a nonzero fallback."""
    import numpy as np

    positive = np.asarray(values)[np.asarray(values) > 0]
    if positive.size == 0:
        return 1.0
    vmax = float(np.quantile(positive, quantile))
    return vmax if np.isfinite(vmax) and vmax > 0 else float(positive.max())


def gene_vector(matrix: Any, present_genes: list[str], gene: str) -> Any:
    """Return one gene's one-dimensional expression vector."""
    import numpy as np

    return np.asarray(matrix[:, present_genes.index(gene)]).ravel()


def source_parity_metrics(adata: Any, gene: str, source_values: Any) -> dict[str, Any]:
    """Compare selected normalized expression with `.X` and count detection."""
    import numpy as np

    metrics: dict[str, Any] = {
        "x_gene_present": gene in set(map(str, adata.var_names)),
        "raw_x_max_abs_difference": math.nan,
        "raw_x_detection_mismatch_cells": math.nan,
        "counts_gene_present": False,
        "counts_detection_mismatch_cells": math.nan,
    }
    source_values = np.asarray(source_values).ravel()

    if metrics["x_gene_present"]:
        x_values = to_dense(adata[:, [gene]].X).ravel()
        metrics["raw_x_max_abs_difference"] = float(
            np.max(np.abs(source_values - x_values))
        )
        metrics["raw_x_detection_mismatch_cells"] = int(
            np.count_nonzero((source_values > 0) != (x_values > 0))
        )

    counts, count_genes = selected_count_matrix(adata, [gene])
    if counts is not None and gene in count_genes:
        count_values = np.asarray(counts).ravel()
        metrics["counts_gene_present"] = True
        metrics["counts_detection_mismatch_cells"] = int(
            np.count_nonzero((source_values > 0) != (count_values > 0))
        )

    return metrics


def summarize_gene(values: Any) -> dict[str, Any]:
    """Calculate absolute expression and detection summaries for one gene."""
    import numpy as np

    values = np.asarray(values).ravel()
    detected = values > 0
    n_cells = int(values.size)
    n_expressing = int(detected.sum())
    return {
        "n_cells": n_cells,
        "n_expressing": n_expressing,
        "percent_expressing": 100 * n_expressing / n_cells if n_cells else math.nan,
        "mean_expression": float(values.mean()) if n_cells else math.nan,
        "median_expression": float(np.median(values)) if n_cells else math.nan,
        "max_expression": float(values.max()) if n_cells else math.nan,
    }


def inspect_sample(
    sample_cfg: dict[str, Any],
    cfg: dict[str, Any],
    genes: list[str],
    categories: dict[str, str],
    input_root: Path,
    require_coordinates: bool,
) -> SampleInspection:
    """Audit one clean object and calculate sample and cluster summaries."""
    import anndata as ad
    import numpy as np

    sample = str(sample_cfg["sample"])
    adata_value = sample_cfg.get("adata_path")
    if not adata_value:
        raise ValueError(f"{sample}: sample entry is missing adata_path")
    adata_path = resolve_path(adata_value, input_root)
    if not adata_path.is_file():
        raise FileNotFoundError(f"{sample}: clean AnnData not found: {adata_path}")

    cluster_column = sample_cfg.get("cluster_column", cfg.get("cluster_column"))
    if not cluster_column:
        raise ValueError(f"{sample}: cluster_column is not configured")

    resolution = str(sample_cfg.get("resolution", cfg.get("resolution", "")))
    expression_source = cfg.get("expression_source", "raw")
    allow_x_fallback = bool(cfg.get("allow_x_fallback", False))
    missing_cluster_label = str(cfg.get("missing_cluster_label", "Missing"))
    cluster_label_map = {
        str(key): str(value)
        for key, value in sample_cfg.get("cluster_labels", {}).items()
    }
    min_cluster_cells = int(cfg.get("min_cluster_cells", 20))
    quantile = float(cfg.get("plot", {}).get("within_sample_quantile", 0.99))

    print(f"Reading {sample}: {adata_path}")
    adata = ad.read_h5ad(adata_path)
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError(f"{sample}: clean AnnData has shape {adata.shape}")
    if cluster_column not in adata.obs:
        raise KeyError(f"{sample}: cluster column not found: {cluster_column}")
    if not adata.obs_names.is_unique:
        raise ValueError(f"{sample}: observation names are not unique")
    if not adata.var_names.is_unique:
        raise ValueError(f"{sample}: gene names are not unique")

    coordinates_source = "not_checked"
    if require_coordinates:
        _, coordinates_source = coordinate_matrix(adata, sample_cfg, cfg)

    matrix, present_genes, source_used = expression_matrix(
        adata,
        genes,
        expression_source,
        allow_x_fallback,
    )
    matrix = np.asarray(matrix)
    if matrix.ndim == 1:
        matrix = matrix[:, None]
    if matrix.shape != (adata.n_obs, len(present_genes)):
        raise ValueError(
            f"{sample}: expression shape {matrix.shape} does not match "
            f"{adata.n_obs} cells and {len(present_genes)} genes"
        )
    if not present_genes:
        raise ValueError(f"{sample}: none of the requested genes are present")
    if not np.isfinite(matrix).all():
        raise ValueError(f"{sample}: selected expression contains non-finite values")
    if np.min(matrix) < 0:
        raise ValueError(f"{sample}: selected normalized expression contains negatives")

    cluster_values = normalized_cluster_labels(
        adata, cluster_column, missing_cluster_label
    )
    cluster_ids = sorted(set(cluster_values), key=cluster_sort_key)
    log1p_metadata_present = "log1p" in adata.uns
    clean_name_matches = adata_path.name.startswith("adata_expression_clean_")

    warnings: list[str] = []
    missing_genes = [gene for gene in genes if gene not in present_genes]
    if missing_genes:
        warnings.append(f"{sample}: missing requested genes: {missing_genes}")
    if not log1p_metadata_present:
        warnings.append(f"{sample}: adata.uns does not contain log1p metadata")
    if not clean_name_matches:
        warnings.append(
            f"{sample}: input filename does not start with adata_expression_clean_"
        )
    if source_used.endswith("fallback"):
        warnings.append(f"{sample}: expression source fell back from .raw to .X")

    common = {
        "sample": sample,
        "resolution": resolution,
        "cluster_column": cluster_column,
        "expression_source_requested": expression_source,
        "expression_source_used": source_used,
        "expression_scale_expected": "normalized_log1p",
        "adata_path": str(adata_path),
    }
    audit_rows: list[dict[str, Any]] = []
    sample_rows: list[dict[str, Any]] = []
    cluster_rows: list[dict[str, Any]] = []
    scale_rows: list[dict[str, Any]] = []

    for gene in genes:
        present = gene in present_genes
        audit = {
            **common,
            "gene": gene,
            "category": categories[gene],
            "status": "ok" if present else "warning",
            "status_detail": "present" if present else "gene_missing",
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "n_clusters": len(cluster_ids),
            "coordinates_source": coordinates_source,
            "raw_present": adata.raw is not None,
            "counts_layer_present": "counts" in adata.layers,
            "log1p_metadata_present": log1p_metadata_present,
            "clean_filename_matches": clean_name_matches,
            "gene_present": present,
        }
        if not present:
            audit_rows.append(audit)
            continue

        values = gene_vector(matrix, present_genes, gene)
        parity = source_parity_metrics(adata, gene, values)
        audit.update(parity)
        if parity["counts_gene_present"] and parity["counts_detection_mismatch_cells"]:
            audit["status"] = "warning"
            audit["status_detail"] = "counts_detection_mismatch"
            warnings.append(
                f"{sample} {gene}: normalized/count detection mismatch in "
                f"{parity['counts_detection_mismatch_cells']} cells"
            )
        audit_rows.append(audit)

        sample_stats = summarize_gene(values)
        sample_rows.append(
            {
                **common,
                "gene": gene,
                "category": categories[gene],
                **sample_stats,
            }
        )
        within_vmax = robust_positive_vmax(values, quantile)
        scale_rows.append(
            {
                "sample": sample,
                "gene": gene,
                "quantile": quantile,
                "within_sample_vmin": 0.0,
                "within_sample_vmax": within_vmax,
                "shared_vmin": 0.0,
                "shared_vmax": math.nan,
                "shared_scale_method": "max_of_sample_positive_quantiles",
            }
        )

        sample_n_expressing = sample_stats["n_expressing"]
        for cluster_id in cluster_ids:
            mask = cluster_values == cluster_id
            cluster_stats = summarize_gene(values[mask])
            cluster_n_expressing = cluster_stats["n_expressing"]
            cluster_rows.append(
                {
                    **common,
                    "cluster_id": cluster_id,
                    "cluster_label": cluster_label_map.get(cluster_id, cluster_id),
                    "gene": gene,
                    "category": categories[gene],
                    **cluster_stats,
                    "fraction_of_sample_expressing_cells": (
                        cluster_n_expressing / sample_n_expressing
                        if sample_n_expressing
                        else 0.0
                    ),
                    "small_cluster": cluster_stats["n_cells"] < min_cluster_cells,
                }
            )

    del adata
    return SampleInspection(
        sample_cfg=sample_cfg,
        audit_rows=audit_rows,
        sample_rows=sample_rows,
        cluster_rows=cluster_rows,
        scale_rows=scale_rows,
        warnings=warnings,
    )


def derive_shared_scales(scale_frame: Any) -> Any:
    """Assign one robust cross-sample colour maximum to each gene."""
    frame = scale_frame.copy()
    shared = frame.groupby("gene", sort=False)["within_sample_vmax"].max()
    frame["shared_vmax"] = frame["gene"].map(shared)
    return frame


def build_review_scaffold(
    sample_frame: Any,
    genes: list[str],
    categories: dict[str, str],
    n_samples_requested: int,
) -> Any:
    """Build a review table without making automatic inclusion decisions."""
    import numpy as np
    import pandas as pd

    rows = []
    for gene in genes:
        subset = sample_frame[sample_frame["gene"] == gene]
        rows.append(
            {
                "gene": gene,
                "category": categories[gene],
                "n_samples_requested": n_samples_requested,
                "n_samples_present": int(subset["sample"].nunique()),
                "n_samples_detected": int((subset["n_expressing"] > 0).sum()),
                "median_percent_expressing": (
                    float(subset["percent_expressing"].median())
                    if not subset.empty
                    else np.nan
                ),
                "maximum_percent_expressing": (
                    float(subset["percent_expressing"].max())
                    if not subset.empty
                    else np.nan
                ),
                "median_mean_expression": (
                    float(subset["mean_expression"].median())
                    if not subset.empty
                    else np.nan
                ),
                "include_for_dotplot": "",
                "review_status": "Pending",
                "spatial_pattern": "",
                "qc_notes": "",
                "reviewer_notes": "",
            }
        )
    return pd.DataFrame(rows)


def output_settings(
    cfg: dict[str, Any], repo_root: Path
) -> tuple[str, Path, Path, dict[str, Any]]:
    """Resolve output labels, directories, and plotting settings."""
    project = str(cfg.get("project", "ptmt"))
    analysis_name = str(cfg.get("analysis_name", "ptmt_marker_expression_spatial_qc"))
    prefix = safe_label(cfg.get("output_prefix", analysis_name))
    output_root = resolve_path(
        cfg.get(
            "output_root",
            f"data/xenium/output/{project}/QC/{analysis_name}",
        ),
        repo_root,
    )
    data_dir = (
        resolve_path(cfg["output_data_dir"], repo_root)
        if "output_data_dir" in cfg
        else output_root / "tables"
    )
    figure_dir = (
        resolve_path(cfg["output_figure_dir"], repo_root)
        if "output_figure_dir" in cfg
        else output_root / "figures"
    )
    return prefix, data_dir, figure_dir, cfg.get("plot", {})


def check_write_target(path: Path, overwrite: bool) -> None:
    """Reject an existing output unless replacement was explicitly allowed."""
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists: {path}. Use --overwrite to replace it."
        )


def write_csv(frame: Any, path: Path, overwrite: bool) -> None:
    """Write one CSV after enforcing the output replacement policy."""
    check_write_target(path, overwrite)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)
    print(f"Wrote {path}")


def save_figure(
    fig: Any,
    output_stem: Path,
    formats: list[str],
    dpi: int,
    overwrite: bool,
) -> None:
    """Save and close a Matplotlib figure in configured formats."""
    import matplotlib.pyplot as plt

    output_stem.parent.mkdir(parents=True, exist_ok=True)
    # Append rather than replace a suffix so decimal labels in the stem cannot
    # silently truncate filenames (for example, r0.70 becoming r0.png).
    targets = [Path(f"{output_stem}.{suffix}") for suffix in formats]
    for target in targets:
        check_write_target(target, overwrite)
    for target in targets:
        kwargs = {"dpi": dpi} if target.suffix == ".png" else {}
        fig.savefig(target, bbox_inches="tight", **kwargs)
        print(f"Wrote {target}")
    plt.close(fig)


def dot_size_legend(ax: Any, max_size: float) -> None:
    """Add a percentage-expressing size legend to a dotplot axis."""
    handles = []
    labels = []
    for percent in (10, 25, 50, 75, 100):
        handles.append(
            ax.scatter([], [], s=max_size * percent / 100, c="#777777", alpha=0.75)
        )
        labels.append(f"{percent}%")
    ax.legend(
        handles,
        labels,
        title="Expressing",
        frameon=False,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
    )


def plot_cluster_map(
    coords: Any,
    clusters: Any,
    sample: str,
    resolution: str,
    output_stem: Path,
    plot_cfg: dict[str, Any],
    overwrite: bool,
) -> None:
    """Plot BANKSY cluster assignments as the tissue reference map."""
    import matplotlib.pyplot as plt
    import numpy as np

    point_size = float(plot_cfg.get("point_size", 1.0))
    dpi = int(plot_cfg.get("dpi", 200))
    formats = list(plot_cfg.get("formats", ["png"]))
    labels = sorted(set(clusters), key=cluster_sort_key)
    cmap = plt.get_cmap(plot_cfg.get("cluster_cmap", "tab20"), max(1, len(labels)))

    fig, ax = plt.subplots(figsize=tuple(plot_cfg.get("spatial_figsize", [10, 8])))
    for index, label in enumerate(labels):
        mask = np.asarray(clusters == label)
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=point_size,
            color=cmap(index),
            label=label,
            linewidths=0,
            rasterized=True,
        )
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title(f"{sample}: BANKSY clusters, r{resolution}")
    ax.legend(
        title="Cluster",
        frameon=False,
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        markerscale=max(2, 5 / max(point_size, 0.2)),
        ncol=1 if len(labels) <= 25 else 2,
    )
    save_figure(fig, output_stem, formats, dpi, overwrite)


def plot_expression_panel(
    coords: Any,
    matrix: Any,
    present_genes: list[str],
    vmax_by_gene: dict[str, float],
    sample: str,
    scale_label: str,
    output_stem: Path,
    plot_cfg: dict[str, Any],
    overwrite: bool,
) -> None:
    """Plot normalized expression maps with highly expressing cells on top."""
    import matplotlib.pyplot as plt
    import numpy as np

    ncols = min(int(plot_cfg.get("gene_panel_ncols", 2)), len(present_genes))
    nrows = math.ceil(len(present_genes) / ncols)
    panel_width = float(plot_cfg.get("gene_panel_width", 5.0))
    panel_height = float(plot_cfg.get("gene_panel_height", 4.5))
    point_size = float(plot_cfg.get("point_size", 1.0))
    cmap = plot_cfg.get("expression_cmap", "magma")
    dpi = int(plot_cfg.get("dpi", 200))
    formats = list(plot_cfg.get("formats", ["png"]))

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(panel_width * ncols, panel_height * nrows),
        squeeze=False,
    )
    for axis, gene in zip(axes.flat, present_genes):
        values = gene_vector(matrix, present_genes, gene)
        order = np.argsort(values, kind="stable")
        scatter = axis.scatter(
            coords[order, 0],
            coords[order, 1],
            c=values[order],
            s=point_size,
            cmap=cmap,
            vmin=0,
            vmax=vmax_by_gene[gene],
            linewidths=0,
            rasterized=True,
        )
        axis.set_aspect("equal")
        axis.set_axis_off()
        axis.set_title(gene)
        colorbar = fig.colorbar(scatter, ax=axis, fraction=0.046, pad=0.02)
        colorbar.set_label("Normalized log1p expression")
    for axis in axes.flat[len(present_genes) :]:
        axis.set_visible(False)
    fig.suptitle(f"{sample}: marker expression ({scale_label})")
    fig.tight_layout()
    save_figure(fig, output_stem, formats, dpi, overwrite)


def plot_detection_panel(
    coords: Any,
    matrix: Any,
    present_genes: list[str],
    sample: str,
    output_stem: Path,
    plot_cfg: dict[str, Any],
    overwrite: bool,
) -> None:
    """Plot binary gene detection to expose sparse or isolated positive cells."""
    import matplotlib.pyplot as plt
    import numpy as np

    ncols = min(int(plot_cfg.get("gene_panel_ncols", 2)), len(present_genes))
    nrows = math.ceil(len(present_genes) / ncols)
    panel_width = float(plot_cfg.get("gene_panel_width", 5.0))
    panel_height = float(plot_cfg.get("gene_panel_height", 4.5))
    point_size = float(plot_cfg.get("point_size", 1.0))
    dpi = int(plot_cfg.get("dpi", 200))
    formats = list(plot_cfg.get("formats", ["png"]))
    negative_color = plot_cfg.get("not_detected_color", "#d0d0d0")
    positive_color = plot_cfg.get("detected_color", "#0072B2")

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(panel_width * ncols, panel_height * nrows),
        squeeze=False,
    )
    for axis, gene in zip(axes.flat, present_genes):
        values = gene_vector(matrix, present_genes, gene)
        detected = values > 0
        for mask, color, label, zorder in [
            (~detected, negative_color, "Not detected", 1),
            (detected, positive_color, "Detected", 2),
        ]:
            if mask.any():
                axis.scatter(
                    coords[mask, 0],
                    coords[mask, 1],
                    s=point_size,
                    color=color,
                    label=label,
                    linewidths=0,
                    rasterized=True,
                    zorder=zorder,
                )
        axis.set_aspect("equal")
        axis.set_axis_off()
        axis.set_title(f"{gene}: {100 * detected.mean():.2f}% detected")
    for axis in axes.flat[len(present_genes) :]:
        axis.set_visible(False)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, loc="upper right")
    fig.suptitle(f"{sample}: marker detection")
    fig.tight_layout()
    save_figure(fig, output_stem, formats, dpi, overwrite)


def plot_cluster_dotplot(
    cluster_frame: Any,
    genes: list[str],
    sample: str,
    resolution: str,
    output_stem: Path,
    plot_cfg: dict[str, Any],
    overwrite: bool,
    zscore: bool = False,
) -> None:
    """Plot per-cluster mean expression and detection for one sample."""
    import matplotlib.pyplot as plt
    import numpy as np

    frame = cluster_frame.copy()
    if zscore:
        grouped = frame.groupby("gene")["mean_expression"]
        means = grouped.transform("mean")
        standard_deviations = grouped.transform("std").replace(0, np.nan)
        frame["plot_expression"] = (
            (frame["mean_expression"] - means) / standard_deviations
        ).fillna(0)
        cmap = plot_cfg.get("zscore_cmap", "RdBu_r")
        color_limit = float(plot_cfg.get("zscore_color_limit", 2.5))
        vmin, vmax = -color_limit, color_limit
        color_label = "Mean expression z-score within gene"
        suffix_label = "gene-wise z-score"
    else:
        frame["plot_expression"] = frame["mean_expression"]
        cmap = plot_cfg.get("dotplot_cmap", "viridis")
        vmax = robust_positive_vmax(
            frame["plot_expression"].to_numpy(),
            float(plot_cfg.get("cluster_color_quantile", 0.99)),
        )
        vmin = 0
        color_label = "Mean normalized log1p expression"
        suffix_label = "absolute expression"

    cluster_order = sorted(frame["cluster_id"].unique(), key=cluster_sort_key)
    label_lookup = (
        frame[["cluster_id", "cluster_label"]]
        .drop_duplicates("cluster_id")
        .set_index("cluster_id")["cluster_label"]
        .to_dict()
    )
    row_labels = [
        str(cluster)
        if str(label_lookup[cluster]) == str(cluster)
        else f"{cluster}: {label_lookup[cluster]}"
        for cluster in cluster_order
    ]
    x_lookup = {gene: index for index, gene in enumerate(genes)}
    y_lookup = {cluster: index for index, cluster in enumerate(cluster_order)}
    max_size = float(plot_cfg.get("max_dot_size", 180))
    width = max(6.0, 1.1 * len(genes) + 3.5)
    height = max(4.5, 0.34 * len(cluster_order) + 2.0)

    fig, ax = plt.subplots(figsize=(width, height))
    scatter = ax.scatter(
        frame["gene"].map(x_lookup),
        frame["cluster_id"].map(y_lookup),
        s=max_size * frame["percent_expressing"] / 100,
        c=frame["plot_expression"],
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        edgecolors="#333333",
        linewidths=0.25,
    )
    ax.set_xticks(range(len(genes)), genes, rotation=45, ha="right")
    ax.set_yticks(range(len(cluster_order)), row_labels)
    ax.set_xlim(-0.6, len(genes) - 0.4)
    ax.set_ylim(len(cluster_order) - 0.4, -0.6)
    ax.grid(color="#e5e5e5", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.set_xlabel("Gene")
    ax.set_ylabel("BANKSY cluster")
    ax.set_title(f"{sample}, r{resolution}: marker QC ({suffix_label})")
    colorbar = fig.colorbar(scatter, ax=ax, pad=0.16)
    colorbar.set_label(color_label)
    dot_size_legend(ax, max_size)
    fig.tight_layout()
    save_figure(
        fig,
        output_stem,
        list(plot_cfg.get("formats", ["png"])),
        int(plot_cfg.get("dpi", 200)),
        overwrite,
    )


def plot_sample_gene_heatmap(
    sample_frame: Any,
    genes: list[str],
    sample_order: list[str],
    value_column: str,
    colorbar_label: str,
    output_stem: Path,
    plot_cfg: dict[str, Any],
    overwrite: bool,
    vmin: float | None = None,
    vmax: float | None = None,
) -> None:
    """Plot a compact sample-by-gene QC heatmap."""
    import matplotlib.pyplot as plt
    import numpy as np

    pivot = sample_frame.pivot(index="sample", columns="gene", values=value_column)
    pivot = pivot.reindex(index=sample_order, columns=genes)
    fig_width = max(6.0, 1.0 * len(genes) + 3.0)
    fig_height = max(4.0, 0.32 * len(sample_order) + 2.0)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    values = pivot.to_numpy()
    image = ax.imshow(
        values,
        aspect="auto",
        cmap="viridis",
        vmin=vmin,
        vmax=vmax,
    )
    ax.set_xticks(range(len(genes)), genes, rotation=45, ha="right")
    ax.set_yticks(range(len(sample_order)), sample_order)
    ax.set_xlabel("Gene")
    ax.set_ylabel("Sample")
    for row in range(len(sample_order)):
        for column in range(len(genes)):
            value = pivot.iloc[row, column]
            if np.isfinite(value):
                label = f"{value:.1f}" if abs(value) >= 1 else f"{value:.2f}"
                normalized = image.norm(value)
                text_color = "white" if normalized < 0.35 else "black"
                ax.text(
                    column,
                    row,
                    label,
                    ha="center",
                    va="center",
                    fontsize=7,
                    color=text_color,
                )
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label(colorbar_label)
    fig.tight_layout()
    save_figure(
        fig,
        output_stem,
        list(plot_cfg.get("formats", ["png"])),
        int(plot_cfg.get("dpi", 200)),
        overwrite,
    )


def load_plot_arrays(
    sample_cfg: dict[str, Any],
    cfg: dict[str, Any],
    genes: list[str],
    input_root: Path,
) -> tuple[Any, Any, list[str], Any, str, str]:
    """Reload only the arrays required for one sample's QC figures."""
    import anndata as ad
    import numpy as np

    sample = str(sample_cfg["sample"])
    adata_path = resolve_path(sample_cfg["adata_path"], input_root)
    adata = ad.read_h5ad(adata_path)
    matrix, present_genes, _ = expression_matrix(
        adata,
        genes,
        cfg.get("expression_source", "raw"),
        bool(cfg.get("allow_x_fallback", False)),
    )
    matrix = np.asarray(matrix)
    if matrix.ndim == 1:
        matrix = matrix[:, None]
    coords, _ = coordinate_matrix(adata, sample_cfg, cfg)
    cluster_column = sample_cfg.get("cluster_column", cfg.get("cluster_column"))
    clusters = normalized_cluster_labels(
        adata,
        cluster_column,
        str(cfg.get("missing_cluster_label", "Missing")),
    )
    resolution = str(sample_cfg.get("resolution", cfg.get("resolution", "")))
    del adata
    return coords, matrix, present_genes, clusters, sample, resolution


def write_manifest(
    path: Path,
    config_path: Path,
    args: argparse.Namespace,
    genes: list[str],
    sample_names: list[str],
    overwrite: bool,
) -> None:
    """Write a small provenance manifest for the completed QC run."""
    check_write_target(path, overwrite)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "script": Path(__file__).name,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "config": str(config_path),
        "repo_root": str(args.repo_root.resolve()),
        "input_root": str(
            args.input_root.resolve() if args.input_root else args.repo_root.resolve()
        ),
        "mode": args.mode,
        "genes": genes,
        "samples": sample_names,
        "source_objects_modified": False,
        "automatic_marker_exclusion": False,
    }
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")
    print(f"Wrote {path}")


def print_validation_summary(audit_frame: Any, warnings: list[str]) -> None:
    """Print concise audit status for review and batch logs."""
    status_counts = audit_frame["status"].value_counts(dropna=False).to_dict()
    samples = audit_frame["sample"].nunique()
    genes_present = int(audit_frame["gene_present"].fillna(False).sum())
    print(
        f"Validated {samples} samples; {genes_present} sample-gene inputs present; "
        f"status counts: {status_counts}"
    )
    for warning in warnings:
        print(f"WARNING: {warning}")


def main() -> None:
    """Run pre-dotplot marker-expression QC from clean normalized objects."""
    import pandas as pd

    args = parse_args()
    repo_root = args.repo_root.resolve()
    input_root = args.input_root.resolve() if args.input_root else repo_root
    config_path = resolve_path(args.config, repo_root)
    cfg = load_config(config_path)
    validate_config(cfg)
    genes, categories = parse_gene_settings(cfg)
    sample_cfgs = select_sample_configs(cfg, args.samples)
    require_coordinates = args.mode in {"all", "spatial"}

    inspections: list[SampleInspection] = []
    errors: list[str] = []
    for sample_cfg in sample_cfgs:
        try:
            inspections.append(
                inspect_sample(
                    sample_cfg,
                    cfg,
                    genes,
                    categories,
                    input_root,
                    require_coordinates,
                )
            )
        except Exception as exc:  # Continue so one validation run audits all samples.
            sample = str(sample_cfg.get("sample", "unknown"))
            message = f"{sample}: {type(exc).__name__}: {exc}"
            errors.append(message)
            print(f"ERROR: {message}")

    if errors:
        raise RuntimeError(
            f"Input validation failed for {len(errors)} sample(s); no outputs were written"
        )

    audit_frame = pd.DataFrame(
        row for inspection in inspections for row in inspection.audit_rows
    )
    sample_frame = pd.DataFrame(
        row for inspection in inspections for row in inspection.sample_rows
    )
    cluster_frame = pd.DataFrame(
        row for inspection in inspections for row in inspection.cluster_rows
    )
    scale_frame = derive_shared_scales(
        pd.DataFrame(row for inspection in inspections for row in inspection.scale_rows)
    )
    warnings = [warning for inspection in inspections for warning in inspection.warnings]
    print_validation_summary(audit_frame, warnings)

    if args.validate_only:
        print("Validation-only mode: no outputs written")
        return

    prefix, data_dir, figure_dir, plot_cfg = output_settings(cfg, repo_root)
    review_frame = build_review_scaffold(
        sample_frame,
        genes,
        categories,
        len(sample_cfgs),
    )
    write_csv(audit_frame, data_dir / f"{prefix}_input_audit.csv", args.overwrite)
    write_csv(
        sample_frame,
        data_dir / f"{prefix}_sample_gene_summary.csv",
        args.overwrite,
    )
    write_csv(
        cluster_frame,
        data_dir / f"{prefix}_cluster_gene_summary.csv",
        args.overwrite,
    )
    write_csv(
        scale_frame,
        data_dir / f"{prefix}_spatial_scale_limits.csv",
        args.overwrite,
    )
    write_csv(
        review_frame,
        data_dir / f"{prefix}_marker_review.csv",
        args.overwrite,
    )

    sample_order = [str(item["sample"]) for item in sample_cfgs]
    if args.mode in {"all", "cluster"}:
        for sample_cfg in sample_cfgs:
            sample = str(sample_cfg["sample"])
            resolution = str(sample_cfg.get("resolution", cfg.get("resolution", "")))
            subset = cluster_frame[cluster_frame["sample"] == sample]
            sample_dir = figure_dir / "cluster_qc" / sample
            plot_cluster_dotplot(
                subset,
                genes,
                sample,
                resolution,
                sample_dir / f"{sample}_cluster_marker_qc_absolute_r{safe_label(resolution)}",
                plot_cfg,
                args.overwrite,
                zscore=False,
            )
            if bool(plot_cfg.get("write_zscore_dotplot", True)):
                plot_cluster_dotplot(
                    subset,
                    genes,
                    sample,
                    resolution,
                    sample_dir
                    / f"{sample}_cluster_marker_qc_zscore_r{safe_label(resolution)}",
                    plot_cfg,
                    args.overwrite,
                    zscore=True,
                )

        cross_dir = figure_dir / "cross_sample_qc"
        plot_sample_gene_heatmap(
            sample_frame,
            genes,
            sample_order,
            "percent_expressing",
            "Cells expressing (%)",
            cross_dir / f"{prefix}_sample_gene_detection_heatmap",
            plot_cfg,
            args.overwrite,
            vmin=0,
            vmax=100,
        )
        plot_sample_gene_heatmap(
            sample_frame,
            genes,
            sample_order,
            "mean_expression",
            "Mean normalized log1p expression",
            cross_dir / f"{prefix}_sample_gene_mean_expression_heatmap",
            plot_cfg,
            args.overwrite,
            vmin=0,
        )

    if args.mode in {"all", "spatial"}:
        for sample_cfg in sample_cfgs:
            coords, matrix, present_genes, clusters, sample, resolution = load_plot_arrays(
                sample_cfg,
                cfg,
                genes,
                input_root,
            )
            sample_scales = scale_frame[scale_frame["sample"] == sample].set_index("gene")
            within_vmax = sample_scales["within_sample_vmax"].to_dict()
            shared_vmax = sample_scales["shared_vmax"].to_dict()
            sample_dir = figure_dir / "spatial_qc" / sample
            plot_cluster_map(
                coords,
                clusters,
                sample,
                resolution,
                sample_dir / f"{sample}_banksy_clusters_r{safe_label(resolution)}",
                plot_cfg,
                args.overwrite,
            )
            plot_expression_panel(
                coords,
                matrix,
                present_genes,
                within_vmax,
                sample,
                "within-sample positive-cell quantile",
                sample_dir / f"{sample}_marker_expression_within_sample_scale",
                plot_cfg,
                args.overwrite,
            )
            plot_expression_panel(
                coords,
                matrix,
                present_genes,
                shared_vmax,
                sample,
                "shared per-gene scale",
                sample_dir / f"{sample}_marker_expression_shared_scale",
                plot_cfg,
                args.overwrite,
            )
            plot_detection_panel(
                coords,
                matrix,
                present_genes,
                sample,
                sample_dir / f"{sample}_marker_detection",
                plot_cfg,
                args.overwrite,
            )

    write_manifest(
        data_dir / f"{prefix}_run_manifest.json",
        config_path,
        args,
        genes,
        sample_order,
        args.overwrite,
    )
    print("Marker-expression spatial QC completed")


if __name__ == "__main__":
    main()
