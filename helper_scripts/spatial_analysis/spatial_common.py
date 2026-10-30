#!/usr/bin/env python
# coding: utf-8

"""
Title: Shared Inputs And Styling For Cross-Sample Spatial Helpers
Date: 2026-10-26
Summary: Follow the helper -> script 10 -> script 08 config chain, describe
each sample's inputs, read only cell-type labels and coordinates from H5AD
files, and provide the shared response colours and figure export used by the
cross-sample spatial helpers in this directory.
"""

import json
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rc_context

try:
    from anndata.io import read_elem
except ImportError:
    from anndata.experimental import read_elem


# Arial when installed (as in script 10); otherwise fall back quietly.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]

RESPONSE_ORDER = ["res", "non_res"]
RESPONSE_COLORS = {"res": "#2a78d6", "non_res": "#eb6834"}
TISSUE_MARKERS = ["o", "s", "^", "D", "v", "P", "X"]
TEXT_DARK = "#0b0b0b"
TEXT_MUTED = "#52514e"
SURFACE = "#fcfcfb"


def resolve_path(root, value):
    """Resolve a configured absolute path or a path relative to its project root."""
    path = Path(value)
    return path if path.is_absolute() else root / path


def load_sample_records(cfg, config_path, selected=None):
    """Follow the script 10 -> script 08 config chain and describe each sample.

    `cfg` is a cross-sample comparison config with `spatial_config` and a
    `samples` mapping of response/tissue labels. Returns a context dict (project
    root, palette order and colours, coordinate keys) and one record per sample.
    """
    root = config_path.parents[int(cfg.get("project_root_parent_depth", 4))]
    spatial_path = resolve_path(root, cfg["spatial_config"])
    spatial = json.loads(spatial_path.read_text())
    spatial_root = spatial_path.parents[int(spatial.get("project_root_parent_depth", 4))]
    annotation_path = resolve_path(spatial_root, spatial["annotation_config"])
    annotation = json.loads(annotation_path.read_text())
    annotation_project = annotation_path.parents[int(annotation.get("project_root_parent_depth", 3))]
    annotation_root = Path(annotation.get("output_root", annotation_project))
    palette = json.loads(resolve_path(annotation_project, annotation["palette_json"]).read_text())
    order = palette["cell_type_order"]

    # Script 10 results may live on another machine (e.g. the HPC mount).
    result_root = Path(cfg.get("result_root", spatial_root))
    result_dir = resolve_path(result_root, spatial["result_dir"])

    annotated = {sample["sample"]: sample for sample in annotation["samples"]}
    unknown = set(cfg["samples"]) - set(annotated)
    if unknown:
        raise ValueError(f"Samples absent from script 08 config: {sorted(unknown)}")
    if selected and set(selected) - set(cfg["samples"]):
        raise ValueError("Requested samples are absent from the comparison config")
    for meta in cfg["samples"].values():
        if meta["response"] not in RESPONSE_ORDER:
            raise ValueError(f"response must be one of {RESPONSE_ORDER}: {meta}")

    records = []
    for name, meta in cfg["samples"].items():
        if selected and name not in selected:
            continue
        sample = annotated[name]
        source = resolve_path(annotation_root, sample["output_h5ad"])
        records.append({
            "sample": name,
            "response": meta["response"],
            "tissue": meta.get("tissue", "unknown"),
            "branch": sample["branch"],
            "resolution": sample["resolution"],
            "annotated_h5ad": source,
            "result_h5ad": result_dir / name / f"{source.stem}_squidpy.h5ad",
        })
    context = {
        "root": root,
        "order": order,
        "colors": palette["colors"],
        "coord_keys": spatial.get("coord_keys", ["x", "y"]),
    }
    return context, records


def read_labels_and_coords(annotated_h5ad, order, coord_keys):
    """Read only `obs['cell_type']` and x/y coordinates from an annotated H5AD.

    Returns string labels, an (n, 2) coordinate array, and the present cell
    types in palette order (absent types dropped, as in script 10).
    """
    with h5py.File(annotated_h5ad, "r") as handle:
        labels = read_elem(handle["obs/cell_type"])
        coords = np.column_stack([read_elem(handle[f"obs/{key}"]) for key in coord_keys[:2]])
    if pd.isna(labels).any() or not np.isfinite(coords).all():
        raise ValueError(f"{annotated_h5ad}: missing labels or nonfinite coordinates")
    labels = pd.Series(labels).astype(str).to_numpy()
    if set(labels) - set(order):
        raise ValueError(f"{annotated_h5ad}: unknown cell types")
    present = [label for label in order if label in set(labels)]
    return labels, coords, present


def tissue_markers(records):
    """Assign one fixed marker per tissue in config order."""
    tissues = list(dict.fromkeys(r["tissue"] for r in records))
    return {tissue: TISSUE_MARKERS[i % len(TISSUE_MARKERS)] for i, tissue in enumerate(tissues)}


def sample_title(record):
    """Short panel title with sample, response group, and resolution."""
    return f"{record['sample']}\n{record['response']} · r{record['resolution']}"


def spread_label_positions(values, min_gap):
    """Nudge end-of-line label y-positions apart, moving each label as little as possible.

    Labels are placed bottom-up with a minimum gap; each run of touching labels
    is then centred on its own lines, so isolated labels stay level with them.
    """
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    targets = values[order]
    clusters = []  # each: [start, end) index range into targets
    for k in range(len(targets)):
        clusters.append([k, k + 1])
        # Merge while the latest cluster, centred on its targets, collides with the previous one.
        while len(clusters) > 1:
            (s0, e0), (s1, e1) = clusters[-2], clusters[-1]
            top_prev = targets[s0:e0].mean() + (e0 - s0 - 1) * min_gap / 2
            bottom_last = targets[s1:e1].mean() - (e1 - s1 - 1) * min_gap / 2
            if bottom_last - top_prev >= min_gap:
                break
            clusters[-2:] = [[s0, e1]]
    placed = np.empty_like(targets)
    for start, end in clusters:
        centre = targets[start:end].mean()
        placed[start:end] = centre + (np.arange(end - start) - (end - start - 1) / 2) * min_gap
    result = np.empty_like(placed)
    result[order] = placed
    return result


def save_figure(fig, output_stem, cfg):
    """Save all requested formats with editable text, then close the figure."""
    try:
        with rc_context({"svg.fonttype": "none", "pdf.fonttype": 42}):
            for fmt in cfg.get("figure_formats", ["png", "svg", "pdf"]):
                path = f"{output_stem}.{fmt}"
                fig.savefig(path, dpi=cfg.get("figure_dpi", 300), bbox_inches="tight")
                print(f"Saved figure: {path}")
    finally:
        plt.close(fig)
