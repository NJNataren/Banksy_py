#!/usr/bin/env python
# coding: utf-8

"""
Title: Shared Inputs And Styling For Cross-Sample Spatial Helpers
Date: 2026-10-26
Summary: Describe each sample's inputs from the cohort config (via
xenium_workflow/cohorts.py), read only cell-type labels and coordinates from H5AD
files, and provide the shared response colours and figure export used by the
cross-sample spatial helpers in this directory. Updated 2026-11-19: the
helper -> script 10 -> script 08 config chain is replaced by one cohort config.
Updated 2026-11-22: helpers name an annotation set and label level
(xenium_workflow/annotation_sets.py); labels are read from `obs/cell_type_<tier>`, inputs and
outputs use derived paths, and each annotated object's stamp is checked.
"""

import json
import sys
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

# The shared xenium_workflow package lives at the repository root, two levels up.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from xenium_workflow.annotation_sets import (  # noqa: E402
    add_annotation_arguments, annotated_h5ad_path, check_stamp, label_column, level_vocabulary,
    require_released, resolve_step, step_dir,
)
from xenium_workflow.spatial_tables import enrichment_table_path  # noqa: E402
from xenium_workflow.cohorts import (  # noqa: E402
    add_root_arguments, resolve_path, select_samples,
)

STEP = "spatial_analysis"


def add_step_arguments(parser):
    """Root, annotation-set, and draft options shared by the spatial helpers."""
    add_root_arguments(parser)
    add_annotation_arguments(parser)
    parser.add_argument("--allow-draft", action="store_true", help="Use a draft annotation set")


# Arial when installed (as in script 10); otherwise fall back quietly.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]

RESPONSE_ORDER = ["res", "non_res"]
RESPONSE_COLORS = {"res": "#2a78d6", "non_res": "#eb6834"}
TISSUE_MARKERS = ["o", "s", "^", "D", "v", "P", "X"]
TEXT_DARK = "#0b0b0b"
TEXT_MUTED = "#52514e"
SURFACE = "#fcfcfb"


def load_sample_records(cfg, args):
    """Describe each sample a cross-sample helper uses, from its annotation set.

    `cfg` names an `annotation` set, an optional `label_level` (default tier_1),
    an optional `output_subdir`, and may list `samples` (names, in plot order).
    Each record carries response, tissue, branch, resolution, the script 08
    object under `--input-root`, and the script 10 neighbourhood enrichment table
    under `--input-root` (written by script 10 or backfilled by script 10a).
    Returns a context dict (label column, level vocabulary, coordinate keys,
    roots, derived figure and table folders) and the records.
    """
    annotation, level = resolve_step(cfg, args)
    require_released(annotation, getattr(args, "allow_draft", False), "This helper")
    order, colors = level_vocabulary(annotation, level)
    input_root, output_root = args.input_root.resolve(), args.output_root.resolve()
    result_dir = resolve_path(input_root, step_dir(annotation, "processed", STEP, level))
    subdir = cfg.get("output_subdir")

    records = []
    for sample in select_samples(annotation["cohort"], cfg.get("samples"), args.sample):
        if sample.get("response") not in RESPONSE_ORDER:
            raise ValueError(f"{sample['sample']}: cohort response must be one of {RESPONSE_ORDER}")
        relative = annotated_h5ad_path(annotation, sample)
        records.append({
            "sample": sample["sample"],
            "response": sample["response"],
            "tissue": sample.get("tissue", "unknown"),
            "branch": sample["branch"],
            "resolution": sample["resolution"],
            "annotated_h5ad": resolve_path(input_root, relative),
            "enrichment_table": enrichment_table_path(result_dir, sample["sample"], relative),
        })
    print(f"Annotation set {annotation['version']} ({annotation['status']}), level {level}")
    context = {
        "annotation": annotation,
        "level": level,
        "column": label_column(level),
        "order": order,
        "colors": colors,
        "coord_keys": annotation["cohort"]["coord_keys"],
        "input_root": input_root,
        "output_root": output_root,
        "figure_dir": resolve_path(output_root, step_dir(annotation, "output", STEP, level, subdir)),
        "table_dir": resolve_path(output_root, step_dir(annotation, "processed", STEP, level, subdir)),
    }
    return context, records


def read_labels_and_coords(annotated_h5ad, context):
    """Read only the level's label column and x/y coordinates from an annotated H5AD.

    Checks the object's annotation stamp against the set in `context`. Returns
    string labels, an (n, 2) coordinate array, and the present labels in
    vocabulary order (absent labels dropped, as in script 10).
    """
    column, order = context["column"], context["order"]
    with h5py.File(annotated_h5ad, "r") as handle:
        found = read_elem(handle["uns/cell_type_annotation"]) if "uns/cell_type_annotation" in handle else {}
        check_stamp(found, context["annotation"], str(annotated_h5ad))
        labels = read_elem(handle[f"obs/{column}"])
        coords = np.column_stack([read_elem(handle[f"obs/{key}"]) for key in context["coord_keys"][:2]])
    if pd.isna(labels).any() or not np.isfinite(coords).all():
        raise ValueError(f"{annotated_h5ad}: missing labels or nonfinite coordinates")
    labels = pd.Series(labels).astype(str).to_numpy()
    if set(labels) - set(order):
        raise ValueError(f"{annotated_h5ad}: labels outside the {context['level']} vocabulary")
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
