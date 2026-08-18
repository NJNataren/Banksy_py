#!/usr/bin/env python
# coding: utf-8

"""
Title: Apply QC Filters For Xenium Reclustering
Date: 2026-08-18
Summary: Read a script 01 QC-annotated Xenium AnnData object, apply the
reviewed filtered_qc_v1 cell filter masks, and write a provenance-preserving
AnnData object with all cells retained for downstream reclustering.
"""

import argparse
import json
import os

import anndata as ad
import pandas as pd


def parse_args():
    """Parse command-line arguments for the QC filtering step."""
    parser = argparse.ArgumentParser(
        prog="apply reviewed QC filters to a QC-annotated Xenium AnnData object"
    )
    parser.add_argument(
        "--config",
        type=str,
        help="JSON config with project, dataset_name, and optional input/output paths.",
        required=True,
    )
    return parser.parse_args()


args = parse_args()

# Each sample has a small JSON config. Keeping paths and labels in config files
# makes this script reusable across VBCT, PTMT, and PC55 runs.
with open(args.config) as f:
    cfg = json.load(f)

# These three values define the sample being annotated and the QC decision label
# that will be carried into script 06.
project = cfg["project"]
dataset_name = cfg["dataset_name"]
output_label = cfg.get("output_label", "filtered_qc_v1")

# `base_dir` defaults to the project Xenium data root, but can be overridden
# for local smoke tests or temporary output checks.
base_dir = cfg.get("base_dir", "data/xenium")
processed_path = os.path.join(base_dir, "processed", project, dataset_name)
output_dir = cfg.get(
    "output_dir",
    os.path.join(base_dir, "output", project, "QC_filtering", dataset_name),
)

# Input is the script 01 QC-annotated clean expression object. It should still
# contain every non-zero-count cell plus the QC columns produced during review.
input_h5ad = cfg.get(
    "input_h5ad",
    os.path.join(
        processed_path,
        f"adata_expression_clean_{dataset_name}_qc_annotated.h5ad",
    ),
)

# Output is also a full-cell clean expression object. The name says
# `qc_annotated` deliberately: no hard-filtered AnnData is written here.
annotated_output_h5ad = cfg.get(
    "annotated_output_h5ad",
    os.path.join(
        processed_path,
        f"adata_expression_clean_{dataset_name}_qc_annotated_{output_label}.h5ad",
    ),
)

# Older drafts of script 05 wrote a hard-filtered object. Block that now so a
# stale config cannot silently drop cells and break provenance.
if "filtered_output_h5ad" in cfg:
    raise ValueError(
        "Script 05 no longer writes hard-filtered AnnData outputs. Remove "
        "`filtered_output_h5ad` from the config and use the annotated output "
        "with `qc_keep_for_reclustering` for downstream reclustering."
    )

summary_csv = os.path.join(
    output_dir,
    f"{dataset_name}_qc_filter_summary_{output_label}.csv",
)

os.makedirs(processed_path, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

print(f"Reading QC-annotated AnnData from: {input_h5ad}")
adata = ad.read_h5ad(input_h5ad)

# These are the reviewed QC masks produced upstream. The script refuses to
# continue if any are missing or contain NA values, because ambiguous QC state
# would make the reclustering subset non-auditable.
required_cols = [
    "min_trans_passed",
    "max_trans_threshold_passed",
    "negative_control_probe_ge2",
    "max_area_threshold_99_by_cluster",
]

missing_cols = [col for col in required_cols if col not in adata.obs.columns]
if missing_cols:
    raise KeyError(f"Missing required QC columns in adata.obs: {missing_cols}")

null_cols = [col for col in required_cols if adata.obs[col].isna().any()]
if null_cols:
    raise ValueError(f"Required QC columns contain missing values: {null_cols}")

# Mask semantics are explicit:
# - min_trans_passed is a pass mask, so True means keep.
# - max_trans_threshold_passed, negative_control_probe_ge2, and
#   max_area_threshold_99_by_cluster are fail masks, so True means remove.
keep_mask = (
    adata.obs["min_trans_passed"].astype(bool)
    & ~adata.obs["max_trans_threshold_passed"].astype(bool)
    & ~adata.obs["negative_control_probe_ge2"].astype(bool)
    & ~adata.obs["max_area_threshold_99_by_cluster"].astype(bool)
)

# This is the single column script 06 should use for its temporary in-memory
# BANKSY subset. Cells with False remain in the saved object.
adata.obs["qc_keep_for_reclustering"] = keep_mask
adata.obs["qc_filter_status"] = (
    keep_mask.map({True: "Pass", False: "Fail"}).astype("category")
)

# Store every individual fail mask with a common prefix so later plotting or
# audits can ask which exact rule excluded a cell.
filter_fail_masks = {
    "min_trans_passed": ~adata.obs["min_trans_passed"].astype(bool),
    "max_trans_threshold_passed": adata.obs["max_trans_threshold_passed"].astype(bool),
    "negative_control_probe_ge2": adata.obs["negative_control_probe_ge2"].astype(bool),
    "max_area_threshold_99_by_cluster": adata.obs[
        "max_area_threshold_99_by_cluster"
    ].astype(bool),
}

# Preserve all cells in the script 05 output. Script 06 can subset in memory
# using `qc_keep_for_reclustering`, then copy reclustering labels back onto the
# full object with failed cells marked as excluded_by_qc.
for filter_name, fail_mask in filter_fail_masks.items():
    adata.obs[f"qc_fail_{filter_name}"] = fail_mask

# `qc_fail_reason` stores the first failed rule for simple colour plots.
# `qc_fail_reason_set` stores all failed rules for cells that fail more than one
# criterion.
fail_reasons = []
fail_reason_sets = []
for idx in adata.obs.index:
    reasons = [
        filter_name
        for filter_name, fail_mask in filter_fail_masks.items()
        if bool(fail_mask.loc[idx])
    ]
    fail_reasons.append(reasons[0] if reasons else "none")
    fail_reason_sets.append(";".join(reasons) if reasons else "none")

adata.obs["qc_fail_reason"] = pd.Categorical(fail_reasons)
adata.obs["qc_fail_reason_set"] = pd.Categorical(fail_reason_sets)

# The summary CSV is an audit table: per-filter failed counts plus the combined
# mask count. Importantly, `n_cells_retained_in_output` should equal the input
# cell count because script 05 no longer removes cells.
summary = pd.DataFrame(
    [
        {
            "sample": dataset_name,
            "output_label": output_label,
            "filter_name": "min_trans_passed",
            "filter_semantics": "pass_filter_keep_true",
            "n_failed": int(filter_fail_masks["min_trans_passed"].sum()),
        },
        {
            "sample": dataset_name,
            "output_label": output_label,
            "filter_name": "max_trans_threshold_passed",
            "filter_semantics": "fail_filter_remove_true",
            "n_failed": int(filter_fail_masks["max_trans_threshold_passed"].sum()),
        },
        {
            "sample": dataset_name,
            "output_label": output_label,
            "filter_name": "negative_control_probe_ge2",
            "filter_semantics": "fail_filter_remove_true",
            "n_failed": int(filter_fail_masks["negative_control_probe_ge2"].sum()),
        },
        {
            "sample": dataset_name,
            "output_label": output_label,
            "filter_name": "max_area_threshold_99_by_cluster",
            "filter_semantics": "fail_filter_remove_true",
            "n_failed": int(
                filter_fail_masks["max_area_threshold_99_by_cluster"].sum()
            ),
        },
        {
            "sample": dataset_name,
            "output_label": output_label,
            "filter_name": "combined_qc_filter",
            "filter_semantics": "combined_keep_mask",
            "n_failed": int((~keep_mask).sum()),
        },
    ]
)

summary["n_cells_before"] = int(adata.n_obs)
summary["n_cells_retained_in_output"] = int(adata.n_obs)
summary["n_cells_marked_keep_for_reclustering"] = int(keep_mask.sum())
summary["n_cells_marked_excluded_by_qc"] = int((~keep_mask).sum())
summary["percent_marked_excluded_by_qc"] = (
    summary["n_cells_marked_excluded_by_qc"] / summary["n_cells_before"] * 100
)
summary["percent_failed"] = summary["n_failed"] / summary["n_cells_before"] * 100

summary.to_csv(summary_csv, index=False)

print(summary)
print(
    f"Marked {int(keep_mask.sum()):,} / {adata.n_obs:,} cells "
    f"({keep_mask.sum() / adata.n_obs * 100:.2f}%) as keep-for-reclustering."
)
print(
    f"Retaining all {adata.n_obs:,} cells in the script 05 output for provenance."
)

# Save the full provenance-preserving AnnData. Script 06 will subset this object
# temporarily, but this file itself keeps all cells and all QC annotations.
print(f"Saving QC-filter annotated AnnData to: {annotated_output_h5ad}")
adata.write_h5ad(annotated_output_h5ad)


print(f"Saved QC filter summary to: {summary_csv}")
