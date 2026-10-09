#!/usr/bin/env python
# coding: utf-8

"""
Title: Apply QC Filters For Xenium Reclustering
Date: 2026-10-09
Summary: Read a script 01 QC-annotated Xenium AnnData object, apply the
reviewed filtered_qc_v1 cell filter masks, and write a provenance-preserving
AnnData object with all cells retained for downstream reclustering. Configurable
area-percentile filters run after the minimum-transcript mask, use
percentile-neutral output columns, and write embedded and sidecar provenance
describing the resolved threshold settings.
"""

import argparse
import json
import os

import anndata as ad
import numpy as np
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


AREA_FAIL_COL = "area_percentile_filter_fail"
AREA_STATUS_COL = "area_percentile_filter_status"
AREA_GROUP_COL = "area_percentile_group"
AREA_THRESHOLD_COL = "area_percentile_threshold"


def make_json_serializable(value):
    """Return `value` converted to objects accepted by `json.dump`."""
    if isinstance(value, dict):
        return {str(key): make_json_serializable(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [make_json_serializable(item) for item in value]
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def validate_upstream_qc_columns(adata):
    """Validate script 01 QC masks required before script 05 filtering."""
    required_cols = [
        "min_trans_passed",
        "max_trans_threshold_passed",
        "negative_control_probe_ge2",
    ]
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        raise KeyError(
            f"Missing required upstream QC columns in adata.obs: {missing_cols}"
        )

    null_cols = [col for col in required_cols if adata.obs[col].isna().any()]
    if null_cols:
        raise ValueError(
            f"Required upstream QC columns contain missing values: {null_cols}"
        )


def build_min_transcript_mask(adata):
    """Return the script 01 minimum-transcript pass mask as a boolean Series."""
    return adata.obs["min_trans_passed"].astype(bool)


def resolve_area_filter_config(cfg):
    """Return a normalized area-filter config dictionary."""
    area_cfg = cfg.get("area_percentile_filter")
    if area_cfg is None:
        source_col = cfg.get(
            "existing_area_percentile_filter_column",
            "max_area_threshold_99_by_cluster",
        )
        return {
            "mode": "existing_script01_mask",
            "source_column": source_col,
            "percentile": 0.99,
            "upper_tail_fraction": 0.01,
            "groupby": "not_recorded_by_script01",
            "group_map": {},
        }

    mode = area_cfg.get("mode", "cluster")
    if mode not in {"cluster", "cell_type_group"}:
        raise ValueError(
            "area_percentile_filter.mode must be either 'cluster' or "
            "'cell_type_group'."
        )

    percentile = float(area_cfg.get("percentile", 0.99))
    if not 0 < percentile < 1:
        raise ValueError("area_percentile_filter.percentile must be between 0 and 1.")

    groupby = area_cfg["groupby"]
    group_map = {
        str(key): str(value)
        for key, value in area_cfg.get("group_map", {}).items()
    }
    if mode == "cell_type_group" and not group_map:
        raise ValueError(
            "area_percentile_filter.group_map is required when mode is "
            "'cell_type_group'."
        )

    return {
        "mode": mode,
        "source_column": AREA_FAIL_COL,
        "percentile": percentile,
        "upper_tail_fraction": 1.0 - percentile,
        "groupby": groupby,
        "group_map": group_map,
    }


def build_area_groups(adata, area_cfg):
    """Return raw cluster labels and resolved area-filter groups per cell."""
    groupby = area_cfg["groupby"]
    if groupby not in adata.obs.columns:
        raise KeyError(f"Missing groupby column for area-percentile filter: {groupby}")

    raw_group_values = adata.obs[groupby]
    if raw_group_values.isna().any():
        n_bad = int(raw_group_values.isna().sum())
        raise ValueError(f"{groupby} contains {n_bad} missing group labels.")

    raw_group = raw_group_values.astype(str)
    if area_cfg["mode"] == "cluster":
        return raw_group, raw_group

    # Unmapped annotations remain separate rather than being pooled into an
    # ambiguous fallback cell type.
    area_group = raw_group.map(area_cfg["group_map"]).fillna(raw_group)
    return raw_group, area_group


def calculate_area_thresholds(area_values, area_group, eligible_mask, percentile):
    """Calculate per-group area thresholds from eligible cells only."""
    threshold_area_values = area_values.loc[eligible_mask]
    threshold_area_group = area_group.loc[eligible_mask]
    if threshold_area_values.empty:
        raise ValueError(
            "No cells passed the upstream filters required for area percentile "
            "threshold calculation."
        )

    thresholds_by_group = threshold_area_values.groupby(
        threshold_area_group, observed=True
    ).quantile(percentile)
    group_thresholds = area_group.map(thresholds_by_group).astype(float)
    return thresholds_by_group, group_thresholds


def build_area_fail_mask(area_values, group_thresholds, eligible_mask):
    """Return the area fail mask after upstream eligibility is applied."""
    return eligible_mask & (area_values >= group_thresholds)


def write_area_filter_obs(adata, fail_mask, area_group, thresholds):
    """Write standardized area-filter columns to `adata.obs`."""
    adata.obs[AREA_FAIL_COL] = fail_mask
    adata.obs[AREA_STATUS_COL] = (
        fail_mask.map({True: "Fail", False: "Pass"}).astype("category")
    )
    adata.obs[AREA_GROUP_COL] = area_group.astype("category")
    adata.obs[AREA_THRESHOLD_COL] = thresholds


def build_area_group_audit(
    area_group,
    raw_group,
    eligible_mask,
    fail_mask,
    thresholds_by_group,
    percentile,
    mode,
    groupby,
):
    """Build the per-area-group threshold and failure audit table."""
    group_audit_rows = []
    for group_name, group_index in area_group.groupby(area_group).groups.items():
        group_fail = fail_mask.loc[group_index]
        group_eligible = eligible_mask.loc[group_index]
        source_clusters = sorted(raw_group.loc[group_index].unique())
        threshold = thresholds_by_group.get(group_name, float("nan"))
        n_eligible = int(group_eligible.sum())
        group_audit_rows.append(
            {
                "group": str(group_name),
                "source_clusters": ";".join(source_clusters),
                "n_cells": n_eligible,
                "n_cells_all": int(len(group_index)),
                "n_cells_area_filter_eligible": n_eligible,
                "percentile": percentile,
                "upper_tail_fraction": 1.0 - percentile,
                "area_threshold": float(threshold),
                "n_failed": int(group_fail.sum()),
                "percent_failed": (
                    float(group_fail.sum() / n_eligible * 100)
                    if n_eligible > 0
                    else 0.0
                ),
                "group_mode": mode,
                "groupby": groupby,
            }
        )
    return pd.DataFrame(group_audit_rows).sort_values("group")


def import_existing_area_filter(adata, area_cfg):
    """Import a historical script 01 area mask into standardized columns."""
    source_col = area_cfg["source_column"]
    if source_col not in adata.obs.columns:
        raise KeyError(
            "No area_percentile_filter config was provided and the existing "
            f"script 01 mask {source_col!r} is missing from adata.obs."
        )
    if adata.obs[source_col].isna().any():
        raise ValueError(
            f"Existing area filter column {source_col!r} contains NA values."
        )

    fail_mask = adata.obs[source_col].astype(bool)
    area_group = pd.Series("not_recorded_by_script01", index=adata.obs_names)
    thresholds = pd.Series(float("nan"), index=adata.obs_names)
    write_area_filter_obs(adata, fail_mask, area_group, thresholds)

    audit = {
        "area_percentile_filter_mode": area_cfg["mode"],
        "area_percentile_filter_source_column": source_col,
        "area_percentile_filter_percentile": area_cfg["percentile"],
        "area_percentile_filter_upper_tail_fraction": area_cfg["upper_tail_fraction"],
        "area_percentile_filter_groupby": area_cfg["groupby"],
        "area_percentile_filter_n_groups": "not_recorded_by_script01",
        "area_percentile_filter_n_mapped_clusters": 0,
    }
    print(f"Using existing script 01 area filter mask: {source_col}")
    return audit, pd.DataFrame(), area_cfg


def apply_area_percentile_filter(adata, cfg, area_filter_eligible_mask=None):
    """Resolve, calculate, write, and audit the configured area filter."""
    area_cfg = resolve_area_filter_config(cfg)
    if area_cfg["mode"] == "existing_script01_mask":
        return import_existing_area_filter(adata, area_cfg)

    if "cell_area" not in adata.obs.columns:
        raise KeyError("Missing column for area-percentile filter: cell_area")
    area_values = pd.to_numeric(adata.obs["cell_area"], errors="coerce")
    if area_values.isna().any():
        n_bad = int(area_values.isna().sum())
        raise ValueError(f"cell_area contains {n_bad} missing or non-numeric values.")

    raw_group, area_group = build_area_groups(adata, area_cfg)
    if area_filter_eligible_mask is None:
        eligible_mask = pd.Series(True, index=adata.obs_names)
    else:
        eligible_mask = pd.Series(area_filter_eligible_mask, index=adata.obs_names)
    if eligible_mask.isna().any():
        raise ValueError("area_filter_eligible_mask contains missing values.")
    eligible_mask = eligible_mask.astype(bool)

    # Minimum-transcript filtering is applied before area filtering: only cells
    # that already pass that upstream mask define percentile thresholds and only
    # those cells can fail the area filter.
    thresholds_by_group, group_thresholds = calculate_area_thresholds(
        area_values,
        area_group,
        eligible_mask,
        area_cfg["percentile"],
    )
    fail_mask = build_area_fail_mask(area_values, group_thresholds, eligible_mask)
    write_area_filter_obs(adata, fail_mask, area_group, group_thresholds)

    group_audit = build_area_group_audit(
        area_group,
        raw_group,
        eligible_mask,
        fail_mask,
        thresholds_by_group,
        area_cfg["percentile"],
        area_cfg["mode"],
        area_cfg["groupby"],
    )
    audit = {
        "area_percentile_filter_mode": area_cfg["mode"],
        "area_percentile_filter_source_column": AREA_FAIL_COL,
        "area_percentile_filter_percentile": area_cfg["percentile"],
        "area_percentile_filter_upper_tail_fraction": area_cfg["upper_tail_fraction"],
        "area_percentile_filter_groupby": area_cfg["groupby"],
        "area_percentile_filter_n_groups": int(area_group.nunique(dropna=False)),
        "area_percentile_filter_n_mapped_clusters": len(area_cfg["group_map"]),
        "area_percentile_filter_threshold_population": "min_trans_passed_cells",
        "area_percentile_filter_n_threshold_cells": int(eligible_mask.sum()),
    }
    provenance = dict(area_cfg)
    provenance["threshold_population"] = "min_trans_passed_cells"
    provenance["n_threshold_cells"] = int(eligible_mask.sum())

    print(
        "Calculated area percentile filter using "
        f"mode={area_cfg['mode']}, groupby={area_cfg['groupby']}, "
        f"percentile={area_cfg['percentile']}."
    )
    print(
        f"Area-filter groups: {audit['area_percentile_filter_n_groups']} "
        f"({audit['area_percentile_filter_n_mapped_clusters']} mapped cluster labels)."
    )
    return audit, group_audit, provenance

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
input_qc_label = cfg.get("input_qc_label", "")
input_qc_filename = f"adata_expression_clean_{dataset_name}_qc_annotated"
if input_qc_label:
    input_qc_filename = f"{input_qc_filename}_{input_qc_label}"
input_h5ad = cfg.get(
    "input_h5ad",
    os.path.join(
        processed_path,
        f"{input_qc_filename}.h5ad",
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
area_group_audit_csv = os.path.join(
    output_dir,
    f"{dataset_name}_area_percentile_filter_groups_{output_label}.csv",
)
resolved_config_json = os.path.join(
    output_dir,
    f"{dataset_name}_qc_filter_resolved_config_{output_label}.json",
)

os.makedirs(processed_path, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

print(f"Reading QC-annotated AnnData from: {input_h5ad}")
adata = ad.read_h5ad(input_h5ad)

# Validate script 01 masks before computing configured area filters so
# low-transcript cells do not set area thresholds.
validate_upstream_qc_columns(adata)
upstream_required_cols = [
    "min_trans_passed",
    "max_trans_threshold_passed",
    "negative_control_probe_ge2",
]
min_trans_pass_mask = build_min_transcript_mask(adata)
area_filter_audit, area_group_audit, area_filter_provenance = (
    apply_area_percentile_filter(
        adata,
        cfg,
        area_filter_eligible_mask=min_trans_pass_mask,
    )
)

# The area filter column is produced above, either from a configured percentile
# calculation or by importing the historical script 01 area mask.
required_cols = upstream_required_cols + [AREA_FAIL_COL]
missing_cols = [col for col in required_cols if col not in adata.obs.columns]
if missing_cols:
    raise KeyError(f"Missing required QC columns in adata.obs: {missing_cols}")

null_cols = [col for col in required_cols if adata.obs[col].isna().any()]
if null_cols:
    raise ValueError(f"Required QC columns contain missing values: {null_cols}")

# ---------------------------------------------------------------------------- #
#                              FILTER APPLICATION                              #
# ---------------------------------------------------------------------------- #

# Mask semantics are explicit:
# - min_trans_passed is a pass mask, so True means keep.
# - max_trans_threshold_passed, negative_control_probe_ge2, and the generic
#   area-percentile result are fail masks, so True means remove.
keep_mask = (
    min_trans_pass_mask
    & ~adata.obs["max_trans_threshold_passed"].astype(bool)
    & ~adata.obs["negative_control_probe_ge2"].astype(bool)
    & ~adata.obs[AREA_FAIL_COL].astype(bool)
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
    "min_trans_passed": ~min_trans_pass_mask,
    "max_trans_threshold_passed": adata.obs["max_trans_threshold_passed"].astype(bool),
    "negative_control_probe_ge2": adata.obs["negative_control_probe_ge2"].astype(bool),
    "area_percentile_filter": adata.obs[AREA_FAIL_COL].astype(bool),
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
            "filter_name": "area_percentile_filter",
            "filter_semantics": "fail_filter_remove_true",
            "n_failed": int(filter_fail_masks["area_percentile_filter"].sum()),
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
for audit_key, audit_value in area_filter_audit.items():
    summary[audit_key] = audit_value

filter_provenance = {
    "output_label": output_label,
    "script01_input_qc_label": input_qc_label or "default",
    "script01_qc_provenance": adata.uns.get("script01_qc_provenance", {}),
    "active_filters": [
        "minimum_transcripts",
        "maximum_transcripts",
        "negative_control_probe",
        "area_percentile",
    ],
    "filter_columns": {
        "minimum_transcripts": "min_trans_passed",
        "maximum_transcripts": "max_trans_threshold_passed",
        "negative_control_probe": "negative_control_probe_ge2",
        "area_percentile": AREA_FAIL_COL,
    },
    "area_percentile_filter": area_filter_provenance,
}
adata.uns["qc_filter_provenance"] = filter_provenance

# Keep a resolved, human-readable snapshot beside the summary so the exact
# settings remain reviewable without opening the AnnData object.
resolved_config = dict(cfg)
resolved_config["resolved_filter_provenance"] = filter_provenance
with open(resolved_config_json, "w") as f:
    json.dump(make_json_serializable(resolved_config), f, indent=2)

summary.to_csv(summary_csv, index=False)
if not area_group_audit.empty:
    area_group_audit.to_csv(area_group_audit_csv, index=False)

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
if not area_group_audit.empty:
    print(f"Saved area-filter group audit to: {area_group_audit_csv}")
print(f"Saved resolved QC filter config to: {resolved_config_json}")
