#!/usr/bin/env python
# coding: utf-8

"""
Title: Compare Xenium Gene Panels
Date: 2026-06-17
Summary: Compare the melanoma Xenium gene panel separately against the hImmuno
and hProstate submitted designer panels. The script creates one overlap
DataFrame for each pairwise comparison, concatenates them, and writes the
combined pairwise overlap table to CSV.
"""

import argparse
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
GENE_MARKER_DIR = REPO_ROOT / "data/xenium/raw_data/gene_markers"
PROSTATE_COMPARISON_DIR = GENE_MARKER_DIR / "prostate_panel_comparison"
DEFAULT_OUTPUT_DIR = PROSTATE_COMPARISON_DIR / "outputs"

DEFAULT_MELANOMA_PANEL = (
    GENE_MARKER_DIR
    / "03_xenium_gene_review_for_dotplot_manual_annotation_2026-06-11.csv"
)
DEFAULT_HIMMUNO_PANEL = (
    PROSTATE_COMPARISON_DIR / "hImmuno_submitted_to_designer_tier1.csv"
)
DEFAULT_HPROSTATE_PANEL = (
    PROSTATE_COMPARISON_DIR / "hProstate_100g_submitted_to_designer.csv"
)
DEFAULT_OUTPUT_CSV = DEFAULT_OUTPUT_DIR / "melanoma_pairwise_panel_overlaps.csv"

MELANOMA_COLUMNS = [
    "Gene",
    "Tier_1_annotation",
    "Tier_2_annotation",
    "Canonical_markers",
    "primary_annotation",
    "secondary_annotation",
    "keep_for_dotplot",
    "drop_reason",
    "manual_review_True_if_keep",
]

COMPARISON_PANEL_CONFIGS = {
    "hImmuno": {
        "path_arg": "himmuno_panel",
        "columns": [
            "Gene",
            "Ensemble ID",
            "Num_Probesets",
            "Codewords",
            "Annotation",
            "Tier_1",
            "present_in_panel_y",
        ],
    },
    "hProstate_100g": {
        "path_arg": "hprostate_panel",
        "columns": ["panel", "Gene", "gene_id", "present_in_panel_y"],
    },
}


def parse_args():
    """Parse command-line arguments for pairwise panel comparison."""
    parser = argparse.ArgumentParser(
        description=(
            "Find pairwise overlaps between the melanoma Xenium panel and the "
            "hImmuno/hProstate designer panels. No three-way comparison is made."
        )
    )
    parser.add_argument(
        "--melanoma-panel",
        default=DEFAULT_MELANOMA_PANEL,
        type=Path,
        help="Melanoma Xenium manual-review gene panel CSV.",
    )
    parser.add_argument(
        "--himmuno-panel",
        default=DEFAULT_HIMMUNO_PANEL,
        type=Path,
        help="hImmuno submitted designer panel CSV.",
    )
    parser.add_argument(
        "--hprostate-panel",
        default=DEFAULT_HPROSTATE_PANEL,
        type=Path,
        help="hProstate 100g submitted designer panel CSV.",
    )
    parser.add_argument(
        "--output-csv",
        default=DEFAULT_OUTPUT_CSV,
        type=Path,
        help="Output CSV for concatenated pairwise overlaps.",
    )
    parser.add_argument(
        "--manual-review-kept-only",
        action="store_true",
        help=(
            "Restrict melanoma genes to manual_review_True_if_keep == TRUE before "
            "calculating overlaps."
        ),
    )
    parser.add_argument(
        "--include-absent-designer-rows",
        action="store_true",
        help=(
            "Keep comparison-panel rows marked present_in_panel_y == FALSE. By "
            "default, those rows are excluded before overlap calculation."
        ),
    )
    return parser.parse_args()


def normalize_gene_symbol(gene_symbol):
    """Return a trimmed uppercase gene symbol for matching across files."""
    if pd.isna(gene_symbol):
        return pd.NA
    normalized = str(gene_symbol).strip().upper()
    return normalized if normalized else pd.NA


def collapse_values(values):
    """Collapse unique non-empty values into a semicolon-delimited string."""
    cleaned = []
    for value in values:
        if pd.isna(value):
            continue
        text = str(value).strip()
        if text and text.lower() != "nan":
            cleaned.append(text)
    return "; ".join(dict.fromkeys(cleaned))


def load_melanoma_panel(path, manual_review_kept_only=False):
    """Load the melanoma panel and keep one row per normalized gene symbol.

    Args:
        path: Melanoma panel CSV path.
        manual_review_kept_only: Whether to keep only genes marked TRUE in
            `manual_review_True_if_keep`.

    Returns:
        DataFrame with melanoma metadata and a `gene_key` match column.
    """
    df = pd.read_csv(path)
    if "Gene" not in df.columns:
        raise ValueError(f"{path} is missing required column 'Gene'")

    if manual_review_kept_only:
        if "manual_review_True_if_keep" not in df.columns:
            raise ValueError(
                "--manual-review-kept-only requires 'manual_review_True_if_keep' "
                "in the melanoma panel."
            )
        keep = (
            df["manual_review_True_if_keep"]
            .astype(str)
            .str.strip()
            .str.upper()
            .eq("TRUE")
        )
        df = df.loc[keep].copy()

    columns = [column for column in MELANOMA_COLUMNS if column in df.columns]
    melanoma = df.loc[:, columns].copy()
    melanoma["gene_key"] = melanoma["Gene"].map(normalize_gene_symbol)
    melanoma = melanoma.dropna(subset=["gene_key"])

    aggregation = {column: collapse_values for column in columns}
    melanoma = melanoma.groupby("gene_key", as_index=False).agg(aggregation)
    return melanoma.rename(columns={column: f"melanoma_{column}" for column in columns})


def load_comparison_panel(path, panel_name, include_absent_designer_rows=False):
    """Load a comparison panel and keep one row per normalized gene symbol.

    Args:
        path: Comparison panel CSV path.
        panel_name: Output label for the comparison panel.
        include_absent_designer_rows: Whether to retain designer rows marked
            absent by `present_in_panel_y`.

    Returns:
        DataFrame with comparison-panel metadata and a `gene_key` match column.
    """
    df = pd.read_csv(path)
    if "Gene" not in df.columns:
        raise ValueError(f"{path} is missing required column 'Gene'")

    if "present_in_panel_y" in df.columns and not include_absent_designer_rows:
        present = df["present_in_panel_y"].astype(str).str.strip().str.upper()
        df = df.loc[present.isin(["TRUE", "1", "YES", "Y"])].copy()

    configured_columns = COMPARISON_PANEL_CONFIGS[panel_name]["columns"]
    columns = [column for column in configured_columns if column in df.columns]
    comparison = df.loc[:, columns].copy()
    comparison["gene_key"] = comparison["Gene"].map(normalize_gene_symbol)
    comparison = comparison.dropna(subset=["gene_key"])

    aggregation = {column: collapse_values for column in columns}
    comparison = comparison.groupby("gene_key", as_index=False).agg(aggregation)
    return comparison.rename(
        columns={column: f"comparison_panel_{column}" for column in columns}
    )


def build_pairwise_overlap(melanoma, comparison_panel, comparison_name):
    """Create one pairwise overlap table for melanoma versus another panel.

    Args:
        melanoma: Standardized melanoma panel DataFrame.
        comparison_panel: Standardized comparison panel DataFrame.
        comparison_name: Label describing the pairwise comparison.

    Returns:
        DataFrame containing genes present in both pairwise inputs.
    """
    overlap = melanoma.merge(comparison_panel, on="gene_key", how="inner")
    overlap.insert(0, "comparison", comparison_name)
    overlap.insert(1, "Gene", overlap["gene_key"])
    return overlap.drop(columns=["gene_key"]).sort_values(["comparison", "Gene"])


def build_summary(pairwise_overlaps):
    """Summarize the concatenated pairwise overlap table."""
    summary = (
        pairwise_overlaps.groupby("comparison", as_index=False)
        .agg(overlap_gene_count=("Gene", "nunique"))
        .sort_values("comparison")
    )
    return summary


def main():
    """Run the pairwise panel overlap workflow."""
    args = parse_args()

    melanoma = load_melanoma_panel(
        args.melanoma_panel,
        manual_review_kept_only=args.manual_review_kept_only,
    )

    overlaps = []
    for panel_name, config in COMPARISON_PANEL_CONFIGS.items():
        panel_path = getattr(args, config["path_arg"])
        comparison_panel = load_comparison_panel(
            panel_path,
            panel_name,
            include_absent_designer_rows=args.include_absent_designer_rows,
        )
        overlaps.append(
            build_pairwise_overlap(
                melanoma,
                comparison_panel,
                comparison_name=f"melanoma_vs_{panel_name}",
            )
        )

    pairwise_overlaps = pd.concat(overlaps, ignore_index=True)
    summary = build_summary(pairwise_overlaps)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    pairwise_overlaps.to_csv(args.output_csv, index=False)

    print(f"Wrote pairwise overlap table: {args.output_csv}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
