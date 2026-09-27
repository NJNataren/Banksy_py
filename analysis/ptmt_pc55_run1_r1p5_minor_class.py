"""
Title: PTMT PC55 run 1 r1.5 comparable minor-class analysis
Date: 2026-09-23
Summary: Reuse the detailed minor-class annotation model and select markers
with the broad workflow's strict/inclusive presence denominators and support
tiers. Outputs are provisional review files; source data and QC configs remain
unchanged.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

import ptmt_pc55_run1_r1p5_marker_selection as detailed


PARENT_LABELS = {"T/NK", "Myeloid", "Stromal", "Epithelial", "B/mast"}
SUPPORT_TIERS = {
    "Pan-cohort",
    "Core conditional",
    "Recurrent conditional",
    "Partial support",
    "Sample-specific",
}


def parse_args() -> argparse.Namespace:
    """Parse repository and Phoenix scratch paths."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument(
        "--scratch-root",
        type=Path,
        default=Path("/home/nnataren/mnt/phoenix_scratch/Banksy_py"),
    )
    return parser.parse_args()


def marker_matches_label(
    gene: str,
    label: str,
    gene_evidence: dict[str, dict[str, float]],
) -> bool:
    """Return whether a reference marker supports a detailed or parent label."""
    allowed_labels = set(gene_evidence.get(str(gene), {}))
    if label in PARENT_LABELS:
        return any(detailed.broad_for(candidate) == label for candidate in allowed_labels)
    return label in allowed_labels


def select_shared_markers(
    review: pd.DataFrame,
    evidence: pd.DataFrame,
    gene_evidence: dict[str, dict[str, float]],
) -> pd.DataFrame:
    """Rank minor-class markers using strict and inclusive presence denominators."""
    assigned = review.query("proposed_label != 'Unresolved'")[
        ["sample", "cluster_id", "proposed_label", "confidence"]
    ]
    positive = evidence.merge(assigned, on=["sample", "cluster_id"], how="inner")
    positive = positive.query("logfoldchanges > 0 and pvals_adj < 0.05").copy()
    rows = []

    for label, clusters in assigned.groupby("proposed_label"):
        target = positive[positive.proposed_label == label]
        other = positive[positive.proposed_label != label]
        inclusive_samples = set(clusters["sample"])
        strict_samples = set(
            clusters.loc[clusters.confidence.isin(["High", "Medium"]), "sample"]
        )
        dmso_samples = {sample for sample in inclusive_samples if sample.endswith("DMSO")}
        treated_samples = inclusive_samples - dmso_samples
        n_clusters = len(clusters)
        candidates = []

        for gene, gene_rows in target.groupby("gene"):
            if (
                gene in detailed.GENERIC_OR_STATE
                or not marker_matches_label(str(gene), str(label), gene_evidence)
            ):
                continue
            marker_samples = set(gene_rows["sample"])
            inclusive_count = len(marker_samples & inclusive_samples)
            strict_count = len(marker_samples & strict_samples)
            cluster_count = len(gene_rows[["sample", "cluster_id"]].drop_duplicates())
            other_count = len(
                other[other.gene == gene][["sample", "cluster_id"]].drop_duplicates()
            )
            inclusive_fraction = inclusive_count / len(inclusive_samples)
            strict_fraction = strict_count / len(strict_samples) if strict_samples else 0.0
            cluster_fraction = cluster_count / n_clusters
            specificity = max(
                0.0,
                cluster_fraction - other_count / max(1, len(assigned) - n_clusters),
            )
            dmso_count = len(marker_samples & dmso_samples)
            treated_count = len(marker_samples & treated_samples)
            dmso_fraction = dmso_count / len(dmso_samples) if dmso_samples else 0.0
            treated_fraction = (
                treated_count / len(treated_samples) if treated_samples else 0.0
            )

            if inclusive_count == review["sample"].nunique():
                tier = "Pan-cohort"
            elif len(strict_samples) >= 3 and strict_count == len(strict_samples):
                tier = "Core conditional"
            elif len(strict_samples) >= 3 and strict_fraction >= 0.75:
                tier = "Recurrent conditional"
            elif inclusive_count >= 2 and inclusive_fraction >= 0.50:
                tier = "Partial support"
            else:
                tier = "Sample-specific"

            selection_score = 1.2 * (
                0.35 * inclusive_fraction
                + 0.25 * cluster_fraction
                + 0.20 * specificity
                + 0.10 * min(1.0, gene_rows.logfoldchanges.median() / 2.0)
                + 0.10 * min(1.0, gene_rows.percent_expressing.median() / 50.0)
            )
            candidates.append({
                "minor_class": label,
                "marker": gene,
                "selection_score": selection_score,
                "support_tier": tier,
                "inclusive_type_sample_n": len(inclusive_samples),
                "inclusive_marker_sample_n": inclusive_count,
                "inclusive_sample_recurrence_fraction": inclusive_fraction,
                "strict_type_sample_n": len(strict_samples),
                "strict_marker_sample_n": strict_count,
                "strict_sample_recurrence_fraction": strict_fraction,
                "dmso_type_sample_n": len(dmso_samples),
                "dmso_marker_sample_n": dmso_count,
                "dmso_sample_recurrence_fraction": dmso_fraction,
                "treated_type_sample_n": len(treated_samples),
                "treated_marker_sample_n": treated_count,
                "treated_sample_recurrence_fraction": treated_fraction,
                "cluster_recurrence_n": cluster_count,
                "cluster_recurrence_fraction": cluster_fraction,
                "median_rank": gene_rows["rank"].median(),
                "median_logfoldchange": gene_rows.logfoldchanges.median(),
                "median_pvals_adj": gene_rows.pvals_adj.median(),
                "median_mean_expression": gene_rows.mean_expression.median(),
                "median_percent_expressing": gene_rows.percent_expressing.median(),
                "specificity_fraction_difference": specificity,
                "present_in_harmonised_reference": True,
                "supporting_samples": "; ".join(sorted(marker_samples)),
                "treatments": "; ".join(sorted({
                    "DMSO" if sample.endswith("DMSO") else "1644"
                    for sample in marker_samples
                })),
                "caveats": (
                    "Treatment-associated; positive in one treatment only."
                    if bool(dmso_count) != bool(treated_count)
                    else ""
                ),
            })

        ranked = pd.DataFrame(candidates)
        if ranked.empty:
            continue
        tier_order = {
            "Pan-cohort": 0,
            "Core conditional": 1,
            "Recurrent conditional": 2,
            "Partial support": 3,
            "Sample-specific": 4,
        }
        ranked["tier_order"] = ranked.support_tier.map(tier_order)
        ranked = ranked.sort_values(
            [
                "tier_order",
                "selection_score",
                "inclusive_marker_sample_n",
                "cluster_recurrence_n",
                "median_rank",
            ],
            ascending=[True, False, False, False, True],
        ).head(5)
        ranked = ranked.drop(columns="tier_order")
        ranked.insert(1, "marker_rank_within_minor_class", range(1, len(ranked) + 1))
        rows.append(ranked)

    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True).sort_values(
        ["minor_class", "marker_rank_within_minor_class"]
    ).reset_index(drop=True)


def build_presence_table(review: pd.DataFrame) -> pd.DataFrame:
    """Classify each minor class as present, possible, or undetected per sample."""
    labels = sorted(set(review.proposed_label) - {"Unresolved"})
    rows = []
    for sample in detailed.SAMPLES:
        sample_rows = review[review["sample"] == sample]
        for label in labels:
            matching = sample_rows[sample_rows.proposed_label == label]
            strict = matching[matching.confidence.isin(["High", "Medium"])]
            low = matching[matching.confidence == "Low"]
            if not strict.empty:
                status = "Present"
            elif not low.empty:
                status = "Possibly present"
            else:
                status = "Not detected"
            rows.append({
                "sample": sample,
                "treatment": "DMSO" if sample.endswith("DMSO") else "1644",
                "minor_class": label,
                "presence_status": status,
                "strict_cluster_n": len(strict),
                "low_confidence_cluster_n": len(low),
                "assigned_cluster_n": len(matching),
                "strict_presence": not strict.empty,
                "inclusive_presence": not matching.empty,
            })
    return pd.DataFrame(rows)


def validate(
    review: pd.DataFrame,
    shared: pd.DataFrame,
    presence: pd.DataFrame,
    panel_genes: set[str],
) -> list[str]:
    """Validate cluster completeness, marker support fields, and presence rows."""
    if len(review) != 230:
        raise ValueError(f"Expected 230 cluster rows, found {len(review)}.")
    if review.duplicated(["sample", "cluster_id"]).any():
        raise ValueError("Duplicate sample/cluster identifiers found.")
    observed = review.groupby("sample").size().to_dict()
    if observed != detailed.SAMPLES:
        raise ValueError(f"Unexpected sample cluster counts: {observed}")
    if not set(shared.marker).issubset(panel_genes):
        raise ValueError("A selected marker is absent from the panel exports.")
    if not set(shared.support_tier).issubset(SUPPORT_TIERS):
        raise ValueError("An unrecognized marker support tier was generated.")
    fraction_columns = [
        column for column in shared.columns if column.endswith("_fraction")
    ]
    if not shared[fraction_columns].apply(lambda values: values.between(0, 1)).all().all():
        raise ValueError("A marker recurrence fraction falls outside [0, 1].")
    expected_presence_rows = len(detailed.SAMPLES) * (
        review.loc[review.proposed_label != "Unresolved", "proposed_label"].nunique()
    )
    if len(presence) != expected_presence_rows:
        raise ValueError("The sample-by-minor-class presence table is incomplete.")
    return [
        "PASS: exactly 230 sample/cluster rows.",
        "PASS: all eight samples have expected cluster counts.",
        "PASS: sample/cluster identifiers are unique.",
        "PASS: selected markers occur in the panel exports.",
        "PASS: marker tiers and recurrence fractions are valid.",
        "PASS: the sample-by-minor-class presence table is complete.",
    ]


def write_summary(
    path: Path,
    review: pd.DataFrame,
    shared: pd.DataFrame,
    presence: pd.DataFrame,
    checks: list[str],
) -> None:
    """Write a concise summary of comparable minor-class review outputs."""
    confidence = review.confidence.value_counts().reindex(
        ["High", "Medium", "Low", "Unresolved"], fill_value=0
    )
    tier_counts = shared.support_tier.value_counts()
    marker_counts = shared.groupby("minor_class").size()
    labels = sorted(set(review.proposed_label) - {"Unresolved"})
    missing = sorted(set(labels) - set(marker_counts.index))
    lines = [
        "# PTMT PC55 Run 1 r1.5 Comparable Minor-Class Analysis",
        "",
        "These assignments and marker selections are provisional. No source marker file or QC configuration was modified.",
        "",
        "## Overview",
        "",
        f"- Clusters: {len(review)} across {review['sample'].nunique()} samples.",
        f"- Labels represented: {len(labels)} plus Unresolved.",
        f"- Confidence: High {confidence['High']}, Medium {confidence['Medium']}, Low {confidence['Low']}, Unresolved {confidence['Unresolved']}.",
        f"- Selected markers: {len(shared)}.",
        f"- Presence rows: {len(presence)}.",
        "",
        "## Marker support tiers",
        "",
        f"- Pan-cohort: {tier_counts.get('Pan-cohort', 0)}.",
        f"- Core conditional: {tier_counts.get('Core conditional', 0)}.",
        f"- Recurrent conditional: {tier_counts.get('Recurrent conditional', 0)}.",
        f"- Partial support: {tier_counts.get('Partial support', 0)}.",
        f"- Sample-specific: {tier_counts.get('Sample-specific', 0)}.",
        "",
        "Strict presence requires a High/Medium assignment. Inclusive presence also accepts Low assignments. At most five markers are retained per minor class.",
        "",
        "## Sparse marker groups",
        "",
    ]
    sparse = marker_counts[marker_counts < 5]
    lines.extend(
        f"- {label}: {count} selected marker(s)." for label, count in sparse.items()
    )
    lines.extend(f"- {label}: no eligible marker." for label in missing)
    if sparse.empty and not missing:
        lines.append("- None.")
    lines.extend([
        "",
        "## Validation",
        "",
        *[f"- {check}" for check in checks],
        "",
        "## Review priorities",
        "",
        "- Review Low, Unresolved, and sensitivity-unstable clusters first.",
        "- Treat sample-specific markers as visual-review candidates, not shared identity markers.",
        "- Closely inspect competing T/NK, myeloid, epithelial, and stromal minor classes.",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    """Run comparable minor-class annotation, marker selection, and export."""
    args = parse_args()
    repo_root = args.repo_root.resolve()
    reference_path = repo_root / (
        "data/xenium/raw_data/gene_markers/ptmt/"
        "2026-08-14_dotplot_marker_genes_harmonised_scRNAseq_analysis_Kat.csv"
    )
    output_dir = repo_root / (
        "data/xenium/processed/ptmt_pc55/cluster_annotation_review/"
        "run_1_r1p5_minor_class"
    )

    _, gene_evidence = detailed.load_reference(reference_path)
    merged, _, panel_genes = detailed.load_inputs(repo_root, args.scratch_root)
    review, evidence = detailed.assign_clusters(merged, gene_evidence)
    shared = select_shared_markers(review, evidence, gene_evidence)
    presence = build_presence_table(review)
    checks = validate(review, shared, presence, panel_genes)

    output_dir.mkdir(parents=True, exist_ok=True)
    review.to_csv(
        output_dir / "ptmt_pc55_run1_r1p5_minor_cluster_label_review.csv",
        index=False,
    )
    shared.to_csv(
        output_dir / "ptmt_pc55_run1_r1p5_minor_shared_markers.csv",
        index=False,
    )
    presence.to_csv(
        output_dir / "ptmt_pc55_run1_r1p5_minor_class_presence.csv",
        index=False,
    )
    write_summary(
        output_dir / "ptmt_pc55_run1_r1p5_minor_analysis_summary.md",
        review,
        shared,
        presence,
        checks,
    )

    print("\n".join(checks))
    print(f"Wrote {len(review)} cluster rows and {len(shared)} marker rows to {output_dir}")


if __name__ == "__main__":
    main()
