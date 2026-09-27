"""
Title: PTMT PC55 run 1 r1.5 minor-class sensitivity analysis
Date: 2026-09-23
Summary: Compare detailed minor-class assignments across the same five marker
scoring variants used for broad-label sensitivity review. Outputs are
review-only and do not modify source marker files or QC configurations.
"""

from __future__ import annotations

import argparse
import math
from collections import defaultdict
from pathlib import Path

import pandas as pd

import ptmt_pc55_run1_r1p5_broad_sensitivity as shared_sensitivity
import ptmt_pc55_run1_r1p5_marker_selection as detailed


VARIANTS = shared_sensitivity.VARIANTS


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


def score_cluster_variant(
    cluster: pd.DataFrame,
    gene_evidence: dict[str, dict[str, float]],
    variant: str,
) -> tuple[pd.DataFrame, dict[str, list[tuple[str, float]]]]:
    """Calculate detailed candidate-label scores under one scoring variant."""
    hits: dict[str, list[tuple[str, float]]] = defaultdict(list)
    for _, row in cluster.iterrows():
        for label, marker_weight in gene_evidence.get(str(row.gene), {}).items():
            value = shared_sensitivity.variant_contribution(
                row, marker_weight, variant
            )
            if value > 0:
                hits[label].append((str(row.gene), value))

    records = []
    for label, label_hits in hits.items():
        top = sorted(label_hits, key=lambda item: item[1], reverse=True)[:5]
        score = sum(value for _, value in top) / math.sqrt(max(1, len(top)))
        records.append({
            "label": label,
            "broad_class": detailed.broad_for(label),
            "score": score,
            "n_hits": len(label_hits),
        })
    scores = pd.DataFrame(records)
    if not scores.empty:
        scores = scores.sort_values(["score", "n_hits"], ascending=False).reset_index(
            drop=True
        )
    return scores, hits


def assign_variant(
    merged: pd.DataFrame,
    gene_evidence: dict[str, dict[str, float]],
    variant: str,
) -> pd.DataFrame:
    """Assign every cluster under one detailed-label scoring variant."""
    rows = []
    for (sample, cluster_id), cluster in merged.groupby(
        ["sample", "cluster_id"], sort=False
    ):
        scores, hits = score_cluster_variant(cluster, gene_evidence, variant)
        if scores.empty:
            best_label, second_label = "Unresolved", ""
            best_score, second_score, support_count = 0.0, 0.0, 0
        else:
            best = scores.iloc[0]
            second = scores.iloc[1] if len(scores) > 1 else None
            best_label = str(best.label)
            second_label = str(second.label) if second is not None else ""
            best_score = float(best.score)
            second_score = float(second.score) if second is not None else 0.0
            support_count = len({gene for gene, _ in hits[best_label]})

            same_broad = (
                second is not None and best.broad_class == second.broad_class
            )
            close = second_score >= best_score * 0.82
            if close and not same_broad:
                best_label = "Unresolved"
            elif (
                close
                and same_broad
                and best.broad_class in {"T/NK", "Myeloid", "Stromal"}
            ):
                # Match the detailed workflow's conservative parent-lineage fallback.
                best_label = str(best.broad_class)

        margin = (best_score - second_score) / best_score if best_score else 0.0
        if best_label == "Unresolved" or best_score < 0.18 or support_count == 0:
            best_label = "Unresolved"
            confidence = "Unresolved"
        elif support_count >= 4 and best_score >= 0.45 and margin >= 0.18:
            confidence = "High"
        elif support_count >= 2 and best_score >= 0.28 and margin >= 0.08:
            confidence = "Medium"
        else:
            confidence = "Low"

        rows.append({
            "sample": sample,
            "cluster_id": int(cluster_id),
            "variant": variant,
            "label": best_label,
            "confidence": confidence,
            "candidate_score": round(best_score, 5),
            "secondary_candidate": second_label,
            "secondary_score": round(second_score, 5),
            "score_margin_fraction": round(margin, 5),
            "supporting_marker_count": support_count,
        })
    return pd.DataFrame(rows)


def build_sensitivity_tables(
    assignments: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Create per-cluster stability and per-variant agreement summaries."""
    current = assignments[assignments.variant == "current_composite"][
        ["sample", "cluster_id", "label", "confidence"]
    ].rename(columns={"label": "current_label", "confidence": "current_confidence"})
    comparison = assignments.merge(current, on=["sample", "cluster_id"], how="left")
    comparison["matches_current_label"] = (
        comparison["label"] == comparison["current_label"]
    )
    comparison["current_resolved"] = comparison.current_label != "Unresolved"
    comparison["variant_resolved"] = comparison.label != "Unresolved"

    stability = comparison.pivot_table(
        index=["sample", "cluster_id"],
        columns="variant",
        values="label",
        aggfunc="first",
    ).reset_index()
    variant_columns = [column for column in VARIANTS if column in stability.columns]
    stability["label_set_across_variants"] = stability[variant_columns].apply(
        lambda row: "; ".join(sorted(set(row.dropna()))),
        axis=1,
    )
    stability["stable_all_variants"] = (
        stability[variant_columns].nunique(axis=1) == 1
    )
    stability["stable_resolved_all_variants"] = (
        stability.stable_all_variants
        & (stability.current_composite != "Unresolved")
    )

    summary_rows = []
    for variant, group in comparison.groupby("variant", sort=False):
        confidence_counts = group.confidence.value_counts()
        summary_rows.append({
            "variant": variant,
            "description": VARIANTS[variant],
            "clusters": len(group),
            "matches_current_label_n": int(group.matches_current_label.sum()),
            "matches_current_label_fraction": group.matches_current_label.mean(),
            "resolved_clusters": int(group.variant_resolved.sum()),
            "high_confidence": int(confidence_counts.get("High", 0)),
            "medium_confidence": int(confidence_counts.get("Medium", 0)),
            "low_confidence": int(confidence_counts.get("Low", 0)),
            "unresolved": int(confidence_counts.get("Unresolved", 0)),
        })
    return stability, pd.DataFrame(summary_rows)


def validate(
    stability: pd.DataFrame,
    summary: pd.DataFrame,
    current_review: pd.DataFrame,
) -> list[str]:
    """Validate completeness and agreement with the main minor-class calls."""
    if len(stability) != 230:
        raise ValueError(f"Expected 230 cluster rows, found {len(stability)}.")
    if set(summary.variant) != set(VARIANTS):
        raise ValueError("Not all five sensitivity variants were summarized.")
    if stability.duplicated(["sample", "cluster_id"]).any():
        raise ValueError("Duplicate sample/cluster identifiers found.")
    expected = current_review[["sample", "cluster_id", "proposed_label"]]
    observed = stability[["sample", "cluster_id", "current_composite"]]
    comparison = expected.merge(
        observed, on=["sample", "cluster_id"], validate="one_to_one"
    )
    if not comparison.proposed_label.equals(comparison.current_composite):
        raise ValueError(
            "Sensitivity current-composite labels differ from the main analysis."
        )
    return [
        "PASS: exactly 230 sample/cluster rows.",
        "PASS: all five scoring variants are summarized.",
        "PASS: sample/cluster identifiers are unique.",
        "PASS: current-composite labels match the main minor-class analysis.",
    ]


def write_summary(
    path: Path,
    stability: pd.DataFrame,
    summary: pd.DataFrame,
    checks: list[str],
) -> None:
    """Write a concise minor-class sensitivity summary."""
    stable_all = int(stability.stable_all_variants.sum())
    stable_resolved = int(stability.stable_resolved_all_variants.sum())
    current_resolved = int((stability.current_composite != "Unresolved").sum())
    lines = [
        "# PTMT PC55 Run 1 r1.5 Minor-Class Sensitivity Analysis",
        "",
        "This review-only analysis applies the same five scoring variants used in the broad workflow to the detailed minor-class vocabulary.",
        "",
        "## Stability overview",
        "",
        f"- Stable labels across all variants: {stable_all} of {len(stability)} clusters.",
        f"- Stable resolved labels: {stable_resolved} of {current_resolved} clusters resolved by the current composite score.",
        f"- At least one variant disagreement: {len(stability) - stable_all} clusters.",
        "",
        "## Variant agreement",
        "",
        "| Variant | Matches current | Resolved | High | Medium | Low | Unresolved |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in summary.itertuples(index=False):
        lines.append(
            f"| {row.variant} | {row.matches_current_label_n}/{row.clusters} "
            f"({row.matches_current_label_fraction:.1%}) | {row.resolved_clusters} | "
            f"{row.high_confidence} | {row.medium_confidence} | "
            f"{row.low_confidence} | {row.unresolved} |"
        )
    lines.extend([
        "",
        "Confidence counts are descriptive because the same thresholds are applied to differently scaled variant scores. Label agreement is the primary sensitivity result.",
        "",
        "## Validation",
        "",
        *[f"- {check}" for check in checks],
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    """Run all minor-class scoring variants and export stability results."""
    args = parse_args()
    repo_root = args.repo_root.resolve()
    reference_path = repo_root / (
        "data/xenium/raw_data/gene_markers/ptmt/"
        "2026-08-14_dotplot_marker_genes_harmonised_scRNAseq_analysis_Kat.csv"
    )
    output_dir = repo_root / (
        "data/xenium/processed/ptmt_pc55/cluster_annotation_review/"
        "run_1_r1p5_minor_class_sensitivity"
    )

    _, gene_evidence = detailed.load_reference(reference_path)
    merged, _, _ = detailed.load_inputs(repo_root, args.scratch_root)
    assignments = pd.concat(
        [
            assign_variant(merged, gene_evidence, variant)
            for variant in VARIANTS
        ],
        ignore_index=True,
    )
    stability, summary = build_sensitivity_tables(assignments)
    current_review, _ = detailed.assign_clusters(merged, gene_evidence)
    checks = validate(stability, summary, current_review)

    output_dir.mkdir(parents=True, exist_ok=True)
    assignments.to_csv(
        output_dir / "ptmt_pc55_run1_r1p5_minor_sensitivity_assignments.csv",
        index=False,
    )
    stability.to_csv(
        output_dir / "ptmt_pc55_run1_r1p5_minor_sensitivity_cluster_stability.csv",
        index=False,
    )
    summary.to_csv(
        output_dir / "ptmt_pc55_run1_r1p5_minor_sensitivity_summary.csv",
        index=False,
    )
    write_summary(
        output_dir / "ptmt_pc55_run1_r1p5_minor_sensitivity_summary.md",
        stability,
        summary,
        checks,
    )

    print("\n".join(checks))
    print(f"Wrote sensitivity outputs for {len(stability)} clusters to {output_dir}")


if __name__ == "__main__":
    main()
