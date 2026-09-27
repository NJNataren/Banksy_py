"""
Title: PTMT PC55 run 1 r1.5 broad annotation sensitivity analysis
Date: 2026-09-23
Summary: Compare the current broad annotation heuristic with simpler scoring
variants using the original Wilcoxon marker tables and expression summaries.
Outputs are review-only CSV and Markdown files; source marker files and QC
configuration files are not modified.
"""

from __future__ import annotations

import argparse
import math
from collections import defaultdict
from pathlib import Path

import pandas as pd

import ptmt_pc55_run1_r1p5_broad_cell_types as broad
import ptmt_pc55_run1_r1p5_marker_selection as detailed


VARIANTS = {
    "current_composite": "Original rank, significance, logFC, and coverage score.",
    "rank_only": "Eligible reference markers weighted only by rank.",
    "unweighted_hits": "Eligible reference markers counted equally.",
    "rank_significance": "Eligible reference markers weighted by rank and adjusted-p significance.",
    "no_coverage": "Original score without percent-expressing coverage weight.",
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


def eligible_marker(row: pd.Series) -> bool:
    """Return whether a gene can contribute evidence under all variants."""
    return (
        row.gene not in detailed.GENERIC_OR_STATE
        and row.logfoldchanges > 0
        and row.scores > 0
        and row.pvals_adj <= 0.1
    )


def variant_contribution(row: pd.Series, marker_weight: float, variant: str) -> float:
    """Calculate one gene's contribution under a named sensitivity variant."""
    if not eligible_marker(row):
        return 0.0
    rank_weight = 1.0 / math.log2(float(row["rank"]) + 1.0)
    significance = min(4.0, max(0.0, -math.log10(max(float(row.pvals_adj), 1e-300)))) / 4.0
    significance_weight = 0.35 + 0.65 * significance
    effect_weight = 0.40 + 0.60 * min(3.0, float(row.logfoldchanges)) / 3.0
    coverage_weight = math.sqrt(min(100.0, max(0.0, float(row.percent_expressing))) / 100.0)

    if variant == "current_composite":
        return detailed.contribution(row, marker_weight)
    if variant == "rank_only":
        return marker_weight * rank_weight
    if variant == "unweighted_hits":
        return marker_weight
    if variant == "rank_significance":
        return marker_weight * rank_weight * significance_weight
    if variant == "no_coverage":
        return marker_weight * rank_weight * significance_weight * effect_weight
    raise ValueError(f"Unknown variant: {variant}")


def score_cluster_variant(
    cluster: pd.DataFrame,
    gene_evidence: dict[str, dict[str, float]],
    variant: str,
) -> tuple[pd.DataFrame, dict[str, list[tuple[str, float]]]]:
    """Calculate broad candidate scores for one cluster under one variant."""
    hits: dict[str, list[tuple[str, float]]] = defaultdict(list)
    for _, row in cluster.iterrows():
        for label, marker_weight in gene_evidence.get(str(row.gene), {}).items():
            value = variant_contribution(row, marker_weight, variant)
            if value > 0:
                hits[label].append((str(row.gene), value))

    records = []
    for label, label_hits in hits.items():
        top = sorted(label_hits, key=lambda item: item[1], reverse=True)[:5]
        score = sum(value for _, value in top) / math.sqrt(max(1, len(top)))
        records.append({
            "label": label,
            "compartment": broad.GROUP_TO_COMPARTMENT[label],
            "score": score,
            "n_hits": len(label_hits),
        })
    scores = pd.DataFrame(records)
    if not scores.empty:
        scores = scores.sort_values(["score", "n_hits"], ascending=False).reset_index(drop=True)
    return scores, hits


def assign_variant(
    merged: pd.DataFrame,
    gene_evidence: dict[str, dict[str, float]],
    variant: str,
) -> pd.DataFrame:
    """Assign broad labels for every cluster under a scoring variant."""
    rows = []
    for (sample, cluster_id), cluster in merged.groupby(["sample", "cluster_id"], sort=False):
        scores, hits = score_cluster_variant(cluster, gene_evidence, variant)
        if scores.empty:
            best_label, second_label = "Unresolved/Mixed", ""
            best_score, second_score, support_count = 0.0, 0.0, 0
        else:
            best = scores.iloc[0]
            second = scores.iloc[1] if len(scores) > 1 else None
            best_label = str(best.label)
            second_label = str(second.label) if second is not None else ""
            best_score = float(best.score)
            second_score = float(second.score) if second is not None else 0.0
            support_count = len({gene for gene, _ in hits[best_label]})

            close = second_score >= best_score * 0.82
            cross_compartment = second is not None and best.compartment != second.compartment
            if close and cross_compartment:
                best_label = "Unresolved/Mixed"

        margin = (best_score - second_score) / best_score if best_score else 0.0
        if best_label == "Unresolved/Mixed" or best_score < 0.18 or support_count == 0:
            best_label = "Unresolved/Mixed"
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


def build_sensitivity_tables(assignments: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Create per-cluster stability and per-variant agreement summaries."""
    current = assignments[assignments.variant == "current_composite"][
        ["sample", "cluster_id", "label", "confidence"]
    ].rename(columns={"label": "current_label", "confidence": "current_confidence"})
    comparison = assignments.merge(current, on=["sample", "cluster_id"], how="left")
    comparison["matches_current_label"] = comparison["label"] == comparison["current_label"]
    comparison["current_resolved"] = comparison["current_label"] != "Unresolved/Mixed"
    comparison["variant_resolved"] = comparison["label"] != "Unresolved/Mixed"

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
    stability["stable_all_variants"] = stability[variant_columns].nunique(axis=1) == 1
    stability["stable_resolved_all_variants"] = (
        stability["stable_all_variants"]
        & (stability["current_composite"] != "Unresolved/Mixed")
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
    summary = pd.DataFrame(summary_rows)
    return stability, summary


def validate(stability: pd.DataFrame, summary: pd.DataFrame) -> list[str]:
    """Validate sensitivity outputs against expected cluster and variant counts."""
    if len(stability) != 230:
        raise ValueError(f"Expected 230 cluster rows, found {len(stability)}.")
    if set(summary.variant) != set(VARIANTS):
        raise ValueError("Not all sensitivity variants were summarized.")
    if stability.duplicated(["sample", "cluster_id"]).any():
        raise ValueError("Duplicate sample/cluster identifiers found.")
    return [
        "PASS: exactly 230 sample/cluster rows in the sensitivity table.",
        "PASS: all five scoring variants were summarized.",
        "PASS: sample/cluster identifiers are unique.",
    ]


def write_summary(path: Path, stability: pd.DataFrame, summary: pd.DataFrame, checks: list[str]) -> None:
    """Write a concise Markdown summary of label stability across variants."""
    stable_all = int(stability.stable_all_variants.sum())
    stable_resolved = int(stability.stable_resolved_all_variants.sum())
    current_resolved = int((stability.current_composite != "Unresolved/Mixed").sum())
    unstable = stability[~stability.stable_all_variants]

    lines = [
        "# PTMT PC55 Run 1 r1.5 Broad Sensitivity Analysis",
        "",
        "This review-only analysis compares the current broad composite score with simpler marker-set scoring variants. It does not modify source marker files, QC configs, or the detailed annotation outputs.",
        "",
        "## Stability overview",
        "",
        f"- Stable labels across all variants: {stable_all} of {len(stability)} clusters.",
        f"- Stable resolved labels across all variants: {stable_resolved} of {current_resolved} clusters resolved by the current composite score.",
        f"- Clusters with at least one variant disagreement: {len(unstable)}.",
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
            f"{row.high_confidence} | {row.medium_confidence} | {row.low_confidence} | {row.unresolved} |"
        )
    lines.extend([
        "",
        "## Review priorities",
        "",
        "- Prioritize clusters where rank-only or unweighted scoring changes the broad label.",
        "- Treat labels that are stable across all variants as stronger candidates for manual review, not as final truth.",
        "- The score constants remain exploratory heuristics; this analysis checks robustness to simpler alternatives rather than calibrating probabilities.",
        "",
        "## Validation",
        "",
        *[f"- {check}" for check in checks],
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    """Run broad-label sensitivity analysis and write review outputs."""
    args = parse_args()
    repo_root = args.repo_root.resolve()
    reference_path = repo_root / "data/xenium/raw_data/gene_markers/ptmt/2026-08-14_dotplot_marker_genes_harmonised_scRNAseq_analysis_Kat.csv"
    output_dir = repo_root / "data/xenium/processed/ptmt_pc55/cluster_annotation_review/run_1_r1p5_broad_sensitivity"

    reference = pd.read_csv(reference_path)
    gene_evidence = broad.build_broad_reference(reference)
    merged, _, _ = detailed.load_inputs(repo_root, args.scratch_root)

    assignments = pd.concat(
        [assign_variant(merged, gene_evidence, variant) for variant in VARIANTS],
        ignore_index=True,
    )
    stability, summary = build_sensitivity_tables(assignments)
    checks = validate(stability, summary)

    output_dir.mkdir(parents=True, exist_ok=True)
    assignments.to_csv(output_dir / "ptmt_pc55_run1_r1p5_broad_sensitivity_assignments.csv", index=False)
    stability.to_csv(output_dir / "ptmt_pc55_run1_r1p5_broad_sensitivity_cluster_stability.csv", index=False)
    summary.to_csv(output_dir / "ptmt_pc55_run1_r1p5_broad_sensitivity_summary.csv", index=False)
    write_summary(
        output_dir / "ptmt_pc55_run1_r1p5_broad_sensitivity_summary.md",
        stability,
        summary,
        checks,
    )

    print("\n".join(checks))
    print(f"Wrote sensitivity outputs for {len(stability)} clusters to {output_dir}")


if __name__ == "__main__":
    main()
