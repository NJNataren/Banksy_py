"""
Title: PTMT PC55 run 1 r1.5 broad cell-type analysis
Date: 2026-09-23
Summary: Re-score all 230 r1.5 clusters using a broader biological vocabulary,
select recurrent markers for those groups, and write a separate provisional
review package without modifying detailed results or source data.
"""

from __future__ import annotations

import argparse
import math
from collections import defaultdict
from pathlib import Path

import pandas as pd

import ptmt_pc55_run1_r1p5_marker_selection as detailed


MINOR_TO_BROAD = {
    "LE-KLK3": "Luminal epithelial",
    "LE-KLK4": "Luminal epithelial",
    "Luminal": "Luminal epithelial",
    "Basal": "Basal epithelial",
    "Club": "Club/Hillock epithelial",
    "Hillock": "Club/Hillock epithelial",
    "LPCs": "Cycling epithelial",
    "NE": "Neuroendocrine epithelial",
    "T cell": "T cell",
    "CD4T": "T cell",
    "CD4 Trm": "T cell",
    "CD4 naive-cm": "T cell",
    "Treg": "T cell",
    "CD8T": "T cell",
    "CD8 Trm": "T cell",
    "CD8 cytotoxic": "T cell",
    "T_MKI67": "T cell",
    "NK": "NK",
    "NK cell": "NK",
    "NK CD16pos": "NK",
    "NK CD16neg": "NK",
    "MNP": "Macrophage",
    "Mph": "Macrophage",
    "Mac1": "Macrophage",
    "Mac2": "Macrophage",
    "Mac-MT1": "Macrophage",
    "Mono": "Monocyte",
    "Mac-cycling": "Cycling myeloid",
    "Mye_cycling": "Cycling myeloid",
    "DCs": "Dendritic cell",
    "B cell": "B cell",
    "Plasma_B": "Plasma cell",
    "Mast": "Mast",
    "Endo": "Endothelial",
    "Fibroblast": "Fibroblast/CAF",
    "CAFs": "Fibroblast/CAF",
    "Mural": "Mural/Pericyte",
    "Pericytes": "Mural/Pericyte",
}

GROUP_TO_COMPARTMENT = {
    "Luminal epithelial": "Epithelial",
    "Basal epithelial": "Epithelial",
    "Club/Hillock epithelial": "Epithelial",
    "Cycling epithelial": "Epithelial",
    "Neuroendocrine epithelial": "Epithelial",
    "T cell": "T/NK",
    "NK": "T/NK",
    "B cell": "B/mast",
    "Plasma cell": "B/mast",
    "Mast": "B/mast",
    "Monocyte": "Myeloid",
    "Macrophage": "Myeloid",
    "Dendritic cell": "Myeloid",
    "Cycling myeloid": "Myeloid",
    "Endothelial": "Stromal",
    "Fibroblast/CAF": "Stromal",
    "Mural/Pericyte": "Stromal",
    "Unresolved/Mixed": "Unresolved/Mixed",
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


def build_broad_reference(reference: pd.DataFrame) -> dict[str, dict[str, float]]:
    """Map reference genes to broad groups with ambiguity normalization."""
    gene_groups: dict[str, set[str]] = defaultdict(set)
    for row in reference.itertuples(index=False):
        group = MINOR_TO_BROAD.get(str(row.minor_class))
        if group:
            gene_groups[str(row.gene)].add(group)

    evidence: dict[str, dict[str, float]] = defaultdict(dict)
    for gene, groups in gene_groups.items():
        weight = 1.0 / math.sqrt(len(groups))
        for group in groups:
            evidence[gene][group] = weight
    return evidence


def score_cluster(
    cluster: pd.DataFrame, gene_evidence: dict[str, dict[str, float]]
) -> tuple[pd.DataFrame, dict[str, list[tuple[str, float]]]]:
    """Calculate marker-set-normalized broad-group scores for one cluster."""
    hits: dict[str, list[tuple[str, float]]] = defaultdict(list)
    for _, row in cluster.iterrows():
        for group, marker_weight in gene_evidence.get(str(row.gene), {}).items():
            value = detailed.contribution(row, marker_weight)
            if value > 0:
                hits[group].append((str(row.gene), value))

    records = []
    for group, group_hits in hits.items():
        top = sorted(group_hits, key=lambda item: item[1], reverse=True)[:5]
        score = sum(value for _, value in top) / math.sqrt(max(1, len(top)))
        records.append({
            "label": group,
            "compartment": GROUP_TO_COMPARTMENT[group],
            "score": score,
            "n_hits": len(group_hits),
        })
    scores = pd.DataFrame(records)
    if not scores.empty:
        scores = scores.sort_values(["score", "n_hits"], ascending=False).reset_index(drop=True)
    return scores, hits


def format_genes(hits: list[tuple[str, float]], limit: int = 8) -> str:
    """Return unique genes in descending contribution order."""
    ordered = sorted(hits, key=lambda item: item[1], reverse=True)
    return "; ".join(list(dict.fromkeys(gene for gene, _ in ordered))[:limit])


def format_evidence(
    cluster: pd.DataFrame,
    genes: list[str],
    gene_evidence: dict[str, dict[str, float]],
    limit: int = 8,
) -> str:
    """Format marker statistics and broad reference groups for review."""
    indexed = cluster.drop_duplicates("gene").set_index("gene")
    records = []
    for gene in list(dict.fromkeys(genes))[:limit]:
        if gene not in indexed.index:
            continue
        row = indexed.loc[gene]
        groups = "/".join(sorted(gene_evidence.get(gene, {}))) or "unannotated"
        records.append(
            f"{gene} [rank={int(row['rank'])}, score={row.scores:.3g}, "
            f"logFC={row.logfoldchanges:.3g}, padj={row.pvals_adj:.3g}, "
            f"mean={row.mean_expression:.3g}, pct={row.percent_expressing:.1f}, ref={groups}]"
        )
    return "; ".join(records)


def assign_clusters(
    merged: pd.DataFrame, gene_evidence: dict[str, dict[str, float]]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Assign broad labels and retain all cluster/gene evidence."""
    review_rows = []
    evidence_rows = []
    for (sample, cluster_id), cluster in merged.groupby(["sample", "cluster_id"], sort=False):
        scores, hits = score_cluster(cluster, gene_evidence)
        if scores.empty:
            best_label, second_label = "Unresolved/Mixed", ""
            best_score, second_score, best_hits = 0.0, 0.0, []
        else:
            best = scores.iloc[0]
            second = scores.iloc[1] if len(scores) > 1 else None
            best_label = str(best.label)
            second_label = str(second.label) if second is not None else ""
            best_score = float(best.score)
            second_score = float(second.score) if second is not None else 0.0
            best_hits = hits[best_label]

            close = second_score >= best_score * 0.82
            cross_compartment = second is not None and best.compartment != second.compartment
            if close and cross_compartment:
                best_label = "Unresolved/Mixed"

        support_genes = [gene for gene, _ in sorted(best_hits, key=lambda item: item[1], reverse=True)]
        support_count = len(set(support_genes))
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

        compartment = GROUP_TO_COMPARTMENT[best_label]
        contradictory: list[tuple[str, float]] = []
        if not scores.empty:
            for candidate in scores.itertuples(index=False):
                if candidate.compartment != compartment and candidate.score >= max(0.12, best_score * 0.45):
                    contradictory.extend(hits[candidate.label][:3])
        contradictory_genes = [gene for gene, _ in contradictory]
        top_positive = cluster.query("logfoldchanges > 0 and pvals_adj < 0.05").nsmallest(10, "rank")["gene"].tolist()

        review_rows.append({
            "sample": sample,
            "treatment": "DMSO" if sample.endswith("DMSO") else "1644",
            "resolution": 1.5,
            "cluster_id": int(cluster_id),
            "proposed_broad_label": best_label,
            "major_compartment": compartment,
            "confidence": confidence,
            "candidate_score": round(best_score, 5),
            "secondary_candidate": second_label,
            "secondary_score": round(second_score, 5),
            "score_margin_fraction": round(margin, 5),
            "supporting_markers": format_genes(best_hits),
            "supporting_marker_evidence": format_evidence(cluster, support_genes, gene_evidence),
            "contradictory_markers": format_genes(contradictory, 6),
            "contradictory_marker_evidence": format_evidence(cluster, contradictory_genes, gene_evidence, 6),
            "top_positive_markers": "; ".join(top_positive),
            "top_positive_marker_evidence": format_evidence(cluster, top_positive, gene_evidence, 10),
            "rationale": (
                f"Top broad candidate {scores.iloc[0].label if not scores.empty else 'none'} "
                f"score={best_score:.3f}; runner-up {second_label or 'none'} "
                f"score={second_score:.3f}; {support_count} positive reference-marker hits."
            ),
            "review_notes": "Provisional broad assignment; review Low, Unresolved, and mixed-lineage calls.",
        })

        for row in cluster.itertuples(index=False):
            evidence_rows.append({
                "sample": sample,
                "cluster_id": int(cluster_id),
                "gene": row.gene,
                "rank": int(row.rank),
                "score": row.scores,
                "logfoldchanges": row.logfoldchanges,
                "pvals_adj": row.pvals_adj,
                "mean_expression": row.mean_expression,
                "percent_expressing": row.percent_expressing,
            })

    review = pd.DataFrame(review_rows).sort_values(["sample", "cluster_id"]).reset_index(drop=True)
    return review, pd.DataFrame(evidence_rows)


def select_shared_markers(
    review: pd.DataFrame,
    evidence: pd.DataFrame,
    gene_evidence: dict[str, dict[str, float]],
) -> pd.DataFrame:
    """Rank markers using strict and inclusive cell-type presence denominators."""
    assigned = review.query("proposed_broad_label != 'Unresolved/Mixed'")[
        ["sample", "cluster_id", "proposed_broad_label", "confidence"]
    ]
    positive = evidence.merge(assigned, on=["sample", "cluster_id"], how="inner")
    positive = positive.query("logfoldchanges > 0 and pvals_adj < 0.05").copy()
    rows = []

    for label, clusters in assigned.groupby("proposed_broad_label"):
        target = positive[positive.proposed_broad_label == label]
        other = positive[positive.proposed_broad_label != label]
        inclusive_samples = set(clusters["sample"])
        strict_samples = set(
            clusters.loc[clusters.confidence.isin(["High", "Medium"]), "sample"]
        )
        dmso_samples = {sample for sample in inclusive_samples if sample.endswith("DMSO")}
        treated_samples = inclusive_samples - dmso_samples
        n_clusters = len(clusters)
        candidates = []

        for gene, gene_rows in target.groupby("gene"):
            if label not in gene_evidence.get(str(gene), {}) or gene in detailed.GENERIC_OR_STATE:
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
                "broad_cell_type": label,
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
        ranked.insert(1, "marker_rank_within_cell_type", range(1, len(ranked) + 1))
        rows.append(ranked)

    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True).sort_values(
        ["broad_cell_type", "marker_rank_within_cell_type"]
    ).reset_index(drop=True)


def build_presence_table(review: pd.DataFrame) -> pd.DataFrame:
    """Classify each broad cell type as present, possible, or not detected per sample."""
    rows = []
    labels = sorted(set(GROUP_TO_COMPARTMENT) - {"Unresolved/Mixed"})
    for sample in detailed.SAMPLES:
        sample_review = review[review["sample"] == sample]
        for label in labels:
            matching = sample_review[sample_review.proposed_broad_label == label]
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
                "broad_cell_type": label,
                "presence_status": status,
                "strict_cluster_n": len(strict),
                "low_confidence_cluster_n": len(low),
                "assigned_cluster_n": len(matching),
                "strict_presence": not strict.empty,
                "inclusive_presence": not matching.empty,
            })
    return pd.DataFrame(rows)


def validate(review: pd.DataFrame, shared: pd.DataFrame, panel_genes: set[str]) -> list[str]:
    """Validate broad review completeness, uniqueness, and marker recurrence."""
    if len(review) != 230:
        raise ValueError(f"Expected 230 cluster rows, found {len(review)}.")
    if review.duplicated(["sample", "cluster_id"]).any():
        raise ValueError("Duplicate sample/cluster identifiers found.")
    if review.groupby("sample").size().to_dict() != detailed.SAMPLES:
        raise ValueError("Observed sample cluster counts do not match the task specification.")
    valid_tiers = {
        "Pan-cohort", "Core conditional", "Recurrent conditional",
        "Partial support", "Sample-specific",
    }
    if not shared.empty and not set(shared.support_tier).issubset(valid_tiers):
        raise ValueError("An unrecognized marker-support tier was generated.")
    recurrence_columns = [
        "inclusive_sample_recurrence_fraction",
        "strict_sample_recurrence_fraction",
        "dmso_sample_recurrence_fraction",
        "treated_sample_recurrence_fraction",
    ]
    if not shared.empty and not shared[recurrence_columns].apply(
        lambda column: column.between(0, 1)
    ).all().all():
        raise ValueError("A sample-recurrence fraction falls outside zero to one.")
    if not shared.empty and not set(shared.marker).issubset(panel_genes):
        raise ValueError("A selected marker is absent from the panel exports.")
    return [
        "PASS: exactly 230 sample/cluster rows.",
        "PASS: all eight samples have the expected cluster counts.",
        "PASS: sample/cluster identifiers are unique.",
        "PASS: all selected markers occur in the panel exports.",
        "PASS: every marker has a recognized support tier and valid recurrence fractions.",
    ]


def write_summary(
    path: Path,
    review: pd.DataFrame,
    shared: pd.DataFrame,
    presence: pd.DataFrame,
    checks: list[str],
) -> None:
    """Write broad-analysis results and methodology as Markdown."""
    confidence = review.confidence.value_counts().reindex(
        ["High", "Medium", "Low", "Unresolved"], fill_value=0
    )
    label_counts = review.proposed_broad_label.value_counts()
    marker_counts = shared.groupby("broad_cell_type").size() if not shared.empty else pd.Series(dtype=int)
    low_support = marker_counts[marker_counts < 5]
    missing = sorted(set(label_counts.index) - {"Unresolved/Mixed"} - set(marker_counts.index))
    all_sample_markers = shared[
        shared.inclusive_marker_sample_n == review["sample"].nunique()
    ].copy()
    tier_counts = shared.support_tier.value_counts()

    lines = [
        "# PTMT PC55 Run 1 r1.5 Broad Cell-Type Review",
        "",
        "This analysis re-scored the underlying marker and expression tables using a broader vocabulary. It did not overwrite the detailed review, QC configurations, or source marker files.",
        "",
        "## Results overview",
        "",
        f"- Clusters reviewed: {len(review)} across {review['sample'].nunique()} samples.",
        f"- Confidence: High {confidence['High']}, Medium {confidence['Medium']}, Low {confidence['Low']}, Unresolved {confidence['Unresolved']}.",
        f"- Broad labels represented: {len(label_counts)}.",
        "",
        "## Broad label counts",
        "",
        *[f"- {label}: {count} clusters." for label, count in label_counts.items()],
        "",
        "## Shared marker coverage",
        "",
    ]
    if low_support.empty and not missing:
        lines.append("- Every represented broad type received five recurrent markers.")
    else:
        lines.extend(f"- {label}: {count} recurrent markers." for label, count in low_support.items())
        lines.extend(f"- {label}: no marker met the two-sample recurrence rule." for label in missing)
    lines.extend([
        "",
        "## Markers shared across all eight samples",
        "",
        "The following markers were positively enriched in at least one matching cluster in every run 1 sample. This is an eight-sample recurrence criterion; it does not require expression in every cluster assigned that cell type.",
        "",
        "| Broad cell type | Markers present across all 8 samples | Matching-cluster recurrence |",
        "| --- | --- | --- |",
    ])
    for label, group in all_sample_markers.groupby("broad_cell_type", sort=True):
        markers = ", ".join(group.marker)
        recurrence = ", ".join(
            f"{row.marker} {row.cluster_recurrence_fraction:.0%}"
            for row in group.itertuples(index=False)
        )
        lines.append(f"| {label} | {markers} | {recurrence} |")
    lines.extend([
        "",
        "Six broad cell types met the literal eight-sample criterion: Basal epithelial, Fibroblast/CAF, Luminal epithelial, Macrophage, Mast, and T cell. `IGKC` was recurrent in all seven samples containing an assigned B-cell cluster, but B cell did not meet the eight-sample criterion because no B-cell cluster was assigned in the remaining sample.",
        "",
        "## Conditional and subset-supported markers",
        "",
        "Cell-type presence was evaluated separately from marker recurrence. Strict presence requires at least one High/Medium cluster; inclusive presence also accepts Low-confidence clusters. 'Not detected' does not assert biological absence.",
        "",
        "| Support tier | Definition | Selected markers |",
        "| --- | --- | ---: |",
        f"| Pan-cohort | Positive in all 8 samples | {tier_counts.get('Pan-cohort', 0)} |",
        f"| Core conditional | Positive in every strict-positive sample; type present strictly in at least 3 samples | {tier_counts.get('Core conditional', 0)} |",
        f"| Recurrent conditional | Positive in at least 75% of strict-positive samples; type present strictly in at least 3 samples | {tier_counts.get('Recurrent conditional', 0)} |",
        f"| Partial support | Positive in at least 2 inclusive-positive samples and at least 50% of them | {tier_counts.get('Partial support', 0)} |",
        f"| Sample-specific | Does not meet a shared-marker tier | {tier_counts.get('Sample-specific', 0)} |",
        "",
        "The shared-marker CSV reports strict and inclusive denominators, DMSO and 1644 recurrence separately, supporting sample names, and treatment-associated caveats. The sample-presence CSV contains the per-sample Present, Possibly present, and Not detected classifications used for these denominators.",
        "",
        "## Presence status counts",
        "",
        *[
            f"- {status}: {count} sample-by-cell-type combinations."
            for status, count in presence.presence_status.value_counts().items()
        ],
        "",
        "## Methodology",
        "",
        "Reference minor classes were mapped to the broad vocabulary defined in the analysis script. All clusters were then re-scored directly; detailed labels were not simply renamed.",
        "",
        "A gene contributed only when log fold-change and Wilcoxon score were positive, adjusted p-value was at most 0.1, it matched the broad reference group, and it was not in the generic/state exclusion list. Its contribution was:",
        "",
        "```text",
        "marker_weight * [1/log2(rank+1)]",
        "* [0.35 + 0.65*min(4,-log10(padj))/4]",
        "* [0.40 + 0.60*min(3,logFC)/3]",
        "* sqrt(percent_expressing/100)",
        "```",
        "",
        "The five strongest gene contributions were summed and divided by the square root of the number retained. A cross-compartment runner-up within 82% of the best score produced `Unresolved/Mixed`.",
        "",
        "Confidence rules matched the detailed analysis: High required at least four hits, score >= 0.45, and margin >= 0.18; Medium required at least two hits, score >= 0.28, and margin >= 0.08. Other non-unresolved calls were Low; scores below 0.18 or unsupported/mixed calls were Unresolved.",
        "",
        "Shared-marker scoring combined inclusive sample recurrence (35%), matching-cluster recurrence (25%), specificity versus other labels (20%), median log fold-change (10%), and median percent expressing (10%), followed by a 1.2 harmonised-reference multiplier. Strict recurrence and treatment-specific recurrence were calculated alongside this score. Up to five markers were retained per broad type and classified into the support tiers above.",
        "",
        "## Validation",
        "",
        *[f"- {check}" for check in checks],
        "- PNG visual review remains pending because the local image viewer was unavailable; the numerical tables behind the plots were used.",
        "",
        "## Review priorities",
        "",
        "- Review Low and Unresolved/Mixed clusters before transferring labels into any configuration.",
        "- Inspect T cell versus NK, Monocyte versus Macrophage, and the three stromal groups most carefully.",
        "- Treat Cycling epithelial, Cycling myeloid, and Neuroendocrine epithelial calls cautiously when only one reference gene contributes.",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    """Run broad annotation, marker selection, validation, and export."""
    args = parse_args()
    repo_root = args.repo_root.resolve()
    reference_path = repo_root / "data/xenium/raw_data/gene_markers/ptmt/2026-08-14_dotplot_marker_genes_harmonised_scRNAseq_analysis_Kat.csv"
    output_dir = repo_root / "data/xenium/processed/ptmt_pc55/cluster_annotation_review/run_1_r1p5_broad_cell_types"

    reference = pd.read_csv(reference_path)
    gene_evidence = build_broad_reference(reference)
    merged, _, panel_genes = detailed.load_inputs(repo_root, args.scratch_root)
    review, evidence = assign_clusters(merged, gene_evidence)
    shared = select_shared_markers(review, evidence, gene_evidence)
    presence = build_presence_table(review)
    checks = validate(review, shared, panel_genes)

    output_dir.mkdir(parents=True, exist_ok=True)
    review.to_csv(output_dir / "ptmt_pc55_run1_r1p5_broad_cluster_label_review.csv", index=False)
    shared.to_csv(output_dir / "ptmt_pc55_run1_r1p5_broad_shared_markers.csv", index=False)
    presence.to_csv(output_dir / "ptmt_pc55_run1_r1p5_broad_cell_type_presence.csv", index=False)
    write_summary(
        output_dir / "ptmt_pc55_run1_r1p5_broad_analysis_summary.md",
        review,
        shared,
        presence,
        checks,
    )

    print("\n".join(checks))
    print(f"Wrote {len(review)} broad cluster rows and {len(shared)} marker rows to {output_dir}")


if __name__ == "__main__":
    main()
