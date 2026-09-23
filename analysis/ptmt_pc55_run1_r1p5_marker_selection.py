"""
Title: PTMT PC55 run 1 r1.5 marker selection
Date: 2026-09-22
Summary: Build provisional cluster annotations and recurrent marker selections for
the eight PTMT PC55 run 1 samples from Wilcoxon ranks, panel expression summaries,
and the harmonised scRNA-seq marker reference. Source files are read-only.
"""

from __future__ import annotations

import argparse
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


SAMPLES = {
    "10329_run_1_DMSO": 19,
    "10331_run_1_1644": 18,
    "10707_run_1_DMSO": 29,
    "10708_run_1_1644": 52,
    "10722_run_1_DMSO": 27,
    "10723_run_1_1644": 31,
    "10733_run_1_DMSO": 26,
    "10734_run_1_1644": 28,
}

LABEL_TO_BROAD = {
    "LE-KLK3": "Epithelial", "LE-KLK4": "Epithelial", "Luminal": "Epithelial",
    "Basal": "Epithelial", "Club": "Epithelial", "Hillock": "Epithelial",
    "LPCs": "Epithelial", "NE": "Epithelial",
    "T cell": "T/NK", "CD4T": "T/NK", "CD4 Trm": "T/NK",
    "CD4 naive-cm": "T/NK", "Treg": "T/NK", "CD8T": "T/NK",
    "CD8 Trm": "T/NK", "CD8 cytotoxic": "T/NK", "T_MKI67": "T/NK",
    "NK": "T/NK", "NK cell": "T/NK", "NK CD16pos": "T/NK",
    "NK CD16neg": "T/NK",
    "MNP": "Myeloid", "Mono": "Myeloid", "Mph": "Myeloid",
    "Mac1": "Myeloid", "Mac2": "Myeloid", "Mac-MT1": "Myeloid",
    "Mac-cycling": "Myeloid", "Mye_cycling": "Myeloid", "DCs": "Myeloid",
    "B cell": "B/mast", "Plasma_B": "B/mast", "Mast": "B/mast",
    "Endo": "Stromal", "Fibroblast": "Stromal", "CAFs": "Stromal",
    "Mural": "Stromal", "Pericytes": "Stromal",
}

# Closely overlapping labels are collapsed when the panel does not separate them.
LABEL_ALIASES = {
    "NK cell": "NK", "CD4T": "T cell", "CD8T": "T cell", "MNP": "Myeloid",
    "Mph": "Macrophage", "Mac1": "Macrophage", "Mac2": "Macrophage",
    "Mural": "Pericytes", "Fibroblast": "Fibroblast", "CAFs": "CAFs",
}

GENERIC_OR_STATE = {
    "ACTB", "B2M", "FOS", "FOSB", "JUN", "JUNB", "DUSP1", "IER2", "IER3",
    "HLA-A", "HLA-B", "HLA-C", "IFITM1", "IFITM2", "IFITM3", "ISG15",
    "MKI67", "TOP2A", "STMN1", "TUBA1B", "TUBB", "EPCAM", "KRT8", "KRT18",
    "PTPRC", "VIM", "MALAT1",
}


def parse_args() -> argparse.Namespace:
    """Parse command-line paths for the reproducible analysis."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument(
        "--scratch-root",
        type=Path,
        default=Path("/home/nnataren/mnt/phoenix_scratch/Banksy_py"),
    )
    return parser.parse_args()


def split_labels(value: object) -> list[str]:
    """Split semicolon-delimited reference labels and retain known vocabulary."""
    if pd.isna(value):
        return []
    return [part.strip() for part in str(value).split(";") if part.strip() in LABEL_TO_BROAD]


def canonical_label(label: str) -> str:
    """Return a conservative canonical label for overlapping reference names."""
    return LABEL_ALIASES.get(label, label)


def load_reference(path: Path) -> tuple[pd.DataFrame, dict[str, dict[str, float]]]:
    """Load the marker reference and construct gene-to-label evidence weights."""
    reference = pd.read_csv(path)
    evidence: dict[str, dict[str, float]] = defaultdict(dict)
    for row in reference.itertuples(index=False):
        # The primary minor-class column defines subtype evidence. The broader
        # cell_types field is retained in the source but is too ambiguous for scoring.
        labels = [canonical_label(row.minor_class)] if row.minor_class in LABEL_TO_BROAD else []
        if not labels:
            continue
        # A reference gene shared by several labels contributes less to each label.
        ambiguity_weight = 1.0 / math.sqrt(len(labels))
        for label in labels:
            primary_bonus = 1.0 if canonical_label(row.minor_class) == label else 0.65
            evidence[str(row.gene)][label] = max(
                evidence[str(row.gene)].get(label, 0.0), primary_bonus * ambiguity_weight
            )
    return reference, evidence


def load_inputs(
    repo_root: Path, scratch_root: Path
) -> tuple[pd.DataFrame, pd.DataFrame, set[str]]:
    """Load and align Wilcoxon and r1.5 expression records for all samples."""
    marker_frames = []
    expression_frames = []
    for sample in SAMPLES:
        marker_path = (
            scratch_root / "data/xenium/output/ptmt_pc55" / sample / "top_marker_tables"
            / f"rank_genes_groups_{sample}_pc55_nc0.20_r1p5_wilcoxon.csv"
        )
        markers = pd.read_csv(marker_path)
        markers["sample"] = sample
        markers["cluster_id"] = markers["group"].astype(int)
        markers["rank"] = markers.groupby("cluster_id").cumcount() + 1
        marker_frames.append(markers)

        expression_path = (
            repo_root / "data/xenium/processed/ptmt_pc55/cross_sample_dotplot_exports"
            / f"{sample}_ptmt_panel_dotplot_summary.csv"
        )
        expression = pd.read_csv(expression_path)
        expression = expression[np.isclose(expression["resolution"].astype(float), 1.5)].copy()
        expression["cluster_id"] = expression["cluster_id"].astype(int)
        expression_frames.append(
            expression[["sample", "cluster_id", "gene", "mean_expression", "percent_expressing", "n_cells"]]
        )

    markers = pd.concat(marker_frames, ignore_index=True).rename(columns={"names": "gene"})
    expression = pd.concat(expression_frames, ignore_index=True)
    merged = markers.merge(expression, on=["sample", "cluster_id", "gene"], how="left", validate="one_to_one")
    if merged[["mean_expression", "percent_expressing"]].isna().any().any():
        raise ValueError("Marker and expression tables did not align for every sample/cluster/gene.")
    return merged, expression, set(expression["gene"])


def contribution(row: pd.Series, marker_weight: float) -> float:
    """Calculate one gene's positive annotation contribution."""
    if row.gene in GENERIC_OR_STATE or row.logfoldchanges <= 0 or row.scores <= 0 or row.pvals_adj > 0.1:
        return 0.0
    rank_weight = 1.0 / math.log2(float(row.get('rank')) + 1.0)
    significance = min(4.0, max(0.0, -math.log10(max(float(row.pvals_adj), 1e-300)))) / 4.0
    effect = min(3.0, float(row.logfoldchanges)) / 3.0
    coverage = math.sqrt(min(100.0, max(0.0, float(row.percent_expressing))) / 100.0)
    return marker_weight * rank_weight * (0.35 + 0.65 * significance) * (0.4 + 0.6 * effect) * coverage


def score_cluster(
    cluster: pd.DataFrame, gene_evidence: dict[str, dict[str, float]]
) -> tuple[pd.DataFrame, dict[str, list[tuple[str, float]]]]:
    """Score candidate labels with marker-set-size-normalized evidence."""
    label_hits: dict[str, list[tuple[str, float]]] = defaultdict(list)
    for _, row in cluster.iterrows():
        for label, marker_weight in gene_evidence.get(row.gene, {}).items():
            value = contribution(row, marker_weight)
            if value > 0:
                label_hits[label].append((row.gene, value))

    records = []
    for label, hits in label_hits.items():
        ordered = sorted(hits, key=lambda item: item[1], reverse=True)
        # The best five independent genes drive the score; larger reference sets gain no automatic advantage.
        top = ordered[:5]
        score = sum(value for _, value in top) / math.sqrt(max(1, len(top)))
        records.append({"label": label, "broad_class": broad_for(label), "score": score, "n_hits": len(hits)})
    scores = pd.DataFrame(records)
    if not scores.empty:
        scores = scores.sort_values(["score", "n_hits"], ascending=False).reset_index(drop=True)
    return scores, label_hits


def broad_for(label: str) -> str:
    """Map canonical or broad fallback labels to a broad class."""
    if label in {"Epithelial", "T/NK", "Myeloid", "B/mast", "Stromal"}:
        return label
    if label in {"Macrophage", "Myeloid"}:
        return "Myeloid"
    return LABEL_TO_BROAD.get(label, "Unresolved")


def format_gene_list(items: list[tuple[str, float]], limit: int = 8) -> str:
    """Format unique genes ordered by descending evidence."""
    ordered = sorted(items, key=lambda item: item[1], reverse=True)
    return "; ".join(list(dict.fromkeys(gene for gene, _ in ordered))[:limit])


def format_marker_evidence(
    cluster: pd.DataFrame, genes: list[str], gene_evidence: dict[str, dict[str, float]], limit: int = 8
) -> str:
    """Format marker statistics and reference labels for reviewer inspection."""
    records = []
    indexed = cluster.drop_duplicates("gene").set_index("gene")
    for gene in list(dict.fromkeys(genes))[:limit]:
        if gene not in indexed.index:
            continue
        row = indexed.loc[gene]
        labels = "/".join(sorted(gene_evidence.get(gene, {}))) or "unannotated"
        records.append(
            f"{gene} [rank={int(row.get('rank'))}, score={row.scores:.3g}, "
            f"logFC={row.logfoldchanges:.3g}, padj={row.pvals_adj:.3g}, "
            f"mean={row.mean_expression:.3g}, pct={row.percent_expressing:.1f}, ref={labels}]"
        )
    return "; ".join(records)


def assign_clusters(
    merged: pd.DataFrame, gene_evidence: dict[str, dict[str, float]]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Assign provisional labels and retain per-gene cluster evidence."""
    review_rows = []
    evidence_rows = []
    for (sample, cluster_id), cluster in merged.groupby(["sample", "cluster_id"], sort=False):
        scores, hits = score_cluster(cluster, gene_evidence)
        if scores.empty:
            best_label, second_label, best_score, second_score, best_hits = "Unresolved", "", 0.0, 0.0, []
        else:
            best = scores.iloc[0]
            second = scores.iloc[1] if len(scores) > 1 else None
            best_label, best_score = str(best.label), float(best.score)
            second_label = str(second.label) if second is not None else ""
            second_score = float(second.score) if second is not None else 0.0
            best_hits = hits[best_label]

            same_broad = second is not None and best.broad_class == second.broad_class
            close = second_score >= best_score * 0.82
            if close and not same_broad:
                best_label = "Unresolved"
            elif close and same_broad and best.broad_class in {"T/NK", "Myeloid", "Stromal"}:
                # Shared immune/stromal markers often resolve lineage better than subtype in this panel.
                best_label = {"T/NK": "T/NK", "Myeloid": "Myeloid", "Stromal": "Stromal"}[str(best.broad_class)]

        supporting = format_gene_list(best_hits)
        supporting_genes = [gene for gene, _ in sorted(best_hits, key=lambda item: item[1], reverse=True)]
        supporting_evidence = format_marker_evidence(cluster, supporting_genes, gene_evidence)
        support_count = len(set(gene for gene, _ in best_hits))
        margin = (best_score - second_score) / best_score if best_score else 0.0
        if best_label == "Unresolved" or best_score < 0.18 or support_count == 0:
            confidence = "Unresolved"
            best_label = "Unresolved"
        elif support_count >= 4 and best_score >= 0.45 and margin >= 0.18:
            confidence = "High"
        elif support_count >= 2 and best_score >= 0.28 and margin >= 0.08:
            confidence = "Medium"
        else:
            confidence = "Low"

        broad_class = broad_for(best_label)
        contradictory = []
        if not scores.empty:
            for candidate in scores.itertuples(index=False):
                if candidate.broad_class != broad_class and candidate.score >= max(0.12, best_score * 0.45):
                    contradictory.extend(hits[candidate.label][:3])
        top_ranked = cluster.query("logfoldchanges > 0 and pvals_adj < 0.05").nsmallest(10, "rank")["gene"].tolist()
        top_ranked_evidence = format_marker_evidence(cluster, top_ranked, gene_evidence, limit=10)
        contradictory_genes = [gene for gene, _ in contradictory]
        contradictory_evidence = format_marker_evidence(cluster, contradictory_genes, gene_evidence, limit=6)
        rationale = (
            f"Top candidate {scores.iloc[0].label if not scores.empty else 'none'} score={best_score:.3f}; "
            f"runner-up {second_label or 'none'} score={second_score:.3f}; "
            f"{support_count} positive reference-marker hits."
        )
        review_rows.append({
            "sample": sample,
            "treatment": "DMSO" if sample.endswith("DMSO") else "1644",
            "resolution": 1.5,
            "cluster_id": int(cluster_id),
            "proposed_label": best_label,
            "broad_class": broad_class,
            "confidence": confidence,
            "candidate_score": round(best_score, 5),
            "secondary_candidate": second_label,
            "secondary_score": round(second_score, 5),
            "score_margin_fraction": round(margin, 5),
            "supporting_markers": supporting,
            "supporting_marker_evidence": supporting_evidence,
            "contradictory_markers": format_gene_list(contradictory, 6),
            "contradictory_marker_evidence": contradictory_evidence,
            "top_positive_markers": "; ".join(top_ranked),
            "top_positive_marker_evidence": top_ranked_evidence,
            "rationale": rationale,
            "review_notes": "Provisional automated assignment; review against dot plot before config update.",
        })

        for row in cluster.itertuples(index=False):
            evidence_rows.append({
                "sample": sample, "cluster_id": int(cluster_id), "gene": row.gene,
                "rank": int(row.rank), "score": row.scores, "logfoldchanges": row.logfoldchanges,
                "pvals_adj": row.pvals_adj, "mean_expression": row.mean_expression,
                "percent_expressing": row.percent_expressing,
            })

    review = pd.DataFrame(review_rows).sort_values(["sample", "cluster_id"]).reset_index(drop=True)
    evidence = pd.DataFrame(evidence_rows)
    return review, evidence


def select_shared_markers(
    review: pd.DataFrame,
    evidence: pd.DataFrame,
    reference: pd.DataFrame,
    panel_genes: set[str],
    gene_evidence: dict[str, dict[str, float]],
) -> pd.DataFrame:
    """Rank up to five recurrent and specific markers for each assigned label."""
    assigned = review.query("proposed_label != 'Unresolved'")[["sample", "cluster_id", "proposed_label"]]
    data = evidence.merge(assigned, on=["sample", "cluster_id"], how="inner")
    positive = data.query("logfoldchanges > 0 and pvals_adj < 0.05").copy()
    ref_genes = set(reference.gene)
    output = []

    for label, clusters in assigned.groupby("proposed_label"):
        label_data = positive[positive.proposed_label == label]
        n_clusters = len(clusters)
        n_samples = clusters["sample"].nunique()
        other = positive[positive.proposed_label != label]
        rows = []
        for gene, gene_rows in label_data.groupby("gene"):
            allowed_labels = set(gene_evidence.get(gene, {}))
            if gene in GENERIC_OR_STATE:
                continue
            if label in {"T/NK", "Myeloid", "Stromal", "Epithelial", "B/mast"}:
                if not any(broad_for(candidate) == label for candidate in allowed_labels):
                    continue
            elif label not in allowed_labels:
                continue
            sample_count = gene_rows["sample"].nunique()
            cluster_count = len(gene_rows[["sample", "cluster_id"]].drop_duplicates())
            if sample_count < 2:
                continue
            other_clusters = len(other[other.gene == gene][["sample", "cluster_id"]].drop_duplicates())
            target_fraction = cluster_count / n_clusters
            other_denominator = max(1, len(assigned) - n_clusters)
            specificity = max(0.0, target_fraction - other_clusters / other_denominator)
            reference_bonus = 1.2 if gene in ref_genes else 1.0
            recurrence = sample_count / n_samples
            quality = (
                0.35 * recurrence + 0.25 * target_fraction + 0.20 * specificity
                + 0.10 * min(1.0, gene_rows.logfoldchanges.median() / 2.0)
                + 0.10 * min(1.0, gene_rows.percent_expressing.median() / 50.0)
            ) * reference_bonus
            rows.append({
                "cell_type": label,
                "marker": gene,
                "selection_score": quality,
                "sample_recurrence_n": sample_count,
                "sample_recurrence_fraction": recurrence,
                "cluster_recurrence_n": cluster_count,
                "cluster_recurrence_fraction": target_fraction,
                "median_rank": gene_rows["rank"].median(),
                "median_logfoldchange": gene_rows.logfoldchanges.median(),
                "median_pvals_adj": gene_rows.pvals_adj.median(),
                "median_mean_expression": gene_rows.mean_expression.median(),
                "median_percent_expressing": gene_rows.percent_expressing.median(),
                "specificity_fraction_difference": specificity,
                "present_in_harmonised_reference": gene in ref_genes,
                "treatments": "; ".join(sorted({"DMSO" if s.endswith("DMSO") else "1644" for s in gene_rows["sample"]})),
                "caveats": "Broad-lineage marker." if label in {"T/NK", "Myeloid", "Stromal", "Epithelial", "B/mast"} else "",
            })
        ranked = pd.DataFrame(rows)
        if ranked.empty:
            continue
        ranked = ranked.sort_values(
            ["selection_score", "sample_recurrence_n", "cluster_recurrence_n", "median_rank"],
            ascending=[False, False, False, True],
        ).head(5)
        ranked.insert(1, "marker_rank_within_cell_type", range(1, len(ranked) + 1))
        output.append(ranked)

    if not output:
        return pd.DataFrame()
    result = pd.concat(output, ignore_index=True)
    if not set(result.marker).issubset(panel_genes):
        raise ValueError("A selected marker is absent from the panel expression exports.")
    return result.sort_values(["cell_type", "marker_rank_within_cell_type"]).reset_index(drop=True)


def validate(review: pd.DataFrame, shared: pd.DataFrame, panel_genes: set[str]) -> list[str]:
    """Run the task's required completeness and marker-support checks."""
    checks = []
    if len(review) != 230:
        raise ValueError(f"Expected 230 cluster rows, found {len(review)}.")
    if review.duplicated(["sample", "cluster_id"]).any():
        raise ValueError("Duplicate sample/cluster identifiers found.")
    observed = review.groupby("sample").size().to_dict()
    if observed != SAMPLES:
        raise ValueError(f"Cluster counts differ from specification: {observed}")
    if not shared.empty and (shared.sample_recurrence_n < 2).any():
        raise ValueError("A shared marker is supported in fewer than two samples.")
    if not shared.empty and not set(shared.marker).issubset(panel_genes):
        raise ValueError("A proposed shared marker is absent from the panel.")
    checks.extend([
        "PASS: exactly 230 sample/cluster rows.",
        "PASS: all eight samples have the expected cluster counts.",
        "PASS: sample/cluster identifiers are unique.",
        "PASS: all selected markers occur in the panel exports.",
        "PASS: every selected shared marker recurs in at least two samples.",
    ])
    return checks


def write_summary(path: Path, review: pd.DataFrame, shared: pd.DataFrame, checks: list[str]) -> None:
    """Write a concise reviewer-facing Markdown summary."""
    confidence = review.confidence.value_counts().reindex(["High", "Medium", "Low", "Unresolved"], fill_value=0)
    label_counts = review.proposed_label.value_counts()
    strong = review[review.confidence.isin(["High", "Medium"])].proposed_label.value_counts().head(12)
    ambiguous = review[review.confidence.isin(["Low", "Unresolved"])][
        ["sample", "cluster_id", "proposed_label", "secondary_candidate", "supporting_markers"]
    ]
    marker_counts = shared.groupby("cell_type").size() if not shared.empty else pd.Series(dtype=int)
    fewer_than_five = marker_counts[marker_counts < 5]

    lines = [
        "# PTMT PC55 Run 1 r1.5 Cluster Annotation Review",
        "",
        "These assignments and marker selections are provisional. No QC configuration or source marker file was modified.",
        "",
        "## Assignment overview",
        "",
        f"- Clusters reviewed: {len(review)} across {review['sample'].nunique()} samples.",
        f"- Confidence: High {confidence['High']}, Medium {confidence['Medium']}, Low {confidence['Low']}, Unresolved {confidence['Unresolved']}.",
        f"- Proposed labels represented: {len(label_counts)}.",
        "",
        "## Stronger assignments",
        "",
    ]
    lines.extend(f"- {label}: {count} High/Medium clusters." for label, count in strong.items())
    lines.extend(["", "## Ambiguous groups", ""])
    if ambiguous.empty:
        lines.append("- None under the current thresholds.")
    else:
        for row in ambiguous.head(30).itertuples(index=False):
            lines.append(
                f"- `{row.sample}` cluster {row.cluster_id}: {row.proposed_label}; "
                f"secondary {row.secondary_candidate or 'none'}; markers {row.supporting_markers or 'none'}."
            )
        if len(ambiguous) > 30:
            lines.append(f"- {len(ambiguous) - 30} additional Low/Unresolved clusters are listed in the review CSV.")
    lines.extend(["", "## Shared marker coverage", ""])
    if fewer_than_five.empty:
        lines.append("- Every represented cell type with recurrent candidates received five markers.")
    else:
        lines.extend(f"- {label}: {count} supported recurrent markers." for label, count in fewer_than_five.items())
    missing = sorted(set(review.proposed_label) - {"Unresolved"} - set(marker_counts.index))
    lines.extend(f"- {label}: no marker met the two-sample recurrence rule." for label in missing)
    lines.extend([
        "",
        "## Validation",
        "",
        *[f"- {check}" for check in checks],
        "- PNG visual spot-check remains pending because the local viewer was unavailable; the exact Wilcoxon tables and plot-summary exports were used as primary evidence.",
        "",
        "## Review priorities",
        "",
        "- Inspect Low and Unresolved clusters first, especially cross-lineage runner-up calls.",
        "- Review T/NK, myeloid, and stromal subtype calls conservatively where shared lineage markers dominate.",
        "- Approve labels before transferring any assignment into script 01 QC configuration.",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    """Run annotation, shared-marker selection, validation, and export."""
    args = parse_args()
    repo_root = args.repo_root.resolve()
    output_dir = repo_root / "data/xenium/processed/ptmt_pc55/cluster_annotation_review/run_1_r1p5"
    reference_path = repo_root / "data/xenium/raw_data/gene_markers/ptmt/2026-08-14_dotplot_marker_genes_harmonised_scRNAseq_analysis_Kat.csv"

    reference, gene_evidence = load_reference(reference_path)
    merged, expression, panel_genes = load_inputs(repo_root, args.scratch_root)
    review, evidence = assign_clusters(merged, gene_evidence)
    shared = select_shared_markers(review, evidence, reference, panel_genes, gene_evidence)
    checks = validate(review, shared, panel_genes)

    output_dir.mkdir(parents=True, exist_ok=True)
    review.to_csv(output_dir / "ptmt_pc55_run1_r1p5_cluster_label_review.csv", index=False)
    shared.to_csv(output_dir / "ptmt_pc55_run1_r1p5_shared_markers_by_cell_type.csv", index=False)
    write_summary(output_dir / "ptmt_pc55_run1_r1p5_annotation_summary.md", review, shared, checks)

    print("\n".join(checks))
    print(f"Wrote {len(review)} cluster rows and {len(shared)} shared-marker rows to {output_dir}")


if __name__ == "__main__":
    main()
