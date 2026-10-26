"""
Title: Plot reviewed cell-type proportions across samples
Date: 2026-10-26
Summary: Validate an existing cell-type count table against configured samples,
branches, resolutions, and a shared palette. Export a zero-filled proportion
table and optionally tissue-faceted 100% stacked bars as PNG, SVG, and PDF
without loading AnnData objects.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd


def resolve_path(root: Path, value: str) -> Path:
    """Resolve configuration paths relative to the project root."""
    path = Path(value)
    return path if path.is_absolute() else root / path


def prepare_proportions(cfg: dict, root: Path) -> tuple[pd.DataFrame, dict]:
    """Validate counts and provenance, then fill missing types and normalize.

    Each denominator contains all annotated post-filtered cells in one sample.
    Reject unexpected samples or conditions rather than silently pooling them.
    """
    palette = json.loads(resolve_path(root, cfg["palette_json"]).read_text())
    order, colors = palette["cell_type_order"], palette["colors"]
    if not order or len(order) != len(set(order)) or set(order) != set(colors):
        raise ValueError("Palette must contain unique ordered labels and matching colors")
    if not all(matplotlib.colors.is_color_like(color) for color in colors.values()):
        raise ValueError("Palette contains invalid colors")
    counts = pd.read_csv(resolve_path(root, cfg["counts_csv"]), dtype={"resolution": str})
    required = ["sample", "branch", "resolution", "cell_type", "n_cells"]
    if not set(required).issubset(counts.columns) or counts[required].isna().any().any():
        raise ValueError("Count table is missing required columns or values")
    samples = cfg["samples"]
    names = [sample["sample"] for sample in samples]
    if not names or len(names) != len(set(names)) or set(counts["sample"]) != set(names):
        raise ValueError("Count table must match the unique configured samples exactly")
    if not set(counts["cell_type"]).issubset(order):
        raise ValueError("Count table contains cell types absent from the shared palette")
    if counts.duplicated(["sample", "cell_type"]).any():
        raise ValueError("Duplicate sample/cell-type rows; conditions must not be pooled")
    values = pd.to_numeric(counts["n_cells"], errors="raise")
    if not np.isfinite(values).all() or (values < 0).any() or (values % 1 != 0).any():
        raise ValueError("Cell counts must be finite, nonnegative integers")
    counts["n_cells"] = values.astype("int64")
    if int(values.sum()) != cfg["expected_total_cells"]:
        raise ValueError("Input total does not match expected_total_cells")
    tables = []
    for sample in samples:
        subset = counts.loc[counts["sample"] == sample["sample"]]
        if not (subset["branch"] == sample["branch"]).all() or not np.isclose(
            pd.to_numeric(subset["resolution"]), float(sample["resolution"])
        ).all():
            raise ValueError(f"Unexpected branch or resolution for {sample['sample']}")
        # Reindex by the canonical palette, not by the types present in this sample.
        table = subset.set_index("cell_type")[["n_cells"]].reindex(order, fill_value=0)
        total = int(table["n_cells"].sum())
        if total <= 0:
            raise ValueError(f"No cells available for {sample['sample']}")
        table["sample"] = sample["sample"]
        table["branch"] = sample["branch"]
        table["resolution"] = sample["resolution"]
        table["sample_total_cells"] = total
        table["proportion"] = table["n_cells"] / total
        if not np.isclose(table["proportion"].sum(), 1):
            raise ValueError("Sample proportions do not sum to one")
        tables.append(table.reset_index())
    return pd.concat(tables, ignore_index=True), palette


def get_panels(cfg: dict) -> list[tuple[str, list[dict]]]:
    """Group samples into configured facets, preserving order within each facet."""
    order = cfg["plot"].get("facet_order")
    if order is None:
        return [("", cfg["samples"])]
    tissues = [sample.get("tissue") for sample in cfg["samples"]]
    if not order or len(order) != len(set(order)) or set(tissues) != set(order):
        raise ValueError("Each sample must have a tissue matching the unique facet_order")
    return [(tissue, [s for s in cfg["samples"] if s["tissue"] == tissue]) for tissue in order]


def plot_proportions(table: pd.DataFrame, palette: dict, cfg: dict, stem: Path) -> None:
    """Render equal-width bars in optional tissue facets with shared scale and legend."""
    settings = cfg["plot"]
    panels = get_panels(cfg)
    faceted = "facet_order" in settings
    plt.rcParams.update(
        {
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": float(settings.get("font_size", 11)),
        }
    )
    fig, axes = plt.subplots(
        1, len(panels), figsize=settings["figsize"], sharey=True, squeeze=False,
        gridspec_kw={"width_ratios": [len(samples) for _, samples in panels]},
    )
    for ax, (title, samples) in zip(axes[0], panels):
        names = [sample["sample"] for sample in samples]
        x = np.arange(len(names))
        bottom = np.zeros(len(names))
        for cell_type in palette["cell_type_order"]:
            heights = table.loc[table["cell_type"] == cell_type].set_index("sample").loc[names, "proportion"].to_numpy()
            ax.bar(x, heights, bottom=bottom, width=0.72,
                   color=palette["colors"][cell_type], label=cell_type)
            bottom += heights
        labels = [sample.get("display_label", f"{sample['sample']}\nr{float(sample['resolution']):.2f}") for sample in samples]
        ax.set_xticks(x, labels, rotation=35, ha="right", rotation_mode="anchor",
                      fontsize=float(settings.get("xtick_fontsize", settings.get("font_size", 11))))
        # Explicit limits plus proportional facet widths keep physical bar widths equal.
        ax.set_xlim(-0.5, len(names) - 0.5)
        ax.set_ylim(0, 1)
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_axisbelow(True)
        ax.grid(axis="y", alpha=0.2)
        if faceted:
            ax.set_title(
                title, pad=14,
                fontsize=float(settings.get("facet_title_fontsize", 12)),
                fontweight="bold",
            )
    for ax in axes[0]:
        ax.tick_params(axis="y", labelsize=float(settings.get("ytick_fontsize", settings.get("font_size", 11))))
    axes[0, 0].set_ylabel(
        "Cell proportion", fontsize=float(settings.get("ylabel_fontsize", settings.get("font_size", 11)))
    )
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        title="Cell type",
        bbox_to_anchor=(0.82, 0.82),
        loc="upper left",
        frameon=False,
        fontsize=float(settings.get("legend_fontsize", settings.get("font_size", 11))),
        title_fontsize=float(settings.get("legend_title_fontsize", settings.get("font_size", 11))),
    )
    show_figure_title = settings.get("show_figure_title", True)
    show_caption = settings.get("show_caption", True)
    if show_figure_title:
        fig.suptitle(
            settings["title"], x=0.08, ha="left", y=0.96,
            fontsize=float(settings.get("title_fontsize", 16)),
        )
    if show_caption:
        fig.text(
            0.08, 0.025, settings["caption"],
            fontsize=float(settings.get("caption_fontsize", 9)), color="#444444",
        )
    fig.subplots_adjust(
        left=0.08,
        right=0.80,
        bottom=0.32 if show_caption else 0.24,
        top=0.82 if show_figure_title else 0.92,
        wspace=0.16,
    )
    for extension in ("png", "svg", "pdf"):
        fig.savefig(stem.with_suffix(f".{extension}"), dpi=settings["dpi"])
    plt.close(fig)


def main() -> None:
    """Read a run config, validate counts, and optionally write plot artifacts."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="JSON run configuration")
    parser.add_argument("--validate-only", action="store_true", help="Check input without writing outputs")
    parser.add_argument("--overwrite", action="store_true", help="Replace existing output files")
    args = parser.parse_args()
    config_path = Path(args.config).resolve()
    cfg = json.loads(config_path.read_text())
    root = config_path.parents[cfg.get("project_root_parent_depth", 3)]
    get_panels(cfg)
    table, palette = prepare_proportions(cfg, root)
    print(f"Validated {table['sample'].nunique()} samples; {table['n_cells'].sum():,} cells; all proportions sum to 1")
    if args.validate_only:
        return
    output = resolve_path(root, cfg["output_dir"])
    stem = output / cfg["output_stem"]
    paths = [stem.with_suffix(f".{extension}") for extension in ("png", "svg", "pdf", "csv")]
    if not args.overwrite and any(path.exists() for path in paths):
        raise FileExistsError("Outputs already exist; use --overwrite to replace them")
    os.makedirs(output, exist_ok=True)
    plot_proportions(table, palette, cfg, stem)
    table.to_csv(stem.with_suffix(".csv"), index=False)
    print(f"Wrote PNG, SVG, PDF, and proportion CSV: {stem}")


if __name__ == "__main__":
    main()
