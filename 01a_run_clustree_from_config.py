#!/usr/bin/env python3
"""
Title: Config-driven clustree batch runner
Date: 2026-08-19
Summary: Read a JSON batch config, resolve one or more BANKSY clustree inputs,
and call 01a_clustree_cluster_resolution_qc.R for each sample. This is intended
for local post-clustering or post-reclustering clustree review without long
hand-written Rscript commands.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run 01a_clustree_cluster_resolution_qc.R from a batch JSON config."
    )
    parser.add_argument("--config", required=True, help="Path to a clustree batch JSON config.")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running Rscript.")
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="Continue after an Rscript failure and report all failures at the end.",
    )
    return parser.parse_args()


def load_json(path: Path) -> dict[str, Any]:
    """Load a JSON file as a dictionary."""
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return data


def repo_root() -> Path:
    """Return the directory containing this runner."""
    return Path(__file__).resolve().parent


def as_list(value: Any) -> list[Any]:
    """Return `value` as a list, treating missing values as empty."""
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def discover_source_configs(batch_config: dict[str, Any], config_path: Path) -> list[Path]:
    """Resolve sample source configs listed directly or discovered from a directory."""
    root = repo_root()
    paths: list[Path] = []

    for item in as_list(batch_config.get("source_configs")):
        item_path = Path(item)
        paths.append(item_path if item_path.is_absolute() else root / item_path)

    source_config_dir = batch_config.get("source_config_dir")
    if source_config_dir:
        paths.extend(sorted((root / source_config_dir).glob("*.json")))

    if not paths:
        raise ValueError(f"No source_configs or source_config_dir configured in {config_path}")

    return paths


def format_context(sample_config: dict[str, Any]) -> dict[str, str]:
    """Build string substitutions shared by path templates."""
    context = {key: str(value) for key, value in sample_config.items() if not isinstance(value, list)}
    res_label = sample_config.get("res_label")
    if not isinstance(res_label, list) or not res_label:
        raise ValueError(
            f"Sample config for {sample_config.get('dataset_name', '<unknown>')} lacks res_label list"
        )
    context["resolutions"] = "_".join(str(value) for value in res_label)
    return context


def resolve_path(template: str, context: dict[str, str]) -> Path:
    """Format a path template against a sample context."""
    value = template.format(**context)
    path = Path(value)
    return path if path.is_absolute() else repo_root() / path


def build_command(batch_config: dict[str, Any], source_config_path: Path) -> tuple[list[str], Path, str]:
    """Build one Rscript command from a source sample config."""
    sample_config = load_json(source_config_path)
    context = format_context(sample_config)
    dataset_name = context.get("dataset_name")
    if not dataset_name:
        raise ValueError(f"Missing dataset_name in {source_config_path}")

    cluster_csv = resolve_path(batch_config["cluster_csv_template"], context)
    output_dir = resolve_path(batch_config["output_dir_template"], context)
    r_script = resolve_path(batch_config.get("r_script", "01a_clustree_cluster_resolution_qc.R"), context)
    cluster_prefix = batch_config["cluster_prefix_template"].format(**context)
    cluster_suffix = batch_config.get("cluster_suffix_template", "").format(**context)

    command = [
        batch_config.get("rscript_executable", "Rscript"),
        str(r_script),
        "--cluster_csv",
        str(cluster_csv),
        "--dataset_name",
        dataset_name,
        "--cluster_prefix",
        cluster_prefix,
        "--cluster_suffix",
        cluster_suffix,
        "--output_dir",
        str(output_dir),
        "--width",
        str(batch_config.get("width", 12)),
        "--height",
        str(batch_config.get("height", 8)),
        "--dpi",
        str(batch_config.get("dpi", 300)),
    ]

    if batch_config.get("include_qc_config", False):
        qc_template = batch_config.get("qc_config_template")
        if not qc_template:
            raise ValueError("include_qc_config is true but qc_config_template is missing")
        qc_config = resolve_path(qc_template, context)
        command.extend(["--qc_config", str(qc_config)])

    return command, cluster_csv, dataset_name


def main() -> None:
    """Run all configured clustree jobs."""
    args = parse_args()
    config_path = Path(args.config).resolve()
    batch_config = load_json(config_path)
    source_configs = discover_source_configs(batch_config, config_path)

    failures: list[str] = []
    for source_config in source_configs:
        command, cluster_csv, dataset_name = build_command(batch_config, source_config)
        if not cluster_csv.exists():
            print(f"Missing input for {dataset_name}: {cluster_csv}")
            if not args.dry_run:
                print(f"Skipping {dataset_name}")
                continue

        action = "Would run" if args.dry_run else "Running"
        print(f"{action}:", " ".join(command))
        if args.dry_run:
            continue

        result = subprocess.run(command, cwd=repo_root(), check=False)
        if result.returncode != 0:
            message = f"{dataset_name} failed with exit code {result.returncode}"
            if not args.keep_going:
                raise SystemExit(message)
            failures.append(message)

    if failures:
        raise SystemExit("Clustree failures:\n" + "\n".join(failures))


if __name__ == "__main__":
    main()
