"""
Title: Shared cohort config and data roots for post-08 workflow steps
Date: 2026-11-19
Updated: 2026-11-22 - cohorts hold clustering facts only. Labels (mapping,
palette) and annotated-object paths moved to versioned annotation sets
(annotation_sets.py, config/annotations/).
Summary: One place for the clustering facts every step from script 08 on
needs: which samples belong to a cohort (response, tissue, branch, resolution,
script 06 object, cluster and UMAP keys) and the coordinate keys. Annotation
sets name a cohort; step configs name an annotation set and may list the
samples they use, in their own order, with plot-only fields. Every step takes
the same `--input-root` and `--output-root` options, both defaulting to the
repository root, so configs never contain machine-specific paths.

Path rules:
- Config files (cohort, annotation set) resolve against the repository root.
- Data read by a step resolves against `--input-root`.
- Data written by a step resolves against `--output-root`.
- Absolute paths are always used as given.
"""

import json
from pathlib import Path

# This module lives in xenium_workflow/, one level below the repository root.
REPO_ROOT = Path(__file__).resolve().parent.parent

REQUIRED_SAMPLE_KEYS = ("sample", "resolution", "recluster_h5ad", "groupby", "umap_key")
# Removed from cohorts on 2026-11-22; they now belong to annotation sets.
RETIRED_COHORT_KEYS = ("palette_json", "mapping_csv")
SHARED_SAMPLE_KEYS = ("project", "branch", "branch_label")


def add_root_arguments(parser):
    """Add the same input/output root options to a step's argument parser."""
    parser.add_argument(
        "--input-root", type=Path, default=REPO_ROOT,
        help="Root for relative input data paths (default: repository root)",
    )
    parser.add_argument(
        "--output-root", type=Path, default=REPO_ROOT,
        help="Root for relative output paths (default: repository root)",
    )


def resolve_path(root, value):
    """Keep absolute paths; join relative paths to `root`."""
    path = Path(value)
    return path if path.is_absolute() else Path(root) / path


def load_json(path):
    """Read one JSON object."""
    return json.loads(Path(path).read_text())


def load_cohort(value):
    """Load and validate a cohort config named relative to the repository root.

    Cohort-level `project`, `branch`, and `branch_label` are copied into every
    sample record so steps see one flat record per sample.
    """
    path = resolve_path(REPO_ROOT, value)
    cohort = load_json(path)
    retired = [key for key in RETIRED_COHORT_KEYS if key in cohort]
    if any("annotated_h5ad" in sample for sample in cohort["samples"]):
        retired.append("samples[].annotated_h5ad")
    if retired:
        raise ValueError(f"{path}: {retired} moved to annotation sets (config/annotations/); remove them")
    if "project" not in cohort or "branch" not in cohort:
        raise ValueError(f"{path}: cohort needs project and branch")
    names = [sample["sample"] for sample in cohort["samples"]]
    if not names or len(names) != len(set(names)):
        raise ValueError(f"{path}: samples must be a nonempty list of unique names")
    shared = {key: cohort[key] for key in SHARED_SAMPLE_KEYS if key in cohort}
    records = []
    for sample in cohort["samples"]:
        missing = [key for key in REQUIRED_SAMPLE_KEYS if key not in sample]
        if missing:
            raise ValueError(f"{path}: {sample.get('sample')} is missing {missing}")
        overlap = sorted(set(shared) & set(sample))
        if overlap:
            raise ValueError(f"{path}: {sample['sample']} repeats cohort-level fields {overlap}")
        records.append({**shared, **sample})
    return {
        "path": path,
        "name": path.stem,
        "project": cohort["project"],
        "branch": cohort["branch"],
        "samples": records,
        "coord_keys": cohort.get("coord_keys", ["x", "y"]),
    }


def select_samples(cohort, entries=None, requested=None):
    """Return sample records in a step's order, merged with its plot-only fields.

    `entries` is the step config's optional `samples` list: names, or objects
    with `sample` plus step-specific fields such as `title` or `facet`. It
    defaults to every cohort sample in cohort order. A step may add fields but
    may not redefine a cohort fact. `requested` (from `--sample`) then filters
    the list without changing its order.
    """
    by_name = {record["sample"]: record for record in cohort["samples"]}
    if entries is None:
        entries = list(by_name)
    records = []
    for entry in entries:
        extra = {"sample": entry} if isinstance(entry, str) else dict(entry)
        name = extra.pop("sample")
        if name not in by_name:
            raise ValueError(f"Sample {name!r} is not in cohort {cohort['path']}")
        fact = by_name[name]
        clash = sorted(key for key in extra if key in fact and extra[key] != fact[key])
        if clash:
            raise ValueError(f"{name}: step config redefines cohort fields {clash}")
        records.append({"title": f"{name} · r{fact['resolution']}", **fact, **extra})
    names = [record["sample"] for record in records]
    if len(names) != len(set(names)):
        raise ValueError("Step config lists a sample more than once")
    if requested:
        unknown = sorted(set(requested) - set(names))
        if unknown:
            raise ValueError(f"Requested samples are not configured for this step: {unknown}")
        records = [record for record in records if record["sample"] in set(requested)]
    return records
