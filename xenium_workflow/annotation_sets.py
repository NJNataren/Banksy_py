"""
Title: Versioned cell-type annotation sets for script 08 and later steps
Date: 2026-11-22
Summary: An annotation set is one reviewed labelling of one cohort's clusters,
kept in `config/annotations/<project>/<cohort>/<version>/`:

    annotation.json  version, cohort, parent, status (draft | released), release hash
    mapping.csv      one row per source cluster: sample, branch, resolution,
                     groupby, cluster_id, cell_type_tier_1[, cell_type_tier_2 ...],
                     annotation_status, notes
    palette.json     per tier: label order and colours; tiers after the first
                     also give each label's parent in the tier above
    CHANGES.md       what the version represents and why (written by hand)

Tiers: `tier_1` is the broad cell type, `tier_2` a finer level, and so on. The
obs column for a tier is `cell_type_<tier>` (e.g. `cell_type_tier_1`). A blank
finer-tier label inherits the label from the tier above (and its colour); the
inheritance is recorded so it is not mistaken for a real split.

Released sets are locked: `release` stores a SHA-256 of mapping.csv and
palette.json, and loading a released set whose files no longer match fails.
Edits after release therefore mean a new version (`new`).

This module also derives every output path from the set, so step configs hold
no folders:

    annotated H5AD  data/xenium/processed/<project>/cell_type_annotation/<cohort>/<version>/<sample>/
                      adata_expression_clean_<sample>_<branch>_cell_type_annotated_<version>_r<res>.h5ad
    step outputs    data/xenium/{output,processed}/<project>/<step>/<cohort>/<version>/[<level>/][<subdir>/]

Command line (run from the repository root with the banksy environment;
moved into the xenium_workflow package on 2026-11-22):

    python -m xenium_workflow.annotation_sets check   config/annotations/vbct/<cohort>/v2
    python -m xenium_workflow.annotation_sets diff    config/annotations/vbct/<cohort>/v2 [--against v1]
    python -m xenium_workflow.annotation_sets release config/annotations/vbct/<cohort>/v2
    python -m xenium_workflow.annotation_sets new     config/annotations/vbct/<cohort>/v1

`check` and `diff` only read. `release` rewrites annotation.json of a draft.
`new` creates the next version folder as a draft copy of a released parent.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import sys
from datetime import date
from pathlib import Path

import pandas as pd
from matplotlib.colors import is_color_like

from .cohorts import REPO_ROOT, load_cohort, load_json, resolve_path

ANNOTATION_ROOT = Path("config/annotations")
VERSION_PATTERN = re.compile(r"^v([1-9][0-9]*)$")
LEVEL_PATTERN = re.compile(r"^tier_([1-9])$")
STATUSES = ("draft", "released")
# Row-level review state in mapping.csv. A released set must be fully reviewed.
ROW_STATUSES = ("reviewed", "proposed")
KEY_COLUMNS = ["sample", "branch", "resolution", "groupby", "cluster_id"]
HASHED_FILES = ("mapping.csv", "palette.json")
STEP_NAMES = ("cell_type_annotation", "cell_type_proportions", "spatial_analysis", "poster")
# Step-config keys replaced by the annotation set and derived paths (2026-11-22).
RETIRED_STEP_KEYS = ("cohort", "figure_dir", "summary_dir", "result_dir", "table_dir",
                     "output_dir", "counts_csv", "script10_result_dir")


def label_column(level: str) -> str:
    """obs column holding the labels of one tier, e.g. `cell_type_tier_1`."""
    if not LEVEL_PATTERN.match(level):
        raise ValueError(f"Label level must look like tier_1, tier_2, ...; got {level!r}")
    return f"cell_type_{level}"


def inherited_column(level: str) -> str:
    """obs/mapping column flagging labels inherited from the tier above."""
    return f"{label_column(level)}_inherited"


def level_title(level: str) -> str:
    """Human-readable tier name for legends, e.g. `tier 1`."""
    return level.replace("_", " ")


# --------------------------------------------------------------------------- palette


def load_palette(path: Path) -> dict:
    """Load and validate a tiered palette.

    Returns `{"levels": [...], "tiers": {level: {"order", "colors", "parent"}}}`
    where `order`, `colors` and `parent` are the *effective* vocabulary: finer
    tiers include the labels inherited from the tier above, placed after that
    parent's own children, with the parent's colour.
    """
    payload = load_json(path)
    levels = sorted(payload, key=lambda name: int(LEVEL_PATTERN.match(name).group(1))
                    if LEVEL_PATTERN.match(name) else 0)
    if not levels or any(not LEVEL_PATTERN.match(level) for level in levels):
        raise ValueError(f"{path}: top-level keys must be tier_1, tier_2, ...")
    if levels != [f"tier_{i}" for i in range(1, len(levels) + 1)]:
        raise ValueError(f"{path}: tiers must be contiguous from tier_1; found {levels}")
    tiers = {}
    for index, level in enumerate(levels):
        spec = payload[level]
        order, colors = spec.get("order"), spec.get("colors")
        if not order or len(order) != len(set(order)) or set(order) != set(colors or {}):
            raise ValueError(f"{path}: {level} needs a unique `order` matching `colors` exactly")
        bad = {label: color for label, color in colors.items() if not is_color_like(color)}
        if bad:
            raise ValueError(f"{path}: {level} has invalid colours {bad}")
        if index == 0:
            if "parent" in spec:
                raise ValueError(f"{path}: tier_1 has no parent tier")
            tiers[level] = {"order": list(order), "colors": dict(colors), "parent": {}}
            continue
        above = tiers[levels[index - 1]]
        parent = spec.get("parent") or {}
        if set(parent) != set(order):
            raise ValueError(f"{path}: {level}.parent must name a parent for every {level} label")
        unknown = sorted(set(parent.values()) - set(above["order"]))
        if unknown:
            raise ValueError(f"{path}: {level} parents are not {levels[index - 1]} labels: {unknown}")
        clash = sorted(label for label in order if label in above["order"] and parent[label] != label)
        if clash:
            raise ValueError(f"{path}: {level} reuses {levels[index - 1]} labels under another parent: {clash}")
        # Effective vocabulary: each parent's children, then the parent itself
        # (for clusters left blank at this tier), unless a child has its name.
        eff_order, eff_colors, eff_parent = [], {}, {}
        for upper in above["order"]:
            children = [label for label in order if parent[label] == upper]
            for label in children:
                eff_order.append(label)
                eff_colors[label] = colors[label]
                eff_parent[label] = upper
            if upper not in children:
                eff_order.append(upper)
                eff_colors[upper] = above["colors"][upper]
                eff_parent[upper] = upper
        tiers[level] = {"order": eff_order, "colors": eff_colors, "parent": eff_parent,
                        "declared": list(order)}
    return {"levels": levels, "tiers": tiers}


# --------------------------------------------------------------------------- set


def set_hash(directory: Path) -> str:
    """SHA-256 over the hashed files' names and bytes, in fixed order."""
    digest = hashlib.sha256()
    for name in HASHED_FILES:
        digest.update(name.encode() + b"\0")
        digest.update((directory / name).read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def find_set_dir(value) -> Path:
    """Resolve a set folder (or its annotation.json) relative to the repository root."""
    path = resolve_path(REPO_ROOT, value)
    if path.name == "annotation.json":
        path = path.parent
    if not (path / "annotation.json").is_file():
        raise FileNotFoundError(f"No annotation.json in {path}")
    return path


def load_annotation(value, verify_release: bool = True) -> dict:
    """Load one annotation set with its cohort and palette.

    The folder must sit at config/annotations/<project>/<cohort name>/<version>/
    so its location always matches what it labels. With `verify_release`, a
    released set whose files changed since release raises an error.
    """
    directory = find_set_dir(value)
    meta = load_json(directory / "annotation.json")
    for key in ("version", "cohort", "status"):
        if key not in meta:
            raise ValueError(f"{directory}/annotation.json is missing {key!r}")
    version = meta["version"]
    if not VERSION_PATTERN.match(version) or version != directory.name:
        raise ValueError(f"{directory}: version {version!r} must be vN and match the folder name")
    if meta["status"] not in STATUSES:
        raise ValueError(f"{directory}: status must be one of {STATUSES}")
    parent = meta.get("parent")
    if parent is not None and (not VERSION_PATTERN.match(parent) or parent == version):
        raise ValueError(f"{directory}: parent must be another vN or null")
    cohort = load_cohort(meta["cohort"])
    expected = REPO_ROOT / ANNOTATION_ROOT / cohort["project"] / cohort["name"] / version
    if directory.resolve() != expected.resolve():
        raise ValueError(f"{directory}: set for cohort {cohort['name']} must live at {expected}")
    palette = load_palette(directory / "palette.json")
    current = set_hash(directory)
    if verify_release and meta["status"] == "released" and meta.get("sha256") != current:
        raise ValueError(
            f"{directory}: released set was edited after release (hash mismatch). "
            "Restore it, or create a new version with `python -m xenium_workflow.annotation_sets new`."
        )
    return {
        "dir": directory,
        "version": version,
        "parent": parent,
        "status": meta["status"],
        "meta": meta,
        "cohort": cohort,
        "project": cohort["project"],
        "cohort_name": cohort["name"],
        "palette": palette,
        "levels": palette["levels"],
        "mapping_path": directory / "mapping.csv",
        "sha256": current,
    }


def sibling_dir(annotation: dict, version: str) -> Path:
    """Folder of another version of the same cohort."""
    return annotation["dir"].parent / version


def stamp(annotation: dict) -> dict:
    """Provenance written to `uns['cell_type_annotation']` by script 08."""
    return {
        "annotation_version": annotation["version"],
        "annotation_parent": annotation["parent"] or "",
        "annotation_status": annotation["status"],
        "annotation_sha256": annotation["sha256"],
        "annotation_dir": str(annotation["dir"].relative_to(REPO_ROOT)),
        "cohort": annotation["cohort_name"],
        "levels": list(annotation["levels"]),
    }


def check_stamp(found: dict, annotation: dict, where: str) -> None:
    """Fail unless an annotated object was written from exactly this set."""
    if not found:
        raise ValueError(f"{where}: no uns['cell_type_annotation'] stamp; rerun script 08 for this set")
    wanted = {"annotation_version": annotation["version"], "cohort": annotation["cohort_name"],
              "annotation_sha256": annotation["sha256"]}
    differ = {key: (found.get(key), value) for key, value in wanted.items() if str(found.get(key)) != value}
    if differ:
        raise ValueError(f"{where}: written from a different annotation set (found, expected): {differ}. "
                         "Rerun script 08 for this set.")


def require_released(annotation: dict, allow_draft: bool, step: str) -> None:
    """Expensive or reported steps need a released set unless explicitly overridden."""
    if annotation["status"] != "released" and not allow_draft:
        raise ValueError(f"{step} needs a released annotation set; {annotation['dir']} is a draft. "
                         "Release it, or pass --allow-draft for a local test.")


# --------------------------------------------------------------------------- mapping


def load_mapping(annotation: dict) -> pd.DataFrame:
    """Load and validate mapping.csv against the set's cohort and palette.

    Returns one row per source cluster with every tier filled in (inherited
    labels included) and a `<column>_inherited` flag for each tier after the first.
    """
    path = annotation["mapping_path"]
    mapping = pd.read_csv(path, dtype=str, keep_default_na=False)
    levels = annotation["levels"]
    tiers = annotation["palette"]["tiers"]
    required = KEY_COLUMNS + [label_column(levels[0]), "annotation_status"]
    missing = sorted(set(required) - set(mapping.columns))
    if missing:
        raise KeyError(f"{path}: missing columns {missing}")
    tier_columns = [column for column in mapping.columns if column.startswith("cell_type_tier_")
                    and not column.endswith("_inherited")]
    extra = sorted(set(tier_columns) - {label_column(level) for level in levels})
    if extra:
        raise ValueError(f"{path}: tier columns {extra} have no palette section")
    for column in KEY_COLUMNS + [label_column(levels[0]), "annotation_status"]:
        if (mapping[column].str.strip() == "").any():
            raise ValueError(f"{path}: blank values in {column}")
    if mapping.duplicated(KEY_COLUMNS).any():
        rows = mapping.loc[mapping.duplicated(KEY_COLUMNS, keep=False), KEY_COLUMNS]
        raise ValueError(f"{path}: duplicate cluster rows:\n{rows.to_string(index=False)}")
    bad_status = sorted(set(mapping["annotation_status"]) - set(ROW_STATUSES))
    if bad_status:
        raise ValueError(f"{path}: annotation_status must be one of {ROW_STATUSES}; found {bad_status}")

    # Rows must describe this cohort's clustering exactly.
    cohort = {record["sample"]: record for record in annotation["cohort"]["samples"]}
    unknown = sorted(set(mapping["sample"]) - set(cohort))
    if unknown:
        raise ValueError(f"{path}: samples not in cohort {annotation['cohort_name']}: {unknown}")
    for name, record in cohort.items():
        rows = mapping[mapping["sample"] == name]
        if rows.empty:
            raise ValueError(f"{path}: no rows for cohort sample {name}")
        wrong = rows[(rows["branch"] != record["branch"]) | (rows["groupby"] != record["groupby"])
                     | (rows["resolution"].astype(float) != float(record["resolution"]))]
        if not wrong.empty:
            raise ValueError(f"{path}: {name} rows do not match the cohort branch/resolution/groupby")
        # Spreadsheet editors rewrite 1.00 as 1 and 0.90 as 0.9; use the cohort's
        # spelling so comparisons between sets and output names stay stable.
        mapping.loc[rows.index, "resolution"] = str(record["resolution"])

    first = label_column(levels[0])
    unknown = sorted(set(mapping[first]) - set(tiers[levels[0]]["order"]))
    if unknown:
        raise ValueError(f"{path}: {first} labels missing from palette: {unknown}")
    for above, level in zip(levels, levels[1:]):
        column, upper = label_column(level), label_column(above)
        if column not in mapping:
            mapping[column] = ""
        blank = mapping[column].str.strip() == ""
        mapping[inherited_column(level)] = blank
        mapping.loc[blank, column] = mapping.loc[blank, upper]
        declared = tiers[level]["declared"]
        unknown = sorted(set(mapping.loc[~blank, column]) - set(declared))
        if unknown:
            raise ValueError(f"{path}: {column} labels missing from palette: {unknown}")
        parent = tiers[level]["parent"]
        bad = mapping[~blank & (mapping[column].map(parent) != mapping[upper])]
        if not bad.empty:
            raise ValueError(f"{path}: {column} labels sit under the wrong {upper}:\n"
                             f"{bad[KEY_COLUMNS + [upper, column]].to_string(index=False)}")
    return mapping


def mapping_for_sample(mapping: pd.DataFrame, record: dict) -> pd.DataFrame:
    """Rows of the mapping for one cohort sample."""
    return mapping[mapping["sample"] == record["sample"]].copy()


# --------------------------------------------------------------------------- changes


def cluster_changes(child: dict, parent: dict) -> pd.DataFrame:
    """Clusters whose label differs between two sets of the same cohort, at any tier.

    A tier the parent lacks is compared against the parent's finest tier
    (inheritance), so splitting a type in tier 2 shows as a change.
    """
    if child["cohort_name"] != parent["cohort_name"]:
        raise ValueError("Sets label different cohorts; cluster-level changes are undefined")
    new, old = load_mapping(child), load_mapping(parent)
    levels = list(dict.fromkeys(parent["levels"] + child["levels"]))
    for frame, own in ((new, child["levels"]), (old, parent["levels"])):
        for level in levels:
            if level not in own:
                frame[label_column(level)] = frame[label_column(own[-1])]
    columns = [label_column(level) for level in levels]
    merged = old[KEY_COLUMNS + columns].merge(
        new[KEY_COLUMNS + columns], on=KEY_COLUMNS, how="outer", suffixes=("_old", "_new"), indicator=True)
    if (merged["_merge"] != "both").any():
        raise ValueError("Sets cover different clusters; compare only sets of the same clustering")
    changed = pd.Series(False, index=merged.index)
    for column in columns:
        changed |= merged[f"{column}_old"] != merged[f"{column}_new"]
    out = merged.loc[changed, KEY_COLUMNS + [f"{c}_{s}" for c in columns for s in ("old", "new")]]
    out.insert(0, "to_version", child["version"])
    out.insert(0, "from_version", parent["version"])
    return out.reset_index(drop=True)


# --------------------------------------------------------------------------- paths


def annotated_h5ad_path(annotation: dict, record: dict) -> Path:
    """Repository-relative path of the script 08 annotated object for one sample."""
    name = (f"adata_expression_clean_{record['sample']}_{record['branch']}"
            f"_cell_type_annotated_{annotation['version']}_r{record['resolution']}.h5ad")
    return (Path("data/xenium/processed") / annotation["project"] / "cell_type_annotation"
            / annotation["cohort_name"] / annotation["version"] / record["sample"] / name)


def step_dir(annotation: dict, kind: str, step: str, level: str | None = None,
             subdir: str | None = None) -> Path:
    """Repository-relative output folder of a step for this set (and level)."""
    if kind not in ("output", "processed"):
        raise ValueError("kind must be output (figures) or processed (tables, H5ADs)")
    if step not in STEP_NAMES:
        raise ValueError(f"Unknown step {step!r}; expected one of {STEP_NAMES}")
    path = (Path("data/xenium") / kind / annotation["project"] / step
            / annotation["cohort_name"] / annotation["version"])
    if level is not None:
        label_column(level)
        path /= level
    if subdir:
        if Path(subdir).is_absolute() or ".." in Path(subdir).parts:
            raise ValueError(f"output_subdir must be a plain relative folder: {subdir!r}")
        path /= subdir
    return path


# --------------------------------------------------------------------------- step configs


def add_annotation_arguments(parser, label_level: bool = True) -> None:
    """Add --annotation (and --label-level) overrides to a step's parser."""
    parser.add_argument("--annotation", help="Annotation set folder; overrides the config's `annotation`")
    if label_level:
        parser.add_argument("--label-level", help="tier_1, tier_2, ...; overrides the config's `label_level`")


def resolve_step(cfg: dict, args, use_level: bool = True) -> tuple[dict, str | None]:
    """Return the annotation set and label level a step config asks for."""
    retired = [key for key in RETIRED_STEP_KEYS if key in cfg]
    if retired:
        raise ValueError(f"Step config keys {retired} were replaced on 2026-11-22 by `annotation` "
                         "and derived paths (xenium_workflow/annotation_sets.py); remove them")
    value = getattr(args, "annotation", None) or cfg.get("annotation")
    if not value:
        raise KeyError("Step config needs `annotation` (or pass --annotation)")
    annotation = load_annotation(value)
    if not use_level:
        return annotation, None
    level = getattr(args, "label_level", None) or cfg.get("label_level", "tier_1")
    if level not in annotation["levels"]:
        raise ValueError(f"Label level {level!r} is not in set {annotation['version']} "
                         f"(levels: {annotation['levels']})")
    return annotation, level


def level_vocabulary(annotation: dict, level: str) -> tuple[list[str], dict[str, str]]:
    """Effective label order and colours for one tier."""
    tier = annotation["palette"]["tiers"][level]
    return list(tier["order"]), dict(tier["colors"])


# --------------------------------------------------------------------------- CLI


def _print_summary(annotation: dict, mapping: pd.DataFrame) -> None:
    print(f"Set      {annotation['dir'].relative_to(REPO_ROOT)}")
    print(f"Cohort   {annotation['cohort_name']} ({len(annotation['cohort']['samples'])} samples)")
    print(f"Status   {annotation['status']}   parent {annotation['parent'] or '-'}")
    print(f"Levels   {', '.join(annotation['levels'])}")
    print(f"Clusters {len(mapping)}   statuses {mapping['annotation_status'].value_counts().to_dict()}")
    for level in annotation["levels"]:
        column = label_column(level)
        used = mapping[column].nunique()
        note = ""
        if level != annotation["levels"][0]:
            note = f"; {int(mapping[inherited_column(level)].sum())} clusters inherit from the tier above"
        print(f"  {column}: {used} labels used{note}")
    print(f"SHA-256  {annotation['sha256']}")


def cmd_check(args) -> None:
    annotation = load_annotation(args.set)
    _print_summary(annotation, load_mapping(annotation))
    print("OK")


def cmd_diff(args) -> None:
    child = load_annotation(args.set)
    against = args.against or child["parent"]
    if not against:
        sys.exit(f"{child['version']} has no parent; pass --against vN")
    parent = load_annotation(sibling_dir(child, against))
    changes = cluster_changes(child, parent)
    if changes.empty:
        print(f"No cluster labels differ between {against} and {child['version']}.")
        return
    with pd.option_context("display.width", 200, "display.max_rows", None, "display.max_colwidth", 40):
        print(changes.drop(columns=["branch", "resolution", "groupby"]).to_string(index=False))
    print(f"\n{len(changes)} clusters changed ({against} -> {child['version']}).")


def cmd_release(args) -> None:
    annotation = load_annotation(args.set)
    if annotation["status"] == "released":
        sys.exit(f"{annotation['version']} is already released.")
    mapping = load_mapping(annotation)
    if set(mapping["annotation_status"]) != {"reviewed"}:
        sys.exit("Release needs every mapping row to have annotation_status 'reviewed'.")
    changes = annotation["dir"] / "CHANGES.md"
    if not changes.is_file() or "TODO" in changes.read_text():
        sys.exit("Write CHANGES.md (no TODO left) before releasing.")
    if annotation["parent"]:
        parent = load_annotation(sibling_dir(annotation, annotation["parent"]))
        for level in parent["levels"]:
            if level not in annotation["levels"]:
                continue
            old = parent["palette"]["tiers"][level]["colors"]
            new = annotation["palette"]["tiers"][level]["colors"]
            moved = {label: (old[label], new[label]) for label in set(old) & set(new) if old[label] != new[label]}
            if moved:
                print(f"note: {level} colours changed from {parent['version']}: {moved}")
    _print_summary(annotation, mapping)
    meta = dict(annotation["meta"])
    meta.update(status="released", released=date.today().isoformat(), sha256=annotation["sha256"])
    (annotation["dir"] / "annotation.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"Released {annotation['version']}.")


def cmd_new(args) -> None:
    parent = load_annotation(args.parent)
    if parent["status"] != "released":
        sys.exit(f"Release {parent['version']} before deriving a new version from it.")
    existing = [int(VERSION_PATTERN.match(p.name).group(1)) for p in parent["dir"].parent.iterdir()
                if p.is_dir() and VERSION_PATTERN.match(p.name)]
    version = f"v{max(existing) + 1}"
    target = parent["dir"].parent / version
    target.mkdir()
    for name in HASHED_FILES:
        shutil.copy2(parent["dir"] / name, target / name)
    meta = {"version": version, "cohort": parent["meta"]["cohort"], "parent": parent["version"],
            "status": "draft", "created": date.today().isoformat()}
    (target / "annotation.json").write_text(json.dumps(meta, indent=2) + "\n")
    (target / "CHANGES.md").write_text(
        f"# Annotation {version}\n\nParent: {parent['version']}\n\n"
        "## What this version represents\n\nTODO\n\n## Changes and reasons\n\nTODO\n"
    )
    print(f"Created draft {target.relative_to(REPO_ROOT)} from {parent['version']}. "
          "Edit mapping.csv/palette.json and CHANGES.md, then run script 08 and release.")


def main() -> None:
    parser = argparse.ArgumentParser(description="Manage versioned cell-type annotation sets.")
    sub = parser.add_subparsers(dest="command", required=True)
    for name, func, helptext in (("check", cmd_check, "validate a set (read-only)"),
                                 ("release", cmd_release, "lock a draft set")):
        p = sub.add_parser(name, help=helptext)
        p.add_argument("set", help="Annotation set folder")
        p.set_defaults(func=func)
    p = sub.add_parser("diff", help="list clusters whose labels changed (read-only)")
    p.add_argument("set", help="Annotation set folder")
    p.add_argument("--against", help="Version to compare with (default: the set's parent)")
    p.set_defaults(func=cmd_diff)
    p = sub.add_parser("new", help="create the next version as a draft copy of a released set")
    p.add_argument("parent", help="Released annotation set folder")
    p.set_defaults(func=cmd_new)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
