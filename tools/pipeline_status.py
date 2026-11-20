"""
Title: Pipeline output status derived from configs
Created: 2026-11-20
Last modified: 2026-11-20
Summary: For one QC/reclustering branch, work out from the step configs which
outputs each sample should have (scripts 01, 05, 06, 08, 10) and report which
exist, which are missing, and which are older than their upstream input. This
replaces trusting dated run notes in handoffs with a check of the filesystem.
Filename rules mirror scripts 05, 06 and 10 and the cohort config read by
scripts 08 and 10; update them here if those scripts change their output names.
Inputs: Script 06 config directory (recursive), matching script 05 configs,
cohort configs under config/cohorts/, script 10 configs that name the cohort;
a data root to check against (repository root by default, or a mount).
Outputs: Status table on stdout (or JSON with --json). Writes nothing.
Usage: python tools/pipeline_status.py --config-dir config/06_recluster_qc_annotated/vbct/grouped_area95_min_transcript_gte15
"""

import argparse
import json
import os
import sys
from collections import OrderedDict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))
from workflow_config import load_cohort, load_json, resolve_path  # noqa: E402

# Steps in pipeline order, with the upstream step whose output each one reads.
# The upstream link drives the "stale" check (output older than its input).
STEPS = OrderedDict([
    ("01-qc", None),
    ("05-annot", "01-qc"),
    ("06-full", "05-annot"),
    ("06-clean", "05-annot"),
    ("06-csv", "05-annot"),
    ("08-annot", "06-clean"),
    ("10-h5ad", "08-annot"),
    ("10-nhood", "08-annot"),
])


def relpath(path):
    """Path as a repository-relative string where possible (for matching configs)."""
    path = Path(path)
    try:
        return str(path.resolve().relative_to(REPO_ROOT)) if path.is_absolute() else os.path.normpath(str(path))
    except ValueError:
        return str(path)


def iter_configs(directory):
    """Yield (path, dict) for every JSON object config under `directory`, sorted."""
    for path in sorted(Path(directory).rglob("*.json")):
        try:
            cfg = load_json(path)
        except (json.JSONDecodeError, UnicodeDecodeError) as err:
            print(f"warning: cannot parse {path}: {err}", file=sys.stderr)
            continue
        if isinstance(cfg, dict):
            yield path, cfg


def script05_lookup():
    """Index script 05 configs by (project, dataset_name, output_label)."""
    index = {}
    for path, cfg in iter_configs(REPO_ROOT / "config" / "05_apply_qc_filters"):
        if "dataset_name" not in cfg:
            continue
        key = (cfg.get("project"), cfg["dataset_name"], cfg.get("output_label", "filtered_qc_v1"))
        index.setdefault(key, []).append((path, cfg))
    return index


def expected_from_script06(cfg, s05_index):
    """Expected output paths for one sample, following scripts 05 and 06.

    Mirrors the defaults and config overrides in
    05_apply_qc_filters_for_reclustering.py and 06_recluster_qc_annotated_with_banksy.py.
    Paths are repository-relative, so they can be checked against any data root.
    """
    project, ds = cfg["project"], cfg["dataset_name"]
    input_label = cfg.get("input_label", cfg.get("output_label", "filtered_qc_v1"))
    run_label = cfg.get("run_label", f"{input_label}_{cfg.get('recluster_inclusion', 'qc_pass_only')}")
    res = cfg["res_label"] if isinstance(cfg["res_label"], list) else [cfg["res_label"]]
    res_str = "_".join(res)
    processed = Path(cfg.get("base_dir", "data/xenium")) / "processed" / project / ds
    paths = OrderedDict()
    notes = []

    # Script 01 output is the script 05 input; only the matching script 05
    # config knows its `input_qc_label`, so link through output_label.
    matches = s05_index.get((project, ds, input_label), [])
    if matches:
        s05_path, s05 = matches[0]
        if len(matches) > 1:
            notes.append(f"{len(matches)} script 05 configs give output_label {input_label}; used {relpath(s05_path)}")
        s05_processed = Path(s05.get("base_dir", "data/xenium")) / "processed" / project / ds
        name = f"adata_expression_clean_{ds}_qc_annotated"
        if s05.get("input_qc_label"):
            name += f"_{s05['input_qc_label']}"
        paths["01-qc"] = s05.get("input_h5ad", s05_processed / f"{name}.h5ad")
    else:
        paths["01-qc"] = None
        notes.append(f"no script 05 config with output_label {input_label}")

    paths["05-annot"] = cfg.get("input_h5ad", processed / f"adata_expression_clean_{ds}_qc_annotated_{input_label}.h5ad")
    paths["06-full"] = processed / (f"adata_expression_clean_{ds}_qc_annotated_{input_label}"
                                    f"_with_banksy_reclusters_{run_label}_{res_str}.h5ad")
    paths["06-clean"] = processed / (f"adata_expression_clean_{ds}_recluster_{run_label}"
                                     f"_cells_used_for_banksy_with_clusters_{res_str}.h5ad")
    paths["06-csv"] = processed / f"{ds}_recluster_{run_label}_cell_cluster_id_across_clustering_res_{res_str}.csv"
    return paths, notes


def find_cohorts():
    """Load every cohort config under config/cohorts/."""
    cohorts = []
    for path, _ in iter_configs(REPO_ROOT / "config" / "cohorts"):
        try:
            cohorts.append(load_cohort(path))
        except (KeyError, ValueError) as err:
            print(f"warning: skipping cohort {relpath(path)}: {err}", file=sys.stderr)
    return cohorts


def script10_configs(cohort_path):
    """Script 10 configs whose `cohort` points at this cohort (they hold result_dir)."""
    found = []
    for directory in sorted((REPO_ROOT / "config").glob("10_*")):
        for path, cfg in iter_configs(directory):
            if "result_dir" in cfg and "cohort" in cfg and \
                    resolve_path(REPO_ROOT, cfg["cohort"]).resolve() == Path(cohort_path).resolve():
                found.append((path, cfg))
    return found


def add_post08(rows, cohort):
    """Add script 08 and 10 expected outputs for samples in `cohort`.

    Script 10 writes <result_dir>/<sample>/<annotated stem>_squidpy.h5ad and
    <result_dir>/<sample>/<annotated stem>_nhood_enrichment.csv
    (see 10_cell_type_spatial_analysis.py and spatial_tables.py).
    """
    s10 = script10_configs(cohort["path"])
    for record in cohort["samples"]:
        row = rows.setdefault(record["sample"], {"paths": OrderedDict(), "notes": []})
        if "06-clean" in row["paths"] and "recluster_h5ad" in record and \
                relpath(row["paths"]["06-clean"]) != relpath(record["recluster_h5ad"]):
            row["notes"].append("cohort recluster_h5ad differs from the script 06 config output")
        row["paths"].setdefault("06-clean", record.get("recluster_h5ad"))
        row["paths"]["08-annot"] = record["annotated_h5ad"]
        if s10:
            cfg = s10[0][1]
            stem = Path(record["annotated_h5ad"]).stem
            base = Path(cfg["result_dir"]) / record["sample"]
            row["paths"]["10-h5ad"] = base / f"{stem}_squidpy.h5ad"
            row["paths"]["10-nhood"] = base / f"{stem}_nhood_enrichment.csv"
        row["resolution"] = record.get("resolution")
    return s10


def check(rows, root):
    """Fill each row with ok / MISSING / stale states checked under `root`."""
    for row in rows.values():
        row["state"] = OrderedDict()
        mtimes = {}
        for step, path in row["paths"].items():
            if path is None:
                row["state"][step] = "?"
                continue
            full = resolve_path(root, path)
            if full.exists():
                mtimes[step] = full.stat().st_mtime
                row["state"][step] = "ok"
            else:
                row["state"][step] = "MISSING"
        # Mark outputs older than the file they were built from. A rerun of an
        # upstream step after a downstream one means the downstream file may
        # not reflect current inputs.
        for step, upstream in STEPS.items():
            if step in mtimes and upstream in mtimes and mtimes[step] < mtimes[upstream]:
                row["state"][step] = "stale"


def print_table(rows, root, cohort, s10, show_missing):
    """Print the status table and per-step summary."""
    columns = [step for step in STEPS if any(step in row["paths"] for row in rows.values())]
    width = max([len("sample")] + [len(name) for name in rows]) + 2
    print(f"Root    {root}")
    print(f"Cohort  {relpath(cohort['path']) if cohort else '(none matched; steps 08/10 not checked)'}")
    if cohort:
        print(f"Step 10 {relpath(s10[0][0]) if s10 else '(no script 10 config names this cohort)'}")
    print()
    print("sample".ljust(width) + "".join(c.ljust(10) for c in columns) + "res")
    for name, row in rows.items():
        cells = "".join(row["state"].get(c, "-").ljust(10) for c in columns)
        print(name.ljust(width) + cells + str(row.get("resolution") or ""))
    print()
    summary = []
    for c in columns:
        states = [row["state"].get(c) for row in rows.values() if c in row["state"]]
        summary.append(f"{c} {states.count('ok')}/{len(states)}")
    print("Summary: " + ", ".join(summary))
    print("States: ok | MISSING | stale = older than its upstream input | ? = cannot derive path | - = not applicable")
    for name, row in rows.items():
        for note in row["notes"]:
            print(f"note: {name}: {note}")
        if show_missing:
            for step, state in row["state"].items():
                if state in ("MISSING", "stale"):
                    print(f"{state.lower()}: {name} {step}: {row['paths'][step]}")


def parse_args():
    """Parse command-line options."""
    parser = argparse.ArgumentParser(
        description="Report which pipeline outputs exist for a branch, derived from its configs (read-only).")
    parser.add_argument("--config-dir", type=Path,
                        help="Script 06 config directory for the branch (searched recursively)")
    parser.add_argument("--cohort", type=Path,
                        help="Cohort config for steps 08/10 (default: auto-detect from script 06 outputs)")
    parser.add_argument("--root", type=Path, default=REPO_ROOT,
                        help="Root that repository-relative data paths resolve against "
                             "(default: repository root; pass a mounted HPC checkout to check there)")
    parser.add_argument("--show-missing", action="store_true", help="List full paths of missing or stale outputs")
    parser.add_argument("--json", action="store_true", help="Print results as JSON")
    args = parser.parse_args()
    if not args.config_dir and not args.cohort:
        parser.error("give --config-dir, --cohort, or both")
    return args


def main():
    """Build expected outputs from configs, check them, and report."""
    args = parse_args()
    root = args.root.expanduser()
    # An unmounted mount point is usually an empty directory; say so rather
    # than reporting every output as missing.
    if not root.is_dir() or not any(root.iterdir()):
        sys.exit(f"Root {root} is missing or empty (mount not active?). Nothing checked.")
    if not (root / "data" / "xenium").is_dir():
        print(f"warning: {root}/data/xenium does not exist; is --root the repository checkout?", file=sys.stderr)

    rows = OrderedDict()
    if args.config_dir:
        config_dir = resolve_path(REPO_ROOT, args.config_dir)
        s05_index = script05_lookup()
        for path, cfg in iter_configs(config_dir):
            if not {"project", "dataset_name", "res_label"} <= set(cfg):
                continue
            paths, notes = expected_from_script06(cfg, s05_index)
            if cfg["dataset_name"] in rows:
                notes.append(f"duplicate config {relpath(path)} ignored")
                rows[cfg["dataset_name"]]["notes"].extend(notes)
                continue
            rows[cfg["dataset_name"]] = {"paths": paths, "notes": notes, "config": relpath(path)}
        if not rows:
            sys.exit(f"No script 06 sample configs found under {config_dir}")

    cohort = None
    if args.cohort:
        cohort = load_cohort(args.cohort)
    else:
        # Auto-detect: the cohort whose recluster_h5ad paths are these script 06 outputs.
        wanted = {relpath(row["paths"]["06-clean"]) for row in rows.values()}
        for candidate in find_cohorts():
            if wanted & {relpath(s.get("recluster_h5ad", "")) for s in candidate["samples"]}:
                cohort = candidate
                break
    s10 = add_post08(rows, cohort) if cohort else []
    check(rows, root)

    if args.json:
        out = {"root": str(root), "cohort": relpath(cohort["path"]) if cohort else None,
               "samples": {name: {"state": row["state"], "notes": row["notes"],
                                  "paths": {k: (str(v) if v is not None else None) for k, v in row["paths"].items()}}
                           for name, row in rows.items()}}
        print(json.dumps(out, indent=2))
    else:
        print_table(rows, root, cohort, s10, args.show_missing)


if __name__ == "__main__":
    main()
