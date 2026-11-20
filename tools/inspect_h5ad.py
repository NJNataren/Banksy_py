"""
Title: Read-only AnnData (.h5ad) inspector
Created: 2026-11-20
Last modified: 2026-11-20
Summary: Describe an .h5ad file without loading it into memory: shape, matrix
type and value range, layers, obsm/uns keys, and grouped .obs columns. Gives a
labelled verdict on whether `.X` holds clean log-normalised expression, raw
counts, or BANKSY-expanded neighbour features, so marker summaries never use
the wrong matrix by accident.
Inputs: One .h5ad file (local or on a read-only mount).
Outputs: Text report on stdout, or JSON with --json. Writes nothing.
Usage: python tools/inspect_h5ad.py path/to/file.h5ad [--obs-values cell_type] [--json]
"""

import argparse
import json
import re
import sys
from collections import OrderedDict
from pathlib import Path

import h5py
import numpy as np

# How many stored values (sparse) or rows (dense) to sample when judging `.X`.
# Large enough to be representative, small enough to stay fast over a mount.
SPARSE_SAMPLE_VALUES = 1_000_000
DENSE_SAMPLE_ROWS = 2_000

# Clean log1p(normalised) Xenium values sit well below this; raw counts and
# BANKSY-scaled features can exceed it.
CLEAN_MAX_VALUE = 20.0

# Filename patterns that identify the pipeline step which wrote an object.
# Order matters: the first match wins, so more specific patterns come first.
KIND_PATTERNS = [
    ("_squidpy.h5ad", "script 10 spatial results object"),
    ("_cell_type_annotated_", "script 08 cell-type annotated clean object"),
    ("_cells_used_for_banksy_with_clusters_", "script 06 clean recluster object (cells used by BANKSY)"),
    ("_with_banksy_reclusters_", "script 06 full provenance object (all cells)"),
    ("adata_spatial_", "BANKSY spatial object (neighbour-expanded features)"),
    ("_qc_annotated", "script 01/05 QC-annotated clean object"),
    ("_with_banksy_clusters_", "script 00 clean clustered object"),
    ("_raw", "raw Xenium object"),
]

# BANKSY label columns look like labels_<decay>_pc<N>_nc<L>_r<res><suffix>.
LABEL_RE = re.compile(r"^(labels_.+?_r)(\d+\.\d+)(.*)$")

# Keywords for script 01 QC metric columns (script 05 decisions start with qc_).
QC_METRIC_HINTS = (
    "threshold", "passed", "negative_control", "neg_control", "neg_probe", "area_percentile",
    "pct_counts", "n_genes_by_counts", "log1p_", "transcripts_per_um2",
    "small_area", "_cat",
)


def decode(value):
    """Return `value` as str if it is bytes, otherwise unchanged."""
    return value.decode() if isinstance(value, bytes) else value


def describe_kind(path):
    """Guess which pipeline step wrote the file, from its name only."""
    name = Path(path).name
    for pattern, kind in KIND_PATTERNS:
        if pattern in name:
            return kind
    return "unknown (filename matches no known pipeline pattern)"


def matrix_info(node):
    """Summarise a dense or sparse AnnData matrix stored in HDF5.

    Args:
        node: h5py Dataset (dense) or Group (csr/csc sparse) for `.X` or a layer.

    Returns:
        Dict with shape, storage, dtype and sampled value statistics.
    """
    info = OrderedDict()
    if isinstance(node, h5py.Group):
        shape = tuple(int(x) for x in node.attrs["shape"])
        data = node["data"]
        nnz = int(data.shape[0])
        sample = np.asarray(data[: min(nnz, SPARSE_SAMPLE_VALUES)])
        info["storage"] = decode(node.attrs.get("encoding-type", "sparse"))
        info["dtype"] = str(data.dtype)
        # Sparse matrices store only non-zeros, so the zero fraction comes
        # from the stored-value count rather than from the sample.
        zero_fraction = 1 - nnz / max(shape[0] * shape[1], 1)
    else:
        shape = tuple(int(x) for x in node.shape)
        sample = np.asarray(node[: min(shape[0], DENSE_SAMPLE_ROWS)])
        info["storage"] = "dense"
        info["dtype"] = str(node.dtype)
        zero_fraction = float((sample == 0).mean()) if sample.size else 0.0
    sample = sample.ravel()
    info["shape"] = shape
    info["sampled_values"] = int(sample.size)
    if sample.size:
        info["min"] = float(sample.min())
        info["max"] = float(sample.max())
        info["negative_fraction"] = float((sample < 0).mean())
        info["integer_valued"] = bool(np.all(np.mod(sample, 1) == 0))
    info["zero_fraction"] = float(zero_fraction)
    return info


def read_column_names(group):
    """Return a dataframe group's column order (obs or var)."""
    return [decode(c) for c in group.attrs.get("column-order", [])]


def read_column(group, name):
    """Read one obs/var column as a NumPy array of Python values.

    Categorical columns are expanded from codes to category labels; missing
    values (code -1) become None.
    """
    node = group[name]
    if isinstance(node, h5py.Group) and "codes" in node:
        categories = np.array([decode(c) for c in node["categories"][:]], dtype=object)
        codes = node["codes"][:]
        values = np.empty(codes.shape, dtype=object)
        valid = codes >= 0
        values[valid] = categories[codes[valid]]
        return values
    if isinstance(node, h5py.Group):
        # Nullable integer/boolean arrays store values and a mask separately.
        values = node["values"][:].astype(object)
        values[node["mask"][:]] = None
        return values
    values = node[:]
    if values.dtype.kind in ("S", "O"):
        values = np.array([decode(v) for v in values], dtype=object)
    return values


def group_obs_columns(columns):
    """Sort obs columns into workflow groups and collapse resolution series.

    Returns:
        OrderedDict mapping group name to a list of display strings.
    """
    groups = OrderedDict((key, []) for key in (
        "cell type annotation (script 08)",
        "recluster labels (script 06)",
        "original BANKSY labels (script 00)",
        "manual annotation labels",
        "QC decisions (script 05)",
        "QC metrics (script 01)",
        "Xenium / other",
    ))
    series = OrderedDict()
    for col in columns:
        match = LABEL_RE.match(col)
        if col.startswith("cell_type"):
            groups["cell type annotation (script 08)"].append(col)
        elif match and match.group(3).endswith("_ann"):
            groups["manual annotation labels"].append(col)
        elif match:
            # Collapse r0.50 ... r1.50 into one line per prefix/suffix pair.
            series.setdefault((match.group(1), match.group(3)), []).append(match.group(2))
        elif col.startswith("qc_"):
            groups["QC decisions (script 05)"].append(col)
        elif any(hint in col for hint in QC_METRIC_HINTS) or col.startswith(("min_", "max_", "cell_area_")):
            groups["QC metrics (script 01)"].append(col)
        else:
            groups["Xenium / other"].append(col)
    for (prefix, suffix), resolutions in series.items():
        key = "recluster labels (script 06)" if "_recluster_" in suffix else "original BANKSY labels (script 00)"
        res = resolutions[0] if len(resolutions) == 1 else f"{{{resolutions[0]}..{resolutions[-1]}}}"
        groups[key].append(f"{prefix}{res}{suffix}  ({len(resolutions)} resolution(s))")
    return OrderedDict((k, v) for k, v in groups.items() if v)


def judge_matrix(x_info, var_columns, banksy_base_genes, panel_size):
    """Classify `.X` as clean expression, raw counts, BANKSY features or unknown.

    This is a heuristic. The verdict is always returned with the evidence used,
    so a reader can check the reasoning rather than trust a label.
    """
    evidence = []
    n_vars = x_info["shape"][1]
    banksy = False
    if "is_nbr" in var_columns:
        banksy = True
        evidence.append(f"var has 'is_nbr' (BANKSY neighbour features); {banksy_base_genes} own-cell genes")
    if x_info.get("negative_fraction", 0) > 0:
        banksy = True
        evidence.append(f"{x_info['negative_fraction']:.0%} of sampled values are negative (scaled features)")
    if panel_size and n_vars != panel_size and n_vars % panel_size == 0:
        banksy = True
        evidence.append(f"n_vars {n_vars} = {n_vars // panel_size} x panel size {panel_size}")
    if banksy:
        return "BANKSY-EXPANDED: do not use for marker expression", evidence
    if panel_size:
        evidence.append(f"n_vars {n_vars} {'matches' if n_vars == panel_size else 'differs from'} panel size {panel_size}")
    if "max" not in x_info:
        return "UNKNOWN: matrix is empty", evidence
    if x_info["integer_valued"]:
        evidence.append(f"all sampled values are whole numbers (max {x_info['max']:.0f})")
        return "RAW COUNTS: not normalised", evidence
    evidence.append(f"non-negative, non-integer, max {x_info['max']:.2f}")
    if x_info["max"] <= CLEAN_MAX_VALUE:
        return "CLEAN EXPRESSION (log-normalised range)", evidence
    return "UNKNOWN: non-negative but above the log1p range", evidence


def read_uns_scalars(node):
    """Return the scalar fields of an uns group (e.g. uns['cell_type_annotation'])."""
    out = OrderedDict()
    for key, child in node.items():
        if isinstance(child, h5py.Dataset) and child.shape == ():
            out[key] = decode(child[()])
        elif isinstance(child, h5py.Dataset) and child.size <= 20:
            out[key] = [decode(v) for v in child[:].tolist()]
    return out


def summarise_values(values, top=30):
    """Value counts for categorical-like data, or a numeric summary."""
    values = np.asarray(values, dtype=object)
    missing = sum(v is None for v in values)
    present = np.array([v for v in values if v is not None], dtype=object)
    numeric = present.size and all(isinstance(v, (int, float, np.number)) and not isinstance(v, (bool, np.bool_))
                                   for v in present[:1000])
    if numeric and len(set(present[:5000].tolist())) > top:
        arr = present.astype(float)
        return {"type": "numeric", "n": int(arr.size), "missing": missing,
                "min": float(arr.min()), "median": float(np.median(arr)), "max": float(arr.max())}
    labels, counts = np.unique(present.astype(str), return_counts=True)
    order = np.argsort(-counts)
    return {"type": "counts", "n": int(present.size), "missing": missing, "n_unique": int(labels.size),
            "counts": OrderedDict((labels[i], int(counts[i])) for i in order[:top])}


def inspect(path, obs_value_columns, panel_size):
    """Collect the full report for one file as a dict."""
    report = OrderedDict()
    report["file"] = str(path)
    report["size_gb"] = round(path.stat().st_size / 1e9, 3)
    report["kind_from_filename"] = describe_kind(path)
    with h5py.File(path, "r") as f:
        if isinstance(f["obs"], h5py.Dataset):
            raise SystemExit("Legacy (pre-0.7 AnnData) layout: open it with anndata instead.")
        var_columns = read_column_names(f["var"])
        banksy_base_genes = None
        if "is_nbr" in var_columns:
            banksy_base_genes = int((~f["var"]["is_nbr"][:].astype(bool)).sum())
        report["X"] = matrix_info(f["X"])
        verdict, evidence = judge_matrix(report["X"], var_columns, banksy_base_genes, panel_size)
        report["verdict"] = verdict
        report["evidence"] = evidence
        report["layers"] = OrderedDict()
        for name in f.get("layers", {}):
            layer = matrix_info(f["layers"][name])
            report["layers"][name] = {k: layer[k] for k in ("storage", "dtype", "max", "integer_valued") if k in layer}
        report["raw_present"] = "raw" in f
        report["obsm"] = list(f.get("obsm", {}).keys())
        uns_keys = list(f.get("uns", {}).keys())
        report["uns_n_keys"] = len(uns_keys)
        report["uns_notable"] = [k for k in uns_keys if not k.endswith(("_colors", "_markers_wilcoxon"))
                                 and not k.startswith("dendrogram_")]
        if "cell_type_annotation" in uns_keys and isinstance(f["uns"]["cell_type_annotation"], h5py.Group):
            report["cell_type_annotation"] = read_uns_scalars(f["uns"]["cell_type_annotation"])
        obs_columns = read_column_names(f["obs"])
        report["obs_n_columns"] = len(obs_columns)
        report["obs_groups"] = group_obs_columns(obs_columns)
        report["var_columns"] = var_columns
        report["obs_values"] = OrderedDict()
        for col in obs_value_columns:
            if col not in obs_columns:
                report["obs_values"][col] = {"error": "column not found"}
                continue
            report["obs_values"][col] = summarise_values(read_column(f["obs"], col))
    return report


def print_report(r):
    """Print the report as readable text."""
    x = r["X"]
    print(f"File     {r['file']}")
    print(f"Size     {r['size_gb']} GB    Shape  {x['shape'][0]:,} cells x {x['shape'][1]:,} features")
    print(f"Kind     {r['kind_from_filename']}  (from filename)")
    print()
    rng = f"min {x['min']:.3g}, max {x['max']:.3g}" if "min" in x else "empty"
    print(f".X       {x['storage']} {x['dtype']}, {rng}, {x['zero_fraction']:.0%} zeros "
          f"(sampled {x['sampled_values']:,} values)")
    print(f"Verdict  {r['verdict']}   [heuristic]")
    for item in r["evidence"]:
        print(f"           - {item}")
    print()
    layers = ", ".join(f"{k} ({v['storage']}, {'integer' if v.get('integer_valued') else 'non-integer'})"
                       for k, v in r["layers"].items()) or "(none)"
    print(f"layers   {layers}")
    print(f"raw      {'present' if r['raw_present'] else '(none)'}")
    print(f"obsm     {', '.join(r['obsm']) or '(none)'}")
    print(f"uns      {r['uns_n_keys']} keys; notable: {', '.join(r['uns_notable']) or '(none)'}")
    if "cell_type_annotation" in r:
        for k, v in r["cell_type_annotation"].items():
            print(f"           cell_type_annotation.{k} = {v}")
    print(f"var cols {', '.join(r['var_columns'])}")
    print(f"obs      {r['obs_n_columns']} columns")
    for group, cols in r["obs_groups"].items():
        print(f"  {group}:")
        for col in cols:
            print(f"      {col}")
    for col, summary in r["obs_values"].items():
        print(f"\nobs['{col}']")
        if "error" in summary:
            print(f"  {summary['error']}")
        elif summary["type"] == "numeric":
            print(f"  numeric n={summary['n']:,} missing={summary['missing']:,} "
                  f"min={summary['min']:.4g} median={summary['median']:.4g} max={summary['max']:.4g}")
        else:
            print(f"  n={summary['n']:,} missing={summary['missing']:,} unique={summary['n_unique']}")
            for label, count in summary["counts"].items():
                print(f"  {count:>10,}  {label}")


def parse_args():
    """Parse command-line options."""
    parser = argparse.ArgumentParser(
        description="Describe an .h5ad file read-only, with a verdict on what .X contains.")
    parser.add_argument("h5ad", type=Path, help="Path to the .h5ad file")
    parser.add_argument("--obs-values", action="append", default=[], metavar="COLUMN",
                        help="Print value counts (or a numeric summary) for an obs column; repeatable")
    parser.add_argument("--panel-size", type=int, default=None,
                        help="Expected gene-panel size (e.g. 304 for VBCT); adds an n_vars check to the verdict")
    parser.add_argument("--json", action="store_true", help="Print the report as JSON")
    return parser.parse_args()


def main():
    """Run the inspector."""
    args = parse_args()
    if not args.h5ad.is_file():
        sys.exit(f"File not found: {args.h5ad}")
    report = inspect(args.h5ad, args.obs_values, args.panel_size)
    if args.json:
        print(json.dumps(report, indent=2, default=str))
    else:
        print_report(report)


if __name__ == "__main__":
    main()
