#!/usr/bin/env python
# coding: utf-8

"""
Title: Backfill neighbourhood enrichment tables from existing script 10 results
Date: 2026-11-20
Updated: 2026-11-22 - paths and the label column come from the annotation set
and label level in the script 10 config (xenium_workflow/annotation_sets.py).
Summary: For script 10 result H5ADs whose `<stem>_nhood_enrichment.csv` is not
available locally, read only `uns/cell_type_<tier>_nhood_enrichment` and
`obs/cell_type_<tier>` (h5py, no `.X`) and write the same table script 10 writes.
Typical local use reads the results from the Phoenix mount and writes the
tables into this repository:

    python 10a_export_nhood_enrichment_tables.py \
      --config config/10_cell_type_spatial_analysis/vbct/min15_area95_selected_resolutions.json \
      --input-root ~/mnt/phoenix_scratch/Banksy_py

Inputs: script 10 config (annotation set, label level); result H5ADs under
`--input-root`. Outputs: one CSV per sample under `--output-root` in script 10's
derived result folder, `<...>/spatial_analysis/<cohort>/<version>/<tier>/<sample>/`. Inputs are never modified;
existing tables are kept unless `--overwrite` is given.
"""

import argparse
import json
from pathlib import Path

import h5py
import numpy as np

try:
    from anndata.io import read_elem
except ImportError:
    from anndata.experimental import read_elem

from xenium_workflow.spatial_tables import (
    enrichment_from_labels, enrichment_table_path, read_enrichment_table, write_enrichment_table,
)
from xenium_workflow.annotation_sets import add_annotation_arguments, annotated_h5ad_path, label_column, resolve_step, step_dir
from xenium_workflow.cohorts import add_root_arguments, resolve_path, select_samples


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog="export script 10 neighbourhood enrichment tables")
    parser.add_argument("--config", required=True, help="Script 10 cell-type config")
    parser.add_argument("--sample", action="append", help="Restrict to a sample; repeat as needed")
    parser.add_argument("--validate-only", action="store_true", help="List inputs and outputs without writing")
    parser.add_argument("--overwrite", action="store_true", help="Replace existing tables")
    add_annotation_arguments(parser)
    add_root_arguments(parser)
    return parser.parse_args()


def main():
    """Export one enrichment table per sample from its script 10 result H5AD."""
    args = parse_args()
    cfg = json.loads(Path(args.config).resolve().read_text())
    annotation, level = resolve_step(cfg, args)
    column = label_column(level)
    result_dir = step_dir(annotation, "processed", "spatial_analysis", level)
    result_in = resolve_path(args.input_root.resolve(), result_dir)
    result_out = resolve_path(args.output_root.resolve(), result_dir)
    for sample in select_samples(annotation["cohort"], cfg.get("samples"), args.sample):
        name = sample["sample"]
        annotated = annotated_h5ad_path(annotation, sample)
        source = result_in / name / f"{annotated.stem}_squidpy.h5ad"
        table = enrichment_table_path(result_out, name, annotated)
        if not source.is_file():
            raise FileNotFoundError(f"{name}: no script 10 result at {source}")
        if args.validate_only:
            print(f"{name}: {source} -> {table}{' (exists)' if table.exists() else ''}")
            continue
        if table.exists() and not args.overwrite:
            print(f"Skipping {name}: table exists ({table}); use --overwrite to replace")
            continue
        with h5py.File(source, "r") as handle:
            parts = enrichment_from_labels(
                read_elem(handle[f"uns/{column}_nhood_enrichment"]), read_elem(handle[f"obs/{column}"])
            )
        write_enrichment_table(table, *parts)
        # Confirm the written table reproduces the stored matrices exactly.
        zscore, count, categories, n_cells = read_enrichment_table(table)
        if (categories != parts[2] or not np.array_equal(zscore, parts[0], equal_nan=True)
                or not np.array_equal(count, parts[1]) or not np.array_equal(n_cells, parts[3])):
            raise RuntimeError(f"{name}: table round trip does not match {source}")
        print(f"Wrote {table} ({len(categories)} cell types)")


if __name__ == "__main__":
    main()
