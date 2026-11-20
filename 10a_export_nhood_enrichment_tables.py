#!/usr/bin/env python
# coding: utf-8

"""
Title: Backfill neighbourhood enrichment tables from existing script 10 results
Date: 2026-11-20
Summary: For script 10 result H5ADs written before script 10 exported
`<stem>_nhood_enrichment.csv`, read only `uns/cell_type_nhood_enrichment` and
`obs/cell_type` (h5py, no `.X`) and write the same table script 10 now writes.
Typical local use reads the results from the Phoenix mount and writes the
tables into this repository:

    python 10a_export_nhood_enrichment_tables.py \
      --config config/10_cell_type_spatial_analysis/vbct/min15_area95_selected_resolutions.json \
      --input-root ~/mnt/phoenix_scratch/Banksy_py

Inputs: script 10 cell-type config (cohort and `result_dir`); result H5ADs
under `--input-root`. Outputs: one CSV per sample under `--output-root` at
`<result_dir>/<sample>/<stem>_nhood_enrichment.csv`. Inputs are never modified;
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

from spatial_tables import (
    enrichment_from_labels, enrichment_table_path, read_enrichment_table, write_enrichment_table,
)
from workflow_config import add_root_arguments, load_cohort, resolve_path, select_samples


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog="export script 10 neighbourhood enrichment tables")
    parser.add_argument("--config", required=True, help="Script 10 cell-type config")
    parser.add_argument("--sample", action="append", help="Restrict to a sample; repeat as needed")
    parser.add_argument("--validate-only", action="store_true", help="List inputs and outputs without writing")
    parser.add_argument("--overwrite", action="store_true", help="Replace existing tables")
    add_root_arguments(parser)
    return parser.parse_args()


def main():
    """Export one enrichment table per sample from its script 10 result H5AD."""
    args = parse_args()
    cfg = json.loads(Path(args.config).resolve().read_text())
    cohort = load_cohort(cfg["cohort"])
    result_in = resolve_path(args.input_root.resolve(), cfg["result_dir"])
    result_out = resolve_path(args.output_root.resolve(), cfg["result_dir"])
    for sample in select_samples(cohort, cfg.get("samples"), args.sample):
        name = sample["sample"]
        stem = Path(sample["annotated_h5ad"]).stem
        source = result_in / name / f"{stem}_squidpy.h5ad"
        table = enrichment_table_path(result_out, name, sample["annotated_h5ad"])
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
                read_elem(handle["uns/cell_type_nhood_enrichment"]), read_elem(handle["obs/cell_type"])
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
