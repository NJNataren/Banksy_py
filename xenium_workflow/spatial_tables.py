"""
Title: Neighbourhood enrichment tables shared by script 10 and its readers
Date: 2026-11-20
Summary: One small CSV per sample holds script 10's cell-type neighbourhood
enrichment: Squidpy edge counts and permutation z-scores for every ordered pair
of present cell types, plus each type's cell count. Script 10 writes it beside
its result H5AD; `10a_export_nhood_enrichment_tables.py` backfills it from
existing result H5ADs; the cross-sample comparison helper reads it, so local
analysis no longer needs the H5ADs on the Phoenix mount.

Format: rows in category order (cell_type_a outer, cell_type_b inner), columns
cell_type_a, cell_type_b, n_cells_a, edge_count, zscore. Floats are written
with Python's shortest round-trip representation and read back with
`float_precision="round_trip"`, so values are bit-identical to the H5AD.
"""

from pathlib import Path

import numpy as np
import pandas as pd

COLUMNS = ["cell_type_a", "cell_type_b", "n_cells_a", "edge_count", "zscore"]


def enrichment_table_path(directory, sample, annotated_h5ad):
    """Table path for one sample, named like script 10's `<stem>_squidpy.h5ad`."""
    return Path(directory) / sample / f"{Path(annotated_h5ad).stem}_nhood_enrichment.csv"


def enrichment_from_labels(enrichment, labels):
    """Return z-scores, edge counts, categories, and cell counts from script 10 results.

    `enrichment` is `uns['<label column>_nhood_enrichment']` (keys `zscore`,
    `count`) and `labels` the categorical label column (e.g.
    `obs['cell_type_tier_1']`); categories keep the order script 10 used.
    """
    labels = pd.Categorical(labels)
    categories = [str(category) for category in labels.categories]
    n_cells = pd.Series(labels).value_counts().reindex(labels.categories, fill_value=0).to_numpy()
    zscore, count = np.asarray(enrichment["zscore"]), np.asarray(enrichment["count"])
    if zscore.shape != (len(categories),) * 2 or count.shape != zscore.shape:
        raise ValueError("Enrichment shape does not match the cell_type categories")
    return zscore, count, categories, n_cells


def write_enrichment_table(path, zscore, count, categories, n_cells):
    """Write one sample's enrichment matrices as a long CSV."""
    n = len(categories)
    rows = [
        {"cell_type_a": categories[i], "cell_type_b": categories[j], "n_cells_a": int(n_cells[i]),
         "edge_count": int(count[i, j]), "zscore": float(zscore[i, j])}
        for i in range(n) for j in range(n)
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=COLUMNS).to_csv(path, index=False)


def read_enrichment_table(path):
    """Read a table written by `write_enrichment_table`; return matrices in category order."""
    table = pd.read_csv(
        path, dtype={"cell_type_a": str, "cell_type_b": str, "n_cells_a": np.int64, "edge_count": np.int64},
        float_precision="round_trip",
    )
    if list(table.columns) != COLUMNS:
        raise ValueError(f"{path}: expected columns {COLUMNS}")
    categories = list(dict.fromkeys(table["cell_type_a"]))
    n = len(categories)
    expected_a = np.repeat(categories, n)
    expected_b = np.tile(categories, n)
    if (len(table) != n * n or (table["cell_type_a"].to_numpy() != expected_a).any()
            or (table["cell_type_b"].to_numpy() != expected_b).any()):
        raise ValueError(f"{path}: rows must cover every cell-type pair in category order")
    n_cells_by_row = table["n_cells_a"].to_numpy().reshape(n, n)
    if (n_cells_by_row != n_cells_by_row[:, :1]).any():
        raise ValueError(f"{path}: inconsistent n_cells_a within a cell type")
    zscore = table["zscore"].to_numpy(dtype=float).reshape(n, n)
    count = table["edge_count"].to_numpy().reshape(n, n)
    return zscore, count, categories, n_cells_by_row[:, 0]
