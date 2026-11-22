# Task: Fix Moran's I Significance in Script 10

## Goal

Make script 10's Moran's I output give a usable, reproducible measure of which
genes are spatially patterned in each sample. At present neither FDR-corrected
p-value column is usable, and the "significant genes" CSV and top-Moran gene
plots that depend on it are artefacts. Decide on a method, implement it in
`10_cell_type_spatial_analysis.py`, and rerun script 10.

Not urgent: no findings in `spatial_analysis.md` use Moran's I. Found on
2026-11-20 while validating the post-08 cleanup HPC test (see
`post08_flow_cleanup.md`).

Related context: `spatial_analysis.md` (analysis status),
`config/10_cell_type_spatial_analysis/vbct/README.md` (says the Moran settings
"have not been upgraded for final inference").

## Current Implementation

`10_cell_type_spatial_analysis.py`, `run_squidpy_for_resolution` (Moran block):

1. Subsample 25% of cells (`moran_subsample_fraction: 0.25`,
   `random_state=1234`).
2. Delaunay spatial graph on the subsample.
3. `sq.gr.spatial_autocorr(mode="moran", n_perms=10, n_jobs=1, seed=1234)` over
   all 304 panel genes. Squidpy returns `I`, `pval_norm` (analytic, normality
   assumption), `pval_z_sim`, `pval_sim` (permutation), and BH-corrected
   versions `*_fdr_bh`.
4. Writes `moran_scores_all_*.csv` (all genes) and `moran_scores_fdr0.1_*.csv`,
   which keeps genes with **`pval_sim_fdr_bh <= moran_fdr_threshold` (0.1)**.
5. Plots the top `top_moran_n` (5) genes from that filtered table, sorted by
   `pval_sim`.

Config: `config/10_cell_type_spatial_analysis/vbct/min15_area95_selected_resolutions.json`
(`moran_subsample_fraction`, `moran_n_perms`, `moran_n_jobs`,
`moran_fdr_threshold`, `top_moran_n`).

## What Was Observed (2026-11-20)

Per-sample Moran outputs (Phoenix: new layout under
`data/xenium/processed/vbct/spatial_analysis/min15_area95_selected_resolutions/<sample>/<branch>_r<res>/`,
old layout under `figures/spatial_analysis/vbct_min15_area95_selected_resolutions/`).
"new" = seeded rerun on 2026-11-20 (job 16757390 for CK_skin_res, job 16758519
in progress for the rest); "old" = unseeded October runs.

| Sample | Source | Genes at pval_sim floor (1/11) | pval_sim_fdr_bh <= 0.1 | min pval_sim_fdr_bh | pval_norm_fdr_bh |
|---|---|---|---|---|---|
| BE_brain_non_res | new | 275/304 | 0 | 0.1005 | valid |
| CK_bowel_res | old | 292/304 | 292 | 0.0946 | valid |
| CK_skin_res | new | 225/304 | 0 | 0.1228 | **all NaN** (LY6D constant) |
| EL_skin_res | old | 303/304 | 303 | 0.0912 | valid |
| GR_lung_non_res_roi | new | 262/304 | 0 | 0.1055 | valid |
| HW_brain_res | old | 283/304 | 283 | 0.0977 | valid |
| MF_skin_non_res | old | 291/304 | 291 | 0.0950 | valid |
| MG_gastric_non_res | old | 284/304 | 284 | 0.0973 | valid |

BH on `pval_norm` after dropping NaN genes gives, at FDR <= 0.1: BE 268/304,
CK_bowel 288, CK_skin 219/303, EL 303, GR 258, HW 280, MF 289, MG 282.

HPC vs local reruns of CK_skin_res with identical package versions (scanpy
1.9.5, anndata 0.10.2, numpy 1.24.4, squidpy 1.6.5): `I` and `pval_norm` agree
to ~1e-14; `pval_sim` differs by up to 0.36. Two local runs with the same thread
count are byte-identical; local runs with 8 vs 16 numba threads differ.

## What I Think Is Happening

Three separate problems, in order of importance.

1. **The permutation FDR filter is a threshold artefact.** With 10
   permutations the smallest possible `pval_sim` is 1/11 = 0.0909, and most
   genes sit exactly at that floor. For tied p-values, BH gives every floored
   gene the adjusted value `(1/11) × 304 / k`, where `k` is the number of genes
   at the floor. That is <= 0.1 only when k >= 277. So each sample flips between
   "almost every gene significant" (k >= 277) and "no gene significant"
   (k < 277): BE misses by 0.0005 with k = 275. The filtered CSV therefore
   reflects how many genes hit the floor, not which genes are patterned.
   (This corrects an earlier statement in the 2026-11-20 session that no gene
   could pass FDR 0.1 with 10 permutations; that was true for CK_skin_res only.)
2. **The top-Moran plots are arbitrary or missing.** They are chosen from the
   filtered table sorted by `pval_sim`. Where hundreds of genes pass, they are
   all tied at 1/11 and `sort_values` (default quicksort, not stable) picks five
   essentially arbitrarily. Where none pass (CK_skin, BE, GR), no plots are made.
3. **`pval_norm_fdr_bh` can be entirely NaN.** If one gene is constant in the
   25% subsample (CK_skin_res: LY6D), its `pval_norm` is NaN and Squidpy's BH
   correction (`statsmodels.multipletests` on the whole vector) returns NaN for
   every gene. The same happens to `pval_z_sim_fdr_bh`. This is Squidpy
   behaviour, not ours, but our output inherits it.

Secondary: `pval_sim` is not reproducible across machines or thread counts.
The permutation RNG is seeded and identical; the differences come from
floating-point differences in the computed Moran's I (~1e-14, likely numba and
`fastmath` summation order) flipping near-ties in Squidpy's count of permuted
values >= observed. With only 10 permutations each flip moves the p-value by
0.09. This matters less once (1) is fixed, but it means permutation p-values
should not be expected to match between HPC and local runs.

Underlying point: with hundreds to hundreds of thousands of cells, Moran's I
p-values call nearly every gene with any spatial structure significant
(200–300 of 304 by `pval_norm`). P-values discriminate poorly here; the size of
`I` is the more useful quantity for ranking genes.

Not yet verified: whether the NaN problem occurs in samples other than
CK_skin_res after the seeded reruns (the subsample changed with the seed), and
whether the near-tie explanation fully accounts for the HPC/local `pval_sim`
differences (it fits the evidence; the exact numba code path was not traced).

## Options

- **(a) Analytic p-values with a correct BH correction (recommended).** Drop
  genes with NaN `pval_norm` (constant in the subsample; record them), apply BH
  to the rest ourselves, and filter on that. Deterministic, machine-independent,
  no floor. Caveat: relies on a normality approximation; reasonable at these
  cell numbers but worth a sentence in methods.
- **(b) More permutations.** `n_perms = 1000` gives a floor of ~0.001 and makes
  ties rare, but the Moran step costs roughly 100 times longer, which is
  significant for MF_skin_non_res and CK_bowel_res (already hours per task),
  and results stay machine-dependent at the margins.
- **(c) Treat Moran's I as descriptive (recommended alongside (a)).** Rank
  genes by `I` and choose the top-Moran plots by `I` (deterministic), and
  present p-values only as a screen, not as evidence.
- Related: consider whether a 25% subsample is still needed for the larger
  samples' runtime, and dropping constant genes before the call rather than
  after.

Recommendation: (a) + (c). Keep `n_perms` small or set it to `None` to skip
permutations entirely if `pval_sim` is no longer used (saves runtime).

## Next Steps

1. Decide the method (scientific decision for Nathalie): (a) + (c), (b), or
   another.
2. Implement in `10_cell_type_spatial_analysis.py`:
   - compute BH on `pval_norm` after dropping NaN genes (new column, e.g.
     `pval_norm_fdr_bh_valid`), log the dropped genes;
   - make the filtered CSV and `moran_fdr_threshold` use that column (or the
     chosen one), and name the column in the output file or a README;
   - pick top-Moran genes by `I` among genes passing the filter, with a stable
     sort and a tie-breaker (gene name);
   - optionally set `n_perms` to `None` in the config if permutations are
     dropped.
3. Validate locally on CK_skin_res (6,101 cells; known LY6D NaN case): filtered
   gene count matches the manual check (219/303 at 0.1 with the current
   subsample), top genes are the highest-I significant genes, two runs are
   byte-identical, and a run with different `NUMBA_NUM_THREADS` gives the same
   filtered list.
4. Rerun script 10 on the HPC for all eight samples (Nathalie submits:
   `sbatch hpc/slurm/run_10_cell_type_spatial_analysis.sl`), then locally
   `10a_export_nhood_enrichment_tables.py --overwrite`. Neighbourhood
   enrichment results are unaffected by this change.
5. Update `config/10_cell_type_spatial_analysis/vbct/README.md` (Moran
   section), the concise README if behaviour changes, and `spatial_analysis.md`.

## Pitfalls

- Do not read anything into the current `moran_scores_fdr0.1_*.csv` files or
  the top-Moran gene plots for any sample; see "What I Think Is Happening".
- `pval_sim` will not match between machines; do not use cross-machine
  byte-identity of `moran_scores_all_*.csv` as a validation check. Compare `I`
  and `pval_norm` (tolerance ~1e-12) instead.
- The job 16758519 rerun (2026-11-20) produces outputs with the same issues;
  a further script 10 rerun is needed after the fix.
- Squidpy 1.6.5 is pinned in both environments; check `spatial_autocorr`
  behaviour (NaN handling, `seed`) again if it is upgraded.

## Session History

- **2026-11-20:** Issue found during the post-08 cleanup HPC test. Seed added
  to `spatial_autocorr` the same day (commit `c269407`): reruns on one machine
  are now reproducible, but the problems above remain. Per-sample evidence
  table compiled from stored and rerun outputs on the Phoenix mount
  (read-only).
