# Script 04 Dotplot Configs

`active/` contains current runnable plot configs for the script 00 clean-object workflow.

`archive/` contains older dated runs, archive-label configs, and old HPC/legacy layouts retained for provenance.

## Active Layout

```text
active/
├── vbct_clean_script00/
│   ├── multi_sample/
│   ├── per_sample_all_genes_all_res/
│   └── per_sample_canonical_markers_all_res/
├── ptmt_clean_script00/
│   └── multi_sample/
└── ptmt_pc55_clean_script00/
    ├── multi_sample/
    ├── per_sample_panel_all_res/run_1/
    └── per_sample_all_genes_all_res/
```

Use `04_plot_multi_sample_dotplot_from_config_local.py --config <config.json>` for these configs.

For PTMT PC55 run 1 per-sample panel dotplots across all available resolutions, run:

```bash
for cfg in config/04_plot_dotplot/active/ptmt_pc55_clean_script00/per_sample_panel_all_res/run_1/*.json; do
  conda run -n banksy python 04_plot_multi_sample_dotplot_from_config_local.py \
    --config "$cfg"
done
```
