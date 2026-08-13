#!/bin/bash

set -euo pipefail

CONFIG="config/03_export_summary/archive/testing/legacy/local_test_CK_skin_res_0p5.json"
OUTPUT="data/xenium/processed/cross_sample_dotplot_exports/local_test_CK_skin_res_0p5_dotplot_summary.csv"

conda run -n banksy python 03_export_dotplot_data_from_config.py --config "$CONFIG"

conda run -n banksy python -c "import pandas as pd; p='$OUTPUT'; df=pd.read_csv(p); required={'sample','resolution','cluster_id','sample_cluster','groupby','gene','marker_group','mean_expression','percent_expressing','n_cells','expression_source','adata_path'}; missing=required-set(df.columns); assert not missing, f'Missing columns: {missing}'; assert len(df)>0, 'No rows exported'; assert df['gene'].nunique()>0, 'No genes exported'; assert df['cluster_id'].nunique()>0, 'No clusters exported'; print(f'PASS: {len(df)} rows, {df.gene.nunique()} genes, {df.cluster_id.nunique()} clusters, expression_source={sorted(df.expression_source.unique())}')"
