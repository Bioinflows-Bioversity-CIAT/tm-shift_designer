import pandas as pd
import glob as glob
import os

qtl_map = pd.read_csv('./config/tg_varians_FFARV3.txt', sep='\t')
qtl_map.rename(columns={'variant_id': 'qtn_name'}, inplace=True)
primer_paths = glob.glob('./results/*/summary/*_primers.csv')
combs_paths = glob.glob('./results/*/summary/*_combinations.csv')

primer_dfs = [pd.read_csv(file) for file in primer_paths]
primer_db = pd.concat(primer_dfs, ignore_index=True)

combs_dfs = [pd.read_csv(file) for file in combs_paths]
combs_db = pd.concat(combs_dfs, ignore_index=True)

combs_db = combs_db.merge(qtl_map, on='qtn_name')
combs_db.to_csv('./results/Combinations_DB.csv')

primer_db.to_csv('./results/primer_DB.csv')
