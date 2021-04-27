import pandas as pd
import numpy as np
import h5py
output = 'sample_list_for_eval.txt'
df = pd.read_csv('/home/t.cri.yliang/labshare/PTRS/query_phenotypes_cleaned_up.csv')
df2 = pd.read_parquet('/home/t.cri.yliang/labshare/ukb_idp/data/imagexcan_phenotype_round_1.parquet')
df = pd.merge(df, df2, on='eid')
df = df[ df.isna().sum(axis=1) == 0 ].reset_index(drop=True)
df_pos = df[ (df.parent_AD > 0) | (df.parent_depression > 0) ]
df_neg = df[ ~ df.eid.isin(df_pos.eid) ]
df_neg = df_neg.iloc[ np.random.choice(df_neg.shape[0], 30000, replace=False), : ]
df_ = pd.concat([df_pos, df_neg], axis=0).reset_index(drop=True)
eids = df_.eid.astype(str)
with open(output, 'w') as f:
    for i in eids:
        f.write(i + '\n')
        
