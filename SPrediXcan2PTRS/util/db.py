import re

import pandas as pd
import sqlite3 as sq

from SPrediXcan2PTRS.util.liftover import liftover

class WeightDB:
    def __init__(self, db_file, func_get_snp_meta):
        self.conn = sq.connect(db_file)
        self.func_get_snp_meta = func_get_snp_meta
        self._set_var()
        self._set_weights()
    def _set_weights(self):
        query = f'select * from weights'
        df = pd.read_sql_query(query, self.conn)
        self.df_weight = df
    def _set_var(self):
        query = f'select distinct rsid, varID, ref_allele, eff_allele from weights'
        df = pd.read_sql_query(query, self.conn)
        df2 = self.func_get_snp_meta(df)
        df2['non_effect_allele'] = df.ref_allele
        df2['effect_allele'] = df.eff_allele
        cols = ['chrom', 'position', 'effect_allele', 'non_effect_allele']
        self.df_var = df2
    def get_variants_of_gene(self, gene_list):
        '''
        Return a pd.DataFrame:
        {
            'chrom': chrom,
            'position': pos,
            'effect_allele': effect_allele,
            'non_effect_allele': non_effect_allele, 
            'weight': weight
        }
        '''
        df = self.df_weight[ self.df_weight.gene.isin(gene_list) ].reset_index(drop=True)
        df = pd.merge(df, self.df_var, on='varID')
        cols = ['gene', 'chrom', 'position', 'effect_allele', 'non_effect_allele', 'weight']
        return df[cols].copy()
    def get_gene_info(self):
        query = "select gene from extra"
        df = pd.read_sql_query(query, self.conn)
        return list(df.gene)

# def _has_infos(format, infos):
#     for i in infos:
#         i = '{' + i + '}'
#         if i not in format:
#             return False
#     return True

def get_snp_meta_from_varID(df, sep='_', liftover_chain_file=None):
    
    df0 = df.copy().reset_index(drop=True)
    
    if 'varID' not in df.columns:
        raise ValueError('Require varID when using `get_snp_meta_from_varID`')
    
    o = [ [], [], [], [] ]
    for ss in df.varID:
        chrom, position, a1, a2 = ss.split(sep)[:4]
        position = int(position)
        for i, v in zip(o, [ chrom, position, a1, a2 ]):
            i.append(v)

    df_snp = pd.DataFrame({
        'chrom': o[0], 'position': o[1], 'effect_allele': o[2], 'non_effect_allele': o[3]
    })
    if liftover_chain_file is not None:
        df_lift = liftover(
            df_snp.chrom.to_list(), df_snp.position.values, 
            chainfile=liftover_chain_file
        )  
        df_snp.chrom = df_lift.liftover_chr
        df_snp.position = df_lift.liftover_pos  
    # remove 'chr' from chrom
    df_snp.chrom = [ re.sub('chr', '', i) for i in df_snp.chrom ]
    
    if df0.shape[0] != df_snp.shape[0]:
        raise ValueError('There are different number of SNPs in df and df_snp. Something wrong!')
    
    return pd.concat([df0, df_snp], axis=1)
    