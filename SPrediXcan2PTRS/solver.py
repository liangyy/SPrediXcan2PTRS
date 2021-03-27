from collections import OrderedDict

import pandas as pd
import numpy as np
from tqdm import tqdm

from transethnic_prs.util.misc import intersect_two_lists
from transethnic_prs.util.genotype_io import snpinfo_to_snpid
from transethnic_prs.util.math_jax import mean_center_col_2d_jax
from transethnic_prs.util.math import diag_mul_mat, mat_mul_diag

class Solver:
    def __init__(self, df_pxcan, weight_db, geno_bed, gene_list=None, lazy=False):
        '''
        Inputs:
            1. df_pxcan: pd.DataFrame({'gene': gene, 'zscore': zscore})
            2. weight_db: WeightDB object 
                          (see SPrediXcan2PTRS/util/db.py)
                          # from prediction model weights in PredictDB sqlite3 format
            3. geno_bed: PlinkBedIo object 
                         (see transethnic_prs/util/genotype_io.py)
                         # from genotype PLINK BED
        Outputs:
            self.R = predicted gene correlation matrix
            self.z = PrediXcan z-score
        Formula:
            # do it one chromosome at a time
            Cov(i, j) = w(i)' Cov(X) w(j)
            R = sqrt(diag(Cov))^-1 Cov sqrt(diag(Cov))^-1
        Caution: 
            Only use the non-ambiguious SNPs
        '''
        self.geno_bed = geno_bed 
        self.weight_db = weight_db
        self.snp_cols = ['chrom', 'position', 'effect_allele', 'non_effect_allele']
        genes = self._get_common_genes(
            list(df_pxcan.gene), 
            gene_list=gene_list
        )
        db_df_gene_var = self._get_variants_per_gene(genes)
        self.db_df_gene_var, self.gene_meta = self._check_in_ref_geno(
            db_df_gene_var
        )
        self.gene_w_vars = list(self.gene_meta[ self.gene_meta.nsnp_in_ref > 0 ].gene)
        if lazy is True:
            return 
        self._build_cor(genes=self.gene_w_vars)
    def build_cor(self, genes=None):
        self._build_cor(genes)
    def _get_common_genes(self, pxcan_gene, gene_list=None):
        db_gene = self.weight_db.get_gene_info()
        genes = intersect_two_lists(pxcan_gene, db_gene)
        if gene_list is not None:
            genes = intersect_two_lists(genes, gene_list)
        return genes
    def _get_variants_per_gene(self, genes, disable_progress_bar=False):
        return self.weight_db.get_variants_of_gene(genes)
        # gene_dict = OrderedDict()
        # for g in tqdm(genes, disable=disable_progress_bar):
        #     gene_dict[g] = self.weight_db.get_variants_of_gene(g)
        # return gene_dict
    def _check_in_ref_geno(self, df_gene_var):
        out_meta = []
        
        snp_all = set(self.geno_bed.get_snplist().snpid)
        
        # annotate df_gene_var with snpid and direction (remove ambiguious snps)
        # and intersect with snp_all (snps from genotype)
        df_var = df_gene_var[self.snp_cols].drop_duplicates()
        vars_ = snpinfo_to_snpid(
            df_var.chrom, df_var.position, 
            df_var.effect_allele, df_var.non_effect_allele, 
            return_complete=True
        )
        df_var = df_var.iloc[ vars_.idx, : ].reset_index(drop=True)
        df_var['snpid'] = vars_.snpid
        df_var['direction'] = vars_.direction
        df_var = df_var[ df_var.snpid.isin(snp_all) ].reset_index(drop=True)
        df_gene_var_new = pd.merge(df_gene_var, df_var, on=self.snp_cols)
        df_meta = df_gene_var.value_counts(subset=['gene']).reset_index()
        df_meta2 = df_gene_var_new.value_counts(subset=['gene']).reset_index()
        df_meta = pd.merge(df_meta, df_meta2, on='gene', how='left').fillna(0)
        df_meta.rename(columns={'0_x': 'nsnp_in_db', '0_y': 'nsnp_in_ref'}, inplace=True)
        df_gene_var_new.chrom = df_gene_var_new.chrom.astype(int)
        df_gene_var_new.position = df_gene_var_new.position.astype(int)
        return df_gene_var_new, df_meta
    def _build_cor(self, genes=None):
        if genes is None:
            genes = list(self.gene_w_vars)
        else:
            genes = intersect_two_lists(list(self.gene_w_vars), genes)
        df_gene_var = self.db_df_gene_var.copy()
        df_gene_var = df_gene_var[ df_gene_var.gene.isin(genes) ].reset_index(drop=True)

        self.R = []
        self.genes = []
        self.var_gene = []
        for i in range(21, 23):
            df_gv = df_gene_var[ df_gene_var.chrom == i ].reset_index(drop=True)
            if df_gv.shape[0] == 0:
                continue
            R, var_gene, genes = self._build_cor_within_chr(df_gv)
            self.R.append(R)
            self.var_gene.append(var_gene)
            self.genes.append(genes)
    def _build_cor_within_chr(self, df_gene_var):
        df_var = df_gene_var[self.snp_cols + ['snpid', 'direction']].drop_duplicates().sort_values(by=self.snp_cols + ['snpid', 'direction'])
        df_gene = df_gene_var[['gene']].drop_duplicates().sort_values(by='gene')
        df_gene.reset_index(drop=True, inplace=True)
        weight_mat = np.zeros((df_var.shape[0], df_gene.shape[0]))
        for idx, g in enumerate(df_gene.gene):
            tmp = df_gene_var[ df_gene_var.gene == g ]
            tmp = pd.merge(df_var, tmp, on='snpid', how='left').weight.values
            tmp[ np.isnan(tmp) ] = 0
            tmp = tmp * df_var.direction.values
            weight_mat[ :, idx ] = tmp
                    
        geno = self.geno_bed.load(df_var.snpid)
        geno = mean_center_col_2d_jax(geno)
        GxW = geno @ weight_mat
        cov_pe = GxW.T @ GxW
        var_pe = cov_pe.diagonal()
        sqrt_var_pe = np.sqrt(var_pe)
        R = diag_mul_mat(1 / sqrt_var_pe, mat_mul_diag(cov_pe, 1 / sqrt_var_pe))
        return R, var_pe, list(df_gene.gene)
        