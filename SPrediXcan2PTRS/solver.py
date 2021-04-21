from collections import OrderedDict

import pandas as pd
import numpy as np
from tqdm import tqdm

import transethnic_prs.model1.Model1Blk as m1b
from transethnic_prs.util.misc import intersect_two_lists
from transethnic_prs.util.genotype_io import snpinfo_to_snpid
from transethnic_prs.util.math import diag_mul_mat, mat_mul_diag

class Solver:
    def __init__(self, df_pxcan, sample_size, weight_db, geno_cov, df_gwas_snp=None, gene_list=None, lazy=False, show_progress_bar=False):
        '''
        Inputs:
            1. df_pxcan: pd.DataFrame({'gene': gene, 'zscore': zscore})
            2. weight_db: WeightDB object 
                          (see SPrediXcan2PTRS/util/db.py)
                          # from prediction model weights in PredictDB sqlite3 format
            3. geno_cov: GenoCov object 
                             (see SPrediXcan2PTRS/geno_cov/cov_constructor.py)
            4. df_gwas_snp: pd.DataFrame({'chrom', 'position', 'effect_allele', 'non_effect_allele'}) showing the SNPs being included in GWAS.
        Outputs:
            self.R = predicted gene correlation matrix
            self.z = PrediXcan z-score
        Formula:
            # do it one chromosome at a time
            # for only SNPs that occur in GWAS.
            Cov(i, j) = w(i)' Cov(X) w(j)
            R = sqrt(diag(Cov))^-1 Cov sqrt(diag(Cov))^-1
        
        '''
        self.geno_cov = geno_cov 
        self.weight_db = weight_db
        self.df_pxcan = df_pxcan
        self.sample_size = sample_size
        self.snp_cols = ['chrom', 'position', 'effect_allele', 'non_effect_allele']
        genes = self._get_common_genes(
            list(df_pxcan.gene), 
            gene_list=gene_list
        )
        db_df_gene_var = self._get_variants_per_gene(genes, df_gwas_snp)
        self.db_df_gene_var, self.gene_meta = self._check_in_ref_geno_cov(
            db_df_gene_var
        )
        self.gene_w_vars = list(self.gene_meta[ self.gene_meta.nsnp_in_db_n_gwas_n_geno_cov > 0 ].gene)
        if lazy is True:
            return 
        self._build_cor(genes=self.gene_w_vars, show_progress_bar=show_progress_bar)
        self._init_z()
    def init_w_genes(self, genes=None, show_progress_bar=False):
        self._build_cor(genes, show_progress_bar=show_progress_bar)
        self._init_z()
    def init_model1blk(self):
        '''
        Initiate an transethnic_prs.model1.Model1Blk instance
        Set self.Xlist as 1 x p zero matrix and y as length-1 zero vector. 
        '''
        self.model1_blk = m1b.Model1Blk(
            Alist=self.R,
            blist=[ z / np.sqrt(self.sample_size) for z in self.z_predixcan ],
            Xlist=[ np.zeros((2, z.shape[0])) for z in self.z_predixcan ],
            y=np.zeros((2, ))
        )
    def _init_z(self):
        zlist = []
        for gg in self.genes:
            zlist.append(
                pd.merge(
                    pd.DataFrame({'gene': gg}),
                    self.df_pxcan[['gene', 'zscore']], on='gene'
                ).zscore.values
            )
        self.z_predixcan = zlist
    def _get_common_genes(self, pxcan_gene, gene_list=None):
        db_gene = self.weight_db.get_gene_info()
        genes = intersect_two_lists(pxcan_gene, db_gene)
        if gene_list is not None:
            genes = intersect_two_lists(genes, gene_list)
        return genes
    def _get_variants_per_gene(self, genes, df_gwas_snp=None, disable_progress_bar=False):
        variants_per_gene = self.weight_db.get_variants_of_gene(genes)
        if df_gwas_snp is not None:
            gwas_snp = snpinfo_to_snpid(
                df_gwas_snp.chrom, 
                df_gwas_snp.position, 
                df_gwas_snp.effect_allele, df_gwas_snp.non_effect_allele,
                allow_ambi=True
            )
            db_snp = snpinfo_to_snpid(
                variants_per_gene.chrom,
                variants_per_gene.position, 
                variants_per_gene.effect_allele, variants_per_gene.non_effect_allele,
                allow_ambi=True
            )
            snp_both = intersect_two_lists(gwas_snp.snpid, db_snp.snpid)
            db_snp_in_both_idx = db_snp[ db_snp.snpid.isin(snp_both) ].idx
            variants_per_gene = variants_per_gene.iloc[ db_snp_in_both_idx, : ].reset_index(drop=True)
        return variants_per_gene
    def _check_in_ref_geno_cov(self, df_gene_var):
        out_meta = []
        
        snp_all = set(self.geno_cov.get_snplist().snpid)
        
        # annotate df_gene_var with snpid and direction
        # and intersect with snp_all (snps from geno cov)
        df_var = df_gene_var[self.snp_cols].drop_duplicates()
        vars_ = snpinfo_to_snpid(
            df_var.chrom, df_var.position, 
            df_var.effect_allele, df_var.non_effect_allele, 
            return_complete=True,
            allow_ambi=True
        )
        df_var = df_var.iloc[ vars_.idx, : ].reset_index(drop=True)
        df_var['snpid'] = vars_.snpid
        df_var['direction'] = vars_.direction
        df_var = df_var[ df_var.snpid.isin(snp_all) ].reset_index(drop=True)
        df_gene_var_new = pd.merge(df_gene_var, df_var, on=self.snp_cols)
        df_meta = df_gene_var.value_counts(subset=['gene']).reset_index()
        df_meta2 = df_gene_var_new.value_counts(subset=['gene']).reset_index()
        df_meta = pd.merge(df_meta, df_meta2, on='gene', how='left').fillna(0)
        df_meta.rename(columns={'0_x': 'nsnp_in_db_n_gwas', '0_y': 'nsnp_in_db_n_gwas_n_geno_cov'}, inplace=True)
        df_meta0 = self.weight_db.get_nsnp_per_gene()
        df_meta = pd.merge(df_meta0, df_meta, on='gene', how='left').fillna(0)
        df_meta.rename(columns={'0': 'nsnp_in_db'})
        df_gene_var_new.chrom = df_gene_var_new.chrom.astype(int)
        df_gene_var_new.position = df_gene_var_new.position.astype(int)
        return df_gene_var_new, df_meta
    def _build_cor(self, genes=None, show_progress_bar=False):
        if genes is None:
            genes = list(self.gene_w_vars)
        else:
            genes = intersect_two_lists(list(self.gene_w_vars), genes)
        df_gene_var = self.db_df_gene_var.copy()
        df_gene_var = df_gene_var[ df_gene_var.gene.isin(genes) ].reset_index(drop=True)

        self.R = []
        self.genes = []
        self.var_gene = []
        for i in tqdm(range(1, 23), disable=not show_progress_bar):
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
        df_weight = pd.DataFrame(weight_mat)
        weight_value_cols = list(df_gene.gene)
        df_weight.columns = weight_value_cols
        df_weight = pd.concat([df_var.reset_index(drop=True), df_weight], axis=1)
        df_weight.rename(
            columns={
                'chrom': 'chr', 'position': 'pos', 
                'effect_allele': 'alt', 'non_effect_allele': 'ref'
            }, 
            inplace=True
        )
        
        cov_pe = self.geno_cov.matmul_xt_cov_x(df_weight, weight_value_cols)
        var_pe = cov_pe.diagonal()
        sqrt_var_pe = np.sqrt(var_pe)
        R = np.array(diag_mul_mat(1 / sqrt_var_pe, mat_mul_diag(cov_pe, 1 / sqrt_var_pe)))
        return R, var_pe, weight_value_cols
        
    def fit_ptrs(self, alpha=0.5, offset=0, nlambda=100, ratio_lambda=100, tol=1e-5, maxiter=1000):
        # rescale offset to eff_offset 
        # so that we fit A + eff_offset * I 
        # instead of (1 - offset) A + offset * I
        eff_offset = offset / (1 - offset)
        beta, lambda_, niter, tol, conv = self.model1_blk.solve_path(
            alpha=alpha, offset=eff_offset, tol=tol, maxiter=maxiter, nlambda=nlambda, ratio_lambda=ratio_lambda
        )
        # need to rescale lambda_
        return beta, lambda_ * (1 - offset), niter, tol, conv
    def fit_ptrs_by_blk(self, alpha=0.5, offset=0, nlambda=100, ratio_lambda=100, tol=1e-5, maxiter=1000):
        # rescale offset to eff_offset 
        # so that we fit A + eff_offset * I 
        # instead of (1 - offset) A + offset * I
        eff_offset = offset / (1 - offset)
        beta, lambda_, niter, tol, conv = self.model1_blk.solve_path_by_blk(
            alpha=alpha, offset=eff_offset, tol=tol, maxiter=maxiter, nlambda=nlambda, ratio_lambda=ratio_lambda
        )
        # need to rescale lambda_
        return beta, lambda_ * (1 - offset), niter, tol, conv
    
