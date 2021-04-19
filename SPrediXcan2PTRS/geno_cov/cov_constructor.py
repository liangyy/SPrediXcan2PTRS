from collections import OrderedDict
import re

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, save_npz, load_npz

from transethnic_prs.util.genotype_io import snpinfo_to_snpid


CHRNUM_WILDCARD = '{chr_num}'

class CovConstructor:
    def __init__(self, data, nbatch=10):
        '''
        Important: self.data is ALWAYS centered
        '''
        self.data = data - data.mean(axis=0)
        self.ncol = self.data.shape[1]
        self.nrow = self.data.shape[0]
        self.nbatch = nbatch
        self._set_batch_anchors()
    def _set_batch_anchors(self):
        ncol = self.ncol
        batch_size = ncol // self.nbatch
        if batch_size * self.nbatch < ncol:
            batch_size += 1
        if batch_size < 10 and self.nbatch > 5:
            raise ValueError('Too many batches. Exit.')
        self.batch_anchors = []
        for i in range(self.nbatch):
            start = i * batch_size
            end = min((i + 1) * batch_size, ncol)
            self.batch_anchors.append([start, end])
    def _flatten_2d_mat(self, mat):
        row_index_mat = np.tile(np.arange(mat.shape[0]), reps=(mat.shape[1], 1)).T
        row = row_index_mat.flatten()
        del row_index_mat
        col_index_mat = np.tile(np.arange(mat.shape[1]), reps=(mat.shape[0], 1))
        col = col_index_mat.flatten()
        del col_index_mat
        return row, col, mat.flatten()
    def compute_to_h5(self, fn, dtype='f'):
        import h5py
        f = h5py.File(fn, 'w')
        dset = f.create_dataset('cov', (self.ncol, self.ncol), dtype=dtype)
        for i, (s1, e1) in enumerate(self.batch_anchors):
            for j, (s2, e2) in enumerate(self.batch_anchors):
                if i > j:
                    continue
                dset[s1 : e1, s2 : e2] = self._compute_cov(
                    s1, e1, s2, e2, 
                    flatten=False, triu=(i == j)
                )
                tmp = self._compute_cov(
                    s1, e1, s2, e2, 
                    flatten=False, triu=(i == j)
                )
        f.close()
    def compute_to_disk(self, mode, output_prefix, param=None):
        if mode == 'naive':
            fn = output_prefix + '.naive.h5'
            self.compute_to_h5(fn, dtype=param)
        elif mode == 'cap':
            fn = output_prefix + '.cap.npz'
            self.compute_to_cap_npz(fn, threshold=param)
        elif mode == 'banded':
            fn = output_prefix + '.banded.npz'
            self.compute_to_banded_npz(fn, band_size=param)
    def _compute_cov(self, s1, e1, s2, e2, flatten=True, triu=True):
        '''
        Given submatrix index: 
            matrix1 = [:, s1 : e1], matrix2 = [:, s2 : e2]
        Return: 
            Pairwise covariance between column in matrix1 and column in matrix2.
            Elements are returned in row, col, val format (flatten = True). 
            And only row <= col ones are returned.
            But if flatten = False, triu could be set to False
            to return the full matrix.
        Formula:
            covariance = col1 * col2 / self.nrow (col is centered in __init__) 
        '''
        tmp = np.einsum('ni,nj->ij', self.data[:, s1 : e1], self.data[:, s2 : e2]) / (self.nrow - 1)
        if flatten is False:
            if triu is True:
                return np.triu(tmp)
            else:
                return tmp
        row, col, val = self._flatten_2d_mat(tmp)
        row += s1
        col += s2
        to_keep = row <= col
        row, col, val = row[to_keep], col[to_keep], val[to_keep]
        del tmp
        return row, col, val
    def compute_to_banded_npz(self, fn, band_size=100):
        row_all, col_all, value_all = [], [], []
        for i, (s1, e1) in enumerate(self.batch_anchors):
            for j, (s2, e2) in enumerate(self.batch_anchors):
                if i > j:
                    continue
                if s2 > e1 + band_size - 1:
                    continue
                row, col, value = self._compute_cov(s1, e1, s2, e2)
                to_keep = col - row <= band_size 
                row, col, value = row[to_keep], col[to_keep], value[to_keep]
                row_all.append(row)
                col_all.append(col)
                value_all.append(value)
        row_all = np.concatenate(row_all, axis=0)
        col_all = np.concatenate(col_all, axis=0)
        value_all = np.concatenate(value_all, axis=0)
        cov_coo = coo_matrix(
            (value_all, (row_all, col_all)), 
            shape=(self.ncol, self.ncol)
        )
        save_npz(fn, cov_coo)  
    def compute_to_cap_npz(self, fn, threshold=1e-5):
        row_all, col_all, value_all = [], [], []
        for i, (s1, e1) in enumerate(self.batch_anchors):
            for j, (s2, e2) in enumerate(self.batch_anchors):
                if i > j:
                    continue
                row, col, value = self._compute_cov(s1, e1, s2, e2)
                to_keep = np.abs(value) > threshold
                row, col, value = row[to_keep], col[to_keep], value[to_keep]
                row_all.append(row)
                col_all.append(col)
                value_all.append(value)
        row_all = np.concatenate(row_all, axis=0)
        col_all = np.concatenate(col_all, axis=0)
        value_all = np.concatenate(value_all, axis=0)
        cov_coo = coo_matrix(
            (value_all, (row_all, col_all)), 
            shape=(self.ncol, self.ncol)
        )
        save_npz(fn, cov_coo)    

class GenoCov:
    def __init__(self, fn, chromosomes=None):
        self._set_chromosomes(fn, chromosomes)
        self._load_cov_meta(fn)
        self._load_cov_mat(fn)
    def _set_chromosomes(self, fn, chrs):
        if CHRNUM_WILDCARD in fn:
            if chrs is None:
                self.chromosomes = [ i for i in range(1, 23) ]
            else:
                self.chromosomes = chrs
        else:
            self.chromosomes = None
    @staticmethod
    def _add_snpid(df_snp):
        df_snp_more = snpinfo_to_snpid(
            df_snp.chr, df_snp.pos, df_snp.ref, df_snp.alt, 
            allow_ambi=True
        )
        if df_snp.shape[0] != df_snp_more.shape[0]:
            raise ValueError('df_snp before and after adding SNPID have different number of rows.')
        df_snp['snpid'] = df_snp_more.snpid
        df_snp['direction'] = df_snp_more.direction
        return df_snp
    def _load_snp_meta(self, snp_meta_file):
        snp_tmp = pd.read_parquet(snp_meta_file)
        snp_tmp.chr = [ int(re.sub('chr', '', i)) for i in snp_tmp.chr ]
        snp_tmp = self._add_snpid(snp_tmp)
        return snp_tmp
    @staticmethod
    def _get_snp_meta_fn(fn_cov_mat):
        fn = '.'.join(fn_cov_mat.split('.')[:-2])
        return fn + '.snp_meta.parquet'
    def _load_cov_meta(self, fn):
        if self.chromosomes is None:
            fn = self._get_snp_meta_fn(fn)
            self.snp_meta = self._load_snp_meta(fn)
        else:
            o = OrderedDict()
            fn = self._get_snp_meta_fn(fn)
            for i in self.chromosomes:
                o[i] = self._load_snp_meta(fn.format(chr_num=i))
            self.snp_meta = o
    def _load_cov_mat(self, fn):
        if self.chromosomes is None:
            self.cov_mat = CovMatrix(fn)
        else:
            o = OrderedDict()
            for i in self.chromosomes:
                o[i] = CovMatrix(fn.format(chr_num=i))
            self.cov_mat = o
    def get_snplist(self):
        if isinstance(self.snp_meta, OrderedDict):
            snp = []
            for _, v in self.snp_meta.items():
                snp.append(v)
            return pd.concat(snp, axis=0)
        else:
            return self.snp_meta
    def _rearrage_df_by_target(df, target, value_cols):
        if 'snpid' in df.columns and 'direction' in df.columns:
            pass
        else:
            df = self._add_snpid(df)
        res = pd.merge(target, df[['snpid', 'direction'] + value_cols], on='left', by='snpid').reset_index(drop=True)
        res.fillna(0, inplace=True)
        res[value_cols] = res[value_cols].values * (res.direction_x.values * res.direction_y.values)[:, np.newaxis]
        res = res.drop(columns=['direction_y']).rename(columns={'direction_x': 'direction'})
        return res
    def rearrange_df(self, df_snp, value_cols):
        '''
        re-arrange df_snp so that df_snp has SNPs in the same order of 
        self.snp_meta/self.cov_mat. 
        Match both chr:pos, the alleles, and the allele direction. 
        Fill 0 for the missing SNPs (in snp_meta but not df_snp).
        Drop SNPs if they do not show up in snp_meta.
        Return df_snp in an OrderedDict if self.chromosomes is not None.
        df_snp = pd.DataFrame({
            'chr', 'pos', 'ref', 'alt',
            'snpid' (optional), 'direction' (optional)
        })
        '''
        if self.chromosomes is not None:
            chrs = df_snp.chr.unique()
            o = OrderedDict()
            for cc in chrs:
                if cc not in self.chromosomes:
                    continue
                snp_meta = self.snp_meta[cc]
                df_snp_sub = df_snp[ df_snp.chrom == cc ].reset_index(drop=True)
                df_snp_sub = self._rearrage_df_by_target(
                    df=df_snp_sub, 
                    target=snp_meta, 
                    value_cols=value_cols
                )
                o[cc] = df_snp_sub
            return o
        else:
            df_snp = self._rearrage_df_by_target(
                df=df_snp, 
                target=self.snp_meta, 
                value_cols=value_cols
            )
            return df_snp
    def matmul_xt_cov_x(self, df, value_cols):
        '''
        Evaluate X.T @ Cov @ X, where X = df[value_cols].
        Procedure:
        1. Call rearrange_df(df, value_cols) to make df have the same list of SNPs 
        (in the same order and direction) as cov.
        2. Evaluate matmal. 
        '''
        df_reordered = self.rearrange_df(df, value_cols)
        if self.chromosomes is not None:
            res_eval = np.zeros((len(values_cols), len(values_cols)))
            for cc in self.chromosomes:
                if cc not in df_reordered:
                    continue
                else:
                    x = df_reordered[value_cols].values
                    cov_x = self.cov_mat.eval_matmul_on_left(x)
                    res_eval += x.T @ cov_x
        else:
            x = df_reordered[value_cols].values
            cov_x = self.cov_mat.eval_matmul_on_left(x)
            res_eval = x.T @ cov_x
        return res_eval
        
class CovMatrix:
    def __init__(self, fn):
        self.fn = fn
        self.mode = self._init_mode()
    def _init_mode(self):
        import pathlib
        if not pathlib.Path(self.fn).is_file():
            raise ValueError('Input file does not exist.') 
        tmp = self.fn.split('.')
        return tmp[-2]
    def eval_matmul_on_left(self, left_mat, param=None):
        '''
        Retur cov @ left_mat along with the diag 
        '''
        if self.mode in ['banded', 'cap']:
            return self._eval_matmul_on_left_npz(left_mat)
        elif self.mode == 'naive':
            return self._eval_matmul_on_left_h5(left_mat, batch_size=param)
    def _eval_matmul_on_left_npz(self, mat):
        csr = load_npz(self.fn).tocsr()
        diag_csr = csr.diagonal()
        return csr.dot(mat) + csr.transpose().dot(mat) - diag_csr[:, np.newaxis] * mat, diag_csr
    def _eval_matmul_on_left_h5(self, mat, batch_size=None):
        import h5py
        f = h5py.File(self.fn, 'r')
        nrow = f['cov'].shape[0]
        ncol = mat.shape[1]
        res = np.zeros((nrow, ncol))
        if batch_size is None:
            nbatch = 1
            batch_size = nrow
        else:
            nbatch = nrow // batch_size
            if nbatch * batch_size < nrow:
                nbatch += 1
        s, e = 0, batch_size
        diag_cov = np.zeros((nrow))
        for i in range(nbatch):
            res[s : e, :] = f['cov'][s : e, :] @ mat
            res[s : e, :] += f['cov'][:, s : e].T @ mat
            diag_cov[s : e] = f['cov'][s : e, s : e].diagonal()
        res -= diag_cov[:, np.newaxis] * mat
        f.close()
        return res, diag_cov
        
