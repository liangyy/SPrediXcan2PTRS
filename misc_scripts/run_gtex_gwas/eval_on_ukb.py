import pandas as pd
import numpy as np
import h5py
from tqdm import tqdm

class PredExpr:
    def __init__(self, fn):
        self.fn = fn
        with h5py.File(self.fn, 'r') as f:
            tmp = f['samples'][:].astype(str)
            self.samples = pd.DataFrame({
                'eid': tmp,
                'idx': [ i for i in range(len(tmp)) ]
            })
            tmp = f['genes'][:].astype(str)
            self.genes = pd.DataFrame({
                'gene': tmp,
                'idx': [ i for i in range(len(tmp)) ]
            })
    @staticmethod
    def _get_range(n, chunksize=500):
        tmp = list(np.arange(0, n, chunksize))
        if tmp[-1] != n:
            tmp = list(tmp) + [ n ]
        return tmp[:-1].copy(), tmp[1:].copy()
    def mul_weights(self, df_weight, samples, max_n=None, chunksize=1000):
        df_sample_sub = pd.merge(
            self.samples, 
            pd.DataFrame({'eid': samples}), 
            on='eid'
        )
        if max_n is not None and max_n < df_sample_sub.shape[0]:
            df_sample_sub = df_sample_sub.iloc[:max_n, :].reset_index(drop=True)
        df_weight_sub = pd.merge(
            self.genes[['gene']], 
            df_weight, 
            on='gene', how='left'
        )
        df_weight_sub.fillna(0, inplace=True)
        weight_mat = df_weight_sub.drop(columns=['gene']).values
        header = list(df_weight_sub.drop(columns=['gene']).columns)
        sample_idx = df_sample_sub.idx.values
        starts, ends = self._get_range(sample_idx.shape[0], chunksize=chunksize)
        o = []
        f = h5py.File(self.fn, 'r')
        for s, e in tqdm(zip(starts, ends), total=len(starts)):
            breakpoint()
            mat = f['pred_expr'][:, sample_idx[s:e]]
            mat = mat.T @ weight_mat
            o.append(mat)
        f.close()
        o = np.concatenate(o, axis=0)
        o = pd.DataFrame(o, columns=header)
        o = pd.concat([
            pd.DataFrame({'eid': df_sample_sub.eid}),
            o
        ], axis=1)
        return o

def pxcan2weight(df_spxcan, pval_cutoffs, weight_col='effect_size'):
    pp = np.sort(np.array(pval_cutoffs))
    oo = None
    cols = [ 'gene' ]
    for p in pp[::-1]:
        sub = df_spxcan[ df_spxcan.pvalue <= p ][['gene', weight_col]].copy()
        if oo is None:
            oo = sub
        else:
            oo = pd.merge(oo, sub, on='gene', how='left')
        cols.append(f'pxcan_{p}')
    oo.fillna(0, inplace=True)
    oo.columns = cols
    return oo, pp[::-1]
def ptrs2weight(fn):
    with h5py.File(fn, 'r') as f:
        weights = f['dataset_0']['betahat'][:]
        lams = f['dataset_0']['lambda_seq'][:]
        genes = f['genes'][:].astype(str)
    df = pd.DataFrame(weights, columns=[ f'lambda_{l}' for l in lams ])
    df = pd.concat([pd.DataFrame({'gene': genes}), df], axis=1)
    return df, lams
# def cor_mat_vec(mat, vec):
#     return np.corrcoef(mat.T, vec[:, np.newaxis].T)[-1, :-1]

PXCAN_PVAL_CUTOFFS = np.concatenate([
    10 ** np.arange(-30, -10, 2).astype(float),
    10 ** np.arange(-10, -2, 0.2).astype(float),
    10 ** np.arange(-2, 0, 0.02).astype(float),
    np.array([1])
], axis=0)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='eval_on_ukb.py', description='''
        Evaluate on UKB data.
    ''')
    parser.add_argument('--pred_expr', help='''
        UKB predicted expression HDF5.
    ''')
    parser.add_argument('--sample_list', help='''
        List of samples. Will randomly choose among these.
    ''')
    parser.add_argument('--n_samples', type=int, default=None, help='''
        Number of samples to use. If None then use all.
    ''')
    parser.add_argument('--chunksize', type=int, default=1000, help='''
        Number of samples per chunk.
    ''')
    parser.add_argument('--list_of_ptrs', nargs='+', help='''
        [Tag-name]:[file-name]
    ''')
    parser.add_argument('--list_of_pxcan', nargs='+', help='''
        [Tag-name]:[file-name]
    ''')
    # parser.add_argument('--seed', type=int, default=1, help='''
    #     Numpy random seed.
    # ''')
    parser.add_argument('--output_parquet', help='''
        Output parquet table for PTRS.
    ''')
    args = parser.parse_args()
    
    import logging, time, sys, os
    # configing util
    logging.basicConfig(
        level = logging.INFO, 
        stream = sys.stderr, 
        format = '%(asctime)s  %(message)s',
        datefmt = '%Y-%m-%d %I:%M:%S %p'
    )
    
    # logging.info(f'Numpy random seed = {args.seed}.')
    # np.random.seed(args.seed)
    logging.info('Loading samples.')
    ss = []
    with open(args.sample_list, 'r') as f:
        for i in f:
            i = i.strip()
            ss.append(i)
    samples = np.array(ss)
    
    logging.info('Initializing predicted expression.')
    pred_expr = PredExpr(args.pred_expr)
    
    df_weight = None
    logging.info('Loading PrediXcan naive scores.')
    pval_dict = {}
    for kk in args.list_of_pxcan:
        tag, fn = kk.split(':')
        logging.info(f'-> Loading {tag}.')
        tmp = pd.read_csv(fn)
        tmp, pval_dict[tag] = pxcan2weight(tmp, PXCAN_PVAL_CUTOFFS)
        tmp.columns = tmp.columns[:1].tolist() + [ f'SP_x_{tag}_x_{i}' for i in tmp.columns[1:] ]
        if df_weight is None:
            df_weight = tmp
        else:
            df_weight = pd.merge(df_weight, tmp, how='outer', on='gene')
        df_weight.fillna(0, inplace=True)
    
    logging.info('Loading PTRS scores.')
    lam_dict = {}
    for kk in args.list_of_ptrs:
        tag, fn = kk.split(':')
        logging.info(f'-> Loading {tag}.')
        tmp, lam_dict[tag] = ptrs2weight(fn)
        tmp.columns = tmp.columns[:1].tolist() + [ f'PT_x_{tag}_x_{i}' for i in tmp.columns[1:] ]
        if df_weight is None:
            df_weight = tmp
        else:
            df_weight = pd.merge(df_weight, tmp, how='outer', on='gene')
        df_weight.fillna(0, inplace=True)
    
    logging.info(f'There are {df_weight.shape[1]} scores to work with.')
    
    logging.info('Calculating PTRSs.')
    out = pred_expr.mul_weights(df_weight, samples, chunksize=args.chunksize, max_n=args.n_samples)
    
    logging.info('Reformatting results.')
    osamples = out.eid.values
    omat = out.iloc[:, 1:].values
    tags = []
    values = []
    models = []
    for cc in out.columns[1:]:
        model, tag, val = cc.split('_x_')
        if model == 'SP':
            model = 'naive'
        elif model == 'PT':
            model = 'en'
        models.append(model)
        tags.append(tag)
        values.append(float(val))
    df_out = pd.DataFrame({'model': models, 'tag': tags, 'param': values})
    df_w = pd.DataFrame(omat.T, columns=osamples)
    df_out = pd.concat([df_out, df_w], axis=1)
    df_out.to_parquet(args.output_parquet, index=False)
        
    logging.info('Done.')
    
    
