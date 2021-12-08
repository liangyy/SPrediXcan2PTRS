from collections import OrderedDict

import re

import pandas as pd
import numpy as np
import h5py

import SPrediXcan2PTRS.util.misc as mi
import SPrediXcan2PTRS.util.liftover as lo

DEFAULT_PARAMS = OrderedDict([
    ('alpha', [ 1. ]),
    ('offset', [ 0.01 ]),
    ('nlambda', 100),
    ('ratio_lambda', 100.),
    ('maxiter', 1000),
    ('tol', 1e-5),
    ('r2_cutoff', [ 0.1 ]),
    ('pval_cutoffs', [ 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.005, 0.01, 0.05, 0.1, 0.5, 1 ])
])

GWAS_COL_SEP = '='

def save_result(f, grp_name, value_dict):
    grp = f.create_group(grp_name)
    
    for k, v in value_dict.items():
        if k in ['betahat', 'lambda_seq', 'pval_cutoffs']:
            grp.create_dataset(k, data=v)
        elif k in ['alpha', 'offset', 'r2_cutoff']:
            grp.attrs[k] = v

def load_gwas_snp(args_gwas, args_gwas_cols, liftover_chain=None):
    df = mi.read_table(args_gwas)
    
    cols_dict = OrderedDict()
    expected = [ 'chromosome', 'position', 'effect_allele', 'non_effect_allele' ]
    exp_patterns = [ '{' + v + '}' for v in expected ]
    for col in args_gwas_cols:
        col_in_table, target_col_pattern = mi.try_parse_gwas_col(col, sep=GWAS_COL_SEP)
        if col_in_table not in df.columns:
            raise ValueError(f'Column {col_in_table} is not in GWAS table.')
        if '{' in target_col_pattern and '}':
            mi.parse_and_update_gwas_col(cols_dict, df[col_in_table], target_col_pattern)
        else:
            if target_col_pattern not in cols_dict:
                cols_dict[target_col_pattern] = df[col_in_table]
            else:
                raise ValueError(f'Duplicated column definition: {target_col_pattern}.')
    df = pd.DataFrame(cols_dict)
 
    if liftover_chain is not None:
        tmp = lo.liftover(df.chromosome, df.position.astype(int), liftover_chain)
        df.chromosome = tmp.liftover_chr
        df.position = tmp.liftover_pos
    
    chrs = []
    pos = []
    for cc, pp in zip(df.chromosome, df.position):
        try:
           cc_ = int(re.sub('chr', '', cc))
        except ValueError:
           cc_ = 'NA'
        chrs.append(cc_)
        pos.append(int(pp))
    df.chromosome = chrs
    df.position = pos
    
    df.rename(columns={
        'chromosome': 'chrom',
    }, inplace=True)
    
    return df
    
def maybe_scale_betahat(beta_mat, divide_factor, if_scale):
    if if_scale is True:
        return beta_mat / divide_factor[:, np.newaxis]
    else:
        return beta_mat
    
    

def load_params(fn_yaml=None):
    if fn_yaml is None:
        my_load = OrderedDict()
    else:
        my_load = mi.read_yaml(fn_yaml)
        mi.try_cast_float(
            my_load, 
            keys=['alpha', 'offset', 'ratio_lambda', 'tol', 
                'r2_cutoff', 'pval_cutoffs']
        )
        for k, v in my_load.items():
            my_load[k] = mi.check_param(k, v)
    
    other_params = OrderedDict()
    for k, v in DEFAULT_PARAMS.items():
        if k in ['alpha', 'offset']:
            continue
        else:
            if k in my_load:
                other_params[k] = my_load[k]
            else:
                other_params[k] = v
    alphas = my_load['alpha'] if 'alpha' in my_load else DEFAULT_PARAMS['alpha']
    offsets = my_load['offset'] if 'offset' in my_load else DEFAULT_PARAMS['offset']
    r2_cutoffs = my_load['r2_cutoff'] if 'r2_cutoff' in my_load else DEFAULT_PARAMS['r2_cutoff']
    
    other_keys = ['ratio_lambda', 'tol']
    pt_keys = ['pval_cutoffs']
    return alphas, offsets, r2_cutoffs, { k: other_params[k] for k in other_keys }, { k: other_params[k] for k in pt_keys }

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='run_pxcan2ptrs.py', description='''
        Convert S-PrediXcan results to PTRS weights along a lambda sequence.  
        Need to have transethnic_prs and SPrediXcan2PTRS in PYTHONPATH.
    ''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--predixcan', help='''
        Path to (S)PrediXcan results 
        (CSV format from SPrediXcan.py).
    ''')
    parser.add_argument('--predictdb', help='''
        Path to the predictDB file.
    ''')
    parser.add_argument('--geno_cov', help='''
        Pattern path to the genotype covariance. 
        This should match the predictDB file, 
        otherwise many genes and SNPs will be missed. 
        Use {chr_num} in replace of the chromosome numbering. 
        For instance, en_Whole_Blood.geno_cov.chr{chr_num}.evd.npz
    ''')
    parser.add_argument('--gwas', help='''
        The GWAS file being used for (S)PrediXcan analysis. 
        We need this since we want to use the same SNPs as the ones used in (S)PrediXcan. 
    ''')
    parser.add_argument('--gwas_cols', nargs='+', help=f'''
        Specify the columns in GWAS table. 
        Expect 4 columns: chromosome, position, effect_allele, non_effect_allele 
        For instance: 
            --gwas_cols 
            chr{GWAS_COL_SEP}chromosomechr 
            pos{GWAS_COL_SEP}position:pos 
            a1{GWAS_COL_SEP}effect_allele 
            a2{GWAS_COL_SEP}non_effect_allele 
        If want to do some parsing on column KK, do KK={{chromosome}}:{{position}}. 
        Another example: 
            --gwas_cols 
            KK{GWAS_COL_SEP}{{chromosome}}:{{position}} 
            a1{GWAS_COL_SEP}effect_allele 
            a2{GWAS_COL_SEP}non_effect_allele 
    ''')
    parser.add_argument('--liftover_chain', default=None, help='''
        If specified, we will liftover GWAS. 
        Note that we expect GWAS and predictDB use the same genome assembly version 
        since we match SNPs by position.
    ''')
    parser.add_argument('--gwas_sample_size', type=int, help='''
        GWAS sample size. 
        We use it as the sample size for running PTRS fitting.
    ''')
    parser.add_argument('--hyperparam_yaml', default=None, help='''
        A YAML file containing the hyper-parameters. 
        For instance (values listed are default): 
        alpha: [ 1. ] # any value in (0, 1] 
        nlambda: 100 # any integer > 0 
        offset: [ 0.01 ] # any value in [0, 1) 
        ratio_lambda: 100. # any value > 1 
        maxiter: 1000 # any integer > 0 
        tol: 1e-5 # any value > 0
        Note that alpha and offset can take multiple values (a list).
    ''')
    parser.add_argument('--output_prefix', help='''
        Output file name in HDF5 format.
    ''')
    parser.add_argument('--mode', default='by_chr', choices=['by_chr', 'jointly'], help='''
        DO NOT need to change usually. 
        Fitting PTRS one chromosome at a time or all jointly. 
        Typically they give very similar result (default = by_chr).
    ''')
    parser.add_argument('--clump', action='store_true', help='''
        If specified, will run P+T PTRS rather than EN PTRS.
        Set the hyperparameters in --hyperparam_yaml (values listed are default)
        r2_cutoff: [ 0.1 ]
        pval_cutoffs: [ 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.005, 0.01, 0.05, 0.1, 
            0.5, 1 ]
    ''')
    parser.add_argument('--original_scale', action='store_true', help='''
        Use original scale if specified. Otherwise, use the standardized scale.
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
    import SPrediXcan2PTRS.geno_cov.cov_constructor as cc
    import SPrediXcan2PTRS.util.db as db
    import SPrediXcan2PTRS.solver as so
    
    alphas, offsets, r2_cutoffs, other_params, pt_params = load_params(args.hyperparam_yaml)
    
    logging.info('Loading (S)PrediXcan results.')
    df_pxcan = pd.read_csv(args.predixcan)
    
    logging.info('Loading GWAS.')
    df_gwas = load_gwas_snp(args.gwas, args.gwas_cols, args.liftover_chain)
    
    logging.info('Loading genotype covariances.')
    geno_cov = cc.GenoCov(args.geno_cov)
    
    logging.info('Loading predictDB.')
    func_get_snp_meta = lambda x: db.get_snp_meta_from_varID(x)
    weight = db.WeightDB(args.predictdb, func_get_snp_meta)
    
    logging.info('Initializing the solver.')
    solver = so.Solver(
        df_pxcan=df_pxcan[['gene', 'zscore']].copy(),
        sample_size=args.gwas_sample_size,
        weight_db=weight,
        geno_cov=geno_cov, 
        df_gwas_snp=df_gwas[['chrom', 'position', 'effect_allele', 'non_effect_allele']].copy(),
        lazy=True  
    )
    
    logging.info('Building gene covariances.')
    solver.init_w_genes(show_progress_bar=True)
    genes = solver.gene_meta.copy()
    
    logging.info(
        'Calculating PTRS weights: nalpha = {}, noffset = {}.'.format(
            len(alphas), len(offsets)
        )
    )
    solver.init_model1blk()
    sd_gene = np.sqrt(np.concatenate(solver.var_gene, axis=0))
    output_handle = h5py.File(args.output_prefix + '.results.h5', 'w')
    result_id = 0
    if args.clump is False:
        for ia, alpha in enumerate(alphas):
            for io, offset in enumerate(offsets):
                logging.info(f'-> Working on {ia + 1}th alpha = {alpha}, {io + 1}th offset = {offset}.')
                if args.mode == 'by_chr':
                    betahat, lambda_seq, _, _, conv = solver.fit_ptrs_by_blk(
                        alpha=alpha, offset=offset, **other_params
                    )
                    conv = np.stack(conv)
                    # only lambdas where all chrs converged are kept
                    conv = conv.mean(axis=0)
                elif args.mode == 'jointly':
                    betahat, lambda_seq, _, _, conv = solver.fit_ptrs(
                        alpha=alpha, offset=offset, **other_params
                    )
                to_keep_idxs = list(np.where(conv[1:] == 1)[0] + 1)
                to_keep_idxs = np.array([ 0 ] + to_keep_idxs)
                betahat = betahat[:, to_keep_idxs ]
                betahat = maybe_scale_betahat(betahat, sd_gene, args.original_scale)
                lambda_seq = lambda_seq[ to_keep_idxs ]
                logging.info('-> {} lambda values converged. Saving results.'.format(lambda_seq.shape[0]))
                save_result(
                    output_handle, f'dataset_{result_id}',
                    value_dict=OrderedDict([
                        ('alpha', alpha),
                        ('offset', offset),
                        ('betahat', betahat),
                        ('lambda_seq', lambda_seq)
                    ])
                )
                result_id += 1
    else:
        for ir, r2_cutoff in enumerate(r2_cutoffs):
            betahat = solver.fit_clump_ptrs(r2_cutoff=r2_cutoff, **pt_params)
            betahat = maybe_scale_betahat(betahat, sd_gene, args.original_scale)
            save_result(
                output_handle, f'dataset_{result_id}',
                value_dict=OrderedDict([
                    ('r2_cutoff', r2_cutoff),
                    ('betahat', betahat),
                    ('pval_cutoffs', pt_params['pval_cutoffs'])
                ])
            )
    
    logging.info('Saving other meta information.')
    output_handle.create_dataset('genes', data=np.concatenate(solver.genes).astype('S'))
    solver.gene_meta.to_csv(
        args.output_prefix + '.gene_meta.tsv.gz', 
        index=False, compression='gzip', sep='\t'
    )
    
    # one more clean up
    output_handle.close()
    
    logging.info('Done.')    
            
    
