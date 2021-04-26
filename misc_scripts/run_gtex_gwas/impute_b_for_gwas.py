if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='impute_b_for_gwas.py', description='''
        Impute effect size (bhat) for GWAS from z-score, allele frequency, and 
        sample_size.
    ''')
    parser.add_argument('--input', help='''
        Input GWAS table.
    ''')
    parser.add_argument('--output', help='''
        Output GWAS table.
    ''')
    parser.add_argument('--zscore', help='''
        Column name of z-score.
    ''')
    parser.add_argument('--freq', help='''
        Column name of alelle frequency.
    ''')
    parser.add_argument('--sample_size', help='''
        Column name of sample size.
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
    
    import numpy as np
    import SPrediXcan2PTRS.util.misc as mi
    
    logging.info('Loading GWAS table.')
    df = mi.read_table(args.input)
    
    logging.info('Imputing.')
    af = df[args.freq].values
    zz = df[args.zscore].values
    nn = df[args.sample_size].values
    bhat_se = 1 / np.sqrt(2 * nn * af * (1 - af))
    bhat = zz * bhat_se
    
    logging.info('Saving.')
    df['effect_size'] = bhat
    df['standard_error'] = bhat_se
    df.to_csv(args.output, compression='gzip', sep='\t', index=False)
    
    logging.info('Done.')
    
    
