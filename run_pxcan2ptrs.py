import pandas as pd

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='run_pxcan2ptrs.py', description='''
        Convert S-PrediXcan results to PTRS weights along a lambda sequence. \\
        Need to have transethnic_prs and SPrediXcan2PTRS in PYTHONPATH.
    ''')
    parser.add_argument('--predixcan', help='''
        Path to (S)PrediXcan results \\
        (CSV format from SPrediXcan.py).
    ''')
    parser.add_argument('--predictdb', help='''
        Path to the predictDB file.
    ''')
    parser.add_argument('--geno_cov', help='''
        Pattern path to the genotype covariance. \\
        This should match the predictDB file, \\
        otherwise many genes and SNPs will be missed. \\
        Use {chr_num} in replace of the chromosome numbering. \\
        For instance, en_Whole_Blood.geno_cov.chr{chr_num}.evd.npz
    ''')
    parser.add_argument('--gwas', nargs='+', help='''
        The GWAS file being used for (S)PrediXcan analysis. \\
        We need this since we want to use the same SNPs as the ones used in (S)PrediXcan. \\
        And specify 4 or 5 columns: chromosome, position, effect_allele, non_effect_allele, sample_size (optional) \\
        For instance: \\
            --gwas filename \\
            chromosome:chr \\
            position:pos \\
            effect_allele:a1 \\
            non_effect_allele: a2 \\
            sample_size:gwas_N
    ''')
    parser.add_argument('--gwas_sample_size', type=int, help='''
        If no sample_size being specified in --gwas, \\
        need to set --gwas_sample_size.
    ''')
    parser.add_argument('--hyperparam_yaml', default=None, help='''
        A YAML file containing the hyper-parameters. \\
        For instance (values listed are default): \\
        alpha: [ 1 ] # any value in (0, 1] \\
        nlambda: 100 # any integer > 0 \\
        offset: [ 0.01 ] # any value in [0, 1) \\
        ratio_lambda: 100 # any value > 1 \\
        Note that alpha and offset can take multiple values (a list).
    ''')
    parser.add_argument('--output', help='''
        Output file name.
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
    
    logging.info('Loading (S)PrediXcan results.')
    df_pxcan = pd.read_csv(args.predixcan)
    
    logging.info('Loading GWAS.')
    df_gwas = load_gwas_snp(args.gwas, args.gwas_sample_size)
    
    logging.info('Loading genotype covariances.')
    geno_cov = cc.GenoCov(args.geno_cov)
    
    logging.info('Loading predictDB.')
    func_get_snp_meta = lambda x: db.get_snp_meta_from_varID(x)
    weight = db.WeightDB(args.predictdb, func_get_snp_meta)
    
    logging.info('Initializing the solver.')
    solver = so.Solver(
        df_pxcan=df_pxcan[['gene', 'zscore']].copy(),
        weight_db=weight,
        geno_cov=geno_cov, 
        df_gwas_snp=df_gwas[['chrom', 'position', 'effect_allele', 'non_effect_allele', 'sample_size']].copy(),
        lazy=True,
        gene_list=genes    
    )
    
    