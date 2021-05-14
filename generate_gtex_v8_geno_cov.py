import os

def load_list(fn):
    o = []
    with open(fn, 'r') as f:
        for i in f:
            i = i.strip()
            o.append(i)
    return o

def load_mode(mode_list):
    if mode_list[0] not in ['naive', 'cap', 'banded', 'evd']:
        raise ValueError('Wrong mode.')
    else:
        if len(mode_list) != 2:
            raise ValueError('Wrong number of parameters for mode.')
        if mode_list[0] == 'cap':
            param = float(mode_list[1])
            ext = 'cap.npz'
        elif mode_list[0] == 'banded':
            param = float(mode_list[1])
            ext = 'banded.npz'
        elif mode_list[0] == 'evd':
            param = float(mode_list[1])
            ext = 'evd.npz'
        elif mode_list[0] == 'naive':
            ext = 'naive.h5'
            if mode_list[1] == 'f32':
                param = np.float32
            elif mode_list[1] == 'f64':
                param = np.float64
            else:
                raise ValueError('Wrong parameter in mode = naive: {}'.format(mode_list[1]))
    return mode_list[0], param, ext

def file_exists(fn):
    return os.path.exists(fn)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='generate_gtex_v8_geno_cov.py', description='''
        Prepare genotype covariance for PredictDBs which 
        were built using GTEx V8 data. 
    ''')
    parser.add_argument('--genotype_vcf', help='''
        Path to the genotype VCF. 
        It contains all GTEx V8 individuals across all chromosomes. 
    ''')
    parser.add_argument('--predictdb', help='''
        PredictDB file. 
    ''')
    parser.add_argument('--sample_list', default=None, help='''
        The list of samples to use in genotype.
    ''')
    parser.add_argument('--mode', nargs='+', help='''
        Indicate the mode and parameter of the mode for genotype covariance:
        1. banded [band-size]; 
        2. cap [threshold-of-cap];
        3. naive [f32|f64];
        4. evd [min-max-threshold]. 
    ''')
    parser.add_argument('--output_prefix', help='''
        Output prefix.
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
    
    import SPrediXcan2PTRS.geno_cov.gtex_v8_loader as v8
    import SPrediXcan2PTRS.geno_cov.cov_constructor as cc
    
    mode, param, out_ext = load_mode(args.mode)
    
    logging.info('Constructing Genotype and DB loaders.')
    weight_v8 = v8.GTExV8DBLoader(args.predictdb)
    geno_v8 = v8.GTExV8GenoLoader(args.genotype_vcf)
    
    if args.sample_list is not None:
        logging.info('Loading sample list for genotype.')
        samples = load_list(args.sample_list)
        logging.info('{} samples being loaded.'.format(len(samples)))
    else:
        samples = None
    
    for chr_num in range(1, 23):
        output_prefix = '{}.chr{}'.format(args.output_prefix, chr_num)
        target_file = f'{output_prefix}.{out_ext}'
        snp_meta_file = f'{output_prefix}.snp_meta.parquet'
        if file_exists(target_file) and file_exists(snp_meta_file):
            logging.info(f'Target file exists for chromosome {chr_num}.')
            continue
        logging.info(f'Extracting SNPs on chromosome {chr_num}.')
        weight_sub = weight_v8.get_by_chr(f'chr{chr_num}')
        tmp = geno_v8.load(weight_sub.rename(columns={'ref_allele': 'ref', 'eff_allele': 'alt'}), samples_only=samples)
        nn = tmp[0].shape[1]
        logging.info(f'Extracting SNPs on chromosome {chr_num}: sample size = {nn}')
        cov_constructor = cc.CovConstructor(tmp[0])
        logging.info(f'Saving geno cov for chromosome {chr_num}.')
        cov_constructor.compute_to_disk(
            mode, 
            output_prefix, 
            param
        )
        logging.info(f'Saving SNP meta for chromosome {chr_num}.')
        tmp[1].to_parquet(snp_meta_file, index=False)
    
    logging.info('Done.') 
       
