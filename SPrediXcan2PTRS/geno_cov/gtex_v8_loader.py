from collections import OrderedDict

import cyvcf2 as cy
import sqlite3 as sq
import pandas as pd
import numpy as np
# from tqdm import tqdm

COMPLEMENT_BASE = {
    'A': 'T',
    'G': 'C',
    'T': 'A',
    'C': 'G'
}
VARID_SEP = '_'

def get_complement(str_):
    o = ''
    for s in str_:
        char_ = s.upper()
        if char_ not in COMPLEMENT_BASE:
            raise ValueError(f'Wrong s in str_: s = {s}.')
        o += COMPLEMENT_BASE[char_]
    return o


class TargetSNP:
    def __init__(self, df_snp):
        self._init_snp_dict(df_snp)
    def _init_snp_dict(self, df_snp):
        self.snp_dict = OrderedDict()
        for chr, pos, ref, alt in zip(df_snp.chr, df_snp.pos, df_snp.ref, df_snp.alt):
            kk = self.gen_chrpos(chr, pos)
            if kk not in self.snp_dict:
                self.snp_dict[kk] = []
            self.snp_dict[kk].append(
                pd.DataFrame({ 
                    'chr': [chr], 'pos': [pos], 
                    'ref': [ref], 'alt': [alt]
                })
            )
    @staticmethod
    def _try_match(a1, a2, b1, b2):
        if a1 == b1 and a2 == b2:
            return 1
        if a1 == get_complement(b1) and a2 == get_complement(b2):
            return 1
        if a1 == b2 and a2 == b1:
            return -1
        if a1 == get_complement(b2) and a2 == get_complement(b1):
            return -1
        return 0
    def get_snp(self, key_, ref, alt):
        if key_ not in self.snp_dict:
            return None, None
        snps = self.snp_dict[key_]
        for i in range(len(snps)):
            res = self._try_match(snps[i].ref.values[0], snps[i].alt.values[0], ref, alt)
            if res == 0:
                continue
            else:
                return snps[i].copy(), res
        return None, None
    @staticmethod
    def gen_chrpos(chr, pos):
        return f'{chr}_x_{pos}'
            
        

class GTExV8GenoLoader:
    def __init__(self, geno):
        self.vcf = cy.VCF(geno)
    def load(self, df_snp):
        '''
        df_snp = pd.DataFrame({
            'chr': chromosome,
            'pos': position,
            'ref': reference allele,
            'alt': alternative allele
        })
        Return the genotype matrix (snp x sample)
        with matched alleles 
        '''
        chrs = df_snp.chr.unique()
        geno_mat = []
        snps = []
        for cc in chrs:
            df_snp_sub = df_snp[ df_snp.chr == cc ]
            start = df_snp_sub.pos.min()
            end = df_snp_sub.pos.max()
            snps_sub, geno_mat_sub = self._load_region(f'{cc}:{start - 1}-{end}', target_snps=df_snp_sub)
            geno_mat.append(geno_mat_sub)
            snps.append(snps_sub)
        geno_mat = np.concatenate(geno_mat, axis=1)
        snps = pd.concat(snps, axis=0)
        return geno_mat, snps
    def _load_region(self, region, target_snps):
        snp_target = TargetSNP(target_snps)
        snps = []
        geno_mat = []
        for kk in self.vcf(region):
            ch, pos, ref, alt = kk.CHROM, kk.POS, kk.REF, kk.ALT
            kk_key = snp_target.gen_chrpos(ch, pos)
            snpi, direction = snp_target.get_snp(kk_key, ref, alt)
            if snpi is None:
                continue
            else:
                mat = np.array(kk.genotypes)[:, :2].sum(axis=1)
                if direction == -1:
                    mat = 2 - mat
                geno_mat.append(mat[:, np.newaxis])
                snps.append(snpi)
        snps = pd.concat(snps, axis=0)
        geno_mat = np.concatenate(geno_mat, axis=1)
        return snps, geno_mat
            

class GTExV8DBLoader:
    def __init__(self, db):
        with sq.connect(db) as conn:
            df_weight = pd.read_sql_query('select * from weights', conn)
        df_weight['chr'], df_weight['pos'] = self._parse_varID(df_weight.varID)
        self.df_weight = df_weight
    @staticmethod
    def _parse_varID(slist):
        cc = []
        pp = []
        for s in slist:
            tmp = s.split(VARID_SEP)
            cc.append(tmp[0])
            pp.append(int(tmp[1]))
        return cc, pp
    def get_by_chr(self, chr_):
        return self.df_weight[ self.df_weight.chr == chr_ ].reset_index(drop=True)    
