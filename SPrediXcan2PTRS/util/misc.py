import os.path
import yaml
import re
import warnings

import pandas as pd

GWAS_PARSER = {
    '{chromosome}': '(?P<chromosome>(?:chr|)[0-9]+)',
    '{position}': '(?P<position>[0-9]+)',
    '{effect_allele}': '(?P<effect_allele>[A-Z]+)',
    '{non_effect_allele}': '(?P<non_effect_allele>[A-Z]+)'
}
GWAS_DEFAULT_NA = {
    'chromosome': 'chr0',
    'position': '1',
    'effect_allele': 'A',
    'non_effect_allele': 'T'
}
# kk = '{chr}:{pos}'; kk = re.sub('{chr}', '(?P<chr>(?:chr|)[0-9]+)', kk); kk = re.sub('{pos}', '(?P<pos>[0-9]+)', kk)

def read_table(fn, indiv_col=None):
    _, fn_ext = os.path.splitext(fn)
    compress_args = {}
    if fn_ext == '.gz':
        fn_new = re.sub('.gz$', '', fn)
        compress_args = {'compression': 'gzip'}
        _, fn_ext = os.path.splitext(fn_new)
    if fn_ext == '.parquet':
        df = pd.read_parquet(fn)
    elif fn_ext == '.csv':
        df = pd.read_csv(fn, **compress_args)
    elif fn_ext == '.txt' or fn_ext == '.tsv':
        df = pd.read_csv(fn, sep='\s+', **compress_args)
    for i in range(df.shape[1]):
        if df.columns[i] == indiv_col:
            break
    if indiv_col is None:
        return df
    col_list = df.columns.to_list()
    col_list.pop(i)
    col_list = [ indiv_col ] + col_list
    df = df.reindex(columns=col_list)
    df.rename(columns={indiv_col: 'indiv'}, inplace=True)
    df.indiv = df.indiv.astype(str)
    return df

def read_yaml(yaml_):
    with open(yaml_, 'r') as f:
        o = yaml.safe_load(f)
    return o
def _try_cast_to_float_list(val_):
    if isinstance(val_, float) or isinstance(val_, int):
        val_ = [ float(val_) ]
    else:
        raise ValueError(f'{val_} is unexpected.')
    return val_
def _check_is(val_, target):
    if not isinstance(val_, target):
        raise ValueError(f'{val_} is not an {target}.')
def _raise_valerr(ss=None):
    raise ValueError(ss)
def _check_range(val, start, end, include_start, include_end):
    if start is not None:
        if val > start or (include_start is True and val == start):
            pass
        else:
            _raise_valerr()
    if end is not None:
        if val < end or (include_end is True and val == end):
            pass
        else:
            _raise_valerr()
def check_param(key_, val_):
    if key_ == 'alpha':
        val_ = _try_cast_to_float_list(val_)
        for v in val_:
            _check_range(v, 0, 1, False, True)
    if key_ == 'offset':
        val_ = _try_cast_to_float_list(val_)
        for v in val_:
            _check_range(v, 0, 1, True, False)
    if key_ == 'nlambda':
        _check_is(val_, int)
        _check_range(val_, 0, None, False, None)
    if key_ == 'ratio_lambda':
        _check_is(val_, float)
        _check_range(val_, 1, None, False, None)
    if key_ == 'maxiter':
        _check_is(val_, int)
        _check_range(val_, 0, None, False, None)
    if key_ == 'tol':
        _check_is(val_, float)
        _check_range(val_, 0, None, False, None)
def try_cast_float(dict_, keys):
    for k in keys:
        if k in dict_:
            if isinstance(dict_[k], list):
                o = dict_[k]
            else:
                o = [ dict_[k] ]
            for i, v in enumerate(o):
                if isinstance(v, str):
                    o[i] = float(v)
            if isinstance(dict_[k], list):
                dict_[k] = o
            else:
                dict_[k] = o[0]
def parse_and_update_gwas_col(cols_dict, col, pattern):
    regex_str = pattern
    for k, v in GWAS_PARSER.items():
        if k in pattern and k not in cols_dict:
            regex_str = re.sub(k, v, regex_str)
            cols_dict[k] = []
            
    for i in col:
        tmp = re.match(regex_str, i)
        if tmp is None:
            warnings.warn(f'Fail to parse {i}.')
            res = GWAS_DEFAULT_NA
        else:
            res = tmp.groupdict()
        for k, v in res.items():
            cols_dict[k].append(res[k])      
    
