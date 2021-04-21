import yaml

def read_yaml(yaml_):
    with open(yaml_, 'r') as f:
        o = yaml.safe_load(f)
    return o
def _try_cast_to_float_list(val_):
    if isinstance(val_, float) or isinstance(val_, int):
        val_ = [ float(val_) ]
    else:
        raise ValueError(f'{val_} is unexpected'.)
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
       
            