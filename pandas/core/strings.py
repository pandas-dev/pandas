import numpy as np
from pandas.util.map import mapwrap, auto_map
import re

startswith = mapwrap(str.startswith)
contains = mapwrap(str.__contains__)
upper = mapwrap(str.upper)
lower = mapwrap(str.lower)

def _re_get_groups(pattern, n):
    def inner(s, *groups):
        m = pattern.search(s)
        if m:
            return m.group(*[int(g) for g in groups])
        return np.nan if n == 1 else [np.nan] * n
    
    return inner

def search_re(arr, pattern, groups=(0,)):
    if isinstance(pattern, str):
        pattern = re.compile(pattern)
    
    if isinstance(groups, np.ndarray):
        if groups.ndim == 1:
            n_groups = 1
        else:
            n_groups = groups.shape[1]
    else:
        n_groups = len(groups)
    
    return auto_map(arr, _re_get_groups(pattern, n_groups), (groups,), n_results=n_groups)
