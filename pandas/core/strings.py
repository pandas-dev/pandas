import numpy as np

from functools import wraps
from itertools import izip
from pandas.core.common import isnull
from pandas.core.series import Series
import re
import pandas.lib as lib
import pandas.core.common as com

class repeat(object):
    def __init__(self, obj):
        self.obj = obj

    def __getitem__(self, i):
        return self.obj

class azip(object):
    def __init__(self, *args):
        self.cols = []
        for a in args:
            if np.isscalar(a):
                self.cols.append(repeat(a))
            else:
                self.cols.append(a)

    def __getitem__(self, i):
        return [col[i] for col in self.cols]


def map_iter_args(arr, f, otherargs, n_otherargs, required, n_results):
    '''
    Substitute for np.vectorize with pandas-friendly dtype inference

    Parameters
    ----------
    arr : ndarray
    f : function

    Returns
    -------
    mapped : ndarray
    '''
    notnull = com.notnull

    n = len(arr)
    result = np.empty((n, n_results), dtype=object)
    for i, val in enumerate(arr):
        args = otherargs[i]
        if notnull(val) and all(notnull(args[r]) for r in required):
            result[i] = f(val, *args)
        else:
            result[i] = [np.nan] * n_results

    return [lib.maybe_convert_objects(col, try_float=0) for col in result.T]


def auto_map(arr, f, otherargs, n_results=1, required='all'):
    from pandas.core.series import Series

    if all(np.isscalar(a) for a in otherargs):
        res = lib.map_infer(arr, lambda v: f(v, *otherargs))
        return Series(res, index=arr.index, copy=False)

    n_otherargs = len(otherargs)
    if required == 'all':
        required = list(range(n_otherargs))
    res = map_iter_args(arr, f, azip(*otherargs), n_otherargs,
                        required, n_results)
    res = [Series(col, index=arr.index, copy=False) for col in res]
    if n_results == 1:
        return res[0]
    return res


def mapwrap(f, n_results_default=1, required='all'):
    # @wraps(f)

    def wrapped(arr, n_results=None, *otherargs):
        n_results = n_results or n_results_default
        return auto_map(arr, f, otherargs, n_results, required)

    return wrapped

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

    return auto_map(arr, _re_get_groups(pattern, n_groups),
                    (groups,), n_results=n_groups)


def _get_array_list(arr, others):
    if isinstance(others[0], (list, np.ndarray)):
        arrays = [arr] + list(others)
    else:
        arrays = [arr, others]

    return [np.asarray(x, dtype=object) for x in arrays]


def str_cat(arr, others=None, sep=None, na_rep=None):
    """
    Concatenate arrays of strings with given separator

    Parameters
    ----------
    arr : list or array-like
    others : list or array, or list of arrays
    sep : string or None, default None
    na_rep : string or None, default None
        If None, an NA in any array will propagate

    Returns
    -------
    concat : array
    """
    if sep is None:
        sep = ''

    if others is not None:
        arrays = _get_array_list(arr, others)

        n = _length_check(arrays)
        masks = np.array([isnull(x) for x in arrays])
        cats = None

        if na_rep is None:
            na_mask = np.logical_or.reduce(masks, axis=0)

            result = np.empty(n, dtype=object)
            np.putmask(result, na_mask, np.nan)

            notmask = -na_mask

            if sep is None:
                for x in arrays:
                    x = x[notmask]
                    if cats is None:
                        cats = x
                    else:
                        cats = cats + x[notmask]
            else:
                tuples = izip(*[x[notmask] for x in arrays])
                cats = [sep.join(tup) for tup in tuples]

            result[notmask] = cats
        else:
            for i, x in enumerate(arrays):
                x = np.where(masks[i], na_rep, x)
                if cats is None:
                    cats = x
                else:
                    cats = cats + sep + x

            result = cats

        return result
    else:
        arr = np.asarray(arr, dtype=object)
        mask = isnull(arr)
        if na_rep is None and mask.any():
            return np.nan
        return sep.join(np.where(mask, na_rep, arr))


def _length_check(others):
    n = None
    for x in others:
        if n is None:
            n = len(x)
        elif len(x) != n:
            raise ValueError('All arrays must be same length')

    return n

def _na_map(f, arr):
    def g(x):
        try:
            return f(x)
        except TypeError:
            return np.nan
    return _map(g, arr)

def _map(f, arr):
    if not isinstance(arr, np.ndarray):
        arr = np.asarray(arr, dtype=object)
    return lib.map_infer(arr, f)


def str_count(arr, pat):
    """
    Count occurrences of pattern in each string

    Parameters
    ----------
    arr : list or array-like
    pat : string, valid regular expression

    Returns
    -------
    counts : arrays
    """
    regex = re.compile(pat)
    f = lambda x: len(regex.findall(x))
    return _na_map(f, arr)


def str_contains(arr, pat):
    """
    Check whether given pattern is contained in each string in the array

    Parameters
    ----------

    Returns
    -------

    """
    regex = re.compile(pat)
    f = lambda x: bool(regex.search(x))
    return _na_map(f, arr)

def str_startswith(arr, pat):
    """

    Parameters
    ----------

    Returns
    -------

    """
    pass

def str_endswith(arr, pat):
    """

    Parameters
    ----------

    Returns
    -------

    """
    pass

def str_lower(arr):
    """

    Parameters
    ----------

    Returns
    -------

    """
    return _na_map(str.lower, arr)


def str_upper(arr):
    """

    Parameters
    ----------

    Returns
    -------

    """
    return _na_map(str.upper, arr)


def str_replace(arr, pat, repl, n=None):
    """

    Parameters
    ----------

    Returns
    -------

    """
    pass


def str_repeat(arr, repeats):
    """

    Parameters
    ----------

    Returns
    -------

    """
    pass


def str_match(arr, pat):
    """

    Parameters
    ----------

    Returns
    -------

    """
    pass


def str_matchall(arr, pat):
    """

    Parameters
    ----------

    Returns
    -------

    """
    pass


def str_join(arr, joiner):
    """

    Parameters
    ----------

    Returns
    -------

    """
    pass


def str_len(arr):
    """

    Parameters
    ----------

    Returns
    -------

    """
    pass


def str_find(arr, pat):
    """

    Parameters
    ----------

    Returns
    -------

    """
    pass


def str_findall(arr, pat):
    """

    Parameters
    ----------

    Returns
    -------

    """
    pass


def str_pad(arr, side='left'):
    """

    Parameters
    ----------

    Returns
    -------

    """
    pass


def str_center(arr, width):
    pass


def str_split(arr, pat, n=None):
    """

    Parameters
    ----------

    Returns
    -------

    """
    pass


def str_slice(arr, i=None, j=None):
    """

    Parameters
    ----------

    Returns
    -------

    """
    pass


def str_slice_replace(arr, i=None, j=None, repl=None):
    """

    Parameters
    ----------

    Returns
    -------

    """
    pass


def str_strip(arr):
    """

    Parameters
    ----------

    Returns
    -------

    """
    pass


def str_lstrip(arr):
    """

    Parameters
    ----------

    Returns
    -------

    """
    pass


def str_rstrip(arr):
    """

    Parameters
    ----------

    Returns
    -------

    """
    pass


def str_wrap(arr, width=80, indent=0):
    """

    Parameters
    ----------

    Returns
    -------

    """
    pass


class StringMethods(object):

    def __init__(self, series):
        self.series = series

    def _wrap_result(self, result):
        return Series(result, index=self.series.index,
                      name=self.series.name)

    def cat(self, others=None, sep=None, na_rep=None):
        result = str_cat(self.series, others=others, sep=sep, na_rep=na_rep)
        return self._wrap_result(result)

    def count(self, pat):
        result = str_count(self.series, pat)
        return self._wrap_result(result)

    contains = str_contains
    lower = str_lower
    upper = str_upper
    replace = str_replace
    repeat = str_repeat
    match = str_match
    matchall = str_matchall
    join = str_join
    len = str_len
    find = str_find
    findall = str_findall
    pad = str_pad
    split = str_split

    slice = str_slice
    slice_replace = str_slice_replace

    strip = str_strip
    rstrip = str_rstrip
    lstrip = str_lstrip
