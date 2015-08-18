import operator
import warnings
from pandas.compat import reduce
from pandas.core.index import Index
import numpy as np
from pandas import algos
from pandas.core import common as com


def match(needles, haystack):
    haystack = Index(haystack)
    needles = Index(needles)
    return haystack.get_indexer(needles)


def cartesian_product(X):
    '''
    Numpy version of itertools.product or pandas.compat.product.
    Sometimes faster (for large inputs)...

    Examples
    --------
    >>> cartesian_product([list('ABC'), [1, 2]])
    [array(['A', 'A', 'B', 'B', 'C', 'C'], dtype='|S1'),
    array([1, 2, 1, 2, 1, 2])]

    '''

    lenX = np.fromiter((len(x) for x in X), dtype=int)
    cumprodX = np.cumproduct(lenX)

    a = np.roll(cumprodX, 1)
    a[0] = 1

    b = cumprodX[-1] / cumprodX

    return [np.tile(np.repeat(np.asarray(com._values_from_object(x)), b[i]),
                    np.product(a[i]))
               for i, x in enumerate(X)]


def _compose2(f, g):
    """Compose 2 callables"""
    return lambda *args, **kwargs: f(g(*args, **kwargs))


def compose(*funcs):
    """Compose 2 or more callables"""
    assert len(funcs) > 1, 'At least 2 callables must be passed to compose'
    return reduce(_compose2, funcs)


def to_numeric(arg, errors='raise', box=True, coerce=None):
    """
    Convert argument to a numeric type.

    Parameters
    ----------
    arg : string, datetime, array of strings (with possible NAs)
    errors : {'ignore', 'raise', 'coerce'}, default 'raise'
        - If 'raise', then invalid parsing will raise an exception
        - If 'coerce', then invalid parsing will be set as NaT
        - If 'ignore', then invalid parsing will return the input
    box : boolean, default True
        - If True returns a Series
        - If False returns ndarray of values.

    Returns
    -------
    ret : numeric if parsing succeeded.
        Return type depends on box


    Examples
    --------
    Take separate series and convert to datetime

    >>> import pandas as pd
    >>> df = pd.DataFrame(['1.0', '2', -3])
    >>> pd.to_numeric(df)
    >>> df = pd.DataFrame(['apple', '1.0', '2', -3])
    >>> pd.to_numeric(df, errors='ignore')
    >>> pd.to_numeric(df, errors='coerce')
    """
    #TODO: Fix examples

    coerce_numeric = False if errors in ('ignore', 'raise') else True
    if errors == 'ignore':
        try:
            values = lib.maybe_convert_numeric(arg,
                                               set(),
                                               coerce_numeric=coerce_numeric)
            return values
        except:
            return arg
    return lib.maybe_convert_numeric(arg, set(), coerce_numeric=coerce_numeric)
