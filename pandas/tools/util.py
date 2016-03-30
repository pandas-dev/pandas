import numpy as np
import pandas.lib as lib

import pandas as pd
from pandas.compat import reduce
from pandas.core.index import Index
from pandas.core import common as com


def match(needles, haystack):
    haystack = Index(haystack)
    needles = Index(needles)
    return haystack.get_indexer(needles)


def cartesian_product(X):
    """
    Numpy version of itertools.product or pandas.compat.product.
    Sometimes faster (for large inputs)...

    Examples
    --------
    >>> cartesian_product([list('ABC'), [1, 2]])
    [array(['A', 'A', 'B', 'B', 'C', 'C'], dtype='|S1'),
    array([1, 2, 1, 2, 1, 2])]

    """

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


def to_numeric(arg, errors='raise'):
    """
    Convert argument to a numeric type.

    Parameters
    ----------
    arg : list, tuple, 1-d array, or Series
    errors : {'ignore', 'raise', 'coerce'}, default 'raise'
        - If 'raise', then invalid parsing will raise an exception
        - If 'coerce', then invalid parsing will be set as NaN
        - If 'ignore', then invalid parsing will return the input

    Returns
    -------
    ret : numeric if parsing succeeded.
        Return type depends on input.  Series if Series, otherwise ndarray

    Examples
    --------
    Take separate series and convert to numeric, coercing when told to

    >>> import pandas as pd
    >>> s = pd.Series(['1.0', '2', -3])
    >>> pd.to_numeric(s)
    >>> s = pd.Series(['apple', '1.0', '2', -3])
    >>> pd.to_numeric(s, errors='ignore')
    >>> pd.to_numeric(s, errors='coerce')
    """
    is_series = False
    is_index = False
    is_scalar = False

    if isinstance(arg, pd.Series):
        is_series = True
        values = arg.values
    elif isinstance(arg, pd.Index):
        is_index = True
        values = arg.asi8
        if values is None:
            values = arg.values
    elif isinstance(arg, (list, tuple)):
        values = np.array(arg, dtype='O')
    elif np.isscalar(arg):
        if com.is_number(arg):
            return arg
        is_scalar = True
        values = np.array([arg], dtype='O')
    elif getattr(arg, 'ndim', 1) > 1:
        raise TypeError('arg must be a list, tuple, 1-d array, or Series')
    else:
        values = arg

    if com.is_numeric_dtype(values):
        pass
    elif com.is_datetime_or_timedelta_dtype(values):
        values = values.astype(np.int64)
    else:
        values = com._ensure_object(values)
        coerce_numeric = False if errors in ('ignore', 'raise') else True

        try:
            values = lib.maybe_convert_numeric(values, set(),
                                               coerce_numeric=coerce_numeric)
        except:
            if errors == 'raise':
                raise

    if is_series:
        return pd.Series(values, index=arg.index, name=arg.name)
    elif is_index:
        # because we want to coerce to numeric if possible,
        # do not use _shallow_copy_with_infer
        return Index(values, name=arg.name)
    elif is_scalar:
        return values[0]
    else:
        return values
