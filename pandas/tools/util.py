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


def to_numeric(arg, errors='raise'):
    """
    Convert argument to a numeric type.

    Parameters
    ----------
    arg : list, tuple or array of objects, or Series
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

    index = name = None
    if isinstance(arg, pd.Series):
        index, name = arg.index, arg.name
    elif isinstance(arg, (list, tuple)):
        arg = np.array(arg, dtype='O')

    conv = arg
    arg = com._ensure_object(arg)

    coerce_numeric = False if errors in ('ignore', 'raise') else True

    try:
        conv = lib.maybe_convert_numeric(arg,
                                         set(),
                                         coerce_numeric=coerce_numeric)
    except:
        if errors == 'raise':
            raise

    if index is not None:
        return pd.Series(conv, index=index, name=name)
    else:
        return conv
