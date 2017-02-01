import numpy as np
import pandas.lib as lib

from pandas.types.common import (is_number,
                                 is_numeric_dtype,
                                 is_datetime_or_timedelta_dtype,
                                 is_list_like,
                                 _ensure_object,
                                 is_decimal,
                                 is_scalar as isscalar)

from pandas.types.cast import _possibly_downcast_to_dtype

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

    Parameters
    ----------
    X : list-like of list-likes

    Returns
    -------
    product : list of ndarrays

    Examples
    --------
    >>> cartesian_product([list('ABC'), [1, 2]])
    [array(['A', 'A', 'B', 'B', 'C', 'C'], dtype='|S1'),
    array([1, 2, 1, 2, 1, 2])]

    See also
    --------
    itertools.product : Cartesian product of input iterables.  Equivalent to
        nested for-loops.
    pandas.compat.product : An alias for itertools.product.
    """
    msg = "Input must be a list-like of list-likes"
    if not is_list_like(X):
        raise TypeError(msg)
    for x in X:
        if not is_list_like(x):
            raise TypeError(msg)

    if len(X) == 0:
        return []

    lenX = np.fromiter((len(x) for x in X), dtype=np.intp)
    cumprodX = np.cumproduct(lenX)

    a = np.roll(cumprodX, 1)
    a[0] = 1

    if cumprodX[-1] != 0:
        b = cumprodX[-1] / cumprodX
    else:
        # if any factor is empty, the cartesian product is empty
        b = np.zeros_like(cumprodX)

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


def to_numeric(arg, errors='raise', downcast=None):
    """
    Convert argument to a numeric type.

    Parameters
    ----------
    arg : list, tuple, 1-d array, or Series
    errors : {'ignore', 'raise', 'coerce'}, default 'raise'
        - If 'raise', then invalid parsing will raise an exception
        - If 'coerce', then invalid parsing will be set as NaN
        - If 'ignore', then invalid parsing will return the input
    downcast : {'integer', 'signed', 'unsigned', 'float'} , default None
        If not None, and if the data has been successfully cast to a
        numerical dtype (or if the data was numeric to begin with),
        downcast that resulting data to the smallest numerical dtype
        possible according to the following rules:

        - 'integer' or 'signed': smallest signed int dtype (min.: np.int8)
        - 'unsigned': smallest unsigned int dtype (min.: np.uint8)
        - 'float': smallest float dtype (min.: np.float32)

        As this behaviour is separate from the core conversion to
        numeric values, any errors raised during the downcasting
        will be surfaced regardless of the value of the 'errors' input.

        In addition, downcasting will only occur if the size
        of the resulting data's dtype is strictly larger than
        the dtype it is to be cast to, so if none of the dtypes
        checked satisfy that specification, no downcasting will be
        performed on the data.

        .. versionadded:: 0.19.0

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
    0    1.0
    1    2.0
    2   -3.0
    dtype: float64
    >>> pd.to_numeric(s, downcast='float')
    0    1.0
    1    2.0
    2   -3.0
    dtype: float32
    >>> pd.to_numeric(s, downcast='signed')
    0    1
    1    2
    2   -3
    dtype: int8
    >>> s = pd.Series(['apple', '1.0', '2', -3])
    >>> pd.to_numeric(s, errors='ignore')
    0    apple
    1      1.0
    2        2
    3       -3
    dtype: object
    >>> pd.to_numeric(s, errors='coerce')
    0    NaN
    1    1.0
    2    2.0
    3   -3.0
    dtype: float64
    """
    if downcast not in (None, 'integer', 'signed', 'unsigned', 'float'):
        raise ValueError('invalid downcasting method provided')

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
    elif isscalar(arg):
        if is_decimal(arg):
            return float(arg)
        if is_number(arg):
            return arg
        is_scalar = True
        values = np.array([arg], dtype='O')
    elif getattr(arg, 'ndim', 1) > 1:
        raise TypeError('arg must be a list, tuple, 1-d array, or Series')
    else:
        values = arg

    try:
        if is_numeric_dtype(values):
            pass
        elif is_datetime_or_timedelta_dtype(values):
            values = values.astype(np.int64)
        else:
            values = _ensure_object(values)
            coerce_numeric = False if errors in ('ignore', 'raise') else True
            values = lib.maybe_convert_numeric(values, set(),
                                               coerce_numeric=coerce_numeric)

    except Exception:
        if errors == 'raise':
            raise

    # attempt downcast only if the data has been successfully converted
    # to a numerical dtype and if a downcast method has been specified
    if downcast is not None and is_numeric_dtype(values):
        typecodes = None

        if downcast in ('integer', 'signed'):
            typecodes = np.typecodes['Integer']
        elif downcast == 'unsigned' and np.min(values) >= 0:
            typecodes = np.typecodes['UnsignedInteger']
        elif downcast == 'float':
            typecodes = np.typecodes['Float']

            # pandas support goes only to np.float32,
            # as float dtypes smaller than that are
            # extremely rare and not well supported
            float_32_char = np.dtype(np.float32).char
            float_32_ind = typecodes.index(float_32_char)
            typecodes = typecodes[float_32_ind:]

        if typecodes is not None:
            # from smallest to largest
            for dtype in typecodes:
                if np.dtype(dtype).itemsize <= values.dtype.itemsize:
                    values = _possibly_downcast_to_dtype(
                        values, dtype)

                    # successful conversion
                    if values.dtype == dtype:
                        break

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
