import numpy as np

from pandas.core.dtypes.common import is_list_like
from pandas.core.dtypes.generic import ABCCategorical

from pandas.core.indexes.api import Index, IntervalIndex


def cartesian_product(X):
    """
    Numpy version of itertools.product.
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

    See Also
    --------
    itertools.product : Cartesian product of input iterables.  Equivalent to
        nested for-loops.
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

    return [_tile_compat(np.repeat(x, b[i]), np.product(a[i])) for i, x in enumerate(X)]


def _broadcast_tile(arr: np.ndarray, num: int) -> np.ndarray:
    """
    Emulate np.tile but using views instead of copies.
    """
    shape = (len(arr), num)
    middle = arr.reshape(len(arr), 1)
    new_arr = np.broadcast_to(middle, shape)

    # Note: doing `ravel` gives us the wrong order
    return new_arr.reshape(-1, order="F")


def _tile_compat(arr, num: int):
    """
    Index compat for np.tile.

    Notes
    -----
    Does not support multi-dimensional `num`.
    """
    if isinstance(arr, np.ndarray):
        return _broadcast_tile(arr, num)

    # Otherwise we have an Index
    values = arr._data

    if isinstance(values, np.ndarray):
        result = _broadcast_tile(values, num)
        return type(arr)._simple_new(result, name=arr.name)

    elif isinstance(values, ABCCategorical):
        codes = _broadcast_tile(values.codes, num)
        result = type(values).from_codes(codes, dtype=values.dtype)
        return type(arr)._simple_new(result, name=arr.name)

    elif isinstance(arr, IntervalIndex):
        new_left = _tile_compat(values.left, num)
        new_right = _tile_compat(values.right, num)
        result = type(values).from_arrays(new_left, new_right, closed=values.closed)
        return type(arr)._simple_new(result, name=arr.name)

    elif isinstance(values._data, np.ndarray):
        # DatetimeIndex, TimedeltaIndex, PeriodIndex
        data = _broadcast_tile(values._data, num)
        result = type(values)._simple_new(data, dtype=values.dtype)
        return type(arr)._simple_new(result, name=arr.name)

    else:
        # As of now this just leaves RangeIndex, which cannot
        #  use type(self)._simple_new
        result = _broadcast_tile(values, num)
        return Index(result, name=arr.name)
