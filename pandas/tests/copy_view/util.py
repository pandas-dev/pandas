from pandas import Series
from pandas.core.arrays import BaseMaskedArray


def get_array(obj, col=None):
    """
    Helper method to get array for a DataFrame column or a Series.

    Equivalent of df[col].values, but without going through normal getitem,
    which triggers tracking references / CoW (and we might be testing that
    this is done by some other operation).
    """
    if isinstance(obj, Series) and (col is None or obj.name == col):
        arr = obj._values
    else:
        assert col is not None
        icol = obj.columns.get_loc(col)
        assert isinstance(icol, int)
        arr = obj._get_column_array(icol)
    if isinstance(arr, BaseMaskedArray):
        return arr._data
    return arr
