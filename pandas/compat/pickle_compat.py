"""
Support pre-0.12 series pickle compatibility.
"""

import copy
import pickle as pkl
from typing import TYPE_CHECKING, Optional
import warnings

from pandas import Index

if TYPE_CHECKING:
    from pandas import Series, DataFrame


def load_reduce(self):
    stack = self.stack
    args = stack.pop()
    func = stack[-1]

    if len(args) and type(args[0]) is type:
        n = args[0].__name__  # noqa

    try:
        stack[-1] = func(*args)
        return
    except TypeError as err:

        # If we have a deprecated function,
        # try to replace and try again.

        msg = "_reconstruct: First argument must be a sub-type of ndarray"

        if msg in str(err):
            try:
                cls = args[0]
                stack[-1] = object.__new__(cls)
                return
            except TypeError:
                pass

        raise


_sparse_msg = """\

Loading a saved '{cls}' as a {new} with sparse values.
'{cls}' is now removed. You should re-save this dataset in its new format.
"""


class _LoadSparseSeries:
    # To load a SparseSeries as a Series[Sparse]

    # https://github.com/python/mypy/issues/1020
    # error: Incompatible return type for "__new__" (returns "Series", but must return
    # a subtype of "_LoadSparseSeries")
    def __new__(cls) -> "Series":  # type: ignore
        from pandas import Series

        warnings.warn(
            _sparse_msg.format(cls="SparseSeries", new="Series"),
            FutureWarning,
            stacklevel=6,
        )

        return Series(dtype=object)


class _LoadSparseFrame:
    # To load a SparseDataFrame as a DataFrame[Sparse]

    # https://github.com/python/mypy/issues/1020
    # error: Incompatible return type for "__new__" (returns "DataFrame", but must
    # return a subtype of "_LoadSparseFrame")
    def __new__(cls) -> "DataFrame":  # type: ignore
        from pandas import DataFrame

        warnings.warn(
            _sparse_msg.format(cls="SparseDataFrame", new="DataFrame"),
            FutureWarning,
            stacklevel=6,
        )

        return DataFrame()


# If classes are moved, provide compat here.
_class_locations_map = {
    ("pandas.core.sparse.array", "SparseArray"): ("pandas.core.arrays", "SparseArray"),
    # 15477
    ("pandas.core.base", "FrozenNDArray"): ("numpy", "ndarray"),
    ("pandas.core.indexes.frozen", "FrozenNDArray"): ("numpy", "ndarray"),
    ("pandas.core.base", "FrozenList"): ("pandas.core.indexes.frozen", "FrozenList"),
    # 10890
    ("pandas.core.series", "TimeSeries"): ("pandas.core.series", "Series"),
    ("pandas.sparse.series", "SparseTimeSeries"): (
        "pandas.core.sparse.series",
        "SparseSeries",
    ),
    # 12588, extensions moving
    ("pandas._sparse", "BlockIndex"): ("pandas._libs.sparse", "BlockIndex"),
    ("pandas.tslib", "Timestamp"): ("pandas._libs.tslib", "Timestamp"),
    # 18543 moving period
    ("pandas._period", "Period"): ("pandas._libs.tslibs.period", "Period"),
    ("pandas._libs.period", "Period"): ("pandas._libs.tslibs.period", "Period"),
    # 18014 moved __nat_unpickle from _libs.tslib-->_libs.tslibs.nattype
    ("pandas.tslib", "__nat_unpickle"): (
        "pandas._libs.tslibs.nattype",
        "__nat_unpickle",
    ),
    ("pandas._libs.tslib", "__nat_unpickle"): (
        "pandas._libs.tslibs.nattype",
        "__nat_unpickle",
    ),
    # 15998 top-level dirs moving
    ("pandas.sparse.array", "SparseArray"): (
        "pandas.core.arrays.sparse",
        "SparseArray",
    ),
    ("pandas.sparse.series", "SparseSeries"): (
        "pandas.compat.pickle_compat",
        "_LoadSparseSeries",
    ),
    ("pandas.sparse.frame", "SparseDataFrame"): (
        "pandas.core.sparse.frame",
        "_LoadSparseFrame",
    ),
    ("pandas.indexes.base", "_new_Index"): ("pandas.core.indexes.base", "_new_Index"),
    ("pandas.indexes.base", "Index"): ("pandas.core.indexes.base", "Index"),
    ("pandas.indexes.numeric", "Int64Index"): (
        "pandas.core.indexes.numeric",
        "Int64Index",
    ),
    ("pandas.indexes.range", "RangeIndex"): ("pandas.core.indexes.range", "RangeIndex"),
    ("pandas.indexes.multi", "MultiIndex"): ("pandas.core.indexes.multi", "MultiIndex"),
    ("pandas.tseries.index", "_new_DatetimeIndex"): (
        "pandas.core.indexes.datetimes",
        "_new_DatetimeIndex",
    ),
    ("pandas.tseries.index", "DatetimeIndex"): (
        "pandas.core.indexes.datetimes",
        "DatetimeIndex",
    ),
    ("pandas.tseries.period", "PeriodIndex"): (
        "pandas.core.indexes.period",
        "PeriodIndex",
    ),
    # 19269, arrays moving
    ("pandas.core.categorical", "Categorical"): ("pandas.core.arrays", "Categorical"),
    # 19939, add timedeltaindex, float64index compat from 15998 move
    ("pandas.tseries.tdi", "TimedeltaIndex"): (
        "pandas.core.indexes.timedeltas",
        "TimedeltaIndex",
    ),
    ("pandas.indexes.numeric", "Float64Index"): (
        "pandas.core.indexes.numeric",
        "Float64Index",
    ),
    ("pandas.core.sparse.series", "SparseSeries"): (
        "pandas.compat.pickle_compat",
        "_LoadSparseSeries",
    ),
    ("pandas.core.sparse.frame", "SparseDataFrame"): (
        "pandas.compat.pickle_compat",
        "_LoadSparseFrame",
    ),
}


# our Unpickler sub-class to override methods and some dispatcher
# functions for compat and uses a non-public class of the pickle module.

# error: Name 'pkl._Unpickler' is not defined
class Unpickler(pkl._Unpickler):  # type: ignore
    def find_class(self, module, name):
        # override superclass
        key = (module, name)
        module, name = _class_locations_map.get(key, key)
        return super().find_class(module, name)


Unpickler.dispatch = copy.copy(Unpickler.dispatch)
Unpickler.dispatch[pkl.REDUCE[0]] = load_reduce


def load_newobj(self):
    args = self.stack.pop()
    cls = self.stack[-1]

    # compat
    if issubclass(cls, Index):
        obj = object.__new__(cls)
    else:
        obj = cls.__new__(cls, *args)

    self.stack[-1] = obj


Unpickler.dispatch[pkl.NEWOBJ[0]] = load_newobj


def load_newobj_ex(self):
    kwargs = self.stack.pop()
    args = self.stack.pop()
    cls = self.stack.pop()

    # compat
    if issubclass(cls, Index):
        obj = object.__new__(cls)
    else:
        obj = cls.__new__(cls, *args, **kwargs)
    self.append(obj)


try:
    Unpickler.dispatch[pkl.NEWOBJ_EX[0]] = load_newobj_ex
except (AttributeError, KeyError):
    pass


def load(fh, encoding: Optional[str] = None, is_verbose: bool = False):
    """
    Load a pickle, with a provided encoding,

    Parameters
    ----------
    fh : a filelike object
    encoding : an optional encoding
    is_verbose : show exception output
    """
    try:
        fh.seek(0)
        if encoding is not None:
            up = Unpickler(fh, encoding=encoding)
        else:
            up = Unpickler(fh)
        up.is_verbose = is_verbose

        return up.load()
    except (ValueError, TypeError):
        raise
