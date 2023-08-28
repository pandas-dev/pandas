"""
Pickle compatibility to pandas version 1.0
"""
from __future__ import annotations

import contextlib
import copyreg
import io
import pickle
from typing import TYPE_CHECKING

import numpy as np

from pandas._libs.arrays import NDArrayBacked
from pandas._libs.tslibs import BaseOffset

from pandas.core.arrays import (
    DatetimeArray,
    PeriodArray,
    TimedeltaArray,
)
from pandas.core.internals import BlockManager

if TYPE_CHECKING:
    from collections.abc import Generator


# If classes are moved, provide compat here.
_class_locations_map = {
    # 50775, remove Int64Index, UInt64Index & Float64Index from codebase
    ("pandas.core.indexes.numeric", "Int64Index"): (
        "pandas.core.indexes.base",
        "Index",
    ),
    ("pandas.core.indexes.numeric", "UInt64Index"): (
        "pandas.core.indexes.base",
        "Index",
    ),
    ("pandas.core.indexes.numeric", "Float64Index"): (
        "pandas.core.indexes.base",
        "Index",
    ),
    ("pandas.core.arrays.sparse.dtype", "SparseDtype"): (
        "pandas.core.dtypes.dtypes",
        "SparseDtype",
    ),
}


def load_reduce(self):
    stack = self.stack
    args = stack.pop()
    func = stack[-1]

    try:
        stack[-1] = func(*args)
        return
    except TypeError:
        # If we have a deprecated function,
        # try to replace and try again.

        if args and isinstance(args[0], type) and issubclass(args[0], BaseOffset):
            # TypeError: object.__new__(Day) is not safe, use Day.__new__()
            cls = args[0]
            stack[-1] = cls.__new__(*args)
            return
        elif args and issubclass(args[0], PeriodArray):
            cls = args[0]
            stack[-1] = NDArrayBacked.__new__(*args)
            return

        raise


def load_newobj(self) -> None:
    args = self.stack.pop()
    cls = self.stack[-1]

    # compat
    if issubclass(cls, DatetimeArray) and not args:
        arr = np.array([], dtype="M8[ns]")
        obj = cls.__new__(cls, arr, arr.dtype)
    elif issubclass(cls, TimedeltaArray) and not args:
        arr = np.array([], dtype="m8[ns]")
        obj = cls.__new__(cls, arr, arr.dtype)
    elif cls is BlockManager and not args:
        obj = cls.__new__(cls, (), [], False)
    else:
        obj = cls.__new__(cls, *args)

    self.stack[-1] = obj


class Unpickler(pickle.Unpickler):
    dispatch_table = copyreg.dispatch_table.copy()
    dispatch_table[pickle.REDUCE[0]] = load_reduce
    dispatch_table[pickle.NEWOBJ[0]] = load_newobj

    def find_class(self, module, name):
        # override superclass
        key = (module, name)
        module, name = _class_locations_map.get(key, key)
        return super().find_class(module, name)


def load(fh, encoding: str | None = None):
    """
    Load a pickle, with a provided encoding,

    Parameters
    ----------
    fh : a filelike object
    encoding : an optional encoding
    """
    fh.seek(0)
    if encoding is not None:
        up = Unpickler(fh, encoding=encoding)
    else:
        up = Unpickler(fh)
    return up.load()


def loads(
    bytes_object: bytes,
    *,
    fix_imports: bool = True,
    encoding: str = "ASCII",
    errors: str = "strict",
):
    """
    Analogous to pickle._loads.
    """
    fd = io.BytesIO(bytes_object)
    return Unpickler(
        fd, fix_imports=fix_imports, encoding=encoding, errors=errors
    ).load()


@contextlib.contextmanager
def patch_pickle() -> Generator[None, None, None]:
    """
    Temporarily patch pickle to use our unpickler.
    """
    orig_loads = pickle.loads
    try:
        setattr(pickle, "loads", loads)
        yield
    finally:
        setattr(pickle, "loads", orig_loads)
