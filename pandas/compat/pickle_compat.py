"""
Pickle compatibility to pandas version 1.0
"""

from __future__ import annotations

import contextlib
import io
import pickle
from typing import (
    TYPE_CHECKING,
    Any,
)

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
    # Re-routing unpickle block logic to go through _unpickle_block instead
    # for pandas <= 1.3.5
    ("pandas.core.internals.blocks", "new_block"): (
        "pandas._libs.internals",
        "_unpickle_block",
    ),
    # Avoid Cython's warning "contradiction to Python 'class private name' rules"
    ("pandas._libs.tslibs.nattype", "__nat_unpickle"): (
        "pandas._libs.tslibs.nattype",
        "_nat_unpickle",
    ),
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


# our Unpickler sub-class to override methods and some dispatcher
# functions for compat and uses a non-public class of the pickle module.
class Unpickler(pickle._Unpickler):
    def find_class(self, module: str, name: str) -> Any:
        key = (module, name)
        module, name = _class_locations_map.get(key, key)
        return super().find_class(module, name)

    dispatch = pickle._Unpickler.dispatch.copy()

    def load_reduce(self) -> None:
        stack = self.stack  # type: ignore[attr-defined]
        args = stack.pop()
        func = stack[-1]

        try:
            stack[-1] = func(*args)
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

    dispatch[pickle.REDUCE[0]] = load_reduce  # type: ignore[assignment]

    def load_newobj(self) -> None:
        args = self.stack.pop()  # type: ignore[attr-defined]
        cls = self.stack.pop()  # type: ignore[attr-defined]

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
        self.append(obj)  # type: ignore[attr-defined]

    dispatch[pickle.NEWOBJ[0]] = load_newobj  # type: ignore[assignment]


def loads(
    bytes_object: bytes,
    *,
    fix_imports: bool = True,
    encoding: str = "ASCII",
    errors: str = "strict",
) -> Any:
    """
    Analogous to pickle._loads.
    """
    fd = io.BytesIO(bytes_object)
    return Unpickler(
        fd, fix_imports=fix_imports, encoding=encoding, errors=errors
    ).load()


@contextlib.contextmanager
def patch_pickle() -> Generator[None]:
    """
    Temporarily patch pickle to use our unpickler.
    """
    orig_loads = pickle.loads
    try:
        setattr(pickle, "loads", loads)
        yield
    finally:
        setattr(pickle, "loads", orig_loads)
