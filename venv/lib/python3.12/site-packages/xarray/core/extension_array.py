from __future__ import annotations

import copy
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any, Generic, cast

import numpy as np
import pandas as pd
from packaging.version import Version

from xarray.core.types import DTypeLikeSave, T_ExtensionArray
from xarray.core.utils import NDArrayMixin, is_allowed_extension_array

HANDLED_EXTENSION_ARRAY_FUNCTIONS: dict[Callable, Callable] = {}


def implements(numpy_function):
    """Register an __array_function__ implementation for MyArray objects."""

    def decorator(func):
        HANDLED_EXTENSION_ARRAY_FUNCTIONS[numpy_function] = func
        return func

    return decorator


@implements(np.issubdtype)
def __extension_duck_array__issubdtype(
    extension_array_dtype: T_ExtensionArray, other_dtype: DTypeLikeSave
) -> bool:
    return False  # never want a function to think a pandas extension dtype is a subtype of numpy


@implements(np.broadcast_to)
def __extension_duck_array__broadcast(arr: T_ExtensionArray, shape: tuple):
    if shape[0] == len(arr) and len(shape) == 1:
        return arr
    raise NotImplementedError("Cannot broadcast 1d-only pandas extension array.")


@implements(np.stack)
def __extension_duck_array__stack(arr: T_ExtensionArray, axis: int):
    raise NotImplementedError("Cannot stack 1d-only pandas extension array.")


@implements(np.concatenate)
def __extension_duck_array__concatenate(
    arrays: Sequence[T_ExtensionArray], axis: int = 0, out=None
) -> T_ExtensionArray:
    return type(arrays[0])._concat_same_type(arrays)  # type: ignore[attr-defined]


@implements(np.where)
def __extension_duck_array__where(
    condition: np.ndarray, x: T_ExtensionArray, y: T_ExtensionArray
) -> T_ExtensionArray:
    if (
        isinstance(x, pd.Categorical)
        and isinstance(y, pd.Categorical)
        and x.dtype != y.dtype
    ):
        x = x.add_categories(set(y.categories).difference(set(x.categories)))  # type: ignore[assignment]
        y = y.add_categories(set(x.categories).difference(set(y.categories)))  # type: ignore[assignment]
    return cast(T_ExtensionArray, pd.Series(x).where(condition, pd.Series(y)).array)


@implements(np.ndim)
def __extension_duck_array__ndim(x: PandasExtensionArray) -> int:
    return x.ndim


@implements(np.reshape)
def __extension_duck_array__reshape(
    arr: T_ExtensionArray, shape: tuple
) -> T_ExtensionArray:
    if (shape[0] == len(arr) and len(shape) == 1) or shape == (-1,):
        return arr
    raise NotImplementedError(
        f"Cannot reshape 1d-only pandas extension array to: {shape}"
    )


@dataclass(frozen=True)
class PandasExtensionArray(NDArrayMixin, Generic[T_ExtensionArray]):
    """NEP-18 compliant wrapper for pandas extension arrays.

    Parameters
    ----------
    array : T_ExtensionArray
        The array to be wrapped upon e.g,. :py:class:`xarray.Variable` creation.
    ```
    """

    array: T_ExtensionArray

    def __post_init__(self):
        if not isinstance(self.array, pd.api.extensions.ExtensionArray):
            raise TypeError(f"{self.array} is not an pandas ExtensionArray.")
        # This does not use the UNSUPPORTED_EXTENSION_ARRAY_TYPES whitelist because
        # we do support extension arrays from datetime, for example, that need
        # duck array support internally via this class.  These can appear from `DatetimeIndex`
        # wrapped by `PandasIndex` internally, for example.
        if not is_allowed_extension_array(self.array):
            raise TypeError(
                f"{self.array.dtype!r} should be converted to a numpy array in `xarray` internally."
            )

    def __array_function__(self, func, types, args, kwargs):
        def replace_duck_with_extension_array(args) -> list:
            args_as_list = list(args)
            for index, value in enumerate(args_as_list):
                if isinstance(value, PandasExtensionArray):
                    args_as_list[index] = value.array
                elif isinstance(
                    value, tuple
                ):  # should handle more than just tuple? iterable?
                    args_as_list[index] = tuple(
                        replace_duck_with_extension_array(value)
                    )
                elif isinstance(value, list):
                    args_as_list[index] = replace_duck_with_extension_array(value)
            return args_as_list

        args = tuple(replace_duck_with_extension_array(args))
        if func not in HANDLED_EXTENSION_ARRAY_FUNCTIONS:
            raise KeyError("Function not registered for pandas extension arrays.")
        res = HANDLED_EXTENSION_ARRAY_FUNCTIONS[func](*args, **kwargs)
        if is_allowed_extension_array(res):
            return PandasExtensionArray(res)
        return res

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        return ufunc(*inputs, **kwargs)

    def __getitem__(self, key) -> PandasExtensionArray[T_ExtensionArray]:
        item = self.array[key]
        if is_allowed_extension_array(item):
            return PandasExtensionArray(item)
        if np.isscalar(item) or isinstance(key, int):
            return PandasExtensionArray(type(self.array)._from_sequence([item]))  # type: ignore[call-arg,attr-defined,unused-ignore]
        return PandasExtensionArray(item)

    def __setitem__(self, key, val):
        self.array[key] = val

    def __eq__(self, other):
        if isinstance(other, PandasExtensionArray):
            return self.array == other.array
        return self.array == other

    def __ne__(self, other):
        return ~(self == other)

    def __len__(self):
        return len(self.array)

    @property
    def ndim(self) -> int:
        return 1

    def __array__(
        self, dtype: np.typing.DTypeLike = None, /, *, copy: bool | None = None
    ) -> np.ndarray:
        if Version(np.__version__) >= Version("2.0.0"):
            return np.asarray(self.array, dtype=dtype, copy=copy)
        else:
            return np.asarray(self.array, dtype=dtype)

    def __getattr__(self, attr: str) -> Any:
        #  with __deepcopy__ or __copy__, the object is first constructed and then the sub-objects are attached (see https://docs.python.org/3/library/copy.html)
        # Thus, if we didn't have `super().__getattribute__("array")` this method would call `self.array` (i.e., `getattr(self, "array")`) again while looking for `__setstate__`
        # (which is apparently the first thing sought in copy.copy from the under-construction copied object),
        # which would cause a recursion error since `array` is not present on the object when it is being constructed during `__{deep}copy__`.
        # Even though we have defined these two methods now below due to `test_extension_array_copy_arrow_type` (cause unknown)
        # we leave this here as it more robust than self.array
        return getattr(super().__getattribute__("array"), attr)

    def __copy__(self) -> PandasExtensionArray[T_ExtensionArray]:
        return PandasExtensionArray(copy.copy(self.array))

    def __deepcopy__(
        self, memo: dict[int, Any] | None = None
    ) -> PandasExtensionArray[T_ExtensionArray]:
        return PandasExtensionArray(copy.deepcopy(self.array, memo=memo))
