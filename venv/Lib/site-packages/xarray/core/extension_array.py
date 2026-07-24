from __future__ import annotations

import copy
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Generic, cast

import numpy as np
import pandas as pd
from packaging.version import Version
from pandas.api.extensions import ExtensionArray, ExtensionDtype
from pandas.api.types import is_scalar as pd_is_scalar

from xarray.core.types import DTypeLikeSave, T_ExtensionArray
from xarray.core.utils import (
    NDArrayMixin,
    is_allowed_extension_array,
    is_allowed_extension_array_dtype,
)

HANDLED_EXTENSION_ARRAY_FUNCTIONS: dict[Callable, Callable] = {}


if TYPE_CHECKING:
    from typing import Any

    from pandas._typing import DtypeObj, Scalar


def is_scalar(value: object) -> bool:
    """Workaround: pandas is_scalar doesn't recognize Categorical nulls for some reason."""
    return value is pd.CategoricalDtype.na_value or pd_is_scalar(value)


def implements(numpy_function_or_name: Callable | str) -> Callable:
    """Register an __array_function__ implementation.

    Pass a function directly if it's guaranteed to exist in all supported numpy versions, or a
    string to first check for its existence.
    """

    def decorator(func):
        if isinstance(numpy_function_or_name, str):
            numpy_function = getattr(np, numpy_function_or_name, None)
        else:
            numpy_function = numpy_function_or_name

        if numpy_function:
            HANDLED_EXTENSION_ARRAY_FUNCTIONS[numpy_function] = func
        return func

    return decorator


@implements(np.issubdtype)
def __extension_duck_array__issubdtype(
    extension_array_dtype: T_ExtensionArray, other_dtype: DTypeLikeSave
) -> bool:
    return False  # never want a function to think a pandas extension dtype is a subtype of numpy


@implements("astype")  # np.astype was added in 2.1.0, but we only require >=1.24
def __extension_duck_array__astype(
    array_or_scalar: T_ExtensionArray,
    dtype: DTypeLikeSave,
    order: str = "K",
    casting: str = "unsafe",
    subok: bool = True,
    copy: bool = True,
    device: str | None = None,
) -> ExtensionArray:
    if (
        not (
            is_allowed_extension_array(array_or_scalar)
            or is_allowed_extension_array_dtype(dtype)
        )
        or casting != "unsafe"
        or not subok
        or order != "K"
    ):
        return NotImplemented

    return as_extension_array(array_or_scalar, dtype, copy=copy)


@implements(np.asarray)
def __extension_duck_array__asarray(
    array_or_scalar: np.typing.ArrayLike | T_ExtensionArray,
    dtype: DTypeLikeSave | None = None,
) -> ExtensionArray:
    if not is_allowed_extension_array(dtype):
        return NotImplemented

    return as_extension_array(array_or_scalar, dtype)


def as_extension_array(
    array_or_scalar: np.typing.ArrayLike | T_ExtensionArray,
    dtype: ExtensionDtype | DTypeLikeSave | None,
    copy: bool = False,
) -> ExtensionArray:
    if is_scalar(array_or_scalar):
        return dtype.construct_array_type()._from_sequence(  # type: ignore[union-attr]
            [array_or_scalar], dtype=dtype
        )
    else:
        # pandas-stubs is overly strict about astype's dtype parameter and return type;
        # ExtensionArray.astype accepts ExtensionDtype and returns ExtensionArray
        return array_or_scalar.astype(dtype, copy=copy)  # type: ignore[union-attr,return-value,arg-type]


@implements(np.result_type)
def __extension_duck_array__result_type(
    *arrays_and_dtypes: list[
        np.typing.ArrayLike | np.typing.DTypeLike | ExtensionDtype | ExtensionArray
    ],
) -> DtypeObj:
    extension_arrays_and_dtypes: list[ExtensionDtype | ExtensionArray] = [
        cast(ExtensionDtype | ExtensionArray, x)
        for x in arrays_and_dtypes
        if is_allowed_extension_array(x) or is_allowed_extension_array_dtype(x)
    ]
    if not extension_arrays_and_dtypes:
        return NotImplemented

    ea_dtypes: list[ExtensionDtype] = [
        getattr(x, "dtype", cast(ExtensionDtype, x))
        for x in extension_arrays_and_dtypes
    ]
    scalars = [
        x for x in arrays_and_dtypes if is_scalar(x) and x not in {pd.NA, np.nan}
    ]
    # other_stuff could include:
    # - arrays such as pd.ABCSeries, np.ndarray, or other array-api duck arrays
    # - dtypes such as pd.DtypeObj, np.dtype, or other array-api duck dtypes
    other_stuff = [
        x
        for x in arrays_and_dtypes
        if not is_allowed_extension_array_dtype(x) and not is_scalar(x)
    ]
    # We implement one special case: when possible, preserve Categoricals (avoid promoting
    # to object) by merging the categories of all given Categoricals + scalars + NA.
    # Ideally this could be upstreamed into pandas find_result_type / find_common_type.
    if not other_stuff and all(
        isinstance(x, pd.CategoricalDtype) and not x.ordered for x in ea_dtypes
    ):
        return union_unordered_categorical_and_scalar(
            cast(list[pd.CategoricalDtype], ea_dtypes),
            scalars,  # type: ignore[arg-type]
        )
    if not other_stuff and all(
        isinstance(x, type(ea_type := ea_dtypes[0])) for x in ea_dtypes
    ):
        return ea_type
    raise ValueError(
        f"Cannot cast values to shared type, found values: {arrays_and_dtypes}"
    )


def union_unordered_categorical_and_scalar(
    categorical_dtypes: list[pd.CategoricalDtype], scalars: list[Scalar]
) -> pd.CategoricalDtype:
    scalars = [x for x in scalars if x is not pd.CategoricalDtype.na_value]
    all_categories = set().union(*(x.categories for x in categorical_dtypes))
    all_categories = all_categories.union(scalars)
    return pd.CategoricalDtype(categories=list(all_categories))


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
    condition: T_ExtensionArray | np.typing.ArrayLike,
    x: T_ExtensionArray,
    y: T_ExtensionArray | np.typing.ArrayLike,
) -> T_ExtensionArray:
    # pd.where won't broadcast 0-dim arrays across a scalar-like series; scalar y's must be preserved
    if hasattr(y, "shape") and len(y.shape) == 1 and y.shape[0] == 1:
        y = y[0]  # type: ignore[index]
    # pandas-stubs has strict overloads for Series.where that don't cover all valid arg types
    return cast(T_ExtensionArray, pd.Series(x).where(condition, y).array)  # type: ignore[call-overload]


def _replace_duck(args, replacer: Callable[[PandasExtensionArray], Any]) -> list:
    args_as_list = list(args)
    for index, value in enumerate(args_as_list):
        if isinstance(value, PandasExtensionArray):
            args_as_list[index] = replacer(value)
        elif isinstance(value, tuple):  # should handle more than just tuple? iterable?
            args_as_list[index] = tuple(_replace_duck(value, replacer))
        elif isinstance(value, list):
            args_as_list[index] = _replace_duck(value, replacer)
    return args_as_list


def replace_duck_with_extension_array(args) -> tuple:
    return tuple(_replace_duck(args, lambda duck: duck.array))


def replace_duck_with_series(args) -> tuple:
    return tuple(_replace_duck(args, lambda duck: pd.Series(duck.array)))


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
            raise TypeError(f"{self.array} is not a pandas ExtensionArray.")
        # This does not use the UNSUPPORTED_EXTENSION_ARRAY_TYPES whitelist because
        # we do support extension arrays from datetime, for example, that need
        # duck array support internally via this class.  These can appear from `DatetimeIndex`
        # wrapped by `PandasIndex` internally, for example.
        if not is_allowed_extension_array(self.array):
            raise TypeError(
                f"{self.array.dtype!r} should be converted to a numpy array in `xarray` internally."
            )

    def __array_function__(self, func, types, args, kwargs):
        if func not in HANDLED_EXTENSION_ARRAY_FUNCTIONS:
            raise KeyError("Function not registered for pandas extension arrays.")
        args = replace_duck_with_extension_array(args)
        res = HANDLED_EXTENSION_ARRAY_FUNCTIONS[func](*args, **kwargs)
        if isinstance(res, ExtensionArray):
            return PandasExtensionArray(res)
        return res

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        return ufunc(*inputs, **kwargs)

    def __getitem__(self, key) -> PandasExtensionArray[T_ExtensionArray]:
        if (
            isinstance(key, tuple) and len(key) == 1
        ):  # pyarrow type arrays can't handle single-length tuples
            (key,) = key
        item = self.array[key]
        if is_allowed_extension_array(item):
            return PandasExtensionArray(item)
        if is_scalar(item) or isinstance(key, int):
            return PandasExtensionArray(type(self.array)._from_sequence([item]))  # type: ignore[call-arg,attr-defined,unused-ignore]
        return PandasExtensionArray(item)

    def __setitem__(self, key, val):
        self.array[key] = val

    def __len__(self):
        return len(self.array)

    def __eq__(self, other):
        if isinstance(other, PandasExtensionArray):
            return self.array == other.array
        return self.array == other

    def __ne__(self, other):
        return ~(self == other)

    @property
    def ndim(self) -> int:
        return 1

    def __array__(
        self, dtype: np.typing.DTypeLike | None = None, /, *, copy: bool | None = None
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
