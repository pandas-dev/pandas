from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from pandas.compat import HAS_PYARROW
from pandas.util._decorators import set_module

from pandas.core.dtypes.base import (
    ExtensionDtype,
    register_extension_dtype,
)
from pandas.core.dtypes.common import is_string_dtype
from pandas.core.dtypes.dtypes import ArrowDtype

from pandas.core.arrays.arrow.array import ArrowExtensionArray

if TYPE_CHECKING:
    from collections.abc import Sequence
    from pandas._typing import (
        type_t,
        ArrayLike,
        AstypeArg,
        DtypeObj,
        Shape,
    )

import re

import pyarrow as pa


def string_to_pyarrow_type(string: str) -> pa.DataType:
    # TODO: combine this with to_pyarrow_type in pandas.core.arrays.arrow ?
    pater = r"list\[(.*)\]"

    if mtch := re.search(pater, string):
        value_type = mtch.groups()[0]
        match value_type:
            # TODO: is there a pyarrow function get a type from the string?
            case "string" | "large_string":
                return pa.large_list(pa.large_string())
            case "int64":
                return pa.large_list(pa.int64())
            # TODO: need to implement many more here, including nested

    raise ValueError(f"Cannot map {string} to a pyarrow list type")


def transpose_homogeneous_list(
    arrays: Sequence[ListArray],
) -> list[ListArray]:
    # TODO: this is the same as transpose_homogeneous_pyarrow
    # but returns the ListArray instead of an ArrowExtensionArray
    # should consolidate these
    arrays = list(arrays)
    nrows, ncols = len(arrays[0]), len(arrays)
    indices = np.arange(nrows * ncols).reshape(ncols, nrows).T.reshape(-1)
    arr = pa.chunked_array([chunk for arr in arrays for chunk in arr._pa_array.chunks])
    arr = arr.take(indices)
    return [ListArray(arr.slice(i * ncols, ncols)) for i in range(nrows)]


@register_extension_dtype
@set_module("pandas")
class ListDtype(ArrowDtype):
    """
    An ExtensionDtype suitable for storing homogeneous lists of data.
    """

    _is_immutable = True  # TODO(wayd): should we allow mutability?

    def __init__(self, value_dtype: pa.DataType) -> None:
        super().__init__(pa.large_list(value_dtype))

    @classmethod
    def construct_from_string(cls, string: str):
        if not isinstance(string, str):
            raise TypeError(
                f"'construct_from_string' expects a string, got {type(string)}"
            )

        try:
            pa_type = string_to_pyarrow_type(string)
        except ValueError as e:
            raise TypeError(
                f"Cannot construct a '{cls.__name__}' from '{string}'"
            ) from e

        return cls(pa_type)

    @property
    def name(self) -> str:  # type: ignore[override]
        """
        A string identifying the data type.
        """
        # TODO: reshaping tests require the name list to match the large_list
        # implementation; assumedly there are some astype(str(dtype)) casts
        # going on. Should fix so this can just be "list[...]" for end user
        return f"large_list[{self.pyarrow_dtype.value_type!s}]"

    @property
    def kind(self) -> str:
        # TODO(wayd): our extension interface says this field should be the
        # NumPy type character, but no such thing exists for list
        # This uses the Arrow C Data exchange code instead
        return "+L"

    @classmethod
    def construct_array_type(cls) -> type_t[ListArray]:
        """
        Return the array type associated with this dtype.

        Returns
        -------
        type
        """
        return ListArray

    def _get_common_dtype(self, dtypes: list[DtypeObj]) -> DtypeObj | None:
        # TODO(wayd): should we implemented value type support?
        for dtype in dtypes:
            if (
                isinstance(dtype, ListDtype)
                and self.pyarrow_dtype.value_type == dtype.pyarrow_dtype.value_type
            ):
                continue
            else:
                return None

        return ListDtype(self.pyarrow_dtype.value_type)


class ListArray(ArrowExtensionArray):
    __array_priority__ = 1000

    def __init__(
        self, values: pa.Array | pa.ChunkedArray | list | ListArray, value_type=None
    ) -> None:
        if not HAS_PYARROW:
            raise NotImplementedError("ListArray requires pyarrow to be installed")

        if isinstance(values, type(self)):
            self._pa_array = values._pa_array
        else:
            if value_type is None:
                if isinstance(values, (pa.Array, pa.ChunkedArray)):
                    value_type = values.type.value_type
                else:
                    value_type = pa.array(values).type.value_type

                # Internally always use large_string instead of string
                if value_type == pa.string():
                    value_type = pa.large_string()

            if not isinstance(values, pa.ChunkedArray):
                # To support NA, we need to create an Array first :-(
                arr = pa.array(values, type=pa.large_list(value_type), from_pandas=True)
                self._pa_array = pa.chunked_array(arr, type=pa.large_list(value_type))
            else:
                self._pa_array = values

    @property
    def _dtype(self):
        return ListDtype(self._pa_array.type.value_type)

    @classmethod
    def _from_sequence(cls, scalars, *, dtype=None, copy: bool = False):
        if isinstance(scalars, ListArray):
            return cls(scalars)
        elif isinstance(scalars, pa.Scalar):
            scalars = [scalars]
            return cls(scalars)

        try:
            values = pa.array(scalars, from_pandas=True)
        except TypeError:
            # TypeError: object of type 'NoneType' has no len() if you have
            # pa.ListScalar(None). Upstream issue in Arrow - see:
            # https://github.com/apache/arrow/issues/40319
            for i in range(len(scalars)):
                if not scalars[i].is_valid:
                    scalars[i] = None

            values = pa.array(scalars, from_pandas=True)

        if values.type == "null" and dtype is not None:
            # TODO: the sequencing here seems wrong; just making the tests pass for now
            # but this needs a comprehensive review
            pa_type = string_to_pyarrow_type(str(dtype))
            values = pa.array(values, type=pa_type)

        return cls(values)

    def __getitem__(self, item):
        # PyArrow does not support NumPy's selection with an equal length
        # mask, so let's convert those to integral positions if needed
        if isinstance(item, np.ndarray) and item.dtype == bool:
            pos = np.array(range(len(item)))
            mask = pos[item]
            return type(self)(self._pa_array.take(mask))
        elif isinstance(item, int):
            return self._pa_array[item]
        elif isinstance(item, list):
            return type(self)(self._pa_array.take(item))

        return type(self)(self._pa_array[item])

    @classmethod
    def _empty(cls, shape: Shape, dtype: ExtensionDtype):
        """
        Create an ExtensionArray with the given shape and dtype.

        See also
        --------
        ExtensionDtype.empty
            ExtensionDtype.empty is the 'official' public version of this API.
        """
        # Implementer note: while ExtensionDtype.empty is the public way to
        # call this method, it is still required to implement this `_empty`
        # method as well (it is called internally in pandas)
        if isinstance(shape, tuple):
            if len(shape) > 1:
                raise ValueError("ListArray may only be 1-D")
            else:
                length = shape[0]
        else:
            length = shape

        return cls._from_sequence([None] * length, dtype=dtype)

    def astype(self, dtype: AstypeArg, copy: bool = True) -> ArrayLike:
        if is_string_dtype(dtype) and not isinstance(dtype, ExtensionDtype):
            return np.array([str(x) for x in self], dtype=dtype)

        return super().astype(dtype, copy)

    def __eq__(self, other):
        if isinstance(other, (pa.ListScalar, pa.LargeListScalar)):
            from pandas.arrays import BooleanArray

            # TODO: pyarrow.compute does not implement broadcasting equality
            # for an array of lists to a listscalar
            # TODO: pyarrow doesn't compare missing values as missing???
            # arr = pa.array([1, 2, None])
            # pc.equal(arr, arr[2]) returns all nulls but
            # arr[2] == arr[2] returns True
            mask = np.array([False] * len(self))
            values = np.array([x == other for x in self._pa_array])
            return BooleanArray(values, mask)

        return super().__eq__(other)
