from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from pandas.compat import HAS_PYARROW
from pandas.util._decorators import set_module

from pandas.core.dtypes.base import (
    ExtensionDtype,
    register_extension_dtype,
)
from pandas.core.dtypes.common import (
    is_bool_dtype,
    is_integer_dtype,
    is_string_dtype,
)
from pandas.core.dtypes.dtypes import ArrowDtype

from pandas.core.arrays.arrow.array import ArrowExtensionArray
from pandas.core.arrays.base import ExtensionArray

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

    _is_immutable = True

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
        return f"list[{self.pyarrow_dtype.value_type!s}]"

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
                    parent_type = values.type
                    if not isinstance(parent_type, (pa.ListType, pa.LargeListType)):
                        # TODO: maybe implement native casts in pyarrow
                        new_values = [
                            [x.as_py()] if x.is_valid else None for x in values
                        ]
                        values = pa.array(new_values, type=pa.large_list(parent_type))

                    value_type = values.type.value_type
                else:
                    value_type = pa.array(values).type.value_type

                if value_type == pa.string():
                    value_type = pa.large_string()

            if not isinstance(values, pa.ChunkedArray):
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
            values = pa.array(scalars.to_pylist(), from_pandas=True)

        if values.type == "null" and dtype is not None:
            pa_type = string_to_pyarrow_type(str(dtype))
            values = pa.array(values, type=pa_type)

        return cls(values)

    @classmethod
    def _box_pa(
        cls, value, pa_type: pa.DataType | None = None
    ) -> pa.Array | pa.ChunkedArray | pa.Scalar:
        """
        Box value into a pyarrow Array, ChunkedArray or Scalar.

        Parameters
        ----------
        value : any
        pa_type : pa.DataType | None

        Returns
        -------
        pa.Array or pa.ChunkedArray or pa.Scalar
        """
        if (
            isinstance(value, (pa.ListScalar, pa.LargeListScalar))
            or isinstance(value, list)
            or value is None
        ):
            return cls._box_pa_scalar(value, pa_type)
        return cls._box_pa_array(value, pa_type)

    def __getitem__(self, item):
        if isinstance(item, (np.ndarray, ExtensionArray)):
            if is_bool_dtype(item.dtype):
                mask_len = len(item)
                if mask_len != len(self):
                    raise IndexError(
                        f"Boolean index has wrong length: {mask_len} "
                        f"instead of {len(self)}"
                    )
                pos = np.array(range(len(item)))

                if isinstance(item, ExtensionArray):
                    mask = pos[item.fillna(False)]
                else:
                    mask = pos[item]
                return type(self)(self._pa_array.take(mask))
            elif is_integer_dtype(item.dtype):
                if isinstance(item, ExtensionArray) and item.isna().any():
                    msg = "Cannot index with an integer indexer containing NA values"
                    raise ValueError(msg)

                indexer = pa.array(item)
                return type(self)(self._pa_array.take(indexer))
        elif isinstance(item, int):
            value = self._pa_array[item]
            if value.is_valid:
                return value.as_py()
            else:
                return self.dtype.na_value
        elif isinstance(item, list):
            # pyarrow does not support taking yet from an empty list
            # https://github.com/apache/arrow/issues/39917
            if item:
                try:
                    result = self._pa_array.take(item)
                except pa.lib.ArrowInvalid as e:
                    if "Could not convert <NA>" in str(e):
                        msg = (
                            "Cannot index with an integer indexer containing NA values"
                        )
                        raise ValueError(msg) from e
                    raise e
            else:
                result = pa.array([], type=self._pa_array.type)

            return type(self)(result)

        try:
            result = type(self)(self._pa_array[item])
        except TypeError as e:
            msg = (
                "only integers, slices (`:`), ellipsis (`...`), numpy.newaxis "
                "(`None`) and integer or boolean arrays are valid indices"
            )
            raise IndexError(msg) from e

        return result

    def __setitem__(self, key, value) -> None:
        msg = "ListArray does not support item assignment via setitem"
        raise TypeError(msg)

    @classmethod
    def _empty(cls, shape: Shape, dtype: ExtensionDtype):
        """
        Create an ExtensionArray with the given shape and dtype.

        See also
        --------
        ExtensionDtype.empty
            ExtensionDtype.empty is the 'official' public version of this API.
        """
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
        if isinstance(other, list):
            from pandas.arrays import BooleanArray

            mask = np.array([False] * len(self))
            values = np.array([x.as_py() == other for x in self._pa_array])
            return BooleanArray(values, mask)
        elif isinstance(other, (pa.ListScalar, pa.LargeListScalar)):
            from pandas.arrays import BooleanArray

            # TODO: pyarrow.compute does not implement equal for lists
            # https://github.com/apache/arrow/issues/45167
            # TODO: pyarrow doesn't compare missing values in Python as missing???
            # arr = pa.array([1, 2, None])
            # pc.equal(arr, arr[2]) returns all nulls but
            # arr[2] == arr[2] returns True
            mask = np.array([False] * len(self))
            values = np.array([x == other for x in self._pa_array])
            return BooleanArray(values, mask)

        return super().__eq__(other)
