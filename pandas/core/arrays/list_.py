from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    ClassVar,
)

import numpy as np

from pandas._libs import missing as libmissing
from pandas.compat import HAS_PYARROW
from pandas.util._decorators import set_module

from pandas.core.dtypes.base import (
    ExtensionDtype,
    register_extension_dtype,
)
from pandas.core.dtypes.common import (
    is_object_dtype,
    is_string_dtype,
)

from pandas.core.arrays import ExtensionArray

if TYPE_CHECKING:
    from pandas._typing import (
        type_t,
        Shape,
    )

import pyarrow as pa


@register_extension_dtype
@set_module("pandas")
class ListDtype(ExtensionDtype):
    """
    An ExtensionDtype suitable for storing homogeneous lists of data.
    """

    type = list
    name: ClassVar[str] = "list"

    @property
    def na_value(self) -> libmissing.NAType:
        return libmissing.NA

    @property
    def kind(self) -> str:
        # TODO: our extension interface says this field should be the
        # NumPy type character, but no such thing exists for list
        # this assumes a PyArrow large list
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


class ListArray(ExtensionArray):
    dtype = ListDtype()
    __array_priority__ = 1000

    def __init__(self, values: pa.Array | pa.ChunkedArray | list | ListArray) -> None:
        if not HAS_PYARROW:
            raise NotImplementedError("ListArray requires pyarrow to be installed")

        if isinstance(values, type(self)):
            self._pa_array = values._pa_array
        elif not isinstance(values, pa.ChunkedArray):
            # To support NA, we need to create an Array first :-(
            arr = pa.array(values, from_pandas=True)
            self._pa_array = pa.chunked_array(arr)
        else:
            self._pa_array = values

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
        if values.type == "null":
            # TODO(wayd): this is a hack to get the tests to pass, but the overall issue
            # is that our extension types don't support parametrization but the pyarrow
            values = pa.array(values, type=pa.list_(pa.null()))

        return cls(values)

    def __getitem__(self, item):
        # PyArrow does not support NumPy's selection with an equal length
        # mask, so let's convert those to integral positions if needed
        if isinstance(item, np.ndarray) and item.dtype == bool:
            pos = np.array(range(len(item)))
            mask = pos[item]
            return type(self)(self._pa_array.take(mask))
        elif isinstance(item, int):  # scalar case
            return self._pa_array[item]

        return type(self)(self._pa_array[item])

    def __len__(self) -> int:
        return len(self._pa_array)

    def isna(self):
        return np.array(self._pa_array.is_null())

    def take(self, indexer, allow_fill=False, fill_value=None):
        # TODO: what do we need to do with allow_fill and fill_value here?
        return type(self)(self._pa_array.take(indexer))

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
        return cls._from_sequence([None] * length, dtype=pa.list_(pa.null()))

    def copy(self):
        mm = pa.default_cpu_memory_manager()

        # TODO(wayd): ChunkedArray does not implement copy_to so this
        # ends up creating an Array
        copied = self._pa_array.combine_chunks().copy_to(mm.device)
        return type(self)(copied)

    def astype(self, dtype, copy=True):
        if isinstance(dtype, type(self.dtype)) and dtype == self.dtype:
            if copy:
                return self.copy()
            return self
        elif is_string_dtype(dtype) and not is_object_dtype(dtype):
            # numpy has problems with astype(str) for nested elements
            # and pyarrow cannot cast from list[string] to string
            return np.array([str(x) for x in self._pa_array], dtype=dtype)

        if not copy:
            raise TypeError(f"astype from ListArray to {dtype} requires a copy")

        return np.array(self._pa_array.to_pylist(), dtype=dtype, copy=copy)

    @classmethod
    def _concat_same_type(cls, to_concat):
        data = [x._pa_array for x in to_concat]
        return cls(data)
