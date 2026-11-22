from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from pandas.core.arrays.boolean import BooleanDtype
from pandas.core.arrays.floating import (
    Float32Dtype,
    Float64Dtype,
)
from pandas.core.arrays.integer import (
    Int8Dtype,
    Int16Dtype,
    Int32Dtype,
    Int64Dtype,
    UInt8Dtype,
    UInt16Dtype,
    UInt32Dtype,
    UInt64Dtype,
)
from pandas.core.arrays.string_ import StringDtype

if TYPE_CHECKING:
    import pyarrow

    from pandas.core.dtypes.dtypes import BaseMaskedDtype


def _arrow_dtype_mapping() -> dict[pyarrow.DataType, BaseMaskedDtype | StringDtype]:
    import pyarrow as pa

    return {
        pa.int8(): Int8Dtype(),
        pa.int16(): Int16Dtype(),
        pa.int32(): Int32Dtype(),
        pa.int64(): Int64Dtype(),
        pa.uint8(): UInt8Dtype(),
        pa.uint16(): UInt16Dtype(),
        pa.uint32(): UInt32Dtype(),
        pa.uint64(): UInt64Dtype(),
        pa.bool_(): BooleanDtype(),
        pa.string(): StringDtype(),
        pa.float32(): Float32Dtype(),
        pa.float64(): Float64Dtype(),
        pa.large_string(): StringDtype(),
    }


def pyarrow_array_to_numpy_and_mask(
    arr, dtype: np.dtype
) -> tuple[np.ndarray, np.ndarray]:
    """
    Convert a primitive pyarrow.Array to a numpy array and boolean mask based
    on the buffers of the Array.

    At the moment pyarrow.BooleanArray is not supported.

    Parameters
    ----------
    arr : pyarrow.Array
    dtype : numpy.dtype

    Returns
    -------
    (data, mask)
        Tuple of two numpy arrays with the raw data (with specified dtype) and
        a boolean mask (validity mask, so False means missing)
    """
    import pyarrow

    dtype = np.dtype(dtype)

    if pyarrow.types.is_null(arr.type):
        # No initialization of data is needed since everything is null
        data = np.empty(len(arr), dtype=dtype)
        mask = np.zeros(len(arr), dtype=bool)
        return data, mask
    buflist = arr.buffers()
    # Since Arrow buffers might contain padding and the data might be offset,
    # the buffer gets sliced here before handing it to numpy.
    # See also https://github.com/pandas-dev/pandas/issues/40896
    offset = arr.offset * dtype.itemsize
    length = len(arr) * dtype.itemsize
    data_buf = buflist[1][offset : offset + length]
    data = np.frombuffer(data_buf, dtype=dtype)
    bitmask = buflist[0]
    if bitmask is not None:
        mask = pyarrow.BooleanArray.from_buffers(
            pyarrow.bool_(), len(arr), [None, bitmask], offset=arr.offset
        )
        mask = np.asarray(mask)
    else:
        mask = np.ones(len(arr), dtype=bool)
    return data, mask
