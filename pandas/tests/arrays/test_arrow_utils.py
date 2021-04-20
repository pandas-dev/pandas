import numpy as np

import pandas.util._test_decorators as td

import pandas._testing as tm


@td.skip_if_no("pyarrow")
def test_pyarrow_array_to_numpy_and_mask():
    """
    Test conversion from pyarrow array to numpy array.

    Also modifies the pyarrow buffer to contain padding and offset, which are
    considered valid buffers by pyarrow.
    See also https://github.com/pandas-dev/pandas/issues/40896
    """
    import pyarrow as pa

    from pandas.core.arrays._arrow_utils import pyarrow_array_to_numpy_and_mask

    dtype = np.int32
    pa_array = pa.array([0, 1, 2], type=pa.int32())
    np_expected = np.array([0, 1, 2], dtype=np.int32)
    mask_expected = np.array([True, True, True])

    data, mask = pyarrow_array_to_numpy_and_mask(pa_array, dtype)
    tm.assert_numpy_array_equal(data, np_expected)
    tm.assert_numpy_array_equal(mask, mask_expected)

    mask_buffer = pa_array.buffers()[0]
    data_buffer = pa_array.buffers()[1].to_pybytes()

    # Add trailing padding to the buffer.
    data_buffer_trail = pa.py_buffer(data_buffer + b"\x00")
    pa_array_trail = pa.Array.from_buffers(
        type=pa_array.type,
        length=len(pa_array),
        buffers=[mask_buffer, data_buffer_trail],
        offset=pa_array.offset,
    )
    pa_array_trail.validate()
    data, mask = pyarrow_array_to_numpy_and_mask(pa_array_trail, dtype)
    tm.assert_numpy_array_equal(data, np_expected)
    tm.assert_numpy_array_equal(mask, mask_expected)

    # Add offset to the buffer.
    offset = b"\x00" * (pa_array.type.bit_width // 8)
    data_buffer_offset = pa.py_buffer(offset + data_buffer)
    mask_buffer_offset = pa.py_buffer(b"\x0F")
    pa_array_offset = pa.Array.from_buffers(
        type=pa_array.type,
        length=len(pa_array),
        buffers=[mask_buffer_offset, data_buffer_offset],
        offset=pa_array.offset + 1,
    )
    pa_array_offset.validate()
    data, mask = pyarrow_array_to_numpy_and_mask(pa_array_offset, dtype)
    tm.assert_numpy_array_equal(data, np_expected)
    tm.assert_numpy_array_equal(mask, mask_expected)
