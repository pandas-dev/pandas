"""encoding.py - methods for reading parquet encoded data blocks."""
import numpy as np
from fastparquet.cencoding import read_bitpacked1, NumpyIO
from fastparquet.speedups import unpack_byte_array
from fastparquet import parquet_thrift

try:
    import pyarrow as _pa  # noqa: F401
    from fastparquet.speedups import unpack_byte_array_arrow as _unpack_arrow
    _HAVE_ARROW = True
except ImportError:
    _HAVE_ARROW = False

# Use the arrow string path only when pandas itself defaults to arrow strings.
# In pandas 3+, pd.options.future.infer_string is True by default.
# In pandas 2, it is False (object dtype is the default for strings).
# This can also be explicitly enabled in pandas 2 with
#   pd.options.future.infer_string = True
# so we respect that too.
try:
    import pandas as _pd
    _USE_ARROW_STRINGS = bool(getattr(getattr(_pd.options, 'future', None),
                                      'infer_string', False))
except Exception:
    _USE_ARROW_STRINGS = False


def read_plain_boolean(raw_bytes, count, out=None):
    data = np.frombuffer(raw_bytes, dtype='uint8')
    out = out or np.empty(count, dtype=bool)
    read_bitpacked1(NumpyIO(data), count, NumpyIO(out.view('uint8')))
    return out[:count]


DECODE_TYPEMAP = {
    parquet_thrift.Type.INT32: np.int32,
    parquet_thrift.Type.INT64: np.int64,
    parquet_thrift.Type.INT96: np.dtype('S12'),
    parquet_thrift.Type.FLOAT: np.float32,
    parquet_thrift.Type.DOUBLE: np.float64,
}


def read_plain(raw_bytes, type_, count, width=0, utf=False, stat=False):
    if type_ in DECODE_TYPEMAP:
        dtype = DECODE_TYPEMAP[type_]
        return np.frombuffer(memoryview(raw_bytes), dtype=dtype, count=count)
    if type_ == parquet_thrift.Type.FIXED_LEN_BYTE_ARRAY:
        if count == 1:
            width = len(raw_bytes)
        dtype = np.dtype('S%i' % width)
        return np.frombuffer(memoryview(raw_bytes), dtype=dtype, count=count)
    if type_ == parquet_thrift.Type.BOOLEAN:
        return read_plain_boolean(raw_bytes, count)
    if type_ == parquet_thrift.Type.BYTE_ARRAY:
        if stat:
            if utf:
                return np.array([bytes(raw_bytes).decode()], dtype='O')
            else:
                return np.array([bytes(raw_bytes)], dtype='O')
        if utf and _HAVE_ARROW and _USE_ARROW_STRINGS:
            # Build an ArrowStringArray directly from the packed bytes, without
            # creating intermediate Python str objects in the inner loop.
            # Only used when pandas defaults to arrow strings (pandas 3+, or
            # pandas 2 with future.infer_string=True).
            raw_np = np.frombuffer(memoryview(raw_bytes), dtype=np.uint8)
            return _unpack_arrow(raw_np, count)
        return unpack_byte_array(raw_bytes, count, utf=utf)
