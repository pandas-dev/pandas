import numpy as np
import pytest

from pandas._libs import index as libindex


@pytest.fixture(params=[
    (libindex.Int64Engine, np.int64),
    (libindex.Int32Engine, np.int32),
    (libindex.Int16Engine, np.int16),
    (libindex.Int8Engine, np.int8),
    (libindex.UInt64Engine, np.uint64),
    (libindex.UInt32Engine, np.uint32),
    (libindex.UInt16Engine, np.uint16),
    (libindex.UInt8Engine, np.uint8),
    (libindex.Float64Engine, np.float64),
    (libindex.Float32Engine, np.float32),
], ids=lambda x: x[0].__name__)
def numeric_indexing_engine_type_and_dtype(request):
    return request.param
