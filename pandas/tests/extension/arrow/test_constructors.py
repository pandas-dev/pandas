from datetime import (
    date,
    datetime,
    time,
    timedelta,
)

import pytest

import pandas as pd
import pandas._testing as tm
from pandas.tests.extension import base

pa = pytest.importorskip("pyarrow", minversion="1.0.1")

from pandas.core.arrays.arrow.dtype import ArrowDtype  # isort:skip


@pytest.fixture(params=tm.ALL_PYARROW_DTYPES)
def dtype(request):
    return ArrowDtype(pa_dtype=request.param)


@pytest.fixture
def data(dtype):
    pa_dtype = dtype.pa_dtype
    if pa.types.is_boolean(pa_dtype):
        data = [True, None, False, None]
    elif pa.types.is_floating(pa_dtype):
        data = [1.0, None, 0.0, None, -2.0, None, 0.5, None, 99.9, None]
    elif pa.types.is_signed_integer(pa_dtype):
        data = [1, None, 0, None, -2, None, 10]
    elif pa.types.is_unsigned_integer(pa_dtype):
        data = [1, None, 0, None, 2, None, 10]
    elif pa.types.is_date(pa_dtype):
        data = [date(2022, 1, 1), None, date(1999, 12, 31)]
    elif pa.types.is_timestamp(pa_dtype):
        data = [
            datetime(2020, 1, 1, 1, 1, 1, 1),
            None,
            datetime(1999, 1, 1, 1, 1, 1, 1),
        ]
    elif pa.types.is_duration(pa_dtype):
        data = [timedelta(1), None, timedelta(1, 1)]
    elif pa.types.is_time(pa_dtype):
        data = [time(12, 0), None, time(0, 12)]
    else:
        data = []
    return pd.array(data, dtype=dtype)


class TestConstructors(base.BaseConstructorsTests):
    pass
