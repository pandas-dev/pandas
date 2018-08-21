import pytest
import numpy as np

import pandas as pd
from pandas.core.sparse.api import SparseDtype


@pytest.mark.parametrize("dtype, fill_value", [
    ('int', 0),
    ('float', np.nan),
    ('bool', False),
    ('object', np.nan),
    ('datetime64[ns]', pd.NaT),
    ('timedelta64[ns]', pd.NaT),
])
def test_inferred_dtype(dtype, fill_value):
    sparse_dtype = SparseDtype(dtype)
    result = sparse_dtype.fill_value
    if pd.isna(fill_value):
        assert pd.isna(result) and type(result) == type(fill_value)
    else:
        assert result == fill_value


def test_from_sparse_dtype():
    dtype = SparseDtype('float', 0)
    result = SparseDtype(dtype)
    assert result.fill_value == 0


@pytest.mark.parametrize('dtype, fill_value', [
    ('int', None),
    ('float', None),
    ('bool', None),
    ('object', None),
    ('datetime64[ns]', None),
    ('timedelta64[ns]', None),
    ('int', np.nan),
    ('float', 0),
])
def test_equal(dtype, fill_value):
    a = SparseDtype(dtype, fill_value)
    b = SparseDtype(dtype, fill_value)
    assert a == b


@pytest.mark.parametrize('a, b', [
    (SparseDtype('float64'), SparseDtype('float32')),
    (SparseDtype('float64'), SparseDtype('float64', 0)),
    (SparseDtype('float64'), SparseDtype('datetime64[ns]', np.nan)),
    (SparseDtype(int, pd.NaT), SparseDtype(float, pd.NaT)),
    (SparseDtype('float64'), np.dtype('float64')),
])
def test_not_equal(a, b):
    assert a != b
