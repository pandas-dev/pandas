import pytest
import pandas as pd


@pytest.mark.parametrize('test_dtype', [object, 'int64'])
def test_dtypes(test_dtype):
    df = pd.DataFrame({'A': pd.Series([1, 2, 3], dtype=test_dtype), 'B': [1, 2, 3]})
    expected = df.dtypes.values[0].type

    result = df.set_index('A').index.dtype.type
    assert result == expected


@pytest.fixture
def mixed_series():
    return pd.Series([1, 2, 3, 'apple', 'corn'], dtype=object)


@pytest.fixture
def int_series():
    return pd.Series([100, 200, 300, 400, 500])


def test_dtypes_between_queries(mixed_series, int_series)
    df = pd.DataFrame({'item': mixed_series, 'cost': int_series})

    orig_dtypes = df.dtypes
    item_dtype = orig_dtypes.get('item').type
    cost_dtype = orig_dtypes.get('cost').type
    expected = {'item': item_dtype, 'cost': cost_dtype}

    # after applying a query that would remove strings from the 'item' series with dtype: object,
    # that series should remain as dtype: object as it becomes an index, and again as it becomes 
    # a column again after calling reset_index()
    dtypes_transformed = df.query('cost < 400').set_index('item').reset_index().dtypes
    item_dtype_transformed = dtypes_transformed.get('item').type
    cost_dtype_transformed = dtypes_transformed.get('cost').type
    result = {'item': item_dtype_transformed, 'cost': cost_dtype_transformed}

    assert result == expected
