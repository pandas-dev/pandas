import pytest

import numpy as np
import pandas as pd
from pandas import (Series, DataFrame, Categorical, concat)
from .common import ensure_clean_store, ensure_clean_path, _maybe_remove

import pandas.util.testing as tm
from pandas.io.pytables import read_hdf


def test_categorical():
    with ensure_clean_store() as store:
        # Basic
        _maybe_remove(store, 's')
        s = Series(Categorical(['a', 'b', 'b', 'a', 'a', 'c'], categories=[
                   'a', 'b', 'c', 'd'], ordered=False))
        store.append('s', s, format='table')
        result = store.select('s')
        tm.assert_series_equal(s, result)

        _maybe_remove(store, 's_ordered')
        s = Series(Categorical(['a', 'b', 'b', 'a', 'a', 'c'], categories=[
                   'a', 'b', 'c', 'd'], ordered=True))
        store.append('s_ordered', s, format='table')
        result = store.select('s_ordered')
        tm.assert_series_equal(s, result)

        _maybe_remove(store, 'df')

        df = DataFrame({"s": s, "vals": [1, 2, 3, 4, 5, 6]})
        store.append('df', df, format='table')
        result = store.select('df')
        tm.assert_frame_equal(result, df)

        # Dtypes
        s = Series([1, 1, 2, 2, 3, 4, 5]).astype('category')
        store.append('si', s)
        result = store.select('si')
        tm.assert_series_equal(result, s)

        s = Series([1, 1, np.nan, 2, 3, 4, 5]).astype('category')
        store.append('si2', s)
        result = store.select('si2')
        tm.assert_series_equal(result, s)

        # Multiple
        df2 = df.copy()
        df2['s2'] = Series(list('abcdefg')).astype('category')
        store.append('df2', df2)
        result = store.select('df2')
        tm.assert_frame_equal(result, df2)

        # Make sure the metadata is OK
        info = store.info()
        assert '/df2   ' in info
        # assert '/df2/meta/values_block_0/meta' in info
        assert '/df2/meta/values_block_1/meta' in info

        # unordered
        s = Series(Categorical(['a', 'b', 'b', 'a', 'a', 'c'], categories=[
                   'a', 'b', 'c', 'd'], ordered=False))
        store.append('s2', s, format='table')
        result = store.select('s2')
        tm.assert_series_equal(result, s)

        # Query
        store.append('df3', df, data_columns=['s'])
        expected = df[df.s.isin(['b', 'c'])]
        result = store.select('df3', where=['s in ["b","c"]'])
        tm.assert_frame_equal(result, expected)

        expected = df[df.s.isin(['b', 'c'])]
        result = store.select('df3', where=['s = ["b","c"]'])
        tm.assert_frame_equal(result, expected)

        expected = df[df.s.isin(['d'])]
        result = store.select('df3', where=['s in ["d"]'])
        tm.assert_frame_equal(result, expected)

        expected = df[df.s.isin(['f'])]
        result = store.select('df3', where=['s in ["f"]'])
        tm.assert_frame_equal(result, expected)

        # Appending with same categories is ok
        store.append('df3', df)

        df = concat([df, df])
        expected = df[df.s.isin(['b', 'c'])]
        result = store.select('df3', where=['s in ["b","c"]'])
        tm.assert_frame_equal(result, expected)

        # Appending must have the same categories
        df3 = df.copy()
        df3['s'].cat.remove_unused_categories(inplace=True)

        with pytest.raises(ValueError):
            store.append('df3', df3)

        # Remove, and make sure meta data is removed (its a recursive
        # removal so should be).
        result = store.select('df3/meta/s/meta')
        assert result is not None
        store.remove('df3')

        with pytest.raises(KeyError):
            store.select('df3/meta/s/meta')


def test_categorical_conversion():
    # GH13322
    # Check that read_hdf with categorical columns doesn't return rows if
    # where criteria isn't met.
    obsids = ['ESP_012345_6789', 'ESP_987654_3210']
    imgids = ['APF00006np', 'APF0001imm']
    data = [4.3, 9.8]

    # Test without categories
    df = DataFrame(dict(obsids=obsids, imgids=imgids, data=data))

    # We are expecting an empty DataFrame matching types of df
    expected = df.iloc[[], :]
    with ensure_clean_path() as path:
        df.to_hdf(path, 'df', format='table', data_columns=True)
        result = read_hdf(path, 'df', where='obsids=B')
        tm.assert_frame_equal(result, expected)

    # Test with categories
    df.obsids = df.obsids.astype('category')
    df.imgids = df.imgids.astype('category')

    # We are expecting an empty DataFrame matching types of df
    expected = df.iloc[[], :]
    with ensure_clean_path() as path:
        df.to_hdf(path, 'df', format='table', data_columns=True)
        result = read_hdf(path, 'df', where='obsids=B')
        tm.assert_frame_equal(result, expected)


def test_categorical_nan_only_columns():
    # GH18413
    # Check that read_hdf with categorical columns with NaN-only values can
    # be read back.
    df = pd.DataFrame({
        'a': ['a', 'b', 'c', np.nan],
        'b': [np.nan, np.nan, np.nan, np.nan],
        'c': [1, 2, 3, 4],
        'd': pd.Series([None] * 4, dtype=object)
    })
    df['a'] = df.a.astype('category')
    df['b'] = df.b.astype('category')
    df['d'] = df.b.astype('category')
    expected = df
    with ensure_clean_path() as path:
        df.to_hdf(path, 'df', format='table', data_columns=True)
        result = read_hdf(path, 'df')
        tm.assert_frame_equal(result, expected)
