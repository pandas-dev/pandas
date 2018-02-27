import pytest

import numpy as np
from pandas import DataFrame, MultiIndex
from pandas.tests.io.pytables.common import (_check_roundtrip,
                                             ensure_clean_store)

import pandas.util.testing as tm
import pandas.util._test_decorators as td


def test_store_hierarchical():
    index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                               ['one', 'two', 'three']],
                       labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                               [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                       names=['foo', 'bar'])
    frame = DataFrame(np.random.randn(10, 3), index=index,
                      columns=['A', 'B', 'C'])

    _check_roundtrip(frame, tm.assert_frame_equal)
    _check_roundtrip(frame.T, tm.assert_frame_equal)
    _check_roundtrip(frame['A'], tm.assert_series_equal)

    # check that the names are stored
    with ensure_clean_store() as store:
        store['frame'] = frame
        recons = store['frame']
        tm.assert_frame_equal(recons, frame)


def test_store_series_name():
    df = tm.makeDataFrame()
    series = df['A']

    with ensure_clean_store() as store:
        store['series'] = series
        recons = store['series']
        tm.assert_series_equal(recons, series)


@pytest.mark.parametrize("compression", [
    False, pytest.param(True, marks=td.skip_if_windows_python_3)
])
def test_store_mixed(compression):
    def _make_one():
        df = tm.makeDataFrame()
        df['obj1'] = 'foo'
        df['obj2'] = 'bar'
        df['bool1'] = df['A'] > 0
        df['bool2'] = df['B'] > 0
        df['int1'] = 1
        df['int2'] = 2
        return df._consolidate()

    df1 = _make_one()
    df2 = _make_one()

    _check_roundtrip(df1, tm.assert_frame_equal)
    _check_roundtrip(df2, tm.assert_frame_equal)

    with ensure_clean_store() as store:
        store['obj'] = df1
        tm.assert_frame_equal(store['obj'], df1)
        store['obj'] = df2
        tm.assert_frame_equal(store['obj'], df2)

    # check that can store Series of all of these types
    _check_roundtrip(df1['obj1'], tm.assert_series_equal,
                     compression=compression)
    _check_roundtrip(df1['bool1'], tm.assert_series_equal,
                     compression=compression)
    _check_roundtrip(df1['int1'], tm.assert_series_equal,
                     compression=compression)
