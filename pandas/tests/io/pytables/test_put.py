import pytest
from warnings import catch_warnings
import datetime

import numpy as np
from pandas import (Series, DataFrame, Index, Timestamp)
from .common import ensure_clean_store, _maybe_remove, _check_roundtrip
import pandas.util.testing as tm
import pandas.util._test_decorators as td
from pandas.compat import range


def test_put():
    with ensure_clean_store() as store:
        ts = tm.makeTimeSeries()
        df = tm.makeTimeDataFrame()
        store['a'] = ts
        store['b'] = df[:10]
        store['foo/bar/bah'] = df[:10]
        store['foo'] = df[:10]
        store['/foo'] = df[:10]
        store.put('c', df[:10], format='table')

        # not OK, not a table
        pytest.raises(
            ValueError, store.put, 'b', df[10:], append=True)

        # can't put to a table (use append instead)
        pytest.raises(ValueError, store.put, 'c', df[10:], append=True)

        # overwrite table
        store.put('c', df[:10], format='table', append=False)
        tm.assert_frame_equal(df[:10], store['c'])


def test_put_string_index():
    with ensure_clean_store() as store:
        index = Index(
            ["I am a very long string index: %s" % i for i in range(20)])
        s = Series(np.arange(20), index=index)
        df = DataFrame({'A': s, 'B': s})

        store['a'] = s
        tm.assert_series_equal(store['a'], s)

        store['b'] = df
        tm.assert_frame_equal(store['b'], df)

        # mixed length
        index = Index(['abcdefghijklmnopqrstuvwxyz1234567890'] +
                      ["I am a very long string index: %s" % i
                       for i in range(20)])
        s = Series(np.arange(21), index=index)
        df = DataFrame({'A': s, 'B': s})
        store['a'] = s
        tm.assert_series_equal(store['a'], s)

        store['b'] = df
        tm.assert_frame_equal(store['b'], df)


def test_put_compression():
    with ensure_clean_store() as store:
        df = tm.makeTimeDataFrame()

        store.put('c', df, format='table', complib='zlib')
        tm.assert_frame_equal(store['c'], df)

        # can't compress if format='fixed'
        pytest.raises(ValueError, store.put, 'b', df,
                      format='fixed', complib='zlib')


def test_put_integer():
    # non-date, non-string index
    df = DataFrame(np.random.randn(50, 100))
    _check_roundtrip(df, tm.assert_frame_equal)


def test_put_mixed_type():
    df = tm.makeTimeDataFrame()
    df['obj1'] = 'foo'
    df['obj2'] = 'bar'
    df['bool1'] = df['A'] > 0
    df['bool2'] = df['B'] > 0
    df['bool3'] = True
    df['int1'] = 1
    df['int2'] = 2
    df['timestamp1'] = Timestamp('20010102')
    df['timestamp2'] = Timestamp('20010103')
    df['datetime1'] = datetime.datetime(2001, 1, 2, 0, 0)
    df['datetime2'] = datetime.datetime(2001, 1, 3, 0, 0)
    df.loc[3:6, ['obj1']] = np.nan
    df = df._consolidate()._convert(datetime=True)

    with ensure_clean_store() as store:
        _maybe_remove(store, 'df')

        # PerformanceWarning
        with catch_warnings(record=True):
            store.put('df', df)

        expected = store.get('df')
        tm.assert_frame_equal(expected, df)


@td.skip_if_windows_python_3
def test_put_compression_blosc():
    df = tm.makeTimeDataFrame()

    with ensure_clean_store() as store:
        # can't compress if format='fixed'
        pytest.raises(ValueError, store.put, 'b', df,
                      format='fixed', complib='blosc')

        store.put('c', df, format='table', complib='blosc')
        tm.assert_frame_equal(store['c'], df)
