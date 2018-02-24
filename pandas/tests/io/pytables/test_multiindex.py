import pytest
import datetime

import numpy as np
import pandas as pd
from pandas import (Series, DataFrame, date_range, Index,
                    concat, MultiIndex, RangeIndex, Int64Index)
from .common import (ensure_clean_store, ensure_clean_path,
                     _maybe_remove, _check_roundtrip)
import pandas.util.testing as tm
from pandas.compat import range
from pandas.io.pytables import HDFStore, read_hdf


def test_column_multiindex():
    # GH 4710
    # recreate multi-indexes properly
    index = MultiIndex.from_tuples([('A', 'a'), ('A', 'b'),
                                    ('B', 'a'), ('B', 'b')],
                                   names=['first', 'second'])
    df = DataFrame(np.arange(12).reshape(3, 4), columns=index)
    expected = df.copy()
    if isinstance(expected.index, RangeIndex):
        expected.index = Int64Index(expected.index)

    with ensure_clean_store() as store:
        store.put('df', df)
        tm.assert_frame_equal(store['df'], expected,
                              check_index_type=True,
                              check_column_type=True)

        store.put('df1', df, format='table')
        tm.assert_frame_equal(store['df1'], expected,
                              check_index_type=True,
                              check_column_type=True)

        pytest.raises(ValueError, store.put, 'df2', df,
                      format='table', data_columns=['A'])
        pytest.raises(ValueError, store.put, 'df3', df,
                      format='table', data_columns=True)

    # appending multi-column on existing table (see GH 6167)
    with ensure_clean_store() as store:
        store.append('df2', df)
        store.append('df2', df)

        tm.assert_frame_equal(store['df2'], concat((df, df)))

    # non_index_axes name
    df = DataFrame(np.arange(12).reshape(3, 4),
                   columns=Index(list('ABCD'), name='foo'))
    expected = df.copy()
    if isinstance(expected.index, RangeIndex):
        expected.index = Int64Index(expected.index)

    with ensure_clean_store() as store:
        store.put('df1', df, format='table')
        tm.assert_frame_equal(store['df1'], expected,
                              check_index_type=True,
                              check_column_type=True)


def test_store_multiindex():
    # validate multi-index names
    # GH 5527
    with ensure_clean_store() as store:
        def make_index(names=None):
            return MultiIndex.from_tuples([(datetime.datetime(2013, 12, d),
                                            s, t)
                                           for d in range(1, 3)
                                           for s in range(2)
                                           for t in range(3)],
                                          names=names)

        # no names
        _maybe_remove(store, 'df')
        df = DataFrame(np.zeros((12, 2)), columns=[
                       'a', 'b'], index=make_index())
        store.append('df', df)
        tm.assert_frame_equal(store.select('df'), df)

        # partial names
        _maybe_remove(store, 'df')
        df = DataFrame(np.zeros((12, 2)), columns=[
                       'a', 'b'], index=make_index(['date', None, None]))
        store.append('df', df)
        tm.assert_frame_equal(store.select('df'), df)

        # series
        _maybe_remove(store, 's')
        s = Series(np.zeros(12), index=make_index(['date', None, None]))
        store.append('s', s)
        xp = Series(np.zeros(12), index=make_index(
            ['date', 'level_1', 'level_2']))
        tm.assert_series_equal(store.select('s'), xp)

        # dup with column
        _maybe_remove(store, 'df')
        df = DataFrame(np.zeros((12, 2)), columns=[
                       'a', 'b'], index=make_index(['date', 'a', 't']))
        pytest.raises(ValueError, store.append, 'df', df)

        # fully names
        _maybe_remove(store, 'df')
        df = DataFrame(np.zeros((12, 2)), columns=[
                       'a', 'b'], index=make_index(['date', 's', 't']))
        store.append('df', df)
        tm.assert_frame_equal(store.select('df'), df)


def test_mi_data_columns():
    # GH 14435
    idx = pd.MultiIndex.from_arrays([date_range('2000-01-01', periods=5),
                                     range(5)], names=['date', 'id'])
    df = pd.DataFrame({'a': [1.1, 1.2, 1.3, 1.4, 1.5]}, index=idx)

    with ensure_clean_store() as store:
        store.append('df', df, data_columns=True)

        actual = store.select('df', where='id == 1')
        expected = df.iloc[[1], :]
        tm.assert_frame_equal(actual, expected)


def test_columns_multiindex_modified():
    # BUG: 7212
    # read_hdf store.select modified the passed columns parameters
    # when multi-indexed.
    df = DataFrame(np.random.rand(4, 5),
                   index=list('abcd'),
                   columns=list('ABCDE'))
    df.index.name = 'letters'
    df = df.set_index(keys='E', append=True)

    data_columns = df.index.names + df.columns.tolist()
    with ensure_clean_path() as path:
        df.to_hdf(path, 'df',
                  mode='a',
                  append=True,
                  data_columns=data_columns,
                  index=False)
        cols2load = list('BCD')
        cols2load_original = list(cols2load)
        df_loaded = read_hdf(path, 'df', columns=cols2load)  # noqa
        assert cols2load_original == cols2load


def test_frame_select_complex2():
    with ensure_clean_path(['parms.hdf', 'hist.hdf']) as paths:
        pp, hh = paths

        # use non-trivial selection criteria
        parms = DataFrame({'A': [1, 1, 2, 2, 3]})
        parms.to_hdf(pp, 'df', mode='w',
                     format='table', data_columns=['A'])

        selection = read_hdf(pp, 'df', where='A=[2,3]')
        hist = DataFrame(np.random.randn(25, 1),
                         columns=['data'],
                         index=MultiIndex.from_tuples(
                             [(i, j) for i in range(5)
                              for j in range(5)],
                             names=['l1', 'l2']))

        hist.to_hdf(hh, 'df', mode='w', format='table')

        expected = read_hdf(hh, 'df', where='l1=[2, 3, 4]')

        # sccope with list like
        l = selection.index.tolist()  # noqa
        store = HDFStore(hh)
        result = store.select('df', where='l1=l')
        tm.assert_frame_equal(result, expected)
        store.close()

        result = read_hdf(hh, 'df', where='l1=l')
        tm.assert_frame_equal(result, expected)

        # index
        index = selection.index  # noqa
        result = read_hdf(hh, 'df', where='l1=index')
        tm.assert_frame_equal(result, expected)

        result = read_hdf(hh, 'df', where='l1=selection.index')
        tm.assert_frame_equal(result, expected)

        result = read_hdf(hh, 'df', where='l1=selection.index.tolist()')
        tm.assert_frame_equal(result, expected)

        result = read_hdf(hh, 'df', where='l1=list(selection.index)')
        tm.assert_frame_equal(result, expected)

        # sccope with index
        store = HDFStore(hh)

        result = store.select('df', where='l1=index')
        tm.assert_frame_equal(result, expected)

        result = store.select('df', where='l1=selection.index')
        tm.assert_frame_equal(result, expected)

        result = store.select('df', where='l1=selection.index.tolist()')
        tm.assert_frame_equal(result, expected)

        result = store.select('df', where='l1=list(selection.index)')
        tm.assert_frame_equal(result, expected)

        store.close()


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
