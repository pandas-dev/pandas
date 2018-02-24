import pytest
from warnings import catch_warnings
import datetime

import numpy as np
import pandas as pd
from pandas import (Series, DataFrame, Timestamp, compat, date_range)
from .common import (ensure_clean_store, _maybe_remove, ensure_clean_path,
                     _check_roundtrip)

import pandas.util.testing as tm
from pandas.compat import (text_type, lrange, u)
from pandas.io.pytables import read_hdf


def test_index_types():
    with catch_warnings(record=True):
        values = np.random.randn(2)

        func = lambda l, r: tm.assert_series_equal(l, r,
                                                   check_dtype=True,
                                                   check_index_type=True,
                                                   check_series_type=True)

    with catch_warnings(record=True):
        ser = Series(values, [0, 'y'])
        _check_roundtrip(ser, func)

    with catch_warnings(record=True):
        ser = Series(values, [datetime.datetime.today(), 0])
        _check_roundtrip(ser, func)

    with catch_warnings(record=True):
        ser = Series(values, ['y', 0])
        _check_roundtrip(ser, func)

    with catch_warnings(record=True):
        ser = Series(values, [datetime.date.today(), 'a'])
        _check_roundtrip(ser, func)

    with catch_warnings(record=True):

        ser = Series(values, [0, 'y'])
        _check_roundtrip(ser, func)

        ser = Series(values, [datetime.datetime.today(), 0])
        _check_roundtrip(ser, func)

        ser = Series(values, ['y', 0])
        _check_roundtrip(ser, func)

        ser = Series(values, [datetime.date.today(), 'a'])
        _check_roundtrip(ser, func)

        ser = Series(values, [1.23, 'b'])
        _check_roundtrip(ser, func)

        ser = Series(values, [1, 1.53])
        _check_roundtrip(ser, func)

        ser = Series(values, [1, 5])
        _check_roundtrip(ser, func)

        ser = Series(values, [datetime.datetime(
            2012, 1, 1), datetime.datetime(2012, 1, 2)])
        _check_roundtrip(ser, func)


def test_store_index_types():
    # GH5386
    # test storing various index types
    with ensure_clean_store() as store:
        def check(format, index):
            df = DataFrame(np.random.randn(10, 2), columns=list('AB'))
            df.index = index(len(df))

            _maybe_remove(store, 'df')
            store.put('df', df, format=format)
            tm.assert_frame_equal(df, store['df'])

        for index in [tm.makeFloatIndex, tm.makeStringIndex,
                      tm.makeIntIndex, tm.makeDateIndex]:

            check('table', index)
            check('fixed', index)

        # period index currently broken for table
        # seee GH7796 FIXME
        check('fixed', tm.makePeriodIndex)
        # check('table',tm.makePeriodIndex)

        # unicode
        index = tm.makeUnicodeIndex
        if compat.PY3:
            check('table', index)
            check('fixed', index)
        else:

            # only support for fixed types (and they have a perf warning)
            pytest.raises(TypeError, check, 'table', index)

            # PerformanceWarning
            with catch_warnings(record=True):
                check('fixed', index)


def test_table_index_incompatible_dtypes():
    df1 = DataFrame({'a': [1, 2, 3]})
    df2 = DataFrame({'a': [4, 5, 6]},
                    index=date_range('1/1/2000', periods=3))

    with ensure_clean_store() as store:
        store.put('frame', df1, format='table')
        pytest.raises(TypeError, store.put, 'frame', df2,
                      format='table', append=True)


def test_store_index_name():
    df = tm.makeDataFrame()
    df.index.name = 'foo'

    with ensure_clean_store() as store:
        store['frame'] = df
        recons = store['frame']
        tm.assert_frame_equal(recons, df)


def test_store_index_name_with_tz():
    # GH 13884
    df = pd.DataFrame({'A': [1, 2]})
    df.index = pd.DatetimeIndex([1234567890123456787, 1234567890123456788])
    df.index = df.index.tz_localize('UTC')
    df.index.name = 'foo'

    with ensure_clean_store() as store:
        store.put('frame', df, format='table')
        recons = store['frame']
        tm.assert_frame_equal(recons, df)


@pytest.mark.parametrize('table_format', ['table', 'fixed'])
def test_store_index_name_numpy_str(table_format):
    # GH #13492
    idx = pd.Index(pd.to_datetime([datetime.date(2000, 1, 1),
                                   datetime.date(2000, 1, 2)]),
                   name=u('cols\u05d2'))
    idx1 = pd.Index(pd.to_datetime([datetime.date(2010, 1, 1),
                                    datetime.date(2010, 1, 2)]),
                    name=u('rows\u05d0'))
    df = pd.DataFrame(np.arange(4).reshape(2, 2), columns=idx, index=idx1)

    # This used to fail, returning numpy strings instead of python strings.
    with ensure_clean_path() as path:
        df.to_hdf(path, 'df', format=table_format)
        df2 = read_hdf(path, 'df')

        tm.assert_frame_equal(df, df2, check_names=True)

        assert type(df2.index.name) == text_type
        assert type(df2.columns.name) == text_type


def test_retain_index_attributes():
    # GH 3499, losing frequency info on index recreation
    df = DataFrame(dict(
        A=Series(lrange(3),
                 index=date_range('2000-1-1', periods=3, freq='H'))))

    with ensure_clean_store() as store:
        _maybe_remove(store, 'data')
        store.put('data', df, format='table')

        result = store.get('data')
        tm.assert_frame_equal(df, result)

        for attr in ['freq', 'tz', 'name']:
            for idx in ['index', 'columns']:
                assert (getattr(getattr(df, idx), attr, None) ==
                        getattr(getattr(result, idx), attr, None))

        # try to append a table with a different frequency
        with catch_warnings(record=True):
            df2 = DataFrame(dict(
                A=Series(lrange(3),
                         index=date_range('2002-1-1',
                                          periods=3, freq='D'))))
            store.append('data', df2)

        assert store.get_storer('data').info['index']['freq'] is None

        # this is ok
        _maybe_remove(store, 'df2')
        df2 = DataFrame(dict(
            A=Series(lrange(3),
                     index=[Timestamp('20010101'), Timestamp('20010102'),
                            Timestamp('20020101')])))
        store.append('df2', df2)
        df3 = DataFrame(dict(
            A=Series(lrange(3),
                     index=date_range('2002-1-1', periods=3,
                                      freq='D'))))
        store.append('df2', df3)


def test_retain_index_attributes2():
    with ensure_clean_path() as path:
        with catch_warnings(record=True):
            df = DataFrame(dict(
                A=Series(lrange(3),
                         index=date_range('2000-1-1',
                                          periods=3, freq='H'))))
            df.to_hdf(path, 'data', mode='w', append=True)
            df2 = DataFrame(dict(
                A=Series(lrange(3),
                         index=date_range('2002-1-1', periods=3,
                                          freq='D'))))
            df2.to_hdf(path, 'data', append=True)

            idx = date_range('2000-1-1', periods=3, freq='H')
            idx.name = 'foo'
            df = DataFrame(dict(A=Series(lrange(3), index=idx)))
            df.to_hdf(path, 'data', mode='w', append=True)

        assert read_hdf(path, 'data').index.name == 'foo'

        with catch_warnings(record=True):
            idx2 = date_range('2001-1-1', periods=3, freq='H')
            idx2.name = 'bar'
            df2 = DataFrame(dict(A=Series(lrange(3), index=idx2)))
            df2.to_hdf(path, 'data', append=True)

        assert read_hdf(path, 'data').index.name is None


def test_float_index():
    # GH #454
    index = np.random.randn(10)
    s = Series(np.random.randn(10), index=index)
    _check_roundtrip(s, tm.assert_series_equal)


def test_tuple_index():
    # GH #492
    col = np.arange(10)
    idx = [(0., 1.), (2., 3.), (4., 5.)]
    data = np.random.randn(30).reshape((3, 10))
    DF = DataFrame(data, index=idx, columns=col)

    with catch_warnings(record=True):
        _check_roundtrip(DF, tm.assert_frame_equal)


def test_create_table_index():
    with ensure_clean_store() as store:
        with catch_warnings(record=True):
            def col(t, column):
                return getattr(store.get_storer(t).table.cols, column)

            # index=False
            wp = tm.makePanel()
            store.append('p5', wp, index=False)
            store.create_table_index('p5', columns=['major_axis'])
            assert(col('p5', 'major_axis').is_indexed is True)
            assert(col('p5', 'minor_axis').is_indexed is False)

            # index=True
            store.append('p5i', wp, index=True)
            assert(col('p5i', 'major_axis').is_indexed is True)
            assert(col('p5i', 'minor_axis').is_indexed is True)

            # default optlevels
            store.get_storer('p5').create_index()
            assert(col('p5', 'major_axis').index.optlevel == 6)
            assert(col('p5', 'minor_axis').index.kind == 'medium')

            # let's change the indexing scheme
            store.create_table_index('p5')
            assert(col('p5', 'major_axis').index.optlevel == 6)
            assert(col('p5', 'minor_axis').index.kind == 'medium')
            store.create_table_index('p5', optlevel=9)
            assert(col('p5', 'major_axis').index.optlevel == 9)
            assert(col('p5', 'minor_axis').index.kind == 'medium')
            store.create_table_index('p5', kind='full')
            assert(col('p5', 'major_axis').index.optlevel == 9)
            assert(col('p5', 'minor_axis').index.kind == 'full')
            store.create_table_index('p5', optlevel=1, kind='light')
            assert(col('p5', 'major_axis').index.optlevel == 1)
            assert(col('p5', 'minor_axis').index.kind == 'light')

            # data columns
            df = tm.makeTimeDataFrame()
            df['string'] = 'foo'
            df['string2'] = 'bar'
            store.append('f', df, data_columns=['string', 'string2'])
            assert(col('f', 'index').is_indexed is True)
            assert(col('f', 'string').is_indexed is True)
            assert(col('f', 'string2').is_indexed is True)

            # specify index=columns
            store.append(
                'f2', df, index=['string'],
                data_columns=['string', 'string2'])
            assert(col('f2', 'index').is_indexed is False)
            assert(col('f2', 'string').is_indexed is True)
            assert(col('f2', 'string2').is_indexed is False)

            # try to index a non-table
            _maybe_remove(store, 'f2')
            store.put('f2', df)
            pytest.raises(TypeError, store.create_table_index, 'f2')
