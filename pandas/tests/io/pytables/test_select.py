import pytest
from warnings import catch_warnings

import numpy as np
from pandas import (Series, DataFrame, date_range, Index, Timestamp,
                    concat, MultiIndex, bdate_range, Panel)
from pandas.tests.io.pytables.common import (ensure_clean_store,
                                             ensure_clean_path, _maybe_remove)

import pandas.util.testing as tm
from pandas.util.testing import assert_panel_equal, assert_frame_equal
from pandas.compat import range
from pandas.io.pytables import read_hdf, Term


def test_select_columns_in_where():
    # GH 6169
    # recreate multi-indexes when columns is passed
    # in the `where` argument
    index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                               ['one', 'two', 'three']],
                       labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                               [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                       names=['foo_name', 'bar_name'])

    # With a DataFrame
    df = DataFrame(np.random.randn(10, 3), index=index,
                   columns=['A', 'B', 'C'])

    with ensure_clean_store() as store:
        store.put('df', df, format='table')
        expected = df[['A']]

        tm.assert_frame_equal(store.select('df', columns=['A']), expected)

        tm.assert_frame_equal(store.select(
            'df', where="columns=['A']"), expected)

    # With a Series
    s = Series(np.random.randn(10), index=index,
               name='A')
    with ensure_clean_store() as store:
        store.put('s', s, format='table')
        tm.assert_series_equal(store.select('s', where="columns=['A']"), s)


def test_select_with_dups():
    # single dtypes
    df = DataFrame(np.random.randn(10, 4), columns=['A', 'A', 'B', 'B'])
    df.index = date_range('20130101 9:30', periods=10, freq='T')

    with ensure_clean_store() as store:
        store.append('df', df)

        result = store.select('df')
        expected = df
        assert_frame_equal(result, expected, by_blocks=True)

        result = store.select('df', columns=df.columns)
        expected = df
        assert_frame_equal(result, expected, by_blocks=True)

        result = store.select('df', columns=['A'])
        expected = df.loc[:, ['A']]
        assert_frame_equal(result, expected)

    # dups across dtypes
    df = concat([DataFrame(np.random.randn(10, 4),
                           columns=['A', 'A', 'B', 'B']),
                 DataFrame(np.random.randint(0, 10, size=20)
                           .reshape(10, 2),
                           columns=['A', 'C'])],
                axis=1)
    df.index = date_range('20130101 9:30', periods=10, freq='T')

    with ensure_clean_store() as store:
        store.append('df', df)

        result = store.select('df')
        expected = df
        assert_frame_equal(result, expected, by_blocks=True)

        result = store.select('df', columns=df.columns)
        expected = df
        assert_frame_equal(result, expected, by_blocks=True)

        expected = df.loc[:, ['A']]
        result = store.select('df', columns=['A'])
        assert_frame_equal(result, expected, by_blocks=True)

        expected = df.loc[:, ['B', 'A']]
        result = store.select('df', columns=['B', 'A'])
        assert_frame_equal(result, expected, by_blocks=True)

    # duplicates on both index and columns
    with ensure_clean_store() as store:
        store.append('df', df)
        store.append('df', df)

        expected = df.loc[:, ['B', 'A']]
        expected = concat([expected, expected])
        result = store.select('df', columns=['B', 'A'])
        assert_frame_equal(result, expected, by_blocks=True)


def test_select():
    with ensure_clean_store() as store:
        with catch_warnings(record=True):
            wp = tm.makePanel()

            # put/select ok
            _maybe_remove(store, 'wp')
            store.put('wp', wp, format='table')
            store.select('wp')

            # non-table ok (where = None)
            _maybe_remove(store, 'wp')
            store.put('wp2', wp)
            store.select('wp2')

            # selection on the non-indexable with a large number of columns
            wp = Panel(np.random.randn(100, 100, 100),
                       items=['Item%03d' % i for i in range(100)],
                       major_axis=date_range('1/1/2000', periods=100),
                       minor_axis=['E%03d' % i for i in range(100)])

            _maybe_remove(store, 'wp')
            store.append('wp', wp)
            items = ['Item%03d' % i for i in range(80)]
            result = store.select('wp', 'items=items')
            expected = wp.reindex(items=items)
            assert_panel_equal(expected, result)

            # selectin non-table with a where
            # pytest.raises(ValueError, store.select,
            #                  'wp2', ('column', ['A', 'D']))

            # select with columns=
            df = tm.makeTimeDataFrame()
            _maybe_remove(store, 'df')
            store.append('df', df)
            result = store.select('df', columns=['A', 'B'])
            expected = df.reindex(columns=['A', 'B'])
            tm.assert_frame_equal(expected, result)

            # equivalentsly
            result = store.select('df', [("columns=['A', 'B']")])
            expected = df.reindex(columns=['A', 'B'])
            tm.assert_frame_equal(expected, result)

            # with a data column
            _maybe_remove(store, 'df')
            store.append('df', df, data_columns=['A'])
            result = store.select('df', ['A > 0'], columns=['A', 'B'])
            expected = df[df.A > 0].reindex(columns=['A', 'B'])
            tm.assert_frame_equal(expected, result)

            # all a data columns
            _maybe_remove(store, 'df')
            store.append('df', df, data_columns=True)
            result = store.select('df', ['A > 0'], columns=['A', 'B'])
            expected = df[df.A > 0].reindex(columns=['A', 'B'])
            tm.assert_frame_equal(expected, result)

            # with a data column, but different columns
            _maybe_remove(store, 'df')
            store.append('df', df, data_columns=['A'])
            result = store.select('df', ['A > 0'], columns=['C', 'D'])
            expected = df[df.A > 0].reindex(columns=['C', 'D'])
            tm.assert_frame_equal(expected, result)


def test_select_dtypes():
    with ensure_clean_store() as store:
        # with a Timestamp data column (GH #2637)
        df = DataFrame(dict(
            ts=bdate_range('2012-01-01', periods=300),
            A=np.random.randn(300)))
        _maybe_remove(store, 'df')
        store.append('df', df, data_columns=['ts', 'A'])

        result = store.select('df', "ts>=Timestamp('2012-02-01')")
        expected = df[df.ts >= Timestamp('2012-02-01')]
        tm.assert_frame_equal(expected, result)

        # bool columns (GH #2849)
        df = DataFrame(np.random.randn(5, 2), columns=['A', 'B'])
        df['object'] = 'foo'
        df.loc[4:5, 'object'] = 'bar'
        df['boolv'] = df['A'] > 0
        _maybe_remove(store, 'df')
        store.append('df', df, data_columns=True)

        expected = (df[df.boolv == True]  # noqa
                    .reindex(columns=['A', 'boolv']))
        for v in [True, 'true', 1]:
            result = store.select('df', 'boolv == %s' % str(v),
                                  columns=['A', 'boolv'])
            tm.assert_frame_equal(expected, result)

        expected = (df[df.boolv == False]  # noqa
                    .reindex(columns=['A', 'boolv']))
        for v in [False, 'false', 0]:
            result = store.select(
                'df', 'boolv == %s' % str(v), columns=['A', 'boolv'])
            tm.assert_frame_equal(expected, result)

        # integer index
        df = DataFrame(dict(A=np.random.rand(20), B=np.random.rand(20)))
        _maybe_remove(store, 'df_int')
        store.append('df_int', df)
        result = store.select(
            'df_int', "index<10 and columns=['A']")
        expected = df.reindex(index=list(df.index)[0:10], columns=['A'])
        tm.assert_frame_equal(expected, result)

        # float index
        df = DataFrame(dict(A=np.random.rand(
            20), B=np.random.rand(20), index=np.arange(20, dtype='f8')))
        _maybe_remove(store, 'df_float')
        store.append('df_float', df)
        result = store.select(
            'df_float', "index<10.0 and columns=['A']")
        expected = df.reindex(index=list(df.index)[0:10], columns=['A'])
        tm.assert_frame_equal(expected, result)

    with ensure_clean_store() as store:
        # floats w/o NaN
        df = DataFrame(
            dict(cols=range(11), values=range(11)), dtype='float64')
        df['cols'] = (df['cols'] + 10).apply(str)

        store.append('df1', df, data_columns=True)
        result = store.select(
            'df1', where='values>2.0')
        expected = df[df['values'] > 2.0]
        tm.assert_frame_equal(expected, result)

        # floats with NaN
        df.iloc[0] = np.nan
        expected = df[df['values'] > 2.0]

        store.append('df2', df, data_columns=True, index=False)
        result = store.select(
            'df2', where='values>2.0')
        tm.assert_frame_equal(expected, result)

        # https://github.com/PyTables/PyTables/issues/282
        # bug in selection when 0th row has a np.nan and an index
        # store.append('df3',df,data_columns=True)
        # result = store.select(
        #    'df3', where='values>2.0')
        # tm.assert_frame_equal(expected, result)

        # not in first position float with NaN ok too
        df = DataFrame(
            dict(cols=range(11), values=range(11)), dtype='float64')
        df['cols'] = (df['cols'] + 10).apply(str)

        df.iloc[1] = np.nan
        expected = df[df['values'] > 2.0]

        store.append('df4', df, data_columns=True)
        result = store.select(
            'df4', where='values>2.0')
        tm.assert_frame_equal(expected, result)

    # test selection with comparison against numpy scalar
    # GH 11283
    with ensure_clean_store() as store:
        df = tm.makeDataFrame()

        expected = df[df['A'] > 0]

        store.append('df', df, data_columns=True)
        np_zero = np.float64(0)  # noqa
        result = store.select('df', where=["A>np_zero"])
        tm.assert_frame_equal(expected, result)


def test_select_with_many_inputs():
    with ensure_clean_store() as store:
        df = DataFrame(dict(ts=bdate_range('2012-01-01', periods=300),
                            A=np.random.randn(300),
                            B=range(300),
                            users=['a'] * 50 + ['b'] * 50 + ['c'] * 100 +
                            ['a%03d' % i for i in range(100)]))
        _maybe_remove(store, 'df')
        store.append('df', df, data_columns=['ts', 'A', 'B', 'users'])

        # regular select
        result = store.select('df', "ts>=Timestamp('2012-02-01')")
        expected = df[df.ts >= Timestamp('2012-02-01')]
        tm.assert_frame_equal(expected, result)

        # small selector
        result = store.select(
            'df',
            "ts>=Timestamp('2012-02-01') & users=['a','b','c']")
        expected = df[(df.ts >= Timestamp('2012-02-01')) &
                      df.users.isin(['a', 'b', 'c'])]
        tm.assert_frame_equal(expected, result)

        # big selector along the columns
        selector = ['a', 'b', 'c'] + ['a%03d' % i for i in range(60)]
        result = store.select(
            'df',
            "ts>=Timestamp('2012-02-01') and users=selector")
        expected = df[(df.ts >= Timestamp('2012-02-01')) &
                      df.users.isin(selector)]
        tm.assert_frame_equal(expected, result)

        selector = range(100, 200)
        result = store.select('df', 'B=selector')
        expected = df[df.B.isin(selector)]
        tm.assert_frame_equal(expected, result)
        assert len(result) == 100

        # big selector along the index
        selector = Index(df.ts[0:100].values)
        result = store.select('df', 'ts=selector')
        expected = df[df.ts.isin(selector.values)]
        tm.assert_frame_equal(expected, result)
        assert len(result) == 100


def test_select_iterator():
    # single table
    with ensure_clean_store() as store:
        df = tm.makeTimeDataFrame(500)
        _maybe_remove(store, 'df')
        store.append('df', df)

        expected = store.select('df')

        results = [s for s in store.select('df', iterator=True)]
        result = concat(results)
        tm.assert_frame_equal(expected, result)

        results = [s for s in store.select('df', chunksize=100)]
        assert len(results) == 5
        result = concat(results)
        tm.assert_frame_equal(expected, result)

        results = [s for s in store.select('df', chunksize=150)]
        result = concat(results)
        tm.assert_frame_equal(result, expected)

    with ensure_clean_path() as path:
        df = tm.makeTimeDataFrame(500)
        df.to_hdf(path, 'df_non_table')
        pytest.raises(TypeError, read_hdf, path,
                      'df_non_table', chunksize=100)
        pytest.raises(TypeError, read_hdf, path,
                      'df_non_table', iterator=True)

    with ensure_clean_path() as path:
        df = tm.makeTimeDataFrame(500)
        df.to_hdf(path, 'df', format='table')

        results = [s for s in read_hdf(path, 'df', chunksize=100)]
        result = concat(results)

        assert len(results) == 5
        tm.assert_frame_equal(result, df)
        tm.assert_frame_equal(result, read_hdf(path, 'df'))

    # multiple
    with ensure_clean_store() as store:
        df1 = tm.makeTimeDataFrame(500)
        store.append('df1', df1, data_columns=True)
        df2 = tm.makeTimeDataFrame(500).rename(
            columns=lambda x: "%s_2" % x)
        df2['foo'] = 'bar'
        store.append('df2', df2)

        df = concat([df1, df2], axis=1)

        # full selection
        expected = store.select_as_multiple(
            ['df1', 'df2'], selector='df1')
        results = [s for s in store.select_as_multiple(
            ['df1', 'df2'], selector='df1', chunksize=150)]
        result = concat(results)
        tm.assert_frame_equal(expected, result)


def test_select_iterator_complete_8014():
    # GH 8014
    # using iterator and where clause
    chunksize = 1e4

    # no iterator
    with ensure_clean_store() as store:
        expected = tm.makeTimeDataFrame(100064, 'S')
        _maybe_remove(store, 'df')
        store.append('df', expected)

        beg_dt = expected.index[0]
        end_dt = expected.index[-1]

        # select w/o iteration and no where clause works
        result = store.select('df')
        tm.assert_frame_equal(expected, result)

        # select w/o iterator and where clause, single term, begin
        # of range, works
        where = "index >= '%s'" % beg_dt
        result = store.select('df', where=where)
        tm.assert_frame_equal(expected, result)

        # select w/o iterator and where clause, single term, end
        # of range, works
        where = "index <= '%s'" % end_dt
        result = store.select('df', where=where)
        tm.assert_frame_equal(expected, result)

        # select w/o iterator and where clause, inclusive range,
        # works
        where = "index >= '%s' & index <= '%s'" % (beg_dt, end_dt)
        result = store.select('df', where=where)
        tm.assert_frame_equal(expected, result)

    # with iterator, full range
    with ensure_clean_store() as store:
        expected = tm.makeTimeDataFrame(100064, 'S')
        _maybe_remove(store, 'df')
        store.append('df', expected)

        beg_dt = expected.index[0]
        end_dt = expected.index[-1]

        # select w/iterator and no where clause works
        results = [s for s in store.select('df', chunksize=chunksize)]
        result = concat(results)
        tm.assert_frame_equal(expected, result)

        # select w/iterator and where clause, single term, begin of range
        where = "index >= '%s'" % beg_dt
        results = [s for s in store.select(
            'df', where=where, chunksize=chunksize)]
        result = concat(results)
        tm.assert_frame_equal(expected, result)

        # select w/iterator and where clause, single term, end of range
        where = "index <= '%s'" % end_dt
        results = [s for s in store.select(
            'df', where=where, chunksize=chunksize)]
        result = concat(results)
        tm.assert_frame_equal(expected, result)

        # select w/iterator and where clause, inclusive range
        where = "index >= '%s' & index <= '%s'" % (beg_dt, end_dt)
        results = [s for s in store.select(
            'df', where=where, chunksize=chunksize)]
        result = concat(results)
        tm.assert_frame_equal(expected, result)


def test_select_iterator_non_complete_8014():
    # GH 8014
    # using iterator and where clause
    chunksize = 1e4

    # with iterator, non complete range
    with ensure_clean_store() as store:
        expected = tm.makeTimeDataFrame(100064, 'S')
        _maybe_remove(store, 'df')
        store.append('df', expected)

        beg_dt = expected.index[1]
        end_dt = expected.index[-2]

        # select w/iterator and where clause, single term, begin of range
        where = "index >= '%s'" % beg_dt
        results = [s for s in store.select(
            'df', where=where, chunksize=chunksize)]
        result = concat(results)
        rexpected = expected[expected.index >= beg_dt]
        tm.assert_frame_equal(rexpected, result)

        # select w/iterator and where clause, single term, end of range
        where = "index <= '%s'" % end_dt
        results = [s for s in store.select(
            'df', where=where, chunksize=chunksize)]
        result = concat(results)
        rexpected = expected[expected.index <= end_dt]
        tm.assert_frame_equal(rexpected, result)

        # select w/iterator and where clause, inclusive range
        where = "index >= '%s' & index <= '%s'" % (beg_dt, end_dt)
        results = [s for s in store.select(
            'df', where=where, chunksize=chunksize)]
        result = concat(results)
        rexpected = expected[(expected.index >= beg_dt) &
                             (expected.index <= end_dt)]
        tm.assert_frame_equal(rexpected, result)

    # with iterator, empty where
    with ensure_clean_store() as store:
        expected = tm.makeTimeDataFrame(100064, 'S')
        _maybe_remove(store, 'df')
        store.append('df', expected)

        end_dt = expected.index[-1]

        # select w/iterator and where clause, single term, begin of range
        where = "index > '%s'" % end_dt
        results = [s for s in store.select(
            'df', where=where, chunksize=chunksize)]
        assert 0 == len(results)


def test_select_iterator_many_empty_frames():
    # GH 8014
    # using iterator and where clause can return many empty
    # frames.
    chunksize = int(1e4)

    # with iterator, range limited to the first chunk
    with ensure_clean_store() as store:
        expected = tm.makeTimeDataFrame(100000, 'S')
        _maybe_remove(store, 'df')
        store.append('df', expected)

        beg_dt = expected.index[0]
        end_dt = expected.index[chunksize - 1]

        # select w/iterator and where clause, single term, begin of range
        where = "index >= '%s'" % beg_dt
        results = [s for s in store.select(
            'df', where=where, chunksize=chunksize)]
        result = concat(results)
        rexpected = expected[expected.index >= beg_dt]
        tm.assert_frame_equal(rexpected, result)

        # select w/iterator and where clause, single term, end of range
        where = "index <= '%s'" % end_dt
        results = [s for s in store.select(
            'df', where=where, chunksize=chunksize)]

        assert len(results) == 1
        result = concat(results)
        rexpected = expected[expected.index <= end_dt]
        tm.assert_frame_equal(rexpected, result)

        # select w/iterator and where clause, inclusive range
        where = "index >= '%s' & index <= '%s'" % (beg_dt, end_dt)
        results = [s for s in store.select(
            'df', where=where, chunksize=chunksize)]

        # should be 1, is 10
        assert len(results) == 1
        result = concat(results)
        rexpected = expected[(expected.index >= beg_dt) &
                             (expected.index <= end_dt)]
        tm.assert_frame_equal(rexpected, result)

        # select w/iterator and where clause which selects
        # *nothing*.
        #
        # To be consistent with Python idiom I suggest this should
        # return [] e.g. `for e in []: print True` never prints
        # True.

        where = "index <= '%s' & index >= '%s'" % (beg_dt, end_dt)
        results = [s for s in store.select(
            'df', where=where, chunksize=chunksize)]

        # should be []
        assert len(results) == 0


def test_select_as_multiple():
    df1 = tm.makeTimeDataFrame()
    df2 = tm.makeTimeDataFrame().rename(columns=lambda x: "%s_2" % x)
    df2['foo'] = 'bar'

    with ensure_clean_store() as store:
        # no tables stored
        pytest.raises(Exception, store.select_as_multiple,
                      None, where=['A>0', 'B>0'], selector='df1')

        store.append('df1', df1, data_columns=['A', 'B'])
        store.append('df2', df2)

        # exceptions
        pytest.raises(Exception, store.select_as_multiple,
                      None, where=['A>0', 'B>0'], selector='df1')
        pytest.raises(Exception, store.select_as_multiple,
                      [None], where=['A>0', 'B>0'], selector='df1')
        pytest.raises(KeyError, store.select_as_multiple,
                      ['df1', 'df3'], where=['A>0', 'B>0'],
                      selector='df1')
        pytest.raises(KeyError, store.select_as_multiple,
                      ['df3'], where=['A>0', 'B>0'], selector='df1')
        pytest.raises(KeyError, store.select_as_multiple,
                      ['df1', 'df2'], where=['A>0', 'B>0'],
                      selector='df4')

        # default select
        result = store.select('df1', ['A>0', 'B>0'])
        expected = store.select_as_multiple(
            ['df1'], where=['A>0', 'B>0'], selector='df1')
        tm.assert_frame_equal(result, expected)
        expected = store.select_as_multiple(
            'df1', where=['A>0', 'B>0'], selector='df1')
        tm.assert_frame_equal(result, expected)

        # multiple
        result = store.select_as_multiple(
            ['df1', 'df2'], where=['A>0', 'B>0'], selector='df1')
        expected = concat([df1, df2], axis=1)
        expected = expected[(expected.A > 0) & (expected.B > 0)]
        tm.assert_frame_equal(result, expected)

        # multiple (diff selector)
        result = store.select_as_multiple(
            ['df1', 'df2'], where='index>df2.index[4]', selector='df2')
        expected = concat([df1, df2], axis=1)
        expected = expected[5:]
        tm.assert_frame_equal(result, expected)

        # test excpection for diff rows
        store.append('df3', tm.makeTimeDataFrame(nper=50))
        pytest.raises(ValueError, store.select_as_multiple,
                      ['df1', 'df3'], where=['A>0', 'B>0'],
                      selector='df1')


def test_select_filter_corner():
    df = DataFrame(np.random.randn(50, 100))
    df.index = ['%.3d' % c for c in df.index]
    df.columns = ['%.3d' % c for c in df.columns]

    with ensure_clean_store() as store:
        store.put('frame', df, format='table')

        crit = 'columns=df.columns[:75]'
        result = store.select('frame', [crit])
        tm.assert_frame_equal(result, df.loc[:, df.columns[:75]])

        crit = 'columns=df.columns[:75:2]'
        result = store.select('frame', [crit])
        tm.assert_frame_equal(result, df.loc[:, df.columns[:75:2]])


def test_panel_select():
    with ensure_clean_store() as store:
        with catch_warnings(record=True):
            wp = tm.makePanel()

            store.put('wp', wp, format='table')
            date = wp.major_axis[len(wp.major_axis) // 2]

            crit1 = ('major_axis>=date')
            crit2 = ("minor_axis=['A', 'D']")

            result = store.select('wp', [crit1, crit2])
            expected = wp.truncate(before=date).reindex(minor=['A', 'D'])
            assert_panel_equal(result, expected)

            result = store.select(
                'wp', ['major_axis>="20000124"',
                       ("minor_axis=['A', 'B']")])
            expected = wp.truncate(
                before='20000124').reindex(minor=['A', 'B'])
            assert_panel_equal(result, expected)


def test_frame_select():
    df = tm.makeTimeDataFrame()

    with ensure_clean_store() as store:
        store.put('frame', df, format='table')
        date = df.index[len(df) // 2]

        crit1 = Term('index>=date')
        assert crit1.env.scope['date'] == date

        crit2 = ("columns=['A', 'D']")
        crit3 = ('columns=A')

        result = store.select('frame', [crit1, crit2])
        expected = df.loc[date:, ['A', 'D']]
        tm.assert_frame_equal(result, expected)

        result = store.select('frame', [crit3])
        expected = df.loc[:, ['A']]
        tm.assert_frame_equal(result, expected)

        # invalid terms
        df = tm.makeTimeDataFrame()
        store.append('df_time', df)
        pytest.raises(ValueError, store.select, 'df_time', "index>0")


def test_frame_select_complex():
    # select via complex criteria
    df = tm.makeTimeDataFrame()
    df['string'] = 'foo'
    df.loc[df.index[0:4], 'string'] = 'bar'

    with ensure_clean_store() as store:
        store.put('df', df, format='table', data_columns=['string'])

        # empty
        result = store.select('df', 'index>df.index[3] & string="bar"')
        expected = df.loc[(df.index > df.index[3]) & (df.string == 'bar')]
        tm.assert_frame_equal(result, expected)

        result = store.select('df', 'index>df.index[3] & string="foo"')
        expected = df.loc[(df.index > df.index[3]) & (df.string == 'foo')]
        tm.assert_frame_equal(result, expected)

        # or
        result = store.select('df', 'index>df.index[3] | string="bar"')
        expected = df.loc[(df.index > df.index[3]) | (df.string == 'bar')]
        tm.assert_frame_equal(result, expected)

        result = store.select('df', '(index>df.index[3] & '
                              'index<=df.index[6]) | string="bar"')
        expected = df.loc[((df.index > df.index[3]) & (
            df.index <= df.index[6])) | (df.string == 'bar')]
        tm.assert_frame_equal(result, expected)

        # invert
        result = store.select('df', 'string!="bar"')
        expected = df.loc[df.string != 'bar']
        tm.assert_frame_equal(result, expected)

        # invert not implemented in numexpr :(
        pytest.raises(NotImplementedError,
                      store.select, 'df', '~(string="bar")')

        # invert ok for filters
        result = store.select('df', "~(columns=['A','B'])")
        expected = df.loc[:, df.columns.difference(['A', 'B'])]
        tm.assert_frame_equal(result, expected)

        # in
        result = store.select(
            'df', "index>df.index[3] & columns in ['A','B']")
        expected = df.loc[df.index > df.index[3]].reindex(columns=[
                                                          'A', 'B'])
        tm.assert_frame_equal(result, expected)
