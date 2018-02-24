import pytest
from warnings import catch_warnings
import datetime
from datetime import timedelta

import numpy as np
import pandas as pd
from pandas import (Series, DataFrame, date_range, Timestamp,
                    concat, MultiIndex, Panel)
from .common import ensure_clean_store, ensure_clean_path, _maybe_remove

import pandas.util.testing as tm
from pandas.io.pytables import read_hdf


def test_append():
    with ensure_clean_store() as store:
        # this is allowed by almost always don't want to do it
        # tables.NaturalNameWarning):
        with catch_warnings(record=True):

            df = tm.makeTimeDataFrame()
            _maybe_remove(store, 'df1')
            store.append('df1', df[:10])
            store.append('df1', df[10:])
            tm.assert_frame_equal(store['df1'], df)

            _maybe_remove(store, 'df2')
            store.put('df2', df[:10], format='table')
            store.append('df2', df[10:])
            tm.assert_frame_equal(store['df2'], df)

            _maybe_remove(store, 'df3')
            store.append('/df3', df[:10])
            store.append('/df3', df[10:])
            tm.assert_frame_equal(store['df3'], df)

            # this is allowed by almost always don't want to do it
            # tables.NaturalNameWarning
            _maybe_remove(store, '/df3 foo')
            store.append('/df3 foo', df[:10])
            store.append('/df3 foo', df[10:])
            tm.assert_frame_equal(store['df3 foo'], df)

            # panel
            wp = tm.makePanel()
            _maybe_remove(store, 'wp1')
            store.append('wp1', wp.iloc[:, :10, :])
            store.append('wp1', wp.iloc[:, 10:, :])
            tm.assert_panel_equal(store['wp1'], wp)

            # test using differt order of items on the non-index axes
            _maybe_remove(store, 'wp1')
            wp_append1 = wp.iloc[:, :10, :]
            store.append('wp1', wp_append1)
            wp_append2 = wp.iloc[:, 10:, :].reindex(items=wp.items[::-1])
            store.append('wp1', wp_append2)
            tm.assert_panel_equal(store['wp1'], wp)

            # dtype issues - mizxed type in a single object column
            df = DataFrame(data=[[1, 2], [0, 1], [1, 2], [0, 0]])
            df['mixed_column'] = 'testing'
            df.loc[2, 'mixed_column'] = np.nan
            _maybe_remove(store, 'df')
            store.append('df', df)
            tm.assert_frame_equal(store['df'], df)

            # uints - test storage of uints
            uint_data = DataFrame({
                'u08': Series(np.random.randint(0, high=255, size=5),
                              dtype=np.uint8),
                'u16': Series(np.random.randint(0, high=65535, size=5),
                              dtype=np.uint16),
                'u32': Series(np.random.randint(0, high=2**30, size=5),
                              dtype=np.uint32),
                'u64': Series([2**58, 2**59, 2**60, 2**61, 2**62],
                              dtype=np.uint64)}, index=np.arange(5))
            _maybe_remove(store, 'uints')
            store.append('uints', uint_data)
            tm.assert_frame_equal(store['uints'], uint_data)

            # uints - test storage of uints in indexable columns
            _maybe_remove(store, 'uints')
            # 64-bit indices not yet supported
            store.append('uints', uint_data, data_columns=[
                         'u08', 'u16', 'u32'])
            tm.assert_frame_equal(store['uints'], uint_data)


def test_append_series():
    with ensure_clean_store() as store:
        # basic
        ss = tm.makeStringSeries()
        ts = tm.makeTimeSeries()
        ns = Series(np.arange(100))

        store.append('ss', ss)
        result = store['ss']
        tm.assert_series_equal(result, ss)
        assert result.name is None

        store.append('ts', ts)
        result = store['ts']
        tm.assert_series_equal(result, ts)
        assert result.name is None

        ns.name = 'foo'
        store.append('ns', ns)
        result = store['ns']
        tm.assert_series_equal(result, ns)
        assert result.name == ns.name

        # select on the values
        expected = ns[ns > 60]
        result = store.select('ns', 'foo>60')
        tm.assert_series_equal(result, expected)

        # select on the index and values
        expected = ns[(ns > 70) & (ns.index < 90)]
        result = store.select('ns', 'foo>70 and index<90')
        tm.assert_series_equal(result, expected)

        # multi-index
        mi = DataFrame(np.random.randn(5, 1), columns=['A'])
        mi['B'] = np.arange(len(mi))
        mi['C'] = 'foo'
        mi.loc[3:5, 'C'] = 'bar'
        mi.set_index(['C', 'B'], inplace=True)
        s = mi.stack()
        s.index = s.index.droplevel(2)
        store.append('mi', s)
        tm.assert_series_equal(store['mi'], s)


def test_append_some_nans():
    with ensure_clean_store() as store:
        df = DataFrame({'A': Series(np.random.randn(20)).astype('int32'),
                        'A1': np.random.randn(20),
                        'A2': np.random.randn(20),
                        'B': 'foo', 'C': 'bar',
                        'D': Timestamp("20010101"),
                        'E': datetime.datetime(2001, 1, 2, 0, 0)},
                       index=np.arange(20))
        # some nans
        _maybe_remove(store, 'df1')
        df.loc[0:15, ['A1', 'B', 'D', 'E']] = np.nan
        store.append('df1', df[:10])
        store.append('df1', df[10:])
        tm.assert_frame_equal(store['df1'], df)

        # first column
        df1 = df.copy()
        df1.loc[:, 'A1'] = np.nan
        _maybe_remove(store, 'df1')
        store.append('df1', df1[:10])
        store.append('df1', df1[10:])
        tm.assert_frame_equal(store['df1'], df1)

        # 2nd column
        df2 = df.copy()
        df2.loc[:, 'A2'] = np.nan
        _maybe_remove(store, 'df2')
        store.append('df2', df2[:10])
        store.append('df2', df2[10:])
        tm.assert_frame_equal(store['df2'], df2)

        # datetimes
        df3 = df.copy()
        df3.loc[:, 'E'] = np.nan
        _maybe_remove(store, 'df3')
        store.append('df3', df3[:10])
        store.append('df3', df3[10:])
        tm.assert_frame_equal(store['df3'], df3)


def test_append_all_nans():
    with ensure_clean_store() as store:
        df = DataFrame({'A1': np.random.randn(20),
                        'A2': np.random.randn(20)},
                       index=np.arange(20))
        df.loc[0:15, :] = np.nan

        # nan some entire rows (dropna=True)
        _maybe_remove(store, 'df')
        store.append('df', df[:10], dropna=True)
        store.append('df', df[10:], dropna=True)
        tm.assert_frame_equal(store['df'], df[-4:])

        # nan some entire rows (dropna=False)
        _maybe_remove(store, 'df2')
        store.append('df2', df[:10], dropna=False)
        store.append('df2', df[10:], dropna=False)
        tm.assert_frame_equal(store['df2'], df)

        # tests the option io.hdf.dropna_table
        pd.set_option('io.hdf.dropna_table', False)
        _maybe_remove(store, 'df3')
        store.append('df3', df[:10])
        store.append('df3', df[10:])
        tm.assert_frame_equal(store['df3'], df)

        pd.set_option('io.hdf.dropna_table', True)
        _maybe_remove(store, 'df4')
        store.append('df4', df[:10])
        store.append('df4', df[10:])
        tm.assert_frame_equal(store['df4'], df[-4:])

        # nan some entire rows (string are still written!)
        df = DataFrame({'A1': np.random.randn(20),
                        'A2': np.random.randn(20),
                        'B': 'foo', 'C': 'bar'},
                       index=np.arange(20))

        df.loc[0:15, :] = np.nan

        _maybe_remove(store, 'df')
        store.append('df', df[:10], dropna=True)
        store.append('df', df[10:], dropna=True)
        tm.assert_frame_equal(store['df'], df)

        _maybe_remove(store, 'df2')
        store.append('df2', df[:10], dropna=False)
        store.append('df2', df[10:], dropna=False)
        tm.assert_frame_equal(store['df2'], df)

        # nan some entire rows (but since we have dates they are still
        # written!)
        df = DataFrame({'A1': np.random.randn(20),
                        'A2': np.random.randn(20),
                        'B': 'foo', 'C': 'bar',
                        'D': Timestamp("20010101"),
                        'E': datetime.datetime(2001, 1, 2, 0, 0)},
                       index=np.arange(20))

        df.loc[0:15, :] = np.nan

        _maybe_remove(store, 'df')
        store.append('df', df[:10], dropna=True)
        store.append('df', df[10:], dropna=True)
        tm.assert_frame_equal(store['df'], df)

        _maybe_remove(store, 'df2')
        store.append('df2', df[:10], dropna=False)
        store.append('df2', df[10:], dropna=False)
        tm.assert_frame_equal(store['df2'], df)

    # Test to make sure defaults are to not drop.
    # Corresponding to Issue 9382
    df_with_missing = DataFrame(
        {'col1': [0, np.nan, 2], 'col2': [1, np.nan, np.nan]})

    with ensure_clean_path() as path:
        df_with_missing.to_hdf(path, 'df_with_missing', format='table')
        reloaded = read_hdf(path, 'df_with_missing')
        tm.assert_frame_equal(df_with_missing, reloaded)

    matrix = [[[np.nan, np.nan, np.nan], [1, np.nan, np.nan]],
              [[np.nan, np.nan, np.nan], [np.nan, 5, 6]],
              [[np.nan, np.nan, np.nan], [np.nan, 3, np.nan]]]

    with catch_warnings(record=True):
        panel_with_missing = Panel(matrix,
                                   items=['Item1', 'Item2', 'Item3'],
                                   major_axis=[1, 2],
                                   minor_axis=['A', 'B', 'C'])

        with ensure_clean_path() as path:
            panel_with_missing.to_hdf(
                path, 'panel_with_missing', format='table')
            reloaded_panel = read_hdf(path, 'panel_with_missing')
            tm.assert_panel_equal(panel_with_missing, reloaded_panel)


def test_append_frame_column_oriented():
    with ensure_clean_store() as store:
        # column oriented
        df = tm.makeTimeDataFrame()
        _maybe_remove(store, 'df1')
        store.append('df1', df.iloc[:, :2], axes=['columns'])
        store.append('df1', df.iloc[:, 2:])
        tm.assert_frame_equal(store['df1'], df)

        result = store.select('df1', 'columns=A')
        expected = df.reindex(columns=['A'])
        tm.assert_frame_equal(expected, result)

        # selection on the non-indexable
        result = store.select(
            'df1', ('columns=A', 'index=df.index[0:4]'))
        expected = df.reindex(columns=['A'], index=df.index[0:4])
        tm.assert_frame_equal(expected, result)

        # this isn't supported
        with pytest.raises(TypeError):
            store.select('df1',
                         'columns=A and index>df.index[4]')


def test_append_with_different_block_ordering():
    # GH 4096; using same frames, but different block orderings
    with ensure_clean_store() as store:
        for i in range(10):
            df = DataFrame(np.random.randn(10, 2), columns=list('AB'))
            df['index'] = range(10)
            df['index'] += i * 10
            df['int64'] = Series([1] * len(df), dtype='int64')
            df['int16'] = Series([1] * len(df), dtype='int16')

            if i % 2 == 0:
                del df['int64']
                df['int64'] = Series([1] * len(df), dtype='int64')
            if i % 3 == 0:
                a = df.pop('A')
                df['A'] = a

            df.set_index('index', inplace=True)
            store.append('df', df)

    # test a different ordering but with more fields (like invalid
    # combinate)
    with ensure_clean_store() as store:
        df = DataFrame(np.random.randn(10, 2),
                       columns=list('AB'), dtype='float64')
        df['int64'] = Series([1] * len(df), dtype='int64')
        df['int16'] = Series([1] * len(df), dtype='int16')
        store.append('df', df)

        # store additional fields in different blocks
        df['int16_2'] = Series([1] * len(df), dtype='int16')
        pytest.raises(ValueError, store.append, 'df', df)

        # store multile additional fields in different blocks
        df['float_3'] = Series([1.] * len(df), dtype='float64')
        pytest.raises(ValueError, store.append, 'df', df)


def test_append_with_strings():
    with ensure_clean_store() as store:
        with catch_warnings(record=True):
            wp = tm.makePanel()
            wp2 = wp.rename_axis(
                {x: "%s_extra" % x for x in wp.minor_axis}, axis=2)

            def check_col(key, name, size):
                assert getattr(store.get_storer(key)
                               .table.description, name).itemsize == size

            store.append('s1', wp, min_itemsize=20)
            store.append('s1', wp2)
            expected = concat([wp, wp2], axis=2)
            expected = expected.reindex(
                minor_axis=sorted(expected.minor_axis))
            tm.assert_panel_equal(store['s1'], expected)
            check_col('s1', 'minor_axis', 20)

            # test dict format
            store.append('s2', wp, min_itemsize={'minor_axis': 20})
            store.append('s2', wp2)
            expected = concat([wp, wp2], axis=2)
            expected = expected.reindex(
                minor_axis=sorted(expected.minor_axis))
            tm.assert_panel_equal(store['s2'], expected)
            check_col('s2', 'minor_axis', 20)

            # apply the wrong field (similar to #1)
            store.append('s3', wp, min_itemsize={'major_axis': 20})
            pytest.raises(ValueError, store.append, 's3', wp2)

            # test truncation of bigger strings
            store.append('s4', wp)
            pytest.raises(ValueError, store.append, 's4', wp2)

            # avoid truncation on elements
            df = DataFrame([[123, 'asdqwerty'], [345, 'dggnhebbsdfbdfb']])
            store.append('df_big', df)
            tm.assert_frame_equal(store.select('df_big'), df)
            check_col('df_big', 'values_block_1', 15)

            # appending smaller string ok
            df2 = DataFrame([[124, 'asdqy'], [346, 'dggnhefbdfb']])
            store.append('df_big', df2)
            expected = concat([df, df2])
            tm.assert_frame_equal(store.select('df_big'), expected)
            check_col('df_big', 'values_block_1', 15)

            # avoid truncation on elements
            df = DataFrame([[123, 'asdqwerty'], [345, 'dggnhebbsdfbdfb']])
            store.append('df_big2', df, min_itemsize={'values': 50})
            tm.assert_frame_equal(store.select('df_big2'), df)
            check_col('df_big2', 'values_block_1', 50)

            # bigger string on next append
            store.append('df_new', df)
            df_new = DataFrame(
                [[124, 'abcdefqhij'], [346, 'abcdefghijklmnopqrtsuvwxyz']])
            pytest.raises(ValueError, store.append, 'df_new', df_new)

            # min_itemsize on Series index (GH 11412)
            df = tm.makeMixedDataFrame().set_index('C')
            store.append('ss', df['B'], min_itemsize={'index': 4})
            tm.assert_series_equal(store.select('ss'), df['B'])

            # same as above, with data_columns=True
            store.append('ss2', df['B'], data_columns=True,
                         min_itemsize={'index': 4})
            tm.assert_series_equal(store.select('ss2'), df['B'])

            # min_itemsize in index without appending (GH 10381)
            store.put('ss3', df, format='table',
                      min_itemsize={'index': 6})
            # just make sure there is a longer string:
            df2 = df.copy().reset_index().assign(C='longer').set_index('C')
            store.append('ss3', df2)
            tm.assert_frame_equal(store.select('ss3'),
                                  pd.concat([df, df2]))

            # same as above, with a Series
            store.put('ss4', df['B'], format='table',
                      min_itemsize={'index': 6})
            store.append('ss4', df2['B'])
            tm.assert_series_equal(store.select('ss4'),
                                   pd.concat([df['B'], df2['B']]))

            # with nans
            _maybe_remove(store, 'df')
            df = tm.makeTimeDataFrame()
            df['string'] = 'foo'
            df.loc[1:4, 'string'] = np.nan
            df['string2'] = 'bar'
            df.loc[4:8, 'string2'] = np.nan
            df['string3'] = 'bah'
            df.loc[1:, 'string3'] = np.nan
            store.append('df', df)
            result = store.select('df')
            tm.assert_frame_equal(result, df)

    with ensure_clean_store() as store:
        def check_col(key, name, size):
            assert getattr(store.get_storer(key)
                           .table.description, name).itemsize, size

        df = DataFrame(dict(A='foo', B='bar'), index=range(10))

        # a min_itemsize that creates a data_column
        _maybe_remove(store, 'df')
        store.append('df', df, min_itemsize={'A': 200})
        check_col('df', 'A', 200)
        assert store.get_storer('df').data_columns == ['A']

        # a min_itemsize that creates a data_column2
        _maybe_remove(store, 'df')
        store.append('df', df, data_columns=['B'], min_itemsize={'A': 200})
        check_col('df', 'A', 200)
        assert store.get_storer('df').data_columns == ['B', 'A']

        # a min_itemsize that creates a data_column2
        _maybe_remove(store, 'df')
        store.append('df', df, data_columns=[
                     'B'], min_itemsize={'values': 200})
        check_col('df', 'B', 200)
        check_col('df', 'values_block_0', 200)
        assert store.get_storer('df').data_columns == ['B']

        # infer the .typ on subsequent appends
        _maybe_remove(store, 'df')
        store.append('df', df[:5], min_itemsize=200)
        store.append('df', df[5:], min_itemsize=200)
        tm.assert_frame_equal(store['df'], df)

        # invalid min_itemsize keys
        df = DataFrame(['foo', 'foo', 'foo', 'barh',
                        'barh', 'barh'], columns=['A'])
        _maybe_remove(store, 'df')
        pytest.raises(ValueError, store.append, 'df',
                      df, min_itemsize={'foo': 20, 'foobar': 20})


def test_append_with_data_columns():
    with ensure_clean_store() as store:
        df = tm.makeTimeDataFrame()
        df.iloc[0, df.columns.get_loc('B')] = 1.
        _maybe_remove(store, 'df')
        store.append('df', df[:2], data_columns=['B'])
        store.append('df', df[2:])
        tm.assert_frame_equal(store['df'], df)

        # check that we have indicies created
        assert(store._handle.root.df.table.cols.index.is_indexed is True)
        assert(store._handle.root.df.table.cols.B.is_indexed is True)

        # data column searching
        result = store.select('df', 'B>0')
        expected = df[df.B > 0]
        tm.assert_frame_equal(result, expected)

        # data column searching (with an indexable and a data_columns)
        result = store.select(
            'df', 'B>0 and index>df.index[3]')
        df_new = df.reindex(index=df.index[4:])
        expected = df_new[df_new.B > 0]
        tm.assert_frame_equal(result, expected)

        # data column selection with a string data_column
        df_new = df.copy()
        df_new['string'] = 'foo'
        df_new.loc[1:4, 'string'] = np.nan
        df_new.loc[5:6, 'string'] = 'bar'
        _maybe_remove(store, 'df')
        store.append('df', df_new, data_columns=['string'])
        result = store.select('df', "string='foo'")
        expected = df_new[df_new.string == 'foo']
        tm.assert_frame_equal(result, expected)

        # using min_itemsize and a data column
        def check_col(key, name, size):
            assert getattr(store.get_storer(key)
                           .table.description, name).itemsize == size

    with ensure_clean_store() as store:
        _maybe_remove(store, 'df')
        store.append('df', df_new, data_columns=['string'],
                     min_itemsize={'string': 30})
        check_col('df', 'string', 30)
        _maybe_remove(store, 'df')
        store.append(
            'df', df_new, data_columns=['string'], min_itemsize=30)
        check_col('df', 'string', 30)
        _maybe_remove(store, 'df')
        store.append('df', df_new, data_columns=['string'],
                     min_itemsize={'values': 30})
        check_col('df', 'string', 30)

    with ensure_clean_store() as store:
        df_new['string2'] = 'foobarbah'
        df_new['string_block1'] = 'foobarbah1'
        df_new['string_block2'] = 'foobarbah2'
        _maybe_remove(store, 'df')
        store.append('df', df_new, data_columns=['string', 'string2'],
                     min_itemsize={'string': 30, 'string2': 40,
                                   'values': 50})
        check_col('df', 'string', 30)
        check_col('df', 'string2', 40)
        check_col('df', 'values_block_1', 50)

    with ensure_clean_store() as store:
        # multiple data columns
        df_new = df.copy()
        df_new.iloc[0, df_new.columns.get_loc('A')] = 1.
        df_new.iloc[0, df_new.columns.get_loc('B')] = -1.
        df_new['string'] = 'foo'

        sl = df_new.columns.get_loc('string')
        df_new.iloc[1:4, sl] = np.nan
        df_new.iloc[5:6, sl] = 'bar'

        df_new['string2'] = 'foo'
        sl = df_new.columns.get_loc('string2')
        df_new.iloc[2:5, sl] = np.nan
        df_new.iloc[7:8, sl] = 'bar'
        _maybe_remove(store, 'df')
        store.append(
            'df', df_new, data_columns=['A', 'B', 'string', 'string2'])
        result = store.select('df',
                              "string='foo' and string2='foo'"
                              " and A>0 and B<0")
        expected = df_new[(df_new.string == 'foo') & (
            df_new.string2 == 'foo') & (df_new.A > 0) & (df_new.B < 0)]
        tm.assert_frame_equal(result, expected, check_index_type=False)

        # yield an empty frame
        result = store.select('df', "string='foo' and string2='cool'")
        expected = df_new[(df_new.string == 'foo') & (
            df_new.string2 == 'cool')]
        tm.assert_frame_equal(result, expected, check_index_type=False)

    with ensure_clean_store() as store:
        # doc example
        df_dc = df.copy()
        df_dc['string'] = 'foo'
        df_dc.loc[4:6, 'string'] = np.nan
        df_dc.loc[7:9, 'string'] = 'bar'
        df_dc['string2'] = 'cool'
        df_dc['datetime'] = Timestamp('20010102')
        df_dc = df_dc._convert(datetime=True)
        df_dc.loc[3:5, ['A', 'B', 'datetime']] = np.nan

        _maybe_remove(store, 'df_dc')
        store.append('df_dc', df_dc,
                     data_columns=['B', 'C', 'string',
                                   'string2', 'datetime'])
        result = store.select('df_dc', 'B>0')

        expected = df_dc[df_dc.B > 0]
        tm.assert_frame_equal(result, expected, check_index_type=False)

        result = store.select(
            'df_dc', ['B > 0', 'C > 0', 'string == foo'])
        expected = df_dc[(df_dc.B > 0) & (df_dc.C > 0) & (
            df_dc.string == 'foo')]
        tm.assert_frame_equal(result, expected, check_index_type=False)

    with ensure_clean_store() as store:
        # doc example part 2
        np.random.seed(1234)
        index = date_range('1/1/2000', periods=8)
        df_dc = DataFrame(np.random.randn(8, 3), index=index,
                          columns=['A', 'B', 'C'])
        df_dc['string'] = 'foo'
        df_dc.loc[4:6, 'string'] = np.nan
        df_dc.loc[7:9, 'string'] = 'bar'
        df_dc.loc[:, ['B', 'C']] = df_dc.loc[:, ['B', 'C']].abs()
        df_dc['string2'] = 'cool'

        # on-disk operations
        store.append('df_dc', df_dc, data_columns=[
                     'B', 'C', 'string', 'string2'])

        result = store.select('df_dc', 'B>0')
        expected = df_dc[df_dc.B > 0]
        tm.assert_frame_equal(result, expected)

        result = store.select(
            'df_dc', ['B > 0', 'C > 0', 'string == "foo"'])
        expected = df_dc[(df_dc.B > 0) & (df_dc.C > 0) &
                         (df_dc.string == 'foo')]
        tm.assert_frame_equal(result, expected)

    with ensure_clean_store() as store:
        with catch_warnings(record=True):
            # panel
            # GH5717 not handling data_columns
            np.random.seed(1234)
            p = tm.makePanel()

            store.append('p1', p)
            tm.assert_panel_equal(store.select('p1'), p)

            store.append('p2', p, data_columns=True)
            tm.assert_panel_equal(store.select('p2'), p)

            result = store.select('p2', where='ItemA>0')
            expected = p.to_frame()
            expected = expected[expected['ItemA'] > 0]
            tm.assert_frame_equal(result.to_frame(), expected)

            result = store.select(
                'p2', where='ItemA>0 & minor_axis=["A","B"]')
            expected = p.to_frame()
            expected = expected[expected['ItemA'] > 0]
            expected = expected[expected.reset_index(
                level=['major']).index.isin(['A', 'B'])]
            tm.assert_frame_equal(result.to_frame(), expected)


def test_append_diff_item_order():
    with catch_warnings(record=True):
        wp = tm.makePanel()
        wp1 = wp.iloc[:, :10, :]
        wp2 = wp.iloc[wp.items.get_indexer(['ItemC', 'ItemB', 'ItemA']),
                      10:, :]

        with ensure_clean_store() as store:
            store.put('panel', wp1, format='table')
            pytest.raises(ValueError, store.put, 'panel', wp2,
                          append=True)


def test_append_hierarchical():
    index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                               ['one', 'two', 'three']],
                       labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                               [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                       names=['foo', 'bar'])
    df = DataFrame(np.random.randn(10, 3), index=index,
                   columns=['A', 'B', 'C'])

    with ensure_clean_store() as store:
        store.append('mi', df)
        result = store.select('mi')
        tm.assert_frame_equal(result, df)

        # GH 3748
        result = store.select('mi', columns=['A', 'B'])
        expected = df.reindex(columns=['A', 'B'])
        tm.assert_frame_equal(result, expected)

    with ensure_clean_path('test.hdf') as path:
        df.to_hdf(path, 'df', format='table')
        result = read_hdf(path, 'df', columns=['A', 'B'])
        expected = df.reindex(columns=['A', 'B'])
        tm.assert_frame_equal(result, expected)


def test_append_misc():
    with ensure_clean_store() as store:
        df = tm.makeDataFrame()
        store.append('df', df, chunksize=1)
        result = store.select('df')
        tm.assert_frame_equal(result, df)

        store.append('df1', df, expectedrows=10)
        result = store.select('df1')
        tm.assert_frame_equal(result, df)

    # more chunksize in append tests
    def check(obj, comparator):
        for c in [10, 200, 1000]:
            with ensure_clean_store(mode='w') as store:
                store.append('obj', obj, chunksize=c)
                result = store.select('obj')
                comparator(result, obj)

    df = tm.makeDataFrame()
    df['string'] = 'foo'
    df['float322'] = 1.
    df['float322'] = df['float322'].astype('float32')
    df['bool'] = df['float322'] > 0
    df['time1'] = Timestamp('20130101')
    df['time2'] = Timestamp('20130102')
    check(df, tm.assert_frame_equal)

    with catch_warnings(record=True):
        p = tm.makePanel()
        check(p, tm.assert_panel_equal)

    # empty frame, GH4273
    with ensure_clean_store() as store:
        # 0 len
        df_empty = DataFrame(columns=list('ABC'))
        store.append('df', df_empty)
        pytest.raises(KeyError, store.select, 'df')

        # repeated append of 0/non-zero frames
        df = DataFrame(np.random.rand(10, 3), columns=list('ABC'))
        store.append('df', df)
        tm.assert_frame_equal(store.select('df'), df)
        store.append('df', df_empty)
        tm.assert_frame_equal(store.select('df'), df)

        # store
        df = DataFrame(columns=list('ABC'))
        store.put('df2', df)
        tm.assert_frame_equal(store.select('df2'), df)

        with catch_warnings(record=True):
            # 0 len
            p_empty = Panel(items=list('ABC'))
            store.append('p', p_empty)
            pytest.raises(KeyError, store.select, 'p')

            # repeated append of 0/non-zero frames
            p = Panel(np.random.randn(3, 4, 5), items=list('ABC'))
            store.append('p', p)
            tm.assert_panel_equal(store.select('p'), p)
            store.append('p', p_empty)
            tm.assert_panel_equal(store.select('p'), p)

            # store
            store.put('p2', p_empty)
            tm.assert_panel_equal(store.select('p2'), p_empty)


def test_append_raise():
    with ensure_clean_store() as store:
        # test append with invalid input to get good error messages
        # list in column
        df = tm.makeDataFrame()
        df['invalid'] = [['a']] * len(df)
        assert df.dtypes['invalid'] == np.object_
        pytest.raises(TypeError, store.append, 'df', df)

        # multiple invalid columns
        df['invalid2'] = [['a']] * len(df)
        df['invalid3'] = [['a']] * len(df)
        pytest.raises(TypeError, store.append, 'df', df)

        # datetime with embedded nans as object
        df = tm.makeDataFrame()
        s = Series(datetime.datetime(2001, 1, 2), index=df.index)
        s = s.astype(object)
        s[0:5] = np.nan
        df['invalid'] = s
        assert df.dtypes['invalid'] == np.object_
        pytest.raises(TypeError, store.append, 'df', df)

        # directly ndarray
        pytest.raises(TypeError, store.append, 'df', np.arange(10))

        # series directly
        pytest.raises(TypeError, store.append,
                      'df', Series(np.arange(10)))

        # appending an incompatible table
        df = tm.makeDataFrame()
        store.append('df', df)

        df['foo'] = 'foo'
        pytest.raises(ValueError, store.append, 'df', df)


def test_append_with_timedelta():
    # GH 3577
    # append timedelta
    df = DataFrame(dict(A=Timestamp('20130101'), B=[Timestamp(
        '20130101') + timedelta(days=i, seconds=10) for i in range(10)]))
    df['C'] = df['A'] - df['B']
    df.loc[3:5, 'C'] = np.nan

    with ensure_clean_store() as store:
        # table
        _maybe_remove(store, 'df')
        store.append('df', df, data_columns=True)
        result = store.select('df')
        tm.assert_frame_equal(result, df)

        result = store.select('df', where="C<100000")
        tm.assert_frame_equal(result, df)

        result = store.select('df', where="C<pd.Timedelta('-3D')")
        tm.assert_frame_equal(result, df.iloc[3:])

        result = store.select('df', "C<'-3D'")
        tm.assert_frame_equal(result, df.iloc[3:])

        # a bit hacky here as we don't really deal with the NaT properly

        result = store.select('df', "C<'-500000s'")
        result = result.dropna(subset=['C'])
        tm.assert_frame_equal(result, df.iloc[6:])

        result = store.select('df', "C<'-3.5D'")
        result = result.iloc[1:]
        tm.assert_frame_equal(result, df.iloc[4:])

        # fixed
        _maybe_remove(store, 'df2')
        store.put('df2', df)
        result = store.select('df2')
        tm.assert_frame_equal(result, df)


def test_append_with_diff_col_name_types_raises_value_error():
    df = DataFrame(np.random.randn(10, 1))
    df2 = DataFrame({'a': np.random.randn(10)})
    df3 = DataFrame({(1, 2): np.random.randn(10)})
    df4 = DataFrame({('1', 2): np.random.randn(10)})
    df5 = DataFrame({('1', 2, object): np.random.randn(10)})

    with ensure_clean_store() as store:
        name = 'df_%s' % tm.rands(10)
        store.append(name, df)

        for d in (df2, df3, df4, df5):
            with pytest.raises(ValueError):
                store.append(name, d)


def test_append_to_multiple():
    df1 = tm.makeTimeDataFrame()
    df2 = tm.makeTimeDataFrame().rename(columns=lambda x: "%s_2" % x)
    df2['foo'] = 'bar'
    df = concat([df1, df2], axis=1)

    with ensure_clean_store() as store:
        # exceptions
        pytest.raises(ValueError, store.append_to_multiple,
                      {'df1': ['A', 'B'], 'df2': None}, df,
                      selector='df3')
        pytest.raises(ValueError, store.append_to_multiple,
                      {'df1': None, 'df2': None}, df, selector='df3')
        pytest.raises(
            ValueError, store.append_to_multiple, 'df1', df, 'df1')

        # regular operation
        store.append_to_multiple(
            {'df1': ['A', 'B'], 'df2': None}, df, selector='df1')
        result = store.select_as_multiple(
            ['df1', 'df2'], where=['A>0', 'B>0'], selector='df1')
        expected = df[(df.A > 0) & (df.B > 0)]
        tm.assert_frame_equal(result, expected)


def test_append_to_multiple_dropna():
    df1 = tm.makeTimeDataFrame()
    df2 = tm.makeTimeDataFrame().rename(columns=lambda x: "%s_2" % x)
    df1.iloc[1, df1.columns.get_indexer(['A', 'B'])] = np.nan
    df = concat([df1, df2], axis=1)

    with ensure_clean_store() as store:
        # dropna=True should guarantee rows are synchronized
        store.append_to_multiple(
            {'df1': ['A', 'B'], 'df2': None}, df, selector='df1',
            dropna=True)
        result = store.select_as_multiple(['df1', 'df2'])
        expected = df.dropna()
        tm.assert_frame_equal(result, expected)
        tm.assert_index_equal(store.select('df1').index,
                              store.select('df2').index)


@pytest.mark.xfail(run=False,
                   reason="append_to_multiple_dropna_false "
                   "is not raising as failed")
def test_append_to_multiple_dropna_false():
    df1 = tm.makeTimeDataFrame()
    df2 = tm.makeTimeDataFrame().rename(columns=lambda x: "%s_2" % x)
    df1.iloc[1, df1.columns.get_indexer(['A', 'B'])] = np.nan
    df = concat([df1, df2], axis=1)

    with ensure_clean_store() as store:
        # dropna=False shouldn't synchronize row indexes
        store.append_to_multiple(
            {'df1a': ['A', 'B'], 'df2a': None}, df, selector='df1a',
            dropna=False)

        with pytest.raises(ValueError):
            store.select_as_multiple(['df1a', 'df2a'])

        assert not store.select('df1a').index.equals(
            store.select('df2a').index)
