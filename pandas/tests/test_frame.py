# -*- coding: utf-8 -*-

from __future__ import print_function
# pylint: disable-msg=W0612,E1101
from copy import deepcopy
from datetime import datetime, timedelta, time, date
import sys
import re
import nose
import itertools
from distutils.version import LooseVersion

from pandas.compat import(
    map, zip, range, long, lrange, lmap, lzip,
    OrderedDict, u, StringIO, is_platform_windows
)
from pandas import compat

from numpy import random, nan, inf
from numpy.random import randn
import numpy as np

import pandas.core.common as com
import pandas.core.format as fmt
import pandas.core.datetools as datetools
from pandas import (DataFrame, Index, Series, Panel, notnull, isnull,
                    MultiIndex, DatetimeIndex, Timestamp, date_range,
                    read_csv, timedelta_range, Timedelta, option_context,
                    period_range)
import pandas as pd

from pandas.util.testing import (assert_almost_equal,
                                 assert_equal,
                                 assert_series_equal,
                                 assert_frame_equal,
                                 assertRaisesRegexp,
                                 makeCustomDataframe as mkdf,
                                 SubclassedDataFrame)

import pandas.util.testing as tm

from numpy.testing.decorators import slow

from pandas.tests.frame.common import (TestData, _check_mixed_float,
                                       _check_mixed_int)

# ---------------------------------------------------------------------
# DataFrame test cases

JOIN_TYPES = ['inner', 'outer', 'left', 'right']


class SafeForSparse(object):

    _multiprocess_can_split_ = True

    def test_copy_index_name_checking(self):
        # don't want to be able to modify the index stored elsewhere after
        # making a copy
        for attr in ('index', 'columns'):
            ind = getattr(self.frame, attr)
            ind.name = None
            cp = self.frame.copy()
            getattr(cp, attr).name = 'foo'
            self.assertIsNone(getattr(self.frame, attr).name)

    def test_getitem_pop_assign_name(self):
        s = self.frame['A']
        self.assertEqual(s.name, 'A')

        s = self.frame.pop('A')
        self.assertEqual(s.name, 'A')

        s = self.frame.ix[:, 'B']
        self.assertEqual(s.name, 'B')

        s2 = s.ix[:]
        self.assertEqual(s2.name, 'B')

    def test_get_value(self):
        for idx in self.frame.index:
            for col in self.frame.columns:
                result = self.frame.get_value(idx, col)
                expected = self.frame[col][idx]
                assert_almost_equal(result, expected)

    def test_join_index(self):
        # left / right

        f = self.frame.reindex(columns=['A', 'B'])[:10]
        f2 = self.frame.reindex(columns=['C', 'D'])

        joined = f.join(f2)
        self.assertTrue(f.index.equals(joined.index))
        self.assertEqual(len(joined.columns), 4)

        joined = f.join(f2, how='left')
        self.assertTrue(joined.index.equals(f.index))
        self.assertEqual(len(joined.columns), 4)

        joined = f.join(f2, how='right')
        self.assertTrue(joined.index.equals(f2.index))
        self.assertEqual(len(joined.columns), 4)

        # inner

        f = self.frame.reindex(columns=['A', 'B'])[:10]
        f2 = self.frame.reindex(columns=['C', 'D'])

        joined = f.join(f2, how='inner')
        self.assertTrue(joined.index.equals(f.index.intersection(f2.index)))
        self.assertEqual(len(joined.columns), 4)

        # outer

        f = self.frame.reindex(columns=['A', 'B'])[:10]
        f2 = self.frame.reindex(columns=['C', 'D'])

        joined = f.join(f2, how='outer')
        self.assertTrue(tm.equalContents(self.frame.index, joined.index))
        self.assertEqual(len(joined.columns), 4)

        assertRaisesRegexp(ValueError, 'join method', f.join, f2, how='foo')

        # corner case - overlapping columns
        for how in ('outer', 'left', 'inner'):
            with assertRaisesRegexp(ValueError, 'columns overlap but no suffix'):
                self.frame.join(self.frame, how=how)

    def test_join_index_more(self):
        af = self.frame.ix[:, ['A', 'B']]
        bf = self.frame.ix[::2, ['C', 'D']]

        expected = af.copy()
        expected['C'] = self.frame['C'][::2]
        expected['D'] = self.frame['D'][::2]

        result = af.join(bf)
        assert_frame_equal(result, expected)

        result = af.join(bf, how='right')
        assert_frame_equal(result, expected[::2])

        result = bf.join(af, how='right')
        assert_frame_equal(result, expected.ix[:, result.columns])

    def test_join_index_series(self):
        df = self.frame.copy()
        s = df.pop(self.frame.columns[-1])
        joined = df.join(s)

        assert_frame_equal(joined, self.frame, check_names=False) # TODO should this check_names ?

        s.name = None
        assertRaisesRegexp(ValueError, 'must have a name', df.join, s)

    def test_join_overlap(self):
        df1 = self.frame.ix[:, ['A', 'B', 'C']]
        df2 = self.frame.ix[:, ['B', 'C', 'D']]

        joined = df1.join(df2, lsuffix='_df1', rsuffix='_df2')
        df1_suf = df1.ix[:, ['B', 'C']].add_suffix('_df1')
        df2_suf = df2.ix[:, ['B', 'C']].add_suffix('_df2')

        no_overlap = self.frame.ix[:, ['A', 'D']]
        expected = df1_suf.join(df2_suf).join(no_overlap)

        # column order not necessarily sorted
        assert_frame_equal(joined, expected.ix[:, joined.columns])

    def test_add_prefix_suffix(self):
        with_prefix = self.frame.add_prefix('foo#')
        expected = ['foo#%s' % c for c in self.frame.columns]
        self.assert_numpy_array_equal(with_prefix.columns, expected)

        with_suffix = self.frame.add_suffix('#foo')
        expected = ['%s#foo' % c for c in self.frame.columns]
        self.assert_numpy_array_equal(with_suffix.columns, expected)


class TestDataFrame(tm.TestCase, SafeForSparse, TestData):

    klass = DataFrame

    _multiprocess_can_split_ = True

    def test_get_axis(self):
        f = self.frame
        self.assertEqual(f._get_axis_number(0), 0)
        self.assertEqual(f._get_axis_number(1), 1)
        self.assertEqual(f._get_axis_number('index'), 0)
        self.assertEqual(f._get_axis_number('rows'), 0)
        self.assertEqual(f._get_axis_number('columns'), 1)

        self.assertEqual(f._get_axis_name(0), 'index')
        self.assertEqual(f._get_axis_name(1), 'columns')
        self.assertEqual(f._get_axis_name('index'), 'index')
        self.assertEqual(f._get_axis_name('rows'), 'index')
        self.assertEqual(f._get_axis_name('columns'), 'columns')

        self.assertIs(f._get_axis(0), f.index)
        self.assertIs(f._get_axis(1), f.columns)

        assertRaisesRegexp(ValueError, 'No axis named', f._get_axis_number, 2)
        assertRaisesRegexp(ValueError, 'No axis.*foo', f._get_axis_name, 'foo')
        assertRaisesRegexp(ValueError, 'No axis.*None', f._get_axis_name, None)
        assertRaisesRegexp(ValueError, 'No axis named', f._get_axis_number,
                           None)

    def test_set_index(self):
        idx = Index(np.arange(len(self.mixed_frame)))

        # cache it
        _ = self.mixed_frame['foo']  # noqa
        self.mixed_frame.index = idx
        self.assertIs(self.mixed_frame['foo'].index, idx)
        with assertRaisesRegexp(ValueError, 'Length mismatch'):
            self.mixed_frame.index = idx[::2]

    def test_set_index_cast(self):

        # issue casting an index then set_index
        df = DataFrame({'A' : [1.1,2.2,3.3], 'B' : [5.0,6.1,7.2]},
                       index = [2010,2011,2012])
        expected = df.ix[2010]
        new_index = df.index.astype(np.int32)
        df.index = new_index
        result = df.ix[2010]
        assert_series_equal(result,expected)

    def test_set_index2(self):
        df = DataFrame({'A': ['foo', 'foo', 'foo', 'bar', 'bar'],
                        'B': ['one', 'two', 'three', 'one', 'two'],
                        'C': ['a', 'b', 'c', 'd', 'e'],
                        'D': np.random.randn(5),
                        'E': np.random.randn(5)})

        # new object, single-column
        result = df.set_index('C')
        result_nodrop = df.set_index('C', drop=False)

        index = Index(df['C'], name='C')

        expected = df.ix[:, ['A', 'B', 'D', 'E']]
        expected.index = index

        expected_nodrop = df.copy()
        expected_nodrop.index = index

        assert_frame_equal(result, expected)
        assert_frame_equal(result_nodrop, expected_nodrop)
        self.assertEqual(result.index.name, index.name)

        # inplace, single
        df2 = df.copy()

        df2.set_index('C', inplace=True)

        assert_frame_equal(df2, expected)

        df3 = df.copy()
        df3.set_index('C', drop=False, inplace=True)

        assert_frame_equal(df3, expected_nodrop)

        # create new object, multi-column
        result = df.set_index(['A', 'B'])
        result_nodrop = df.set_index(['A', 'B'], drop=False)

        index = MultiIndex.from_arrays([df['A'], df['B']], names=['A', 'B'])

        expected = df.ix[:, ['C', 'D', 'E']]
        expected.index = index

        expected_nodrop = df.copy()
        expected_nodrop.index = index

        assert_frame_equal(result, expected)
        assert_frame_equal(result_nodrop, expected_nodrop)
        self.assertEqual(result.index.names, index.names)

        # inplace
        df2 = df.copy()
        df2.set_index(['A', 'B'], inplace=True)
        assert_frame_equal(df2, expected)

        df3 = df.copy()
        df3.set_index(['A', 'B'], drop=False, inplace=True)
        assert_frame_equal(df3, expected_nodrop)

        # corner case
        with assertRaisesRegexp(ValueError, 'Index has duplicate keys'):
            df.set_index('A', verify_integrity=True)

        # append
        result = df.set_index(['A', 'B'], append=True)
        xp = df.reset_index().set_index(['index', 'A', 'B'])
        xp.index.names = [None, 'A', 'B']
        assert_frame_equal(result, xp)

        # append to existing multiindex
        rdf = df.set_index(['A'], append=True)
        rdf = rdf.set_index(['B', 'C'], append=True)
        expected = df.set_index(['A', 'B', 'C'], append=True)
        assert_frame_equal(rdf, expected)

        # Series
        result = df.set_index(df.C)
        self.assertEqual(result.index.name, 'C')

    def test_set_index_nonuniq(self):
        df = DataFrame({'A': ['foo', 'foo', 'foo', 'bar', 'bar'],
                        'B': ['one', 'two', 'three', 'one', 'two'],
                        'C': ['a', 'b', 'c', 'd', 'e'],
                        'D': np.random.randn(5),
                        'E': np.random.randn(5)})
        with assertRaisesRegexp(ValueError, 'Index has duplicate keys'):
            df.set_index('A', verify_integrity=True, inplace=True)
        self.assertIn('A', df)

    def test_set_index_bug(self):
        # GH1590
        df = DataFrame({'val': [0, 1, 2], 'key': ['a', 'b', 'c']})
        df2 = df.select(lambda indx: indx >= 1)
        rs = df2.set_index('key')
        xp = DataFrame({'val': [1, 2]},
                       Index(['b', 'c'], name='key'))
        assert_frame_equal(rs, xp)

    def test_set_index_pass_arrays(self):
        df = DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                              'foo', 'bar', 'foo', 'foo'],
                        'B': ['one', 'one', 'two', 'three',
                              'two', 'two', 'one', 'three'],
                        'C': np.random.randn(8),
                        'D': np.random.randn(8)})

        # multiple columns
        result = df.set_index(['A', df['B'].values], drop=False)
        expected = df.set_index(['A', 'B'], drop=False)

        # TODO should set_index check_names ?
        assert_frame_equal(result, expected, check_names=False)

    def test_construction_with_categorical_index(self):

        ci = tm.makeCategoricalIndex(10)

        # with Categorical
        df = DataFrame({'A' : np.random.randn(10),
                        'B' : ci.values})
        idf = df.set_index('B')
        str(idf)
        tm.assert_index_equal(idf.index, ci, check_names=False)
        self.assertEqual(idf.index.name, 'B')

        # from a CategoricalIndex
        df = DataFrame({'A' : np.random.randn(10),
                        'B' : ci})
        idf = df.set_index('B')
        str(idf)
        tm.assert_index_equal(idf.index, ci, check_names=False)
        self.assertEqual(idf.index.name, 'B')

        idf = df.set_index('B').reset_index().set_index('B')
        str(idf)
        tm.assert_index_equal(idf.index, ci, check_names=False)
        self.assertEqual(idf.index.name, 'B')

        new_df = idf.reset_index()
        new_df.index = df.B
        tm.assert_index_equal(new_df.index, ci, check_names=False)
        self.assertEqual(idf.index.name, 'B')

    def test_set_index_cast_datetimeindex(self):
        df = DataFrame({'A': [datetime(2000, 1, 1) + timedelta(i)
                              for i in range(1000)],
                        'B': np.random.randn(1000)})

        idf = df.set_index('A')
        tm.assertIsInstance(idf.index, DatetimeIndex)

        # don't cast a DatetimeIndex WITH a tz, leave as object
        # GH 6032
        i = pd.DatetimeIndex(pd.tseries.tools.to_datetime(['2013-1-1 13:00','2013-1-2 14:00'], errors="raise")).tz_localize('US/Pacific')
        df = DataFrame(np.random.randn(2,1),columns=['A'])

        expected = Series(np.array([pd.Timestamp('2013-01-01 13:00:00-0800', tz='US/Pacific'),
                                    pd.Timestamp('2013-01-02 14:00:00-0800', tz='US/Pacific')], dtype="object"))

        # convert index to series
        result = Series(i)
        assert_series_equal(result, expected)

        # assignt to frame
        df['B'] = i
        result = df['B']
        assert_series_equal(result, expected, check_names=False)
        self.assertEqual(result.name, 'B')

        # keep the timezone
        result = i.to_series(keep_tz=True)
        assert_series_equal(result.reset_index(drop=True), expected)

        # convert to utc
        df['C'] = i.to_series().reset_index(drop=True)
        result = df['C']
        comp = DatetimeIndex(expected.values).copy()
        comp.tz = None
        self.assert_numpy_array_equal(result.values, comp.values)

        # list of datetimes with a tz
        df['D'] = i.to_pydatetime()
        result = df['D']
        assert_series_equal(result, expected, check_names=False)
        self.assertEqual(result.name, 'D')

        # GH 6785
        # set the index manually
        import pytz
        df = DataFrame([{'ts':datetime(2014, 4, 1, tzinfo=pytz.utc), 'foo':1}])
        expected = df.set_index('ts')
        df.index = df['ts']
        df.pop('ts')
        assert_frame_equal(df, expected)

        # GH 3950
        # reset_index with single level
        for tz in ['UTC', 'Asia/Tokyo', 'US/Eastern']:
            idx = pd.date_range('1/1/2011', periods=5, freq='D', tz=tz, name='idx')
            df = pd.DataFrame({'a': range(5), 'b': ['A', 'B', 'C', 'D', 'E']}, index=idx)

            expected = pd.DataFrame({'idx': [datetime(2011, 1, 1), datetime(2011, 1, 2),
                                             datetime(2011, 1, 3), datetime(2011, 1, 4),
                                             datetime(2011, 1, 5)],
                                     'a': range(5), 'b': ['A', 'B', 'C', 'D', 'E']},
                                     columns=['idx', 'a', 'b'])
            expected['idx'] = expected['idx'].apply(lambda d: pd.Timestamp(d, tz=tz))
            assert_frame_equal(df.reset_index(), expected)

    def test_set_index_multiindexcolumns(self):
        columns = MultiIndex.from_tuples([('foo', 1), ('foo', 2), ('bar', 1)])
        df = DataFrame(np.random.randn(3, 3), columns=columns)
        rs = df.set_index(df.columns[0])
        xp = df.ix[:, 1:]
        xp.index = df.ix[:, 0].values
        xp.index.names = [df.columns[0]]
        assert_frame_equal(rs, xp)

    def test_set_index_empty_column(self):
        # #1971
        df = DataFrame([
            dict(a=1, p=0),
            dict(a=2, m=10),
            dict(a=3, m=11, p=20),
            dict(a=4, m=12, p=21)
        ], columns=('a', 'm', 'p', 'x'))

        # it works!
        result = df.set_index(['a', 'x'])
        repr(result)

    def test_set_columns(self):
        cols = Index(np.arange(len(self.mixed_frame.columns)))
        self.mixed_frame.columns = cols
        with assertRaisesRegexp(ValueError, 'Length mismatch'):
            self.mixed_frame.columns = cols[::2]

    def test_keys(self):
        getkeys = self.frame.keys
        self.assertIs(getkeys(), self.frame.columns)

    def test_column_contains_typeerror(self):
        try:
            self.frame.columns in self.frame
        except TypeError:
            pass

    def test_insert_error_msmgs(self):

        # GH 7432
        df = DataFrame({'foo':['a', 'b', 'c'], 'bar':[1,2,3], 'baz':['d','e','f']}).set_index('foo')
        s = DataFrame({'foo':['a', 'b', 'c', 'a'], 'fiz':['g','h','i','j']}).set_index('foo')
        msg = 'cannot reindex from a duplicate axis'
        with assertRaisesRegexp(ValueError, msg):
            df['newcol'] = s

        # GH 4107, more descriptive error message
        df = DataFrame(np.random.randint(0,2,(4,4)),
                       columns=['a', 'b', 'c', 'd'])

        msg = 'incompatible index of inserted column with frame index'
        with assertRaisesRegexp(TypeError, msg):
            df['gr'] = df.groupby(['b', 'c']).count()

    def test_frame_subclassing_and_slicing(self):
        # Subclass frame and ensure it returns the right class on slicing it
        # In reference to PR 9632

        class CustomSeries(Series):
            @property
            def _constructor(self):
                return CustomSeries

            def custom_series_function(self):
                return 'OK'

        class CustomDataFrame(DataFrame):
            "Subclasses pandas DF, fills DF with simulation results, adds some custom plotting functions."

            def __init__(self, *args, **kw):
                super(CustomDataFrame, self).__init__(*args, **kw)

            @property
            def _constructor(self):
                return CustomDataFrame

            _constructor_sliced = CustomSeries

            def custom_frame_function(self):
                return 'OK'

        data = {'col1': range(10),
                'col2': range(10)}
        cdf = CustomDataFrame(data)

        # Did we get back our own DF class?
        self.assertTrue(isinstance(cdf, CustomDataFrame))

        # Do we get back our own Series class after selecting a column?
        cdf_series = cdf.col1
        self.assertTrue(isinstance(cdf_series, CustomSeries))
        self.assertEqual(cdf_series.custom_series_function(), 'OK')

        # Do we get back our own DF class after slicing row-wise?
        cdf_rows = cdf[1:5]
        self.assertTrue(isinstance(cdf_rows, CustomDataFrame))
        self.assertEqual(cdf_rows.custom_frame_function(), 'OK')

        # Make sure sliced part of multi-index frame is custom class
        mcol = pd.MultiIndex.from_tuples([('A', 'A'), ('A', 'B')])
        cdf_multi = CustomDataFrame([[0, 1], [2, 3]], columns=mcol)
        self.assertTrue(isinstance(cdf_multi['A'], CustomDataFrame))

        mcol = pd.MultiIndex.from_tuples([('A', ''), ('B', '')])
        cdf_multi2 = CustomDataFrame([[0, 1], [2, 3]], columns=mcol)
        self.assertTrue(isinstance(cdf_multi2['A'], CustomSeries))

    def test_column_dups_operations(self):

        def check(result, expected=None):
            if expected is not None:
                assert_frame_equal(result,expected)
            result.dtypes
            str(result)

        # assignment
        # GH 3687
        arr = np.random.randn(3, 2)
        idx = lrange(2)
        df = DataFrame(arr, columns=['A', 'A'])
        df.columns = idx
        expected = DataFrame(arr,columns=idx)
        check(df,expected)

        idx = date_range('20130101',periods=4,freq='Q-NOV')
        df = DataFrame([[1,1,1,5],[1,1,2,5],[2,1,3,5]],columns=['a','a','a','a'])
        df.columns = idx
        expected = DataFrame([[1,1,1,5],[1,1,2,5],[2,1,3,5]],columns=idx)
        check(df,expected)

        # insert
        df = DataFrame([[1,1,1,5],[1,1,2,5],[2,1,3,5]],columns=['foo','bar','foo','hello'])
        df['string'] = 'bah'
        expected = DataFrame([[1,1,1,5,'bah'],[1,1,2,5,'bah'],[2,1,3,5,'bah']],columns=['foo','bar','foo','hello','string'])
        check(df,expected)
        with assertRaisesRegexp(ValueError, 'Length of value'):
            df.insert(0, 'AnotherColumn', range(len(df.index) - 1))

        # insert same dtype
        df['foo2'] = 3
        expected = DataFrame([[1,1,1,5,'bah',3],[1,1,2,5,'bah',3],[2,1,3,5,'bah',3]],columns=['foo','bar','foo','hello','string','foo2'])
        check(df,expected)

        # set (non-dup)
        df['foo2'] = 4
        expected = DataFrame([[1,1,1,5,'bah',4],[1,1,2,5,'bah',4],[2,1,3,5,'bah',4]],columns=['foo','bar','foo','hello','string','foo2'])
        check(df,expected)
        df['foo2'] = 3

        # delete (non dup)
        del df['bar']
        expected = DataFrame([[1,1,5,'bah',3],[1,2,5,'bah',3],[2,3,5,'bah',3]],columns=['foo','foo','hello','string','foo2'])
        check(df,expected)

        # try to delete again (its not consolidated)
        del df['hello']
        expected = DataFrame([[1,1,'bah',3],[1,2,'bah',3],[2,3,'bah',3]],columns=['foo','foo','string','foo2'])
        check(df,expected)

        # consolidate
        df = df.consolidate()
        expected = DataFrame([[1,1,'bah',3],[1,2,'bah',3],[2,3,'bah',3]],columns=['foo','foo','string','foo2'])
        check(df,expected)

        # insert
        df.insert(2,'new_col',5.)
        expected = DataFrame([[1,1,5.,'bah',3],[1,2,5.,'bah',3],[2,3,5.,'bah',3]],columns=['foo','foo','new_col','string','foo2'])
        check(df,expected)

        # insert a dup
        assertRaisesRegexp(ValueError, 'cannot insert', df.insert, 2, 'new_col', 4.)
        df.insert(2,'new_col',4.,allow_duplicates=True)
        expected = DataFrame([[1,1,4.,5.,'bah',3],[1,2,4.,5.,'bah',3],[2,3,4.,5.,'bah',3]],columns=['foo','foo','new_col','new_col','string','foo2'])
        check(df,expected)

        # delete (dup)
        del df['foo']
        expected = DataFrame([[4.,5.,'bah',3],[4.,5.,'bah',3],[4.,5.,'bah',3]],columns=['new_col','new_col','string','foo2'])
        assert_frame_equal(df,expected)

        # dup across dtypes
        df = DataFrame([[1,1,1.,5],[1,1,2.,5],[2,1,3.,5]],columns=['foo','bar','foo','hello'])
        check(df)

        df['foo2'] = 7.
        expected = DataFrame([[1,1,1.,5,7.],[1,1,2.,5,7.],[2,1,3.,5,7.]],columns=['foo','bar','foo','hello','foo2'])
        check(df,expected)

        result = df['foo']
        expected = DataFrame([[1,1.],[1,2.],[2,3.]],columns=['foo','foo'])
        check(result,expected)

        # multiple replacements
        df['foo'] = 'string'
        expected = DataFrame([['string',1,'string',5,7.],['string',1,'string',5,7.],['string',1,'string',5,7.]],columns=['foo','bar','foo','hello','foo2'])
        check(df,expected)

        del df['foo']
        expected = DataFrame([[1,5,7.],[1,5,7.],[1,5,7.]],columns=['bar','hello','foo2'])
        check(df,expected)

        # values
        df = DataFrame([[1,2.5],[3,4.5]], index=[1,2], columns=['x','x'])
        result = df.values
        expected = np.array([[1,2.5],[3,4.5]])
        self.assertTrue((result == expected).all().all())

        # rename, GH 4403
        df4 = DataFrame({'TClose': [22.02],
                         'RT': [0.0454],
                         'TExg': [0.0422]},
                        index=MultiIndex.from_tuples([(600809, 20130331)], names=['STK_ID', 'RPT_Date']))

        df5 = DataFrame({'STK_ID': [600809] * 3,
                         'RPT_Date': [20120930,20121231,20130331],
                         'STK_Name': [u('饡驦'), u('饡驦'), u('饡驦')],
                         'TClose': [38.05, 41.66, 30.01]},
                        index=MultiIndex.from_tuples([(600809, 20120930), (600809, 20121231),(600809,20130331)], names=['STK_ID', 'RPT_Date']))

        k = pd.merge(df4,df5,how='inner',left_index=True,right_index=True)
        result = k.rename(columns={'TClose_x':'TClose', 'TClose_y':'QT_Close'})
        str(result)
        result.dtypes

        expected = DataFrame([[0.0454, 22.02, 0.0422, 20130331, 600809, u('饡驦'), 30.01 ]],
                             columns=['RT','TClose','TExg','RPT_Date','STK_ID','STK_Name','QT_Close']).set_index(['STK_ID','RPT_Date'],drop=False)
        assert_frame_equal(result,expected)

        # reindex is invalid!
        df = DataFrame([[1,5,7.],[1,5,7.],[1,5,7.]],columns=['bar','a','a'])
        self.assertRaises(ValueError, df.reindex, columns=['bar'])
        self.assertRaises(ValueError, df.reindex, columns=['bar','foo'])

        # drop
        df = DataFrame([[1,5,7.],[1,5,7.],[1,5,7.]],columns=['bar','a','a'])
        result = df.drop(['a'],axis=1)
        expected = DataFrame([[1],[1],[1]],columns=['bar'])
        check(result,expected)
        result = df.drop('a',axis=1)
        check(result,expected)

        # describe
        df = DataFrame([[1,1,1],[2,2,2],[3,3,3]],columns=['bar','a','a'],dtype='float64')
        result = df.describe()
        s = df.iloc[:,0].describe()
        expected = pd.concat([ s, s, s],keys=df.columns,axis=1)
        check(result,expected)

        # check column dups with index equal and not equal to df's index
        df = DataFrame(np.random.randn(5, 3), index=['a', 'b', 'c', 'd', 'e'],
                       columns=['A', 'B', 'A'])
        for index in [df.index, pd.Index(list('edcba'))]:
            this_df = df.copy()
            expected_ser = pd.Series(index.values, index=this_df.index)
            expected_df = DataFrame.from_items([('A', expected_ser),
                                                ('B', this_df['B']),
                                                ('A', expected_ser)])
            this_df['A'] = index
            check(this_df, expected_df)

        # operations
        for op in ['__add__','__mul__','__sub__','__truediv__']:
            df = DataFrame(dict(A = np.arange(10), B = np.random.rand(10)))
            expected = getattr(df,op)(df)
            expected.columns = ['A','A']
            df.columns = ['A','A']
            result = getattr(df,op)(df)
            check(result,expected)

        # multiple assignments that change dtypes
        # the location indexer is a slice
        # GH 6120
        df = DataFrame(np.random.randn(5,2), columns=['that', 'that'])
        expected = DataFrame(1.0, index=range(5), columns=['that', 'that'])

        df['that'] = 1.0
        check(df, expected)

        df = DataFrame(np.random.rand(5,2), columns=['that', 'that'])
        expected = DataFrame(1, index=range(5), columns=['that', 'that'])

        df['that'] = 1
        check(df, expected)

    def test_column_dups2(self):

        # drop buggy GH 6240
        df = DataFrame({'A' : np.random.randn(5),
                        'B' : np.random.randn(5),
                        'C' : np.random.randn(5),
                        'D' : ['a','b','c','d','e'] })

        expected = df.take([0,1,1], axis=1)
        df2 = df.take([2,0,1,2,1], axis=1)
        result = df2.drop('C',axis=1)
        assert_frame_equal(result, expected)

        # dropna
        df = DataFrame({'A' : np.random.randn(5),
                        'B' : np.random.randn(5),
                        'C' : np.random.randn(5),
                        'D' : ['a','b','c','d','e'] })
        df.iloc[2,[0,1,2]] = np.nan
        df.iloc[0,0] = np.nan
        df.iloc[1,1] = np.nan
        df.iloc[:,3] = np.nan
        expected = df.dropna(subset=['A','B','C'],how='all')
        expected.columns = ['A','A','B','C']

        df.columns = ['A','A','B','C']

        result = df.dropna(subset=['A','C'],how='all')
        assert_frame_equal(result, expected)

    def test_column_dups_indexing(self):
        def check(result, expected=None):
            if expected is not None:
                assert_frame_equal(result,expected)
            result.dtypes
            str(result)

        # boolean indexing
        # GH 4879
        dups = ['A', 'A', 'C', 'D']
        df = DataFrame(np.arange(12).reshape(3,4), columns=['A', 'B', 'C', 'D'],dtype='float64')
        expected = df[df.C > 6]
        expected.columns = dups
        df = DataFrame(np.arange(12).reshape(3,4), columns=dups,dtype='float64')
        result = df[df.C > 6]
        check(result,expected)

        # where
        df = DataFrame(np.arange(12).reshape(3,4), columns=['A', 'B', 'C', 'D'],dtype='float64')
        expected = df[df > 6]
        expected.columns = dups
        df = DataFrame(np.arange(12).reshape(3,4), columns=dups,dtype='float64')
        result = df[df > 6]
        check(result,expected)

        # boolean with the duplicate raises
        df = DataFrame(np.arange(12).reshape(3,4), columns=dups,dtype='float64')
        self.assertRaises(ValueError, lambda : df[df.A > 6])

        # dup aligining operations should work
        # GH 5185
        df1 = DataFrame([1, 2, 3, 4, 5], index=[1, 2, 1, 2, 3])
        df2 = DataFrame([1, 2, 3], index=[1, 2, 3])
        expected = DataFrame([0,2,0,2,2],index=[1,1,2,2,3])
        result = df1.sub(df2)
        assert_frame_equal(result,expected)

        # equality
        df1 = DataFrame([[1,2],[2,np.nan],[3,4],[4,4]],columns=['A','B'])
        df2 = DataFrame([[0,1],[2,4],[2,np.nan],[4,5]],columns=['A','A'])

        # not-comparing like-labelled
        self.assertRaises(ValueError, lambda : df1 == df2)

        df1r = df1.reindex_like(df2)
        result = df1r == df2
        expected = DataFrame([[False,True],[True,False],[False,False],[True,False]],columns=['A','A'])
        assert_frame_equal(result,expected)

        # mixed column selection
        # GH 5639
        dfbool = DataFrame({'one' : Series([True, True, False], index=['a', 'b', 'c']),
                            'two' : Series([False, False, True, False], index=['a', 'b', 'c', 'd']),
                            'three': Series([False, True, True, True], index=['a', 'b', 'c', 'd'])})
        expected = pd.concat([dfbool['one'],dfbool['three'],dfbool['one']],axis=1)
        result = dfbool[['one', 'three', 'one']]
        check(result,expected)

        # multi-axis dups
        # GH 6121
        df = DataFrame(np.arange(25.).reshape(5,5),
                       index=['a', 'b', 'c', 'd', 'e'],
                       columns=['A', 'B', 'C', 'D', 'E'])
        z = df[['A', 'C', 'A']].copy()
        expected = z.ix[['a', 'c', 'a']]

        df = DataFrame(np.arange(25.).reshape(5,5),
                       index=['a', 'b', 'c', 'd', 'e'],
                       columns=['A', 'B', 'C', 'D', 'E'])
        z = df[['A', 'C', 'A']]
        result = z.ix[['a', 'c', 'a']]
        check(result,expected)

    def test_column_dups_indexing2(self):

        # GH 8363
        # datetime ops with a non-unique index
        df = DataFrame({'A' : np.arange(5,dtype='int64'),
                        'B' : np.arange(1,6,dtype='int64')},
                       index=[2,2,3,3,4])
        result = df.B-df.A
        expected = Series(1,index=[2,2,3,3,4])
        assert_series_equal(result,expected)

        df = DataFrame({'A' : date_range('20130101',periods=5), 'B' : date_range('20130101 09:00:00', periods=5)},index=[2,2,3,3,4])
        result = df.B-df.A
        expected = Series(Timedelta('9 hours'),index=[2,2,3,3,4])
        assert_series_equal(result,expected)

    def test_insert_benchmark(self):
        # from the vb_suite/frame_methods/frame_insert_columns
        N = 10
        K = 5
        df = DataFrame(index=lrange(N))
        new_col = np.random.randn(N)
        for i in range(K):
            df[i] = new_col
        expected = DataFrame(np.repeat(new_col,K).reshape(N,K),index=lrange(N))
        assert_frame_equal(df,expected)

    def test_not_hashable(self):
        df = pd.DataFrame([1])
        self.assertRaises(TypeError, hash, df)
        self.assertRaises(TypeError, hash, self.empty)

    def test_timedeltas(self):

        df = DataFrame(dict(A = Series(date_range('2012-1-1', periods=3, freq='D')),
                            B = Series([ timedelta(days=i) for i in range(3) ])))
        result = df.get_dtype_counts().sort_values()
        expected = Series({'datetime64[ns]': 1, 'timedelta64[ns]' : 1 }).sort_values()
        assert_series_equal(result, expected)

        df['C'] = df['A'] + df['B']
        expected = Series({'datetime64[ns]': 2, 'timedelta64[ns]' : 1 }).sort_values()
        result = df.get_dtype_counts().sort_values()
        assert_series_equal(result, expected)

        # mixed int types
        df['D'] = 1
        expected = Series({'datetime64[ns]': 2, 'timedelta64[ns]' : 1, 'int64' : 1 }).sort_values()
        result = df.get_dtype_counts().sort_values()
        assert_series_equal(result, expected)

    def test_operators_timedelta64(self):

        from datetime import timedelta
        df = DataFrame(dict(A = date_range('2012-1-1', periods=3, freq='D'),
                            B = date_range('2012-1-2', periods=3, freq='D'),
                            C = Timestamp('20120101')-timedelta(minutes=5,seconds=5)))

        diffs = DataFrame(dict(A = df['A']-df['C'],
                               B = df['A']-df['B']))


        # min
        result = diffs.min()
        self.assertEqual(result[0], diffs.ix[0,'A'])
        self.assertEqual(result[1], diffs.ix[0,'B'])

        result = diffs.min(axis=1)
        self.assertTrue((result == diffs.ix[0,'B']).all() == True)

        # max
        result = diffs.max()
        self.assertEqual(result[0], diffs.ix[2,'A'])
        self.assertEqual(result[1], diffs.ix[2,'B'])

        result = diffs.max(axis=1)
        self.assertTrue((result == diffs['A']).all() == True)

        # abs
        result = diffs.abs()
        result2 = abs(diffs)
        expected = DataFrame(dict(A = df['A']-df['C'],
                                  B = df['B']-df['A']))
        assert_frame_equal(result,expected)
        assert_frame_equal(result2, expected)

        # mixed frame
        mixed = diffs.copy()
        mixed['C'] = 'foo'
        mixed['D'] = 1
        mixed['E'] = 1.
        mixed['F'] = Timestamp('20130101')

        # results in an object array
        from pandas.tseries.timedeltas import _coerce_scalar_to_timedelta_type
        result = mixed.min()
        expected = Series([_coerce_scalar_to_timedelta_type(timedelta(seconds=5*60+5)),
                           _coerce_scalar_to_timedelta_type(timedelta(days=-1)),
                           'foo',
                           1,
                           1.0,
                           Timestamp('20130101')],
                          index=mixed.columns)
        assert_series_equal(result,expected)

        # excludes numeric
        result = mixed.min(axis=1)
        expected = Series([1, 1, 1.],index=[0, 1, 2])
        assert_series_equal(result,expected)

        # works when only those columns are selected
        result = mixed[['A','B']].min(1)
        expected = Series([ timedelta(days=-1) ] * 3)
        assert_series_equal(result,expected)

        result = mixed[['A','B']].min()
        expected = Series([ timedelta(seconds=5*60+5), timedelta(days=-1) ],index=['A','B'])
        assert_series_equal(result,expected)

        # GH 3106
        df = DataFrame({'time' : date_range('20130102',periods=5),
                        'time2' : date_range('20130105',periods=5) })
        df['off1'] = df['time2']-df['time']
        self.assertEqual(df['off1'].dtype, 'timedelta64[ns]')

        df['off2'] = df['time']-df['time2']
        df._consolidate_inplace()
        self.assertTrue(df['off1'].dtype == 'timedelta64[ns]')
        self.assertTrue(df['off2'].dtype == 'timedelta64[ns]')

    def test_datetimelike_setitem_with_inference(self):
        # GH 7592
        # assignment of timedeltas with NaT

        one_hour = timedelta(hours=1)
        df = DataFrame(index=date_range('20130101',periods=4))
        df['A'] = np.array([1*one_hour]*4, dtype='m8[ns]')
        df.loc[:,'B'] = np.array([2*one_hour]*4, dtype='m8[ns]')
        df.loc[:3,'C'] = np.array([3*one_hour]*3, dtype='m8[ns]')
        df.ix[:,'D'] = np.array([4*one_hour]*4, dtype='m8[ns]')
        df.ix[:3,'E'] = np.array([5*one_hour]*3, dtype='m8[ns]')
        df['F'] = np.timedelta64('NaT')
        df.ix[:-1,'F'] = np.array([6*one_hour]*3, dtype='m8[ns]')
        df.ix[-3:,'G'] = date_range('20130101',periods=3)
        df['H'] = np.datetime64('NaT')
        result = df.dtypes
        expected = Series([np.dtype('timedelta64[ns]')]*6+[np.dtype('datetime64[ns]')]*2,index=list('ABCDEFGH'))
        assert_series_equal(result,expected)

    def test_setitem_datetime_coercion(self):
        # GH 1048
        df = pd.DataFrame({'c': [pd.Timestamp('2010-10-01')]*3})
        df.loc[0:1, 'c'] = np.datetime64('2008-08-08')
        self.assertEqual(pd.Timestamp('2008-08-08'), df.loc[0, 'c'])
        self.assertEqual(pd.Timestamp('2008-08-08'), df.loc[1, 'c'])
        df.loc[2, 'c'] = date(2005, 5, 5)
        self.assertEqual(pd.Timestamp('2005-05-05'), df.loc[2, 'c'])

    def test_new_empty_index(self):
        df1 = DataFrame(randn(0, 3))
        df2 = DataFrame(randn(0, 3))
        df1.index.name = 'foo'
        self.assertIsNone(df2.index.name)

    def test_astype(self):
        casted = self.frame.astype(int)
        expected = DataFrame(self.frame.values.astype(int),
                             index=self.frame.index,
                             columns=self.frame.columns)
        assert_frame_equal(casted, expected)

        casted = self.frame.astype(np.int32)
        expected = DataFrame(self.frame.values.astype(np.int32),
                             index=self.frame.index,
                             columns=self.frame.columns)
        assert_frame_equal(casted, expected)

        self.frame['foo'] = '5'
        casted = self.frame.astype(int)
        expected = DataFrame(self.frame.values.astype(int),
                             index=self.frame.index,
                             columns=self.frame.columns)
        assert_frame_equal(casted, expected)

        # mixed casting
        def _check_cast(df, v):
            self.assertEqual(list(set([ s.dtype.name for _, s in compat.iteritems(df) ]))[0], v)

        mn = self.all_mixed._get_numeric_data().copy()
        mn['little_float'] = np.array(12345.,dtype='float16')
        mn['big_float']    = np.array(123456789101112.,dtype='float64')

        casted = mn.astype('float64')
        _check_cast(casted, 'float64')

        casted = mn.astype('int64')
        _check_cast(casted, 'int64')

        casted = self.mixed_float.reindex(columns = ['A','B']).astype('float32')
        _check_cast(casted, 'float32')

        casted = mn.reindex(columns = ['little_float']).astype('float16')
        _check_cast(casted, 'float16')

        casted = self.mixed_float.reindex(columns = ['A','B']).astype('float16')
        _check_cast(casted, 'float16')

        casted = mn.astype('float32')
        _check_cast(casted, 'float32')

        casted = mn.astype('int32')
        _check_cast(casted, 'int32')

        # to object
        casted = mn.astype('O')
        _check_cast(casted, 'object')

    def test_astype_with_exclude_string(self):
        df = self.frame.copy()
        expected = self.frame.astype(int)
        df['string'] = 'foo'
        casted = df.astype(int, raise_on_error = False)

        expected['string'] = 'foo'
        assert_frame_equal(casted, expected)

        df = self.frame.copy()
        expected = self.frame.astype(np.int32)
        df['string'] = 'foo'
        casted = df.astype(np.int32, raise_on_error = False)

        expected['string'] = 'foo'
        assert_frame_equal(casted, expected)

    def test_astype_with_view(self):

        tf = self.mixed_float.reindex(columns = ['A','B','C'])

        casted = tf.astype(np.int64)

        casted = tf.astype(np.float32)

        # this is the only real reason to do it this way
        tf = np.round(self.frame).astype(np.int32)
        casted = tf.astype(np.float32, copy = False)

        tf = self.frame.astype(np.float64)
        casted = tf.astype(np.int64, copy = False)

    def test_astype_cast_nan_int(self):
        df = DataFrame(data={"Values": [1.0, 2.0, 3.0, np.nan]})
        self.assertRaises(ValueError, df.astype, np.int64)

    def test_astype_str(self):
        # GH9757
        a = Series(date_range('2010-01-04', periods=5))
        b = Series(date_range('3/6/2012 00:00', periods=5, tz='US/Eastern'))
        c = Series([Timedelta(x, unit='d') for x in range(5)])
        d = Series(range(5))
        e = Series([0.0, 0.2, 0.4, 0.6, 0.8])

        df = DataFrame({'a' : a, 'b' : b, 'c' : c, 'd' : d, 'e' : e})

        # datetimelike
        # Test str and unicode on python 2.x and just str on python 3.x
        for tt in set([str, compat.text_type]):
            result = df.astype(tt)

            expected = DataFrame({
                'a' : list(map(tt, map(lambda x: Timestamp(x)._date_repr, a._values))),
                'b' : list(map(tt, map(Timestamp, b._values))),
                'c' : list(map(tt, map(lambda x: Timedelta(x)._repr_base(format='all'), c._values))),
                'd' : list(map(tt, d._values)),
                'e' : list(map(tt, e._values)),
                })

            assert_frame_equal(result, expected)

        # float/nan
        # 11302
        # consistency in astype(str)
        for tt in set([str, compat.text_type]):
            result = DataFrame([np.NaN]).astype(tt)
            expected = DataFrame(['nan'])
            assert_frame_equal(result, expected)

            result = DataFrame([1.12345678901234567890]).astype(tt)
            expected = DataFrame(['1.12345678901'])
            assert_frame_equal(result, expected)

    def test_array_interface(self):
        result = np.sqrt(self.frame)
        tm.assertIsInstance(result, type(self.frame))
        self.assertIs(result.index, self.frame.index)
        self.assertIs(result.columns, self.frame.columns)

        assert_frame_equal(result, self.frame.apply(np.sqrt))

    def test_pickle(self):
        unpickled = self.round_trip_pickle(self.mixed_frame)
        assert_frame_equal(self.mixed_frame, unpickled)

        # buglet
        self.mixed_frame._data.ndim

        # empty
        unpickled = self.round_trip_pickle(self.empty)
        repr(unpickled)

        # tz frame
        unpickled = self.round_trip_pickle(self.tzframe)
        assert_frame_equal(self.tzframe, unpickled)

    def test_to_dict(self):
        test_data = {
            'A': {'1': 1, '2': 2},
            'B': {'1': '1', '2': '2', '3': '3'},
        }
        recons_data = DataFrame(test_data).to_dict()

        for k, v in compat.iteritems(test_data):
            for k2, v2 in compat.iteritems(v):
                self.assertEqual(v2, recons_data[k][k2])

        recons_data = DataFrame(test_data).to_dict("l")

        for k, v in compat.iteritems(test_data):
            for k2, v2 in compat.iteritems(v):
                self.assertEqual(v2, recons_data[k][int(k2) - 1])

        recons_data = DataFrame(test_data).to_dict("s")

        for k, v in compat.iteritems(test_data):
            for k2, v2 in compat.iteritems(v):
                self.assertEqual(v2, recons_data[k][k2])

        recons_data = DataFrame(test_data).to_dict("sp")

        expected_split = {'columns': ['A', 'B'], 'index': ['1', '2', '3'],
                          'data': [[1.0, '1'], [2.0, '2'], [nan, '3']]}

        tm.assert_almost_equal(recons_data, expected_split)

        recons_data = DataFrame(test_data).to_dict("r")

        expected_records = [{'A': 1.0, 'B': '1'},
                            {'A': 2.0, 'B': '2'},
                            {'A': nan, 'B': '3'}]

        tm.assert_almost_equal(recons_data, expected_records)

        # GH10844
        recons_data = DataFrame(test_data).to_dict("i")

        for k, v in compat.iteritems(test_data):
            for k2, v2 in compat.iteritems(v):
                self.assertEqual(v2, recons_data[k2][k])

    def test_latex_repr(self):
        result=r"""\begin{tabular}{llll}
\toprule
{} &         0 &  1 &  2 \\
\midrule
0 &  $\alpha$ &  b &  c \\
1 &         1 &  2 &  3 \\
\bottomrule
\end{tabular}
"""
        with option_context("display.latex.escape",False):
            df=DataFrame([[r'$\alpha$','b','c'],[1,2,3]])
            self.assertEqual(result,df._repr_latex_())

    def test_to_dict_timestamp(self):

        # GH11247
        # split/records producing np.datetime64 rather than Timestamps
        # on datetime64[ns] dtypes only

        tsmp = Timestamp('20130101')
        test_data = DataFrame({'A': [tsmp, tsmp], 'B': [tsmp, tsmp]})
        test_data_mixed = DataFrame({'A': [tsmp, tsmp], 'B': [1, 2]})

        expected_records = [{'A': tsmp, 'B': tsmp},
                            {'A': tsmp, 'B': tsmp}]
        expected_records_mixed = [{'A': tsmp, 'B': 1},
                            {'A': tsmp, 'B': 2}]

        tm.assert_almost_equal(test_data.to_dict(
            orient='records'), expected_records)
        tm.assert_almost_equal(test_data_mixed.to_dict(
            orient='records'), expected_records_mixed)

        expected_series = {
            'A': Series([tsmp, tsmp]),
            'B': Series([tsmp, tsmp]),
        }
        expected_series_mixed = {
            'A': Series([tsmp, tsmp]),
            'B': Series([1, 2]),
        }

        tm.assert_almost_equal(test_data.to_dict(
            orient='series'), expected_series)
        tm.assert_almost_equal(test_data_mixed.to_dict(
            orient='series'), expected_series_mixed)

        expected_split = {
            'index': [0, 1],
            'data': [[tsmp, tsmp],
                     [tsmp, tsmp]],
            'columns': ['A', 'B']
        }
        expected_split_mixed = {
            'index': [0, 1],
            'data': [[tsmp, 1],
                     [tsmp, 2]],
            'columns': ['A', 'B']
        }

        tm.assert_almost_equal(test_data.to_dict(
            orient='split'), expected_split)
        tm.assert_almost_equal(test_data_mixed.to_dict(
            orient='split'), expected_split_mixed)

    def test_to_dict_invalid_orient(self):
        df = DataFrame({'A':[0, 1]})
        self.assertRaises(ValueError, df.to_dict, orient='xinvalid')

    def test_to_records_dt64(self):
        df = DataFrame([["one", "two", "three"],
                        ["four", "five", "six"]],
                       index=date_range("2012-01-01", "2012-01-02"))
        self.assertEqual(df.to_records()['index'][0], df.index[0])

        rs = df.to_records(convert_datetime64=False)
        self.assertEqual(rs['index'][0], df.index.values[0])

    def test_to_records_with_multindex(self):
        # GH3189
        index = [['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux'],
                 ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']]
        data = np.zeros((8, 4))
        df = DataFrame(data, index=index)
        r = df.to_records(index=True)['level_0']
        self.assertTrue('bar' in r)
        self.assertTrue('one' not in r)

    def test_to_records_with_Mapping_type(self):
        import email
        from email.parser import Parser
        import collections

        collections.Mapping.register(email.message.Message)

        headers = Parser().parsestr('From: <user@example.com>\n'
                'To: <someone_else@example.com>\n'
                'Subject: Test message\n'
                '\n'
                'Body would go here\n')

        frame = DataFrame.from_records([headers])
        all( x in frame for x in ['Type','Subject','From'])

    def test_to_records_floats(self):
        df = DataFrame(np.random.rand(10, 10))
        df.to_records()

    def test_to_recods_index_name(self):
        df = DataFrame(np.random.randn(3, 3))
        df.index.name = 'X'
        rs = df.to_records()
        self.assertIn('X', rs.dtype.fields)

        df = DataFrame(np.random.randn(3, 3))
        rs = df.to_records()
        self.assertIn('index', rs.dtype.fields)

        df.index = MultiIndex.from_tuples([('a', 'x'), ('a', 'y'), ('b', 'z')])
        df.index.names = ['A', None]
        rs = df.to_records()
        self.assertIn('level_0', rs.dtype.fields)

    def test_join_str_datetime(self):
        str_dates = ['20120209', '20120222']
        dt_dates = [datetime(2012, 2, 9), datetime(2012, 2, 22)]

        A = DataFrame(str_dates, index=lrange(2), columns=['aa'])
        C = DataFrame([[1, 2], [3, 4]], index=str_dates, columns=dt_dates)

        tst = A.join(C, on='aa')

        self.assertEqual(len(tst.columns), 3)

    def test_join_multiindex_leftright(self):
        # GH 10741
        df1 = pd.DataFrame([['a', 'x', 0.471780], ['a','y', 0.774908],
                            ['a', 'z', 0.563634], ['b', 'x', -0.353756],
                            ['b', 'y', 0.368062], ['b', 'z', -1.721840],
                            ['c', 'x', 1], ['c', 'y', 2], ['c', 'z', 3]],
                           columns=['first', 'second', 'value1']).set_index(['first', 'second'])
        df2 = pd.DataFrame([['a', 10], ['b', 20]], columns=['first', 'value2']).set_index(['first'])

        exp = pd.DataFrame([[0.471780, 10], [0.774908, 10], [0.563634, 10],
                            [-0.353756, 20], [0.368062, 20], [-1.721840, 20],
                            [1.000000, np.nan], [2.000000, np.nan], [3.000000, np.nan]],
                           index=df1.index, columns=['value1', 'value2'])

        # these must be the same results (but columns are flipped)
        assert_frame_equal(df1.join(df2, how='left'), exp)
        assert_frame_equal(df2.join(df1, how='right'),
                           exp[['value2', 'value1']])

        exp_idx = pd.MultiIndex.from_product([['a', 'b'], ['x', 'y', 'z']],
                                             names=['first', 'second'])
        exp = pd.DataFrame([[0.471780, 10], [0.774908, 10], [0.563634, 10],
                            [-0.353756, 20], [0.368062, 20], [-1.721840, 20]],
                           index=exp_idx, columns=['value1', 'value2'])

        assert_frame_equal(df1.join(df2, how='right'), exp)
        assert_frame_equal(df2.join(df1, how='left'),
                           exp[['value2', 'value1']])

    def test_get_agg_axis(self):
        cols = self.frame._get_agg_axis(0)
        self.assertIs(cols, self.frame.columns)

        idx = self.frame._get_agg_axis(1)
        self.assertIs(idx, self.frame.index)

        self.assertRaises(ValueError, self.frame._get_agg_axis, 2)

    def test_nonzero(self):
        self.assertTrue(self.empty.empty)

        self.assertFalse(self.frame.empty)
        self.assertFalse(self.mixed_frame.empty)

        # corner case
        df = DataFrame({'A': [1., 2., 3.],
                        'B': ['a', 'b', 'c']},
                       index=np.arange(3))
        del df['A']
        self.assertFalse(df.empty)

    def test_repr_empty(self):
        buf = StringIO()

        # empty
        foo = repr(self.empty)

        # empty with index
        frame = DataFrame(index=np.arange(1000))
        foo = repr(frame)

    def test_repr_mixed(self):
        buf = StringIO()

        # mixed
        foo = repr(self.mixed_frame)
        self.mixed_frame.info(verbose=False, buf=buf)

    @slow
    def test_repr_mixed_big(self):
        # big mixed
        biggie = DataFrame({'A': randn(200),
                            'B': tm.makeStringIndex(200)},
                           index=lrange(200))
        biggie.loc[:20,'A'] = nan
        biggie.loc[:20,'B'] = nan

        foo = repr(biggie)

    def test_repr(self):
        buf = StringIO()

        # small one
        foo = repr(self.frame)
        self.frame.info(verbose=False, buf=buf)

        # even smaller
        self.frame.reindex(columns=['A']).info(verbose=False, buf=buf)
        self.frame.reindex(columns=['A', 'B']).info(verbose=False, buf=buf)

        # exhausting cases in DataFrame.info

        # columns but no index
        no_index = DataFrame(columns=[0, 1, 3])
        foo = repr(no_index)

        # no columns or index
        self.empty.info(buf=buf)

        df = DataFrame(["a\n\r\tb"], columns=["a\n\r\td"], index=["a\n\r\tf"])
        self.assertFalse("\t" in repr(df))
        self.assertFalse("\r" in repr(df))
        self.assertFalse("a\n" in repr(df))

    def test_repr_dimensions(self):
        df = DataFrame([[1, 2,], [3, 4]])
        with option_context('display.show_dimensions', True):
            self.assertTrue("2 rows x 2 columns" in repr(df))

        with option_context('display.show_dimensions', False):
            self.assertFalse("2 rows x 2 columns" in repr(df))

        with option_context('display.show_dimensions', 'truncate'):
            self.assertFalse("2 rows x 2 columns" in repr(df))

    @slow
    def test_repr_big(self):
        buf = StringIO()

        # big one
        biggie = DataFrame(np.zeros((200, 4)), columns=lrange(4),
                           index=lrange(200))
        foo = repr(biggie)

    def test_repr_unsortable(self):
        # columns are not sortable
        import warnings
        warn_filters = warnings.filters
        warnings.filterwarnings('ignore',
                                category=FutureWarning,
                                module=".*format")

        unsortable = DataFrame({'foo': [1] * 50,
                                datetime.today(): [1] * 50,
                                'bar': ['bar'] * 50,
                                datetime.today(
                                ) + timedelta(1): ['bar'] * 50},
                               index=np.arange(50))
        foo = repr(unsortable)

        fmt.set_option('display.precision', 3, 'display.column_space', 10)
        repr(self.frame)

        fmt.set_option('display.max_rows', 10, 'display.max_columns', 2)
        repr(self.frame)

        fmt.set_option('display.max_rows', 1000, 'display.max_columns', 1000)
        repr(self.frame)

        self.reset_display_options()

        warnings.filters = warn_filters

    def test_repr_unicode(self):
        uval = u('\u03c3\u03c3\u03c3\u03c3')
        bval = uval.encode('utf-8')
        df = DataFrame({'A': [uval, uval]})

        result = repr(df)
        ex_top = '      A'
        self.assertEqual(result.split('\n')[0].rstrip(), ex_top)

        df = DataFrame({'A': [uval, uval]})
        result = repr(df)
        self.assertEqual(result.split('\n')[0].rstrip(), ex_top)

    def test_unicode_string_with_unicode(self):
        df = DataFrame({'A': [u("\u05d0")]})

        if compat.PY3:
            str(df)
        else:
            compat.text_type(df)

    def test_bytestring_with_unicode(self):
        df = DataFrame({'A': [u("\u05d0")]})
        if compat.PY3:
            bytes(df)
        else:
            str(df)

    def test_very_wide_info_repr(self):
        df = DataFrame(np.random.randn(10, 20),
                       columns=tm.rands_array(10, 20))
        repr(df)

    def test_repr_column_name_unicode_truncation_bug(self):
        # #1906
        df = DataFrame({'Id': [7117434],
                        'StringCol': ('Is it possible to modify drop plot code'
                                      ' so that the output graph is displayed '
                                      'in iphone simulator, Is it possible to '
                                      'modify drop plot code so that the '
                                      'output graph is \xe2\x80\xa8displayed '
                                      'in iphone simulator.Now we are adding '
                                      'the CSV file externally. I want to Call'
                                      ' the File through the code..')})

        result = repr(df)
        self.assertIn('StringCol', result)

    def test_head_tail(self):
        assert_frame_equal(self.frame.head(), self.frame[:5])
        assert_frame_equal(self.frame.tail(), self.frame[-5:])

        assert_frame_equal(self.frame.head(0), self.frame[0:0])
        assert_frame_equal(self.frame.tail(0), self.frame[0:0])

        assert_frame_equal(self.frame.head(-1), self.frame[:-1])
        assert_frame_equal(self.frame.tail(-1), self.frame[1:])
        assert_frame_equal(self.frame.head(1), self.frame[:1])
        assert_frame_equal(self.frame.tail(1), self.frame[-1:])
        # with a float index
        df = self.frame.copy()
        df.index = np.arange(len(self.frame)) + 0.1
        assert_frame_equal(df.head(), df.iloc[:5])
        assert_frame_equal(df.tail(), df.iloc[-5:])
        assert_frame_equal(df.head(0), df[0:0])
        assert_frame_equal(df.tail(0), df[0:0])
        assert_frame_equal(df.head(-1), df.iloc[:-1])
        assert_frame_equal(df.tail(-1), df.iloc[1:])
        #test empty dataframe
        empty_df = DataFrame()
        assert_frame_equal(empty_df.tail(), empty_df)
        assert_frame_equal(empty_df.head(), empty_df)

    def test_insert(self):
        df = DataFrame(np.random.randn(5, 3), index=np.arange(5),
                       columns=['c', 'b', 'a'])

        df.insert(0, 'foo', df['a'])
        self.assert_numpy_array_equal(df.columns, ['foo', 'c', 'b', 'a'])
        assert_almost_equal(df['a'], df['foo'])

        df.insert(2, 'bar', df['c'])
        self.assert_numpy_array_equal(df.columns, ['foo', 'c', 'bar', 'b', 'a'])
        assert_almost_equal(df['c'], df['bar'])

        # diff dtype

        # new item
        df['x'] = df['a'].astype('float32')
        result = Series(dict(float64 = 5, float32 = 1))
        self.assertTrue((df.get_dtype_counts() == result).all())

        # replacing current (in different block)
        df['a'] = df['a'].astype('float32')
        result = Series(dict(float64 = 4, float32 = 2))
        self.assertTrue((df.get_dtype_counts() == result).all())

        df['y'] = df['a'].astype('int32')
        result = Series(dict(float64 = 4, float32 = 2, int32 = 1))
        self.assertTrue((df.get_dtype_counts() == result).all())

        with assertRaisesRegexp(ValueError, 'already exists'):
            df.insert(1, 'a', df['b'])
        self.assertRaises(ValueError, df.insert, 1, 'c', df['b'])

        df.columns.name = 'some_name'
        # preserve columns name field
        df.insert(0, 'baz', df['c'])
        self.assertEqual(df.columns.name, 'some_name')

    def test_delitem(self):
        del self.frame['A']
        self.assertNotIn('A', self.frame)

    def test_pop(self):
        self.frame.columns.name = 'baz'

        A = self.frame.pop('A')
        self.assertNotIn('A', self.frame)

        self.frame['foo'] = 'bar'
        foo = self.frame.pop('foo')
        self.assertNotIn('foo', self.frame)
        # TODO self.assertEqual(self.frame.columns.name, 'baz')

        # 10912
        # inplace ops cause caching issue
        a = DataFrame([[1,2,3],[4,5,6]], columns=['A','B','C'], index=['X','Y'])
        b = a.pop('B')
        b += 1

        # original frame
        expected = DataFrame([[1,3],[4,6]], columns=['A','C'], index=['X','Y'])
        assert_frame_equal(a, expected)

        # result
        expected = Series([2,5],index=['X','Y'],name='B')+1
        assert_series_equal(b, expected)

    def test_pop_non_unique_cols(self):
        df = DataFrame({0: [0, 1], 1: [0, 1], 2: [4, 5]})
        df.columns = ["a", "b", "a"]

        res = df.pop("a")
        self.assertEqual(type(res), DataFrame)
        self.assertEqual(len(res), 2)
        self.assertEqual(len(df.columns), 1)
        self.assertTrue("b" in df.columns)
        self.assertFalse("a" in df.columns)
        self.assertEqual(len(df.index), 2)

    def test_iter(self):
        self.assertTrue(tm.equalContents(list(self.frame), self.frame.columns))

    def test_iterrows(self):
        for i, (k, v) in enumerate(self.frame.iterrows()):
            exp = self.frame.xs(self.frame.index[i])
            assert_series_equal(v, exp)

        for i, (k, v) in enumerate(self.mixed_frame.iterrows()):
            exp = self.mixed_frame.xs(self.mixed_frame.index[i])
            assert_series_equal(v, exp)

    def test_itertuples(self):
        for i, tup in enumerate(self.frame.itertuples()):
            s = Series(tup[1:])
            s.name = tup[0]
            expected = self.frame.ix[i, :].reset_index(drop=True)
            assert_series_equal(s, expected)

        df = DataFrame({'floats': np.random.randn(5),
                        'ints': lrange(5)}, columns=['floats', 'ints'])

        for tup in df.itertuples(index=False):
            tm.assertIsInstance(tup[1], np.integer)

        df = DataFrame(data={"a": [1, 2, 3], "b": [4, 5, 6]})
        dfaa = df[['a', 'a']]
        self.assertEqual(list(dfaa.itertuples()), [(0, 1, 1), (1, 2, 2), (2, 3, 3)])

        self.assertEqual(repr(list(df.itertuples(name=None))), '[(0, 1, 4), (1, 2, 5), (2, 3, 6)]')

        tup = next(df.itertuples(name='TestName'))

        # no support for field renaming in Python 2.6, regular tuples are returned
        if sys.version >= LooseVersion('2.7'):
            self.assertEqual(tup._fields, ('Index', 'a', 'b'))
            self.assertEqual((tup.Index, tup.a, tup.b), tup)
            self.assertEqual(type(tup).__name__, 'TestName')

        df.columns = ['def', 'return']
        tup2 = next(df.itertuples(name='TestName'))
        self.assertEqual(tup2, (0, 1, 4))

        if sys.version >= LooseVersion('2.7'):
            self.assertEqual(tup2._fields, ('Index', '_1', '_2'))

        df3 = DataFrame(dict(('f'+str(i), [i]) for i in range(1024)))
        # will raise SyntaxError if trying to create namedtuple
        tup3 = next(df3.itertuples())
        self.assertFalse(hasattr(tup3, '_fields'))
        self.assertIsInstance(tup3, tuple)

    def test_len(self):
        self.assertEqual(len(self.frame), len(self.frame.index))

    def test_constructor_lists_to_object_dtype(self):
        # from #1074
        d = DataFrame({'a': [np.nan, False]})
        self.assertEqual(d['a'].dtype, np.object_)
        self.assertFalse(d['a'][1])

    def test_constructor_with_nas(self):
        # GH 5016
        # na's in indicies

        def check(df):
            for i in range(len(df.columns)):
                df.iloc[:,i]

            # allow single nans to succeed
            indexer = np.arange(len(df.columns))[isnull(df.columns)]

            if len(indexer) == 1:
                assert_series_equal(df.iloc[:,indexer[0]],df.loc[:,np.nan])


            # multiple nans should fail
            else:

                def f():
                    df.loc[:,np.nan]
                self.assertRaises(TypeError, f)


        df = DataFrame([[1,2,3],[4,5,6]], index=[1,np.nan])
        check(df)

        df = DataFrame([[1,2,3],[4,5,6]], columns=[1.1,2.2,np.nan])
        check(df)

        df = DataFrame([[0,1,2,3],[4,5,6,7]], columns=[np.nan,1.1,2.2,np.nan])
        check(df)

        df = DataFrame([[0.0,1,2,3.0],[4,5,6,7]], columns=[np.nan,1.1,2.2,np.nan])
        check(df)

    def test_first_last_valid(self):
        N = len(self.frame.index)
        mat = randn(N)
        mat[:5] = nan
        mat[-5:] = nan

        frame = DataFrame({'foo': mat}, index=self.frame.index)
        index = frame.first_valid_index()

        self.assertEqual(index, frame.index[5])

        index = frame.last_valid_index()
        self.assertEqual(index, frame.index[-6])

    def test_equals_different_blocks(self):
        # GH 9330
        df0 = pd.DataFrame({"A": ["x","y"], "B": [1,2],
                            "C": ["w","z"]})
        df1 = df0.reset_index()[["A","B","C"]]
        # this assert verifies that the above operations have
        # induced a block rearrangement
        self.assertTrue(df0._data.blocks[0].dtype !=
                        df1._data.blocks[0].dtype)
        # do the real tests
        assert_frame_equal(df0, df1)
        self.assertTrue(df0.equals(df1))
        self.assertTrue(df1.equals(df0))

    def test_copy_blocks(self):
        # API/ENH 9607
        df = DataFrame(self.frame, copy=True)
        column = df.columns[0]

        # use the default copy=True, change a column
        blocks = df.as_blocks()
        for dtype, _df in blocks.items():
            if column in _df:
                _df.ix[:, column] = _df[column] + 1

        # make sure we did not change the original DataFrame
        self.assertFalse(_df[column].equals(df[column]))

    def test_no_copy_blocks(self):
        # API/ENH 9607
        df = DataFrame(self.frame, copy=True)
        column = df.columns[0]

        # use the copy=False, change a column
        blocks = df.as_blocks(copy=False)
        for dtype, _df in blocks.items():
            if column in _df:
                _df.ix[:, column] = _df[column] + 1

        # make sure we did change the original DataFrame
        self.assertTrue(_df[column].equals(df[column]))

    def test_info(self):
        io = StringIO()
        self.frame.info(buf=io)
        self.tsframe.info(buf=io)

        frame = DataFrame(np.random.randn(5, 3))

        import sys
        sys.stdout = StringIO()
        frame.info()
        frame.info(verbose=False)
        sys.stdout = sys.__stdout__

    def test_info_wide(self):
        from pandas import set_option, reset_option
        io = StringIO()
        df = DataFrame(np.random.randn(5, 101))
        df.info(buf=io)

        io = StringIO()
        df.info(buf=io, max_cols=101)
        rs = io.getvalue()
        self.assertTrue(len(rs.splitlines()) > 100)
        xp = rs

        set_option('display.max_info_columns', 101)
        io = StringIO()
        df.info(buf=io)
        self.assertEqual(rs, xp)
        reset_option('display.max_info_columns')

    def test_info_duplicate_columns(self):
        io = StringIO()

        # it works!
        frame = DataFrame(np.random.randn(1500, 4),
                          columns=['a', 'a', 'b', 'b'])
        frame.info(buf=io)

    def test_info_duplicate_columns_shows_correct_dtypes(self):
        # GH11761
        io = StringIO()

        frame = DataFrame([[1, 2.0]],
                          columns=['a', 'a'])
        frame.info(buf=io)
        io.seek(0)
        lines = io.readlines()
        self.assertEqual('a    1 non-null int64\n', lines[3])
        self.assertEqual('a    1 non-null float64\n', lines[4])

    def test_info_shows_column_dtypes(self):
        dtypes = ['int64', 'float64', 'datetime64[ns]', 'timedelta64[ns]',
                  'complex128', 'object', 'bool']
        data = {}
        n = 10
        for i, dtype in enumerate(dtypes):
            data[i] = np.random.randint(2, size=n).astype(dtype)
        df = DataFrame(data)
        buf = StringIO()
        df.info(buf=buf)
        res = buf.getvalue()
        for i, dtype in enumerate(dtypes):
            name = '%d    %d non-null %s' % (i, n, dtype)
            assert name in res

    def test_info_max_cols(self):
        df = DataFrame(np.random.randn(10, 5))
        for len_, verbose in [(5, None), (5, False), (10, True)]:
        # For verbose always      ^ setting  ^ summarize ^ full output
            with option_context('max_info_columns', 4):
                buf = StringIO()
                df.info(buf=buf, verbose=verbose)
                res = buf.getvalue()
                self.assertEqual(len(res.strip().split('\n')), len_)

        for len_, verbose in [(10, None), (5, False), (10, True)]:

            # max_cols no exceeded
            with option_context('max_info_columns', 5):
                buf = StringIO()
                df.info(buf=buf, verbose=verbose)
                res = buf.getvalue()
                self.assertEqual(len(res.strip().split('\n')), len_)

        for len_, max_cols in [(10, 5), (5, 4)]:
            # setting truncates
            with option_context('max_info_columns', 4):
                buf = StringIO()
                df.info(buf=buf, max_cols=max_cols)
                res = buf.getvalue()
                self.assertEqual(len(res.strip().split('\n')), len_)

            # setting wouldn't truncate
            with option_context('max_info_columns', 5):
                buf = StringIO()
                df.info(buf=buf, max_cols=max_cols)
                res = buf.getvalue()
                self.assertEqual(len(res.strip().split('\n')), len_)

    def test_info_memory_usage(self):
        # Ensure memory usage is displayed, when asserted, on the last line
        dtypes = ['int64', 'float64', 'datetime64[ns]', 'timedelta64[ns]',
                  'complex128', 'object', 'bool']
        data = {}
        n = 10
        for i, dtype in enumerate(dtypes):
            data[i] = np.random.randint(2, size=n).astype(dtype)
        df = DataFrame(data)
        buf = StringIO()
        # display memory usage case
        df.info(buf=buf, memory_usage=True)
        res = buf.getvalue().splitlines()
        self.assertTrue("memory usage: " in res[-1])
        # do not display memory usage cas
        df.info(buf=buf, memory_usage=False)
        res = buf.getvalue().splitlines()
        self.assertTrue("memory usage: " not in res[-1])

        df.info(buf=buf, memory_usage=True)
        res = buf.getvalue().splitlines()
        # memory usage is a lower bound, so print it as XYZ+ MB
        self.assertTrue(re.match(r"memory usage: [^+]+\+", res[-1]))

        df.iloc[:, :5].info(buf=buf, memory_usage=True)
        res = buf.getvalue().splitlines()
        # excluded column with object dtype, so estimate is accurate
        self.assertFalse(re.match(r"memory usage: [^+]+\+", res[-1]))

        df_with_object_index = pd.DataFrame({'a': [1]}, index=['foo'])
        df_with_object_index.info(buf=buf, memory_usage=True)
        res = buf.getvalue().splitlines()
        self.assertTrue(re.match(r"memory usage: [^+]+\+", res[-1]))

        df_with_object_index.info(buf=buf, memory_usage='deep')
        res = buf.getvalue().splitlines()
        self.assertTrue(re.match(r"memory usage: [^+]+$", res[-1]))

        self.assertTrue(df_with_object_index.memory_usage(index=True,
                                                          deep=True).sum()
                        > df_with_object_index.memory_usage(index=True).sum())

        df_object = pd.DataFrame({'a': ['a']})
        self.assertTrue(df_object.memory_usage(deep=True).sum() \
                        > df_object.memory_usage().sum())

        # Test a DataFrame with duplicate columns
        dtypes = ['int64', 'int64', 'int64', 'float64']
        data = {}
        n = 100
        for i, dtype in enumerate(dtypes):
            data[i] = np.random.randint(2, size=n).astype(dtype)
        df = DataFrame(data)
        df.columns = dtypes
        # Ensure df size is as expected
        df_size = df.memory_usage().sum()
        exp_size = (len(dtypes) + 1) * n * 8  # (cols + index) * rows * bytes
        self.assertEqual(df_size, exp_size)
        # Ensure number of cols in memory_usage is the same as df
        size_df = np.size(df.columns.values) + 1  # index=True; default
        self.assertEqual(size_df, np.size(df.memory_usage()))

        # assert deep works only on object
        self.assertEqual(df.memory_usage().sum(),
                         df.memory_usage(deep=True).sum())

        # test for validity
        DataFrame(1, index=['a'], columns=['A']
                  ).memory_usage(index=True)
        DataFrame(1, index=['a'], columns=['A']
                  ).index.nbytes
        df = DataFrame(
            data=1,
            index=pd.MultiIndex.from_product(
                [['a'], range(1000)]),
            columns=['A']
        )
        df.index.nbytes
        df.memory_usage(index=True)
        df.index.values.nbytes

        # sys.getsizeof will call the .memory_usage with
        # deep=True, and add on some GC overhead
        diff = df.memory_usage(deep=True).sum() - sys.getsizeof(df)
        self.assertTrue(abs(diff) < 100)

    def test_dtypes(self):
        self.mixed_frame['bool'] = self.mixed_frame['A'] > 0
        result = self.mixed_frame.dtypes
        expected = Series(dict((k, v.dtype)
                               for k, v in compat.iteritems(self.mixed_frame)),
                          index=result.index)
        assert_series_equal(result, expected)

        # compat, GH 8722
        with option_context('use_inf_as_null',True):
            df = DataFrame([[1]])
            result = df.dtypes
            assert_series_equal(result,Series({0:np.dtype('int64')}))

    def test_convert_objects(self):

        oops = self.mixed_frame.T.T
        converted = oops._convert(datetime=True)
        assert_frame_equal(converted, self.mixed_frame)
        self.assertEqual(converted['A'].dtype, np.float64)

        # force numeric conversion
        self.mixed_frame['H'] = '1.'
        self.mixed_frame['I'] = '1'

        # add in some items that will be nan
        l = len(self.mixed_frame)
        self.mixed_frame['J'] = '1.'
        self.mixed_frame['K'] = '1'
        self.mixed_frame.ix[0:5,['J','K']] = 'garbled'
        converted = self.mixed_frame._convert(datetime=True, numeric=True)
        self.assertEqual(converted['H'].dtype, 'float64')
        self.assertEqual(converted['I'].dtype, 'int64')
        self.assertEqual(converted['J'].dtype, 'float64')
        self.assertEqual(converted['K'].dtype, 'float64')
        self.assertEqual(len(converted['J'].dropna()), l-5)
        self.assertEqual(len(converted['K'].dropna()), l-5)

        # via astype
        converted = self.mixed_frame.copy()
        converted['H'] = converted['H'].astype('float64')
        converted['I'] = converted['I'].astype('int64')
        self.assertEqual(converted['H'].dtype, 'float64')
        self.assertEqual(converted['I'].dtype, 'int64')

        # via astype, but errors
        converted = self.mixed_frame.copy()
        with assertRaisesRegexp(ValueError, 'invalid literal'):
            converted['H'].astype('int32')

        # mixed in a single column
        df = DataFrame(dict(s = Series([1, 'na', 3 ,4])))
        result = df._convert(datetime=True, numeric=True)
        expected = DataFrame(dict(s = Series([1, np.nan, 3 ,4])))
        assert_frame_equal(result, expected)

    def test_convert_objects_no_conversion(self):
        mixed1 = DataFrame(
            {'a': [1, 2, 3], 'b': [4.0, 5, 6], 'c': ['x', 'y', 'z']})
        mixed2 = mixed1._convert(datetime=True)
        assert_frame_equal(mixed1, mixed2)

    def test_append_series_dict(self):
        df = DataFrame(np.random.randn(5, 4),
                       columns=['foo', 'bar', 'baz', 'qux'])

        series = df.ix[4]
        with  assertRaisesRegexp(ValueError, 'Indexes have overlapping values'):
            df.append(series, verify_integrity=True)
        series.name = None
        with assertRaisesRegexp(TypeError, 'Can only append a Series if '
                                'ignore_index=True'):
            df.append(series, verify_integrity=True)

        result = df.append(series[::-1], ignore_index=True)
        expected = df.append(DataFrame({0: series[::-1]}, index=df.columns).T,
                             ignore_index=True)
        assert_frame_equal(result, expected)

        # dict
        result = df.append(series.to_dict(), ignore_index=True)
        assert_frame_equal(result, expected)

        result = df.append(series[::-1][:3], ignore_index=True)
        expected = df.append(DataFrame({0: series[::-1][:3]}).T,
                             ignore_index=True)
        assert_frame_equal(result, expected.ix[:, result.columns])

        # can append when name set
        row = df.ix[4]
        row.name = 5
        result = df.append(row)
        expected = df.append(df[-1:], ignore_index=True)
        assert_frame_equal(result, expected)

    def test_append_list_of_series_dicts(self):
        df = DataFrame(np.random.randn(5, 4),
                       columns=['foo', 'bar', 'baz', 'qux'])

        dicts = [x.to_dict() for idx, x in df.iterrows()]

        result = df.append(dicts, ignore_index=True)
        expected = df.append(df, ignore_index=True)
        assert_frame_equal(result, expected)

        # different columns
        dicts = [{'foo': 1, 'bar': 2, 'baz': 3, 'peekaboo': 4},
                 {'foo': 5, 'bar': 6, 'baz': 7, 'peekaboo': 8}]
        result = df.append(dicts, ignore_index=True)
        expected = df.append(DataFrame(dicts), ignore_index=True)
        assert_frame_equal(result, expected)

    def test_append_empty_dataframe(self):

        # Empty df append empty df
        df1 = DataFrame([])
        df2 = DataFrame([])
        result = df1.append(df2)
        expected = df1.copy()
        assert_frame_equal(result, expected)

        # Non-empty df append empty df
        df1 = DataFrame(np.random.randn(5, 2))
        df2 = DataFrame()
        result = df1.append(df2)
        expected = df1.copy()
        assert_frame_equal(result, expected)

        # Empty df with columns append empty df
        df1 = DataFrame(columns=['bar', 'foo'])
        df2 = DataFrame()
        result = df1.append(df2)
        expected = df1.copy()
        assert_frame_equal(result, expected)

        # Non-Empty df with columns append empty df
        df1 = DataFrame(np.random.randn(5, 2), columns=['bar', 'foo'])
        df2 = DataFrame()
        result = df1.append(df2)
        expected = df1.copy()
        assert_frame_equal(result, expected)

    def test_append_dtypes(self):

        # GH 5754
        # row appends of different dtypes (so need to do by-item)
        # can sometimes infer the correct type

        df1 = DataFrame({ 'bar' : Timestamp('20130101') }, index=lrange(5))
        df2 = DataFrame()
        result = df1.append(df2)
        expected = df1.copy()
        assert_frame_equal(result, expected)

        df1 = DataFrame({ 'bar' : Timestamp('20130101') }, index=lrange(1))
        df2 = DataFrame({ 'bar' : 'foo' }, index=lrange(1,2))
        result = df1.append(df2)
        expected = DataFrame({ 'bar' : [ Timestamp('20130101'), 'foo' ]})
        assert_frame_equal(result, expected)

        df1 = DataFrame({ 'bar' : Timestamp('20130101') }, index=lrange(1))
        df2 = DataFrame({ 'bar' : np.nan }, index=lrange(1,2))
        result = df1.append(df2)
        expected = DataFrame({ 'bar' : Series([ Timestamp('20130101'), np.nan ],dtype='M8[ns]') })
        assert_frame_equal(result, expected)

        df1 = DataFrame({ 'bar' : Timestamp('20130101') }, index=lrange(1))
        df2 = DataFrame({ 'bar' : np.nan }, index=lrange(1,2), dtype=object)
        result = df1.append(df2)
        expected = DataFrame({ 'bar' : Series([ Timestamp('20130101'), np.nan ],dtype='M8[ns]') })
        assert_frame_equal(result, expected)

        df1 = DataFrame({ 'bar' : np.nan }, index=lrange(1))
        df2 = DataFrame({ 'bar' : Timestamp('20130101') }, index=lrange(1,2))
        result = df1.append(df2)
        expected = DataFrame({ 'bar' : Series([ np.nan, Timestamp('20130101')] ,dtype='M8[ns]') })
        assert_frame_equal(result, expected)

        df1 = DataFrame({ 'bar' : Timestamp('20130101') }, index=lrange(1))
        df2 = DataFrame({ 'bar' : 1 }, index=lrange(1,2), dtype=object)
        result = df1.append(df2)
        expected = DataFrame({ 'bar' : Series([ Timestamp('20130101'), 1 ]) })
        assert_frame_equal(result, expected)

    def test_asfreq(self):
        offset_monthly = self.tsframe.asfreq(datetools.bmonthEnd)
        rule_monthly = self.tsframe.asfreq('BM')

        assert_almost_equal(offset_monthly['A'], rule_monthly['A'])

        filled = rule_monthly.asfreq('B', method='pad')
        # TODO: actually check that this worked.

        # don't forget!
        filled_dep = rule_monthly.asfreq('B', method='pad')

        # test does not blow up on length-0 DataFrame
        zero_length = self.tsframe.reindex([])
        result = zero_length.asfreq('BM')
        self.assertIsNot(result, zero_length)

    def test_asfreq_datetimeindex(self):
        df = DataFrame({'A': [1, 2, 3]},
                       index=[datetime(2011, 11, 1), datetime(2011, 11, 2),
                              datetime(2011, 11, 3)])
        df = df.asfreq('B')
        tm.assertIsInstance(df.index, DatetimeIndex)

        ts = df['A'].asfreq('B')
        tm.assertIsInstance(ts.index, DatetimeIndex)

    def test_at_time_between_time_datetimeindex(self):
        index = date_range("2012-01-01", "2012-01-05", freq='30min')
        df = DataFrame(randn(len(index), 5), index=index)
        akey = time(12, 0, 0)
        bkey = slice(time(13, 0, 0), time(14, 0, 0))
        ainds = [24, 72, 120, 168]
        binds = [26, 27, 28, 74, 75, 76, 122, 123, 124, 170, 171, 172]

        result = df.at_time(akey)
        expected = df.ix[akey]
        expected2 = df.ix[ainds]
        assert_frame_equal(result, expected)
        assert_frame_equal(result, expected2)
        self.assertEqual(len(result), 4)

        result = df.between_time(bkey.start, bkey.stop)
        expected = df.ix[bkey]
        expected2 = df.ix[binds]
        assert_frame_equal(result, expected)
        assert_frame_equal(result, expected2)
        self.assertEqual(len(result), 12)

        result = df.copy()
        result.ix[akey] = 0
        result = result.ix[akey]
        expected = df.ix[akey].copy()
        expected.ix[:] = 0
        assert_frame_equal(result, expected)

        result = df.copy()
        result.ix[akey] = 0
        result.ix[akey] = df.ix[ainds]
        assert_frame_equal(result, df)

        result = df.copy()
        result.ix[bkey] = 0
        result = result.ix[bkey]
        expected = df.ix[bkey].copy()
        expected.ix[:] = 0
        assert_frame_equal(result, expected)

        result = df.copy()
        result.ix[bkey] = 0
        result.ix[bkey] = df.ix[binds]
        assert_frame_equal(result, df)

    def test_as_matrix(self):
        frame = self.frame
        mat = frame.as_matrix()

        frameCols = frame.columns
        for i, row in enumerate(mat):
            for j, value in enumerate(row):
                col = frameCols[j]
                if np.isnan(value):
                    self.assertTrue(np.isnan(frame[col][i]))
                else:
                    self.assertEqual(value, frame[col][i])

        # mixed type
        mat = self.mixed_frame.as_matrix(['foo', 'A'])
        self.assertEqual(mat[0, 0], 'bar')

        df = DataFrame({'real': [1, 2, 3], 'complex': [1j, 2j, 3j]})
        mat = df.as_matrix()
        self.assertEqual(mat[0, 0], 1j)

        # single block corner case
        mat = self.frame.as_matrix(['A', 'B'])
        expected = self.frame.reindex(columns=['A', 'B']).values
        assert_almost_equal(mat, expected)

    def test_as_matrix_duplicates(self):
        df = DataFrame([[1, 2, 'a', 'b'],
                        [1, 2, 'a', 'b']],
                       columns=['one', 'one', 'two', 'two'])

        result = df.values
        expected = np.array([[1, 2, 'a', 'b'], [1, 2, 'a', 'b']],
                            dtype=object)

        self.assertTrue(np.array_equal(result, expected))

    def test_ftypes(self):
        frame = self.mixed_float
        expected = Series(dict(A = 'float32:dense',
                               B = 'float32:dense',
                               C = 'float16:dense',
                               D = 'float64:dense')).sort_values()
        result = frame.ftypes.sort_values()
        assert_series_equal(result,expected)

    def test_values(self):
        self.frame.values[:, 0] = 5.
        self.assertTrue((self.frame.values[:, 0] == 5).all())

    def test_deepcopy(self):
        cp = deepcopy(self.frame)
        series = cp['A']
        series[:] = 10
        for idx, value in compat.iteritems(series):
            self.assertNotEqual(self.frame['A'][idx], value)

    def test_copy(self):
        cop = self.frame.copy()
        cop['E'] = cop['A']
        self.assertNotIn('E', self.frame)

        # copy objects
        copy = self.mixed_frame.copy()
        self.assertIsNot(copy._data, self.mixed_frame._data)

    def test_drop_names(self):
        df = DataFrame([[1, 2, 3],[3, 4, 5],[5, 6, 7]], index=['a', 'b', 'c'],
                       columns=['d', 'e', 'f'])
        df.index.name, df.columns.name = 'first', 'second'
        df_dropped_b = df.drop('b')
        df_dropped_e = df.drop('e', axis=1)
        df_inplace_b, df_inplace_e = df.copy(), df.copy()
        df_inplace_b.drop('b', inplace=True)
        df_inplace_e.drop('e', axis=1, inplace=True)
        for obj in (df_dropped_b, df_dropped_e, df_inplace_b, df_inplace_e):
            self.assertEqual(obj.index.name, 'first')
            self.assertEqual(obj.columns.name, 'second')
        self.assertEqual(list(df.columns), ['d', 'e', 'f'])

        self.assertRaises(ValueError, df.drop, ['g'])
        self.assertRaises(ValueError, df.drop, ['g'], 1)

        # errors = 'ignore'
        dropped = df.drop(['g'], errors='ignore')
        expected = Index(['a', 'b', 'c'], name='first')
        self.assert_index_equal(dropped.index, expected)

        dropped = df.drop(['b', 'g'], errors='ignore')
        expected = Index(['a', 'c'], name='first')
        self.assert_index_equal(dropped.index, expected)

        dropped = df.drop(['g'], axis=1, errors='ignore')
        expected = Index(['d', 'e', 'f'], name='second')
        self.assert_index_equal(dropped.columns, expected)

        dropped = df.drop(['d', 'g'], axis=1, errors='ignore')
        expected = Index(['e', 'f'], name='second')
        self.assert_index_equal(dropped.columns, expected)

    def test_dropEmptyRows(self):
        N = len(self.frame.index)
        mat = randn(N)
        mat[:5] = nan

        frame = DataFrame({'foo': mat}, index=self.frame.index)
        original = Series(mat, index=self.frame.index, name='foo')
        expected = original.dropna()
        inplace_frame1, inplace_frame2 = frame.copy(), frame.copy()

        smaller_frame = frame.dropna(how='all')
        # check that original was preserved
        assert_series_equal(frame['foo'], original)
        inplace_frame1.dropna(how='all', inplace=True)
        assert_series_equal(smaller_frame['foo'], expected)
        assert_series_equal(inplace_frame1['foo'], expected)

        smaller_frame = frame.dropna(how='all', subset=['foo'])
        inplace_frame2.dropna(how='all', subset=['foo'], inplace=True)
        assert_series_equal(smaller_frame['foo'], expected)
        assert_series_equal(inplace_frame2['foo'], expected)

    def test_dropIncompleteRows(self):
        N = len(self.frame.index)
        mat = randn(N)
        mat[:5] = nan

        frame = DataFrame({'foo': mat}, index=self.frame.index)
        frame['bar'] = 5
        original = Series(mat, index=self.frame.index, name='foo')
        inp_frame1, inp_frame2 = frame.copy(), frame.copy()

        smaller_frame = frame.dropna()
        assert_series_equal(frame['foo'], original)
        inp_frame1.dropna(inplace=True)
        self.assert_numpy_array_equal(smaller_frame['foo'], mat[5:])
        self.assert_numpy_array_equal(inp_frame1['foo'], mat[5:])

        samesize_frame = frame.dropna(subset=['bar'])
        assert_series_equal(frame['foo'], original)
        self.assertTrue((frame['bar'] == 5).all())
        inp_frame2.dropna(subset=['bar'], inplace=True)
        self.assertTrue(samesize_frame.index.equals(self.frame.index))
        self.assertTrue(inp_frame2.index.equals(self.frame.index))

    def test_dropna(self):
        df = DataFrame(np.random.randn(6, 4))
        df[2][:2] = nan

        dropped = df.dropna(axis=1)
        expected = df.ix[:, [0, 1, 3]]
        inp = df.copy()
        inp.dropna(axis=1, inplace=True)
        assert_frame_equal(dropped, expected)
        assert_frame_equal(inp, expected)

        dropped = df.dropna(axis=0)
        expected = df.ix[lrange(2, 6)]
        inp = df.copy()
        inp.dropna(axis=0, inplace=True)
        assert_frame_equal(dropped, expected)
        assert_frame_equal(inp, expected)

        # threshold
        dropped = df.dropna(axis=1, thresh=5)
        expected = df.ix[:, [0, 1, 3]]
        inp = df.copy()
        inp.dropna(axis=1, thresh=5, inplace=True)
        assert_frame_equal(dropped, expected)
        assert_frame_equal(inp, expected)

        dropped = df.dropna(axis=0, thresh=4)
        expected = df.ix[lrange(2, 6)]
        inp = df.copy()
        inp.dropna(axis=0, thresh=4, inplace=True)
        assert_frame_equal(dropped, expected)
        assert_frame_equal(inp, expected)

        dropped = df.dropna(axis=1, thresh=4)
        assert_frame_equal(dropped, df)

        dropped = df.dropna(axis=1, thresh=3)
        assert_frame_equal(dropped, df)

        # subset
        dropped = df.dropna(axis=0, subset=[0, 1, 3])
        inp = df.copy()
        inp.dropna(axis=0, subset=[0, 1, 3], inplace=True)
        assert_frame_equal(dropped, df)
        assert_frame_equal(inp, df)

        # all
        dropped = df.dropna(axis=1, how='all')
        assert_frame_equal(dropped, df)

        df[2] = nan
        dropped = df.dropna(axis=1, how='all')
        expected = df.ix[:, [0, 1, 3]]
        assert_frame_equal(dropped, expected)

        # bad input
        self.assertRaises(ValueError, df.dropna, axis=3)


    def test_drop_and_dropna_caching(self):
        # tst that cacher updates
        original = Series([1, 2, np.nan], name='A')
        expected = Series([1, 2], dtype=original.dtype, name='A')
        df = pd.DataFrame({'A': original.values.copy()})
        df2 = df.copy()
        df['A'].dropna()
        assert_series_equal(df['A'], original)
        df['A'].dropna(inplace=True)
        assert_series_equal(df['A'], expected)
        df2['A'].drop([1])
        assert_series_equal(df2['A'], original)
        df2['A'].drop([1], inplace=True)
        assert_series_equal(df2['A'], original.drop([1]))

    def test_dropna_corner(self):
        # bad input
        self.assertRaises(ValueError, self.frame.dropna, how='foo')
        self.assertRaises(TypeError, self.frame.dropna, how=None)
        # non-existent column - 8303
        self.assertRaises(KeyError, self.frame.dropna, subset=['A','X'])

    def test_dropna_multiple_axes(self):
        df = DataFrame([[1, np.nan, 2, 3],
                        [4, np.nan, 5, 6],
                        [np.nan, np.nan, np.nan, np.nan],
                        [7, np.nan, 8, 9]])
        cp = df.copy()
        result = df.dropna(how='all', axis=[0, 1])
        result2 = df.dropna(how='all', axis=(0, 1))
        expected = df.dropna(how='all').dropna(how='all', axis=1)

        assert_frame_equal(result, expected)
        assert_frame_equal(result2, expected)
        assert_frame_equal(df, cp)

        inp = df.copy()
        inp.dropna(how='all', axis=(0, 1), inplace=True)
        assert_frame_equal(inp, expected)

    def test_drop_duplicates(self):
        df = DataFrame({'AAA': ['foo', 'bar', 'foo', 'bar',
                                'foo', 'bar', 'bar', 'foo'],
                        'B': ['one', 'one', 'two', 'two',
                              'two', 'two', 'one', 'two'],
                        'C': [1, 1, 2, 2, 2, 2, 1, 2],
                        'D': lrange(8)})

        # single column
        result = df.drop_duplicates('AAA')
        expected = df[:2]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('AAA', keep='last')
        expected = df.ix[[6, 7]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('AAA', keep=False)
        expected = df.ix[[]]
        assert_frame_equal(result, expected)
        self.assertEqual(len(result), 0)

        # deprecate take_last
        with tm.assert_produces_warning(FutureWarning):
            result = df.drop_duplicates('AAA', take_last=True)
            expected = df.ix[[6, 7]]
            assert_frame_equal(result, expected)

        # multi column
        expected = df.ix[[0, 1, 2, 3]]
        result = df.drop_duplicates(np.array(['AAA', 'B']))
        assert_frame_equal(result, expected)
        result = df.drop_duplicates(['AAA', 'B'])
        assert_frame_equal(result, expected)

        result = df.drop_duplicates(('AAA', 'B'), keep='last')
        expected = df.ix[[0, 5, 6, 7]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates(('AAA', 'B'), keep=False)
        expected = df.ix[[0]]
        assert_frame_equal(result, expected)

        # deprecate take_last
        with tm.assert_produces_warning(FutureWarning):
            result = df.drop_duplicates(('AAA', 'B'), take_last=True)
        expected = df.ix[[0, 5, 6, 7]]
        assert_frame_equal(result, expected)

        # consider everything
        df2 = df.ix[:, ['AAA', 'B', 'C']]

        result = df2.drop_duplicates()
        # in this case only
        expected = df2.drop_duplicates(['AAA', 'B'])
        assert_frame_equal(result, expected)

        result = df2.drop_duplicates(keep='last')
        expected = df2.drop_duplicates(['AAA', 'B'], keep='last')
        assert_frame_equal(result, expected)

        result = df2.drop_duplicates(keep=False)
        expected = df2.drop_duplicates(['AAA', 'B'], keep=False)
        assert_frame_equal(result, expected)

        # deprecate take_last
        with tm.assert_produces_warning(FutureWarning):
            result = df2.drop_duplicates(take_last=True)
        with tm.assert_produces_warning(FutureWarning):
            expected = df2.drop_duplicates(['AAA', 'B'], take_last=True)
        assert_frame_equal(result, expected)

        # integers
        result = df.drop_duplicates('C')
        expected = df.iloc[[0,2]]
        assert_frame_equal(result, expected)
        result = df.drop_duplicates('C',keep='last')
        expected = df.iloc[[-2,-1]]
        assert_frame_equal(result, expected)

        df['E'] = df['C'].astype('int8')
        result = df.drop_duplicates('E')
        expected = df.iloc[[0,2]]
        assert_frame_equal(result, expected)
        result = df.drop_duplicates('E',keep='last')
        expected = df.iloc[[-2,-1]]
        assert_frame_equal(result, expected)

        # GH 11376
        df = pd.DataFrame({'x': [7, 6, 3, 3, 4, 8, 0],
                           'y': [0, 6, 5, 5, 9, 1, 2]})
        expected = df.loc[df.index != 3]
        assert_frame_equal(df.drop_duplicates(), expected)

        df = pd.DataFrame([[1 , 0], [0, 2]])
        assert_frame_equal(df.drop_duplicates(), df)

        df = pd.DataFrame([[-2, 0], [0, -4]])
        assert_frame_equal(df.drop_duplicates(), df)

        x = np.iinfo(np.int64).max / 3 * 2
        df = pd.DataFrame([[-x, x], [0, x + 4]])
        assert_frame_equal(df.drop_duplicates(), df)

        df = pd.DataFrame([[-x, x], [x, x + 4]])
        assert_frame_equal(df.drop_duplicates(), df)

        # GH 11864
        df = pd.DataFrame([i] * 9 for i in range(16))
        df = df.append([[1] + [0] * 8], ignore_index=True)

        for keep in ['first', 'last', False]:
            assert_equal(df.duplicated(keep=keep).sum(), 0)

    def test_drop_duplicates_for_take_all(self):
        df = DataFrame({'AAA': ['foo', 'bar', 'baz', 'bar',
                                'foo', 'bar', 'qux', 'foo'],
                        'B': ['one', 'one', 'two', 'two',
                              'two', 'two', 'one', 'two'],
                        'C': [1, 1, 2, 2, 2, 2, 1, 2],
                        'D': lrange(8)})

        # single column
        result = df.drop_duplicates('AAA')
        expected = df.iloc[[0, 1, 2, 6]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('AAA', keep='last')
        expected = df.iloc[[2, 5, 6, 7]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('AAA', keep=False)
        expected = df.iloc[[2, 6]]
        assert_frame_equal(result, expected)

        # multiple columns
        result = df.drop_duplicates(['AAA', 'B'])
        expected = df.iloc[[0, 1, 2, 3, 4, 6]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates(['AAA', 'B'], keep='last')
        expected = df.iloc[[0, 1, 2, 5, 6, 7]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates(['AAA', 'B'], keep=False)
        expected = df.iloc[[0, 1, 2, 6]]
        assert_frame_equal(result, expected)

    def test_drop_duplicates_deprecated_warning(self):
        df = DataFrame({'AAA': ['foo', 'bar', 'foo', 'bar',
                                'foo', 'bar', 'bar', 'foo'],
                        'B': ['one', 'one', 'two', 'two',
                              'two', 'two', 'one', 'two'],
                        'C': [1, 1, 2, 2, 2, 2, 1, 2],
                        'D': lrange(8)})
        expected = df[:2]

        # Raises warning
        with tm.assert_produces_warning(False):
            result = df.drop_duplicates(subset='AAA')
        assert_frame_equal(result, expected)

        with tm.assert_produces_warning(FutureWarning):
            result = df.drop_duplicates(cols='AAA')
        assert_frame_equal(result, expected)

        # Does not allow both subset and cols
        self.assertRaises(TypeError, df.drop_duplicates,
                          kwargs={'cols': 'AAA', 'subset': 'B'})

        # Does not allow unknown kwargs
        self.assertRaises(TypeError, df.drop_duplicates,
                          kwargs={'subset': 'AAA', 'bad_arg': True})

        # deprecate take_last
        # Raises warning
        with tm.assert_produces_warning(FutureWarning):
            result = df.drop_duplicates(take_last=False, subset='AAA')
        assert_frame_equal(result, expected)

        self.assertRaises(ValueError, df.drop_duplicates, keep='invalid_name')

    def test_drop_duplicates_tuple(self):
        df = DataFrame({('AA', 'AB'): ['foo', 'bar', 'foo', 'bar',
                                       'foo', 'bar', 'bar', 'foo'],
                        'B': ['one', 'one', 'two', 'two',
                              'two', 'two', 'one', 'two'],
                        'C': [1, 1, 2, 2, 2, 2, 1, 2],
                        'D': lrange(8)})

        # single column
        result = df.drop_duplicates(('AA', 'AB'))
        expected = df[:2]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates(('AA', 'AB'), keep='last')
        expected = df.ix[[6, 7]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates(('AA', 'AB'), keep=False)
        expected = df.ix[[]] # empty df
        self.assertEqual(len(result), 0)
        assert_frame_equal(result, expected)

        # deprecate take_last
        with tm.assert_produces_warning(FutureWarning):
            result = df.drop_duplicates(('AA', 'AB'), take_last=True)
        expected = df.ix[[6, 7]]
        assert_frame_equal(result, expected)

        # multi column
        expected = df.ix[[0, 1, 2, 3]]
        result = df.drop_duplicates((('AA', 'AB'), 'B'))
        assert_frame_equal(result, expected)

    def test_drop_duplicates_NA(self):
        # none
        df = DataFrame({'A': [None, None, 'foo', 'bar',
                              'foo', 'bar', 'bar', 'foo'],
                        'B': ['one', 'one', 'two', 'two',
                              'two', 'two', 'one', 'two'],
                        'C': [1.0, np.nan, np.nan, np.nan, 1., 1., 1, 1.],
                        'D': lrange(8)})

        # single column
        result = df.drop_duplicates('A')
        expected = df.ix[[0, 2, 3]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('A', keep='last')
        expected = df.ix[[1, 6, 7]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('A', keep=False)
        expected = df.ix[[]] # empty df
        assert_frame_equal(result, expected)
        self.assertEqual(len(result), 0)

        # deprecate take_last
        with tm.assert_produces_warning(FutureWarning):
            result = df.drop_duplicates('A', take_last=True)
        expected = df.ix[[1, 6, 7]]
        assert_frame_equal(result, expected)

        # multi column
        result = df.drop_duplicates(['A', 'B'])
        expected = df.ix[[0, 2, 3, 6]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates(['A', 'B'], keep='last')
        expected = df.ix[[1, 5, 6, 7]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates(['A', 'B'], keep=False)
        expected = df.ix[[6]]
        assert_frame_equal(result, expected)

        # deprecate take_last
        with tm.assert_produces_warning(FutureWarning):
            result = df.drop_duplicates(['A', 'B'], take_last=True)
        expected = df.ix[[1, 5, 6, 7]]
        assert_frame_equal(result, expected)

        # nan
        df = DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                              'foo', 'bar', 'bar', 'foo'],
                        'B': ['one', 'one', 'two', 'two',
                              'two', 'two', 'one', 'two'],
                        'C': [1.0, np.nan, np.nan, np.nan, 1., 1., 1, 1.],
                        'D': lrange(8)})

        # single column
        result = df.drop_duplicates('C')
        expected = df[:2]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('C', keep='last')
        expected = df.ix[[3, 7]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('C', keep=False)
        expected = df.ix[[]] # empty df
        assert_frame_equal(result, expected)
        self.assertEqual(len(result), 0)

        # deprecate take_last
        with tm.assert_produces_warning(FutureWarning):
            result = df.drop_duplicates('C', take_last=True)
        expected = df.ix[[3, 7]]
        assert_frame_equal(result, expected)

        # multi column
        result = df.drop_duplicates(['C', 'B'])
        expected = df.ix[[0, 1, 2, 4]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates(['C', 'B'], keep='last')
        expected = df.ix[[1, 3, 6, 7]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates(['C', 'B'], keep=False)
        expected = df.ix[[1]]
        assert_frame_equal(result, expected)

        # deprecate take_last
        with tm.assert_produces_warning(FutureWarning):
            result = df.drop_duplicates(['C', 'B'], take_last=True)
        expected = df.ix[[1, 3, 6, 7]]
        assert_frame_equal(result, expected)

    def test_drop_duplicates_NA_for_take_all(self):
        # none
        df = DataFrame({'A': [None, None, 'foo', 'bar',
                              'foo', 'baz', 'bar', 'qux'],
                        'C': [1.0, np.nan, np.nan, np.nan, 1., 2., 3, 1.]})

        # single column
        result = df.drop_duplicates('A')
        expected = df.iloc[[0, 2, 3, 5, 7]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('A', keep='last')
        expected = df.iloc[[1, 4, 5, 6, 7]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('A', keep=False)
        expected = df.iloc[[5, 7]]
        assert_frame_equal(result, expected)

        # nan

        # single column
        result = df.drop_duplicates('C')
        expected = df.iloc[[0, 1, 5, 6]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('C', keep='last')
        expected = df.iloc[[3, 5, 6, 7]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('C', keep=False)
        expected = df.iloc[[5, 6]]
        assert_frame_equal(result, expected)

    def test_drop_duplicates_inplace(self):
        orig = DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                                'foo', 'bar', 'bar', 'foo'],
                          'B': ['one', 'one', 'two', 'two',
                                'two', 'two', 'one', 'two'],
                          'C': [1, 1, 2, 2, 2, 2, 1, 2],
                          'D': lrange(8)})

        # single column
        df = orig.copy()
        df.drop_duplicates('A', inplace=True)
        expected = orig[:2]
        result = df
        assert_frame_equal(result, expected)

        df = orig.copy()
        df.drop_duplicates('A', keep='last', inplace=True)
        expected = orig.ix[[6, 7]]
        result = df
        assert_frame_equal(result, expected)

        df = orig.copy()
        df.drop_duplicates('A', keep=False, inplace=True)
        expected = orig.ix[[]]
        result = df
        assert_frame_equal(result, expected)
        self.assertEqual(len(df), 0)

        # deprecate take_last
        df = orig.copy()
        with tm.assert_produces_warning(FutureWarning):
            df.drop_duplicates('A', take_last=True, inplace=True)
        expected = orig.ix[[6, 7]]
        result = df
        assert_frame_equal(result, expected)

        # multi column
        df = orig.copy()
        df.drop_duplicates(['A', 'B'], inplace=True)
        expected = orig.ix[[0, 1, 2, 3]]
        result = df
        assert_frame_equal(result, expected)

        df = orig.copy()
        df.drop_duplicates(['A', 'B'], keep='last', inplace=True)
        expected = orig.ix[[0, 5, 6, 7]]
        result = df
        assert_frame_equal(result, expected)

        df = orig.copy()
        df.drop_duplicates(['A', 'B'], keep=False, inplace=True)
        expected = orig.ix[[0]]
        result = df
        assert_frame_equal(result, expected)

        # deprecate take_last
        df = orig.copy()
        with tm.assert_produces_warning(FutureWarning):
            df.drop_duplicates(['A', 'B'], take_last=True, inplace=True)
        expected = orig.ix[[0, 5, 6, 7]]
        result = df
        assert_frame_equal(result, expected)

        # consider everything
        orig2 = orig.ix[:, ['A', 'B', 'C']].copy()

        df2 = orig2.copy()
        df2.drop_duplicates(inplace=True)
        # in this case only
        expected = orig2.drop_duplicates(['A', 'B'])
        result = df2
        assert_frame_equal(result, expected)

        df2 = orig2.copy()
        df2.drop_duplicates(keep='last', inplace=True)
        expected = orig2.drop_duplicates(['A', 'B'], keep='last')
        result = df2
        assert_frame_equal(result, expected)

        df2 = orig2.copy()
        df2.drop_duplicates(keep=False, inplace=True)
        expected = orig2.drop_duplicates(['A', 'B'], keep=False)
        result = df2
        assert_frame_equal(result, expected)

        # deprecate take_last
        df2 = orig2.copy()
        with tm.assert_produces_warning(FutureWarning):
            df2.drop_duplicates(take_last=True, inplace=True)
        with tm.assert_produces_warning(FutureWarning):
            expected = orig2.drop_duplicates(['A', 'B'], take_last=True)
        result = df2
        assert_frame_equal(result, expected)

    def test_duplicated_deprecated_warning(self):
        df = DataFrame({'AAA': ['foo', 'bar', 'foo', 'bar',
                                'foo', 'bar', 'bar', 'foo'],
                        'B': ['one', 'one', 'two', 'two',
                              'two', 'two', 'one', 'two'],
                        'C': [1, 1, 2, 2, 2, 2, 1, 2],
                        'D': lrange(8)})

        # Raises warning
        with tm.assert_produces_warning(False):
            result = df.duplicated(subset='AAA')

        with tm.assert_produces_warning(FutureWarning):
            result = df.duplicated(cols='AAA')

        # Does not allow both subset and cols
        self.assertRaises(TypeError, df.duplicated,
                          kwargs={'cols': 'AAA', 'subset': 'B'})

        # Does not allow unknown kwargs
        self.assertRaises(TypeError, df.duplicated,
                          kwargs={'subset': 'AAA', 'bad_arg': True})

    def test_drop_col_still_multiindex(self):
        arrays = [['a', 'b', 'c', 'top'],
                  ['', '', '', 'OD'],
                  ['', '', '', 'wx']]

        tuples = sorted(zip(*arrays))
        index = MultiIndex.from_tuples(tuples)

        df = DataFrame(randn(3, 4), columns=index)
        del df[('a', '', '')]
        assert(isinstance(df.columns, MultiIndex))

    def test_drop(self):
        simple = DataFrame({"A": [1, 2, 3, 4], "B": [0, 1, 2, 3]})
        assert_frame_equal(simple.drop("A", axis=1), simple[['B']])
        assert_frame_equal(simple.drop(["A", "B"], axis='columns'),
                           simple[[]])
        assert_frame_equal(simple.drop([0, 1, 3], axis=0), simple.ix[[2], :])
        assert_frame_equal(simple.drop([0, 3], axis='index'), simple.ix[[1, 2], :])

        self.assertRaises(ValueError, simple.drop, 5)
        self.assertRaises(ValueError, simple.drop, 'C', 1)
        self.assertRaises(ValueError, simple.drop, [1, 5])
        self.assertRaises(ValueError, simple.drop, ['A', 'C'], 1)

        # errors = 'ignore'
        assert_frame_equal(simple.drop(5, errors='ignore'), simple)
        assert_frame_equal(simple.drop([0, 5], errors='ignore'),
                           simple.ix[[1, 2, 3], :])
        assert_frame_equal(simple.drop('C', axis=1, errors='ignore'), simple)
        assert_frame_equal(simple.drop(['A', 'C'], axis=1, errors='ignore'),
                           simple[['B']])

        #non-unique - wheee!
        nu_df = DataFrame(lzip(range(3), range(-3, 1), list('abc')),
                          columns=['a', 'a', 'b'])
        assert_frame_equal(nu_df.drop('a', axis=1), nu_df[['b']])
        assert_frame_equal(nu_df.drop('b', axis='columns'), nu_df['a'])

        nu_df = nu_df.set_index(pd.Index(['X', 'Y', 'X']))
        nu_df.columns = list('abc')
        assert_frame_equal(nu_df.drop('X', axis='rows'), nu_df.ix[["Y"], :])
        assert_frame_equal(nu_df.drop(['X', 'Y'], axis=0), nu_df.ix[[], :])

        # inplace cache issue
        # GH 5628
        df = pd.DataFrame(np.random.randn(10,3), columns=list('abc'))
        expected = df[~(df.b>0)]
        df.drop(labels=df[df.b>0].index, inplace=True)
        assert_frame_equal(df,expected)

    def test_fillna(self):
        self.tsframe.ix[:5,'A'] = nan
        self.tsframe.ix[-5:,'A'] = nan

        zero_filled = self.tsframe.fillna(0)
        self.assertTrue((zero_filled.ix[:5,'A'] == 0).all())

        padded = self.tsframe.fillna(method='pad')
        self.assertTrue(np.isnan(padded.ix[:5,'A']).all())
        self.assertTrue((padded.ix[-5:,'A'] == padded.ix[-5,'A']).all())

        # mixed type
        self.mixed_frame.ix[5:20,'foo'] = nan
        self.mixed_frame.ix[-10:,'A'] = nan
        result = self.mixed_frame.fillna(value=0)
        result = self.mixed_frame.fillna(method='pad')

        self.assertRaises(ValueError, self.tsframe.fillna)
        self.assertRaises(ValueError, self.tsframe.fillna, 5, method='ffill')

        # mixed numeric (but no float16)
        mf = self.mixed_float.reindex(columns=['A','B','D'])
        mf.ix[-10:,'A'] = nan
        result = mf.fillna(value=0)
        _check_mixed_float(result, dtype = dict(C = None))

        result = mf.fillna(method='pad')
        _check_mixed_float(result, dtype = dict(C = None))

        # empty frame (GH #2778)
        df = DataFrame(columns=['x'])
        for m in ['pad','backfill']:
            df.x.fillna(method=m,inplace=1)
            df.x.fillna(method=m)

        # with different dtype (GH3386)
        df = DataFrame([['a','a',np.nan,'a'],['b','b',np.nan,'b'],['c','c',np.nan,'c']])

        result = df.fillna({ 2: 'foo' })
        expected = DataFrame([['a','a','foo','a'],['b','b','foo','b'],['c','c','foo','c']])
        assert_frame_equal(result, expected)

        df.fillna({ 2: 'foo' }, inplace=True)
        assert_frame_equal(df, expected)

        # limit and value
        df = DataFrame(np.random.randn(10,3))
        df.iloc[2:7,0] = np.nan
        df.iloc[3:5,2] = np.nan

        expected = df.copy()
        expected.iloc[2,0] = 999
        expected.iloc[3,2] = 999
        result = df.fillna(999,limit=1)
        assert_frame_equal(result, expected)

        # with datelike
        # GH 6344
        df = DataFrame({
            'Date':[pd.NaT, Timestamp("2014-1-1")],
            'Date2':[ Timestamp("2013-1-1"), pd.NaT]
            })

        expected = df.copy()
        expected['Date'] = expected['Date'].fillna(df.ix[0,'Date2'])
        result = df.fillna(value={'Date':df['Date2']})
        assert_frame_equal(result, expected)

    def test_fillna_dtype_conversion(self):
        # make sure that fillna on an empty frame works
        df = DataFrame(index=["A","B","C"], columns = [1,2,3,4,5])
        result = df.get_dtype_counts().sort_values()
        expected = Series({ 'object' : 5 })
        assert_series_equal(result, expected)

        result = df.fillna(1)
        expected = DataFrame(1, index=["A","B","C"], columns = [1,2,3,4,5])
        result = result.get_dtype_counts().sort_values()
        expected = Series({ 'int64' : 5 })
        assert_series_equal(result, expected)

        # empty block
        df = DataFrame(index=lrange(3),columns=['A','B'],dtype='float64')
        result = df.fillna('nan')
        expected = DataFrame('nan',index=lrange(3),columns=['A','B'])
        assert_frame_equal(result, expected)

        # equiv of replace
        df = DataFrame(dict(A = [1,np.nan], B = [1.,2.]))
        for v in ['',1,np.nan,1.0]:
            expected = df.replace(np.nan,v)
            result = df.fillna(v)
            assert_frame_equal(result, expected)

    def test_fillna_datetime_columns(self):
        # GH 7095
        df = pd.DataFrame({'A': [-1, -2, np.nan],
                           'B': date_range('20130101', periods=3),
                           'C': ['foo', 'bar', None],
                           'D': ['foo2', 'bar2', None]},
                          index=date_range('20130110', periods=3))
        result = df.fillna('?')
        expected = pd.DataFrame({'A': [-1, -2, '?'],
                                 'B': date_range('20130101', periods=3),
                                 'C': ['foo', 'bar', '?'],
                                 'D': ['foo2', 'bar2', '?']},
                                index=date_range('20130110', periods=3))
        self.assert_frame_equal(result, expected)

        df = pd.DataFrame({'A': [-1, -2, np.nan],
                           'B': [pd.Timestamp('2013-01-01'), pd.Timestamp('2013-01-02'), pd.NaT],
                           'C': ['foo', 'bar', None],
                           'D': ['foo2', 'bar2', None]},
                          index=date_range('20130110', periods=3))
        result = df.fillna('?')
        expected = pd.DataFrame({'A': [-1, -2, '?'],
                                 'B': [pd.Timestamp('2013-01-01'), pd.Timestamp('2013-01-02'), '?'],
                                 'C': ['foo', 'bar', '?'],
                                 'D': ['foo2', 'bar2', '?']},
                                index=date_range('20130110', periods=3))
        self.assert_frame_equal(result, expected)

    def test_ffill(self):
        self.tsframe['A'][:5] = nan
        self.tsframe['A'][-5:] = nan

        assert_frame_equal(self.tsframe.ffill(),
                           self.tsframe.fillna(method='ffill'))

    def test_bfill(self):
        self.tsframe['A'][:5] = nan
        self.tsframe['A'][-5:] = nan

        assert_frame_equal(self.tsframe.bfill(),
                           self.tsframe.fillna(method='bfill'))

    def test_fillna_skip_certain_blocks(self):
        # don't try to fill boolean, int blocks

        df = DataFrame(np.random.randn(10, 4).astype(int))

        # it works!
        df.fillna(np.nan)

    def test_fillna_inplace(self):
        df = DataFrame(np.random.randn(10, 4))
        df[1][:4] = np.nan
        df[3][-4:] = np.nan

        expected = df.fillna(value=0)
        self.assertIsNot(expected, df)

        df.fillna(value=0, inplace=True)
        assert_frame_equal(df, expected)

        df[1][:4] = np.nan
        df[3][-4:] = np.nan
        expected = df.fillna(method='ffill')
        self.assertIsNot(expected, df)

        df.fillna(method='ffill', inplace=True)
        assert_frame_equal(df, expected)

    def test_fillna_dict_series(self):
        df = DataFrame({'a': [nan, 1, 2, nan, nan],
                        'b': [1, 2, 3, nan, nan],
                        'c': [nan, 1, 2, 3, 4]})

        result = df.fillna({'a': 0, 'b': 5})

        expected = df.copy()
        expected['a'] = expected['a'].fillna(0)
        expected['b'] = expected['b'].fillna(5)
        assert_frame_equal(result, expected)

        # it works
        result = df.fillna({'a': 0, 'b': 5, 'd': 7})

        # Series treated same as dict
        result = df.fillna(df.max())
        expected = df.fillna(df.max().to_dict())
        assert_frame_equal(result, expected)

        # disable this for now
        with assertRaisesRegexp(NotImplementedError, 'column by column'):
            df.fillna(df.max(1), axis=1)

    def test_fillna_dataframe(self):
        # GH 8377
        df = DataFrame({'a': [nan, 1, 2, nan, nan],
                        'b': [1, 2, 3, nan, nan],
                        'c': [nan, 1, 2, 3, 4]},
                       index = list('VWXYZ'))

        # df2 may have different index and columns
        df2 = DataFrame({'a': [nan, 10, 20, 30, 40],
                         'b': [50, 60, 70, 80, 90],
                         'foo': ['bar']*5},
                        index = list('VWXuZ'))

        result = df.fillna(df2)

        # only those columns and indices which are shared get filled
        expected = DataFrame({'a': [nan, 1, 2, nan, 40],
                              'b': [1, 2, 3, nan, 90],
                              'c': [nan, 1, 2, 3, 4]},
                             index = list('VWXYZ'))

        assert_frame_equal(result, expected)

    def test_fillna_columns(self):
        df = DataFrame(np.random.randn(10, 10))
        df.values[:, ::2] = np.nan

        result = df.fillna(method='ffill', axis=1)
        expected = df.T.fillna(method='pad').T
        assert_frame_equal(result, expected)

        df.insert(6, 'foo', 5)
        result = df.fillna(method='ffill', axis=1)
        expected = df.astype(float).fillna(method='ffill', axis=1)
        assert_frame_equal(result, expected)


    def test_fillna_invalid_method(self):
        with assertRaisesRegexp(ValueError, 'ffil'):
            self.frame.fillna(method='ffil')

    def test_fillna_invalid_value(self):
        # list
        self.assertRaises(TypeError, self.frame.fillna, [1, 2])
        # tuple
        self.assertRaises(TypeError, self.frame.fillna, (1, 2))
        # frame with series
        self.assertRaises(ValueError, self.frame.iloc[:,0].fillna, self.frame)

    def test_interpolate(self):
        pass

    def test_combine_multiple_frames_dtypes(self):

        # GH 2759
        A = DataFrame(data=np.ones((10, 2)), columns=['foo', 'bar'], dtype=np.float64)
        B = DataFrame(data=np.ones((10, 2)), dtype=np.float32)
        results = pd.concat((A, B), axis=1).get_dtype_counts()
        expected = Series(dict( float64 = 2, float32 = 2 ))
        assert_series_equal(results,expected)

    def test_ops(self):

        # tst ops and reversed ops in evaluation
        # GH7198

        # smaller hits python, larger hits numexpr
        for n in [ 4, 4000 ]:

            df = DataFrame(1,index=range(n),columns=list('abcd'))
            df.iloc[0] = 2
            m = df.mean()

            for op_str, op, rop in [('+','__add__','__radd__'),
                                    ('-','__sub__','__rsub__'),
                                    ('*','__mul__','__rmul__'),
                                    ('/','__truediv__','__rtruediv__')]:

                base = DataFrame(np.tile(m.values,n).reshape(n,-1),columns=list('abcd'))
                expected = eval("base{op}df".format(op=op_str))

                # ops as strings
                result = eval("m{op}df".format(op=op_str))
                assert_frame_equal(result,expected)

                # these are commutative
                if op in ['+','*']:
                    result = getattr(df,op)(m)
                    assert_frame_equal(result,expected)

                # these are not
                elif op in ['-','/']:
                    result = getattr(df,rop)(m)
                    assert_frame_equal(result,expected)

        # GH7192
        df = DataFrame(dict(A=np.random.randn(25000)))
        df.iloc[0:5] = np.nan
        expected = (1-np.isnan(df.iloc[0:25]))
        result = (1-np.isnan(df)).iloc[0:25]
        assert_frame_equal(result,expected)

    def test_truncate(self):
        offset = datetools.bday

        ts = self.tsframe[::3]

        start, end = self.tsframe.index[3], self.tsframe.index[6]

        start_missing = self.tsframe.index[2]
        end_missing = self.tsframe.index[7]

        # neither specified
        truncated = ts.truncate()
        assert_frame_equal(truncated, ts)

        # both specified
        expected = ts[1:3]

        truncated = ts.truncate(start, end)
        assert_frame_equal(truncated, expected)

        truncated = ts.truncate(start_missing, end_missing)
        assert_frame_equal(truncated, expected)

        # start specified
        expected = ts[1:]

        truncated = ts.truncate(before=start)
        assert_frame_equal(truncated, expected)

        truncated = ts.truncate(before=start_missing)
        assert_frame_equal(truncated, expected)

        # end specified
        expected = ts[:3]

        truncated = ts.truncate(after=end)
        assert_frame_equal(truncated, expected)

        truncated = ts.truncate(after=end_missing)
        assert_frame_equal(truncated, expected)

        self.assertRaises(ValueError, ts.truncate,
                          before=ts.index[-1] - 1,
                          after=ts.index[0] +1)

    def test_truncate_copy(self):
        index = self.tsframe.index
        truncated = self.tsframe.truncate(index[5], index[10])
        truncated.values[:] = 5.
        self.assertFalse((self.tsframe.values[5:11] == 5).any())

    def test_xs(self):
        idx = self.frame.index[5]
        xs = self.frame.xs(idx)
        for item, value in compat.iteritems(xs):
            if np.isnan(value):
                self.assertTrue(np.isnan(self.frame[item][idx]))
            else:
                self.assertEqual(value, self.frame[item][idx])

        # mixed-type xs
        test_data = {
            'A': {'1': 1, '2': 2},
            'B': {'1': '1', '2': '2', '3': '3'},
        }
        frame = DataFrame(test_data)
        xs = frame.xs('1')
        self.assertEqual(xs.dtype, np.object_)
        self.assertEqual(xs['A'], 1)
        self.assertEqual(xs['B'], '1')

        with tm.assertRaises(KeyError):
            self.tsframe.xs(self.tsframe.index[0] - datetools.bday)

        # xs get column
        series = self.frame.xs('A', axis=1)
        expected = self.frame['A']
        assert_series_equal(series, expected)

        # view is returned if possible
        series = self.frame.xs('A', axis=1)
        series[:] = 5
        self.assertTrue((expected == 5).all())

    def test_xs_corner(self):
        # pathological mixed-type reordering case
        df = DataFrame(index=[0])
        df['A'] = 1.
        df['B'] = 'foo'
        df['C'] = 2.
        df['D'] = 'bar'
        df['E'] = 3.

        xs = df.xs(0)
        assert_almost_equal(xs, [1., 'foo', 2., 'bar', 3.])

        # no columns but Index(dtype=object)
        df = DataFrame(index=['a', 'b', 'c'])
        result = df.xs('a')
        expected = Series([], name='a', index=pd.Index([], dtype=object))
        assert_series_equal(result, expected)

    def test_xs_duplicates(self):
        df = DataFrame(randn(5, 2), index=['b', 'b', 'c', 'b', 'a'])

        cross = df.xs('c')
        exp = df.iloc[2]
        assert_series_equal(cross, exp)

    def test_xs_keep_level(self):
        df = DataFrame({'day': {0: 'sat', 1: 'sun'},
                        'flavour': {0: 'strawberry', 1: 'strawberry'},
                        'sales': {0: 10, 1: 12},
                        'year': {0: 2008, 1: 2008}}).set_index(['year','flavour','day'])
        result = df.xs('sat', level='day', drop_level=False)
        expected = df[:1]
        assert_frame_equal(result, expected)

        result = df.xs([2008, 'sat'], level=['year', 'day'], drop_level=False)
        assert_frame_equal(result, expected)

    def test_reindex(self):
        newFrame = self.frame.reindex(self.ts1.index)

        for col in newFrame.columns:
            for idx, val in compat.iteritems(newFrame[col]):
                if idx in self.frame.index:
                    if np.isnan(val):
                        self.assertTrue(np.isnan(self.frame[col][idx]))
                    else:
                        self.assertEqual(val, self.frame[col][idx])
                else:
                    self.assertTrue(np.isnan(val))

        for col, series in compat.iteritems(newFrame):
            self.assertTrue(tm.equalContents(series.index, newFrame.index))
        emptyFrame = self.frame.reindex(Index([]))
        self.assertEqual(len(emptyFrame.index), 0)

        # Cython code should be unit-tested directly
        nonContigFrame = self.frame.reindex(self.ts1.index[::2])

        for col in nonContigFrame.columns:
            for idx, val in compat.iteritems(nonContigFrame[col]):
                if idx in self.frame.index:
                    if np.isnan(val):
                        self.assertTrue(np.isnan(self.frame[col][idx]))
                    else:
                        self.assertEqual(val, self.frame[col][idx])
                else:
                    self.assertTrue(np.isnan(val))

        for col, series in compat.iteritems(nonContigFrame):
            self.assertTrue(tm.equalContents(series.index,
                                          nonContigFrame.index))

        # corner cases

        # Same index, copies values but not index if copy=False
        newFrame = self.frame.reindex(self.frame.index, copy=False)
        self.assertIs(newFrame.index, self.frame.index)

        # length zero
        newFrame = self.frame.reindex([])
        self.assertTrue(newFrame.empty)
        self.assertEqual(len(newFrame.columns), len(self.frame.columns))

        # length zero with columns reindexed with non-empty index
        newFrame = self.frame.reindex([])
        newFrame = newFrame.reindex(self.frame.index)
        self.assertEqual(len(newFrame.index), len(self.frame.index))
        self.assertEqual(len(newFrame.columns), len(self.frame.columns))

        # pass non-Index
        newFrame = self.frame.reindex(list(self.ts1.index))
        self.assertTrue(newFrame.index.equals(self.ts1.index))

        # copy with no axes
        result = self.frame.reindex()
        assert_frame_equal(result,self.frame)
        self.assertFalse(result is self.frame)

    def test_reindex_nan(self):
        df = pd.DataFrame([[1, 2], [3, 5], [7, 11], [9, 23]],
                          index=[2, np.nan, 1, 5],
                          columns=['joe', 'jim'])

        i, j = [np.nan, 5, 5, np.nan, 1, 2, np.nan], [1, 3, 3, 1, 2, 0, 1]
        assert_frame_equal(df.reindex(i), df.iloc[j])

        df.index = df.index.astype('object')
        assert_frame_equal(df.reindex(i), df.iloc[j], check_index_type=False)

        # GH10388
        df = pd.DataFrame({'other': ['a', 'b', np.nan, 'c'],
                           'date': ['2015-03-22', np.nan,
                                    '2012-01-08', np.nan],
                           'amount': [2, 3, 4, 5]})

        df['date'] = pd.to_datetime(df.date)
        df['delta'] = (pd.to_datetime('2015-06-18') - df['date']).shift(1)

        left = df.set_index(['delta', 'other', 'date']).reset_index()
        right = df.reindex(columns=['delta', 'other', 'date', 'amount'])
        assert_frame_equal(left, right)

    def test_reindex_name_remains(self):
        s = Series(random.rand(10))
        df = DataFrame(s, index=np.arange(len(s)))
        i = Series(np.arange(10), name='iname')

        df = df.reindex(i)
        self.assertEqual(df.index.name, 'iname')

        df = df.reindex(Index(np.arange(10), name='tmpname'))
        self.assertEqual(df.index.name, 'tmpname')

        s = Series(random.rand(10))
        df = DataFrame(s.T, index=np.arange(len(s)))
        i = Series(np.arange(10), name='iname')
        df = df.reindex(columns=i)
        self.assertEqual(df.columns.name, 'iname')

    def test_reindex_int(self):
        smaller = self.intframe.reindex(self.intframe.index[::2])

        self.assertEqual(smaller['A'].dtype, np.int64)

        bigger = smaller.reindex(self.intframe.index)
        self.assertEqual(bigger['A'].dtype, np.float64)

        smaller = self.intframe.reindex(columns=['A', 'B'])
        self.assertEqual(smaller['A'].dtype, np.int64)

    def test_reindex_like(self):
        other = self.frame.reindex(index=self.frame.index[:10],
                                   columns=['C', 'B'])

        assert_frame_equal(other, self.frame.reindex_like(other))

    def test_reindex_columns(self):
        newFrame = self.frame.reindex(columns=['A', 'B', 'E'])

        assert_series_equal(newFrame['B'], self.frame['B'])
        self.assertTrue(np.isnan(newFrame['E']).all())
        self.assertNotIn('C', newFrame)

        # length zero
        newFrame = self.frame.reindex(columns=[])
        self.assertTrue(newFrame.empty)

    def test_reindex_axes(self):

        # GH 3317, reindexing by both axes loses freq of the index
        from datetime import datetime
        df = DataFrame(np.ones((3, 3)), index=[datetime(2012, 1, 1), datetime(2012, 1, 2), datetime(2012, 1, 3)], columns=['a', 'b', 'c'])
        time_freq = date_range('2012-01-01', '2012-01-03', freq='d')
        some_cols = ['a', 'b']

        index_freq = df.reindex(index=time_freq).index.freq
        both_freq = df.reindex(index=time_freq, columns=some_cols).index.freq
        seq_freq = df.reindex(index=time_freq).reindex(columns=some_cols).index.freq
        self.assertEqual(index_freq, both_freq)
        self.assertEqual(index_freq, seq_freq)

    def test_reindex_fill_value(self):
        df = DataFrame(np.random.randn(10, 4))

        # axis=0
        result = df.reindex(lrange(15))
        self.assertTrue(np.isnan(result.values[-5:]).all())

        result = df.reindex(lrange(15), fill_value=0)
        expected = df.reindex(lrange(15)).fillna(0)
        assert_frame_equal(result, expected)

        # axis=1
        result = df.reindex(columns=lrange(5), fill_value=0.)
        expected = df.copy()
        expected[4] = 0.
        assert_frame_equal(result, expected)

        result = df.reindex(columns=lrange(5), fill_value=0)
        expected = df.copy()
        expected[4] = 0
        assert_frame_equal(result, expected)

        result = df.reindex(columns=lrange(5), fill_value='foo')
        expected = df.copy()
        expected[4] = 'foo'
        assert_frame_equal(result, expected)

        # reindex_axis
        result = df.reindex_axis(lrange(15), fill_value=0., axis=0)
        expected = df.reindex(lrange(15)).fillna(0)
        assert_frame_equal(result, expected)

        result = df.reindex_axis(lrange(5), fill_value=0., axis=1)
        expected = df.reindex(columns=lrange(5)).fillna(0)
        assert_frame_equal(result, expected)

        # other dtypes
        df['foo'] = 'foo'
        result = df.reindex(lrange(15), fill_value=0)
        expected = df.reindex(lrange(15)).fillna(0)
        assert_frame_equal(result, expected)

    def test_reindex_dups(self):

        # GH4746, reindex on duplicate index error messages
        arr = np.random.randn(10)
        df = DataFrame(arr,index=[1,2,3,4,5,1,2,3,4,5])

        # set index is ok
        result = df.copy()
        result.index = list(range(len(df)))
        expected = DataFrame(arr,index=list(range(len(df))))
        assert_frame_equal(result,expected)

        # reindex fails
        self.assertRaises(ValueError, df.reindex, index=list(range(len(df))))

    def test_align(self):
        af, bf = self.frame.align(self.frame)
        self.assertIsNot(af._data, self.frame._data)

        af, bf = self.frame.align(self.frame, copy=False)
        self.assertIs(af._data, self.frame._data)

        # axis = 0
        other = self.frame.ix[:-5, :3]
        af, bf = self.frame.align(other, axis=0, fill_value=-1)
        self.assertTrue(bf.columns.equals(other.columns))
        # test fill value
        join_idx = self.frame.index.join(other.index)
        diff_a = self.frame.index.difference(join_idx)
        diff_b = other.index.difference(join_idx)
        diff_a_vals = af.reindex(diff_a).values
        diff_b_vals = bf.reindex(diff_b).values
        self.assertTrue((diff_a_vals == -1).all())

        af, bf = self.frame.align(other, join='right', axis=0)
        self.assertTrue(bf.columns.equals(other.columns))
        self.assertTrue(bf.index.equals(other.index))
        self.assertTrue(af.index.equals(other.index))

        # axis = 1
        other = self.frame.ix[:-5, :3].copy()
        af, bf = self.frame.align(other, axis=1)
        self.assertTrue(bf.columns.equals(self.frame.columns))
        self.assertTrue(bf.index.equals(other.index))

        # test fill value
        join_idx = self.frame.index.join(other.index)
        diff_a = self.frame.index.difference(join_idx)
        diff_b = other.index.difference(join_idx)
        diff_a_vals = af.reindex(diff_a).values
        diff_b_vals = bf.reindex(diff_b).values
        self.assertTrue((diff_a_vals == -1).all())

        af, bf = self.frame.align(other, join='inner', axis=1)
        self.assertTrue(bf.columns.equals(other.columns))

        af, bf = self.frame.align(other, join='inner', axis=1, method='pad')
        self.assertTrue(bf.columns.equals(other.columns))

        # test other non-float types
        af, bf = self.intframe.align(other, join='inner', axis=1, method='pad')
        self.assertTrue(bf.columns.equals(other.columns))

        af, bf = self.mixed_frame.align(self.mixed_frame,
                                        join='inner', axis=1, method='pad')
        self.assertTrue(bf.columns.equals(self.mixed_frame.columns))

        af, bf = self.frame.align(other.ix[:, 0], join='inner', axis=1,
                                  method=None, fill_value=None)
        self.assertTrue(bf.index.equals(Index([])))

        af, bf = self.frame.align(other.ix[:, 0], join='inner', axis=1,
                                  method=None, fill_value=0)
        self.assertTrue(bf.index.equals(Index([])))

        # mixed floats/ints
        af, bf = self.mixed_float.align(other.ix[:, 0], join='inner', axis=1,
                                        method=None, fill_value=0)
        self.assertTrue(bf.index.equals(Index([])))

        af, bf = self.mixed_int.align(other.ix[:, 0], join='inner', axis=1,
                                        method=None, fill_value=0)
        self.assertTrue(bf.index.equals(Index([])))

        # try to align dataframe to series along bad axis
        self.assertRaises(ValueError, self.frame.align, af.ix[0, :3],
                          join='inner', axis=2)

        # align dataframe to series with broadcast or not
        idx = self.frame.index
        s = Series(range(len(idx)), index=idx)

        left, right = self.frame.align(s, axis=0)
        tm.assert_index_equal(left.index, self.frame.index)
        tm.assert_index_equal(right.index, self.frame.index)
        self.assertTrue(isinstance(right, Series))

        left, right = self.frame.align(s, broadcast_axis=1)
        tm.assert_index_equal(left.index, self.frame.index)
        expected = {}
        for c in self.frame.columns:
            expected[c] = s
        expected = DataFrame(expected, index=self.frame.index,
                             columns=self.frame.columns)
        assert_frame_equal(right, expected)

        # GH 9558
        df = DataFrame({'a':[1,2,3], 'b':[4,5,6]})
        result = df[df['a'] == 2]
        expected = DataFrame([[2, 5]], index=[1], columns=['a', 'b'])
        assert_frame_equal(result, expected)

        result = df.where(df['a'] == 2, 0)
        expected = DataFrame({'a':[0, 2, 0], 'b':[0, 5, 0]})
        assert_frame_equal(result, expected)

    def _check_align(self, a, b, axis, fill_axis, how, method, limit=None):
        aa, ab = a.align(b, axis=axis, join=how, method=method, limit=limit,
                         fill_axis=fill_axis)

        join_index, join_columns = None, None

        ea, eb = a, b
        if axis is None or axis == 0:
            join_index = a.index.join(b.index, how=how)
            ea = ea.reindex(index=join_index)
            eb = eb.reindex(index=join_index)

        if axis is None or axis == 1:
            join_columns = a.columns.join(b.columns, how=how)
            ea = ea.reindex(columns=join_columns)
            eb = eb.reindex(columns=join_columns)

        ea = ea.fillna(axis=fill_axis, method=method, limit=limit)
        eb = eb.fillna(axis=fill_axis, method=method, limit=limit)

        assert_frame_equal(aa, ea)
        assert_frame_equal(ab, eb)

    def test_align_fill_method_inner(self):
        for meth in ['pad', 'bfill']:
            for ax in [0, 1, None]:
                for fax in [0, 1]:
                    self._check_align_fill('inner', meth, ax, fax)

    def test_align_fill_method_outer(self):
        for meth in ['pad', 'bfill']:
            for ax in [0, 1, None]:
                for fax in [0, 1]:
                    self._check_align_fill('outer', meth, ax, fax)

    def test_align_fill_method_left(self):
        for meth in ['pad', 'bfill']:
            for ax in [0, 1, None]:
                for fax in [0, 1]:
                    self._check_align_fill('left', meth, ax, fax)

    def test_align_fill_method_right(self):
        for meth in ['pad', 'bfill']:
            for ax in [0, 1, None]:
                for fax in [0, 1]:
                    self._check_align_fill('right', meth, ax, fax)

    def _check_align_fill(self, kind, meth, ax, fax):
        left = self.frame.ix[0:4, :10]
        right = self.frame.ix[2:, 6:]
        empty = self.frame.ix[:0, :0]

        self._check_align(left, right, axis=ax, fill_axis=fax,
                          how=kind, method=meth)
        self._check_align(left, right, axis=ax, fill_axis=fax,
                          how=kind, method=meth, limit=1)

        # empty left
        self._check_align(empty, right, axis=ax, fill_axis=fax,
                          how=kind, method=meth)
        self._check_align(empty, right, axis=ax, fill_axis=fax,
                          how=kind, method=meth, limit=1)

        # empty right
        self._check_align(left, empty, axis=ax, fill_axis=fax,
                          how=kind, method=meth)
        self._check_align(left, empty, axis=ax, fill_axis=fax,
                          how=kind, method=meth, limit=1)

        # both empty
        self._check_align(empty, empty, axis=ax, fill_axis=fax,
                          how=kind, method=meth)
        self._check_align(empty, empty, axis=ax, fill_axis=fax,
                          how=kind, method=meth, limit=1)

    def test_align_int_fill_bug(self):
        # GH #910
        X = np.arange(10*10, dtype='float64').reshape(10, 10)
        Y = np.ones((10, 1), dtype=int)

        df1 = DataFrame(X)
        df1['0.X'] = Y.squeeze()

        df2 = df1.astype(float)

        result = df1 - df1.mean()
        expected = df2 - df2.mean()
        assert_frame_equal(result, expected)

    def test_align_multiindex(self):
        # GH 10665
        # same test cases as test_align_multiindex in test_series.py

        midx = pd.MultiIndex.from_product([range(2), range(3), range(2)],
                                 names=('a', 'b', 'c'))
        idx = pd.Index(range(2), name='b')
        df1 = pd.DataFrame(np.arange(12,dtype='int64'), index=midx)
        df2 = pd.DataFrame(np.arange(2,dtype='int64'), index=idx)

        # these must be the same results (but flipped)
        res1l, res1r = df1.align(df2, join='left')
        res2l, res2r = df2.align(df1, join='right')

        expl = df1
        assert_frame_equal(expl, res1l)
        assert_frame_equal(expl, res2r)
        expr = pd.DataFrame([0, 0, 1, 1, np.nan, np.nan] * 2, index=midx)
        assert_frame_equal(expr, res1r)
        assert_frame_equal(expr, res2l)

        res1l, res1r = df1.align(df2, join='right')
        res2l, res2r = df2.align(df1, join='left')

        exp_idx = pd.MultiIndex.from_product([range(2), range(2), range(2)],
                                             names=('a', 'b', 'c'))
        expl = pd.DataFrame([0, 1, 2, 3, 6, 7, 8, 9], index=exp_idx)
        assert_frame_equal(expl, res1l)
        assert_frame_equal(expl, res2r)
        expr = pd.DataFrame([0, 0, 1, 1] * 2, index=exp_idx)
        assert_frame_equal(expr, res1r)
        assert_frame_equal(expr, res2l)

    def test_where(self):
        default_frame = DataFrame(np.random.randn(5, 3),
                                  columns=['A', 'B', 'C'])

        def _safe_add(df):
            # only add to the numeric items
            def is_ok(s):
                return issubclass(s.dtype.type, (np.integer,np.floating)) and s.dtype != 'uint8'
            return DataFrame(dict([ (c,s+1) if is_ok(s) else (c,s) for c, s in compat.iteritems(df) ]))

        def _check_get(df, cond, check_dtypes = True):
            other1 = _safe_add(df)
            rs = df.where(cond, other1)
            rs2 = df.where(cond.values, other1)
            for k, v in rs.iteritems():
                exp = Series(np.where(cond[k], df[k], other1[k]),index=v.index)
                assert_series_equal(v, exp, check_names=False)
            assert_frame_equal(rs, rs2)

            # dtypes
            if check_dtypes:
                self.assertTrue((rs.dtypes == df.dtypes).all() == True)

        # check getting
        for df in [ default_frame, self.mixed_frame, self.mixed_float, self.mixed_int ]:
            cond = df > 0
            _check_get(df, cond)


        # upcasting case (GH # 2794)
        df = DataFrame(dict([ (c,Series([1]*3,dtype=c)) for c in ['int64','int32','float32','float64'] ]))
        df.ix[1,:] = 0
        result = df.where(df>=0).get_dtype_counts()

        #### when we don't preserve boolean casts ####
        #expected = Series({ 'float32' : 1, 'float64' : 3 })

        expected = Series({ 'float32' : 1, 'float64' : 1, 'int32' : 1, 'int64' : 1 })
        assert_series_equal(result, expected)

        # aligning
        def _check_align(df, cond, other, check_dtypes = True):
            rs = df.where(cond, other)
            for i, k in enumerate(rs.columns):
                result = rs[k]
                d = df[k].values
                c = cond[k].reindex(df[k].index).fillna(False).values

                if np.isscalar(other):
                    o = other
                else:
                    if isinstance(other,np.ndarray):
                        o = Series(other[:,i],index=result.index).values
                    else:
                        o = other[k].values

                new_values = d if c.all() else np.where(c, d, o)
                expected = Series(new_values, index=result.index, name=k)

                # since we can't always have the correct numpy dtype
                # as numpy doesn't know how to downcast, don't check
                assert_series_equal(result, expected, check_dtype=False)

            # dtypes
            # can't check dtype when other is an ndarray

            if check_dtypes and not isinstance(other,np.ndarray):
                self.assertTrue((rs.dtypes == df.dtypes).all() == True)

        for df in [ self.mixed_frame, self.mixed_float, self.mixed_int ]:

            # other is a frame
            cond = (df > 0)[1:]
            _check_align(df, cond, _safe_add(df))

            # check other is ndarray
            cond = df > 0
            _check_align(df, cond, (_safe_add(df).values))

            # integers are upcast, so don't check the dtypes
            cond = df > 0
            check_dtypes = all([ not issubclass(s.type,np.integer) for s in df.dtypes ])
            _check_align(df, cond, np.nan, check_dtypes = check_dtypes)

        # invalid conditions
        df = default_frame
        err1 = (df + 1).values[0:2, :]
        self.assertRaises(ValueError, df.where, cond, err1)

        err2 = cond.ix[:2, :].values
        other1 = _safe_add(df)
        self.assertRaises(ValueError, df.where, err2, other1)

        self.assertRaises(ValueError, df.mask, True)
        self.assertRaises(ValueError, df.mask, 0)

        # where inplace
        def _check_set(df, cond, check_dtypes = True):
            dfi = df.copy()
            econd = cond.reindex_like(df).fillna(True)
            expected = dfi.mask(~econd)

            dfi.where(cond, np.nan, inplace=True)
            assert_frame_equal(dfi, expected)

            # dtypes (and confirm upcasts)x
            if check_dtypes:
                for k, v in compat.iteritems(df.dtypes):
                    if issubclass(v.type,np.integer) and not cond[k].all():
                        v = np.dtype('float64')
                    self.assertEqual(dfi[k].dtype, v)

        for df in [ default_frame, self.mixed_frame, self.mixed_float, self.mixed_int ]:

            cond = df > 0
            _check_set(df, cond)

            cond = df >= 0
            _check_set(df, cond)

            # aligining
            cond = (df >= 0)[1:]
            _check_set(df, cond)

        # GH 10218
        # test DataFrame.where with Series slicing
        df = DataFrame({'a': range(3), 'b': range(4, 7)})
        result = df.where(df['a'] == 1)
        expected = df[df['a'] == 1].reindex(df.index)
        assert_frame_equal(result, expected)

    def test_where_bug(self):

        # GH 2793

        df = DataFrame({'a': [1.0, 2.0, 3.0, 4.0], 'b': [4.0, 3.0, 2.0, 1.0]}, dtype = 'float64')
        expected = DataFrame({'a': [np.nan, np.nan, 3.0, 4.0], 'b': [4.0, 3.0, np.nan, np.nan]}, dtype = 'float64')
        result   = df.where(df > 2, np.nan)
        assert_frame_equal(result, expected)

        result = df.copy()
        result.where(result > 2, np.nan, inplace=True)
        assert_frame_equal(result, expected)

        # mixed
        for dtype in ['int16','int8','int32','int64']:
            df = DataFrame({'a': np.array([1, 2, 3, 4],dtype=dtype), 'b': np.array([4.0, 3.0, 2.0, 1.0], dtype = 'float64') })
            expected = DataFrame({'a': [np.nan, np.nan, 3.0, 4.0], 'b': [4.0, 3.0, np.nan, np.nan]}, dtype = 'float64')
            result   = df.where(df > 2, np.nan)
            assert_frame_equal(result, expected)

            result = df.copy()
            result.where(result > 2, np.nan, inplace=True)
            assert_frame_equal(result, expected)

        # transpositional issue
        # GH7506
        a = DataFrame({ 0 : [1,2], 1 : [3,4], 2 : [5,6]})
        b = DataFrame({ 0 : [np.nan,8], 1:[9,np.nan], 2:[np.nan,np.nan]})
        do_not_replace = b.isnull() | (a > b)

        expected = a.copy()
        expected[~do_not_replace] = b

        result = a.where(do_not_replace,b)
        assert_frame_equal(result,expected)

        a = DataFrame({ 0 : [4,6], 1 : [1,0]})
        b = DataFrame({ 0 : [np.nan,3],1:[3,np.nan]})
        do_not_replace = b.isnull() | (a > b)

        expected = a.copy()
        expected[~do_not_replace] = b

        result = a.where(do_not_replace,b)
        assert_frame_equal(result,expected)

    def test_where_datetime(self):

        # GH 3311
        df = DataFrame(dict(A = date_range('20130102',periods=5),
                            B = date_range('20130104',periods=5),
                            C = np.random.randn(5)))

        stamp = datetime(2013,1,3)
        result = df[df>stamp]
        expected = df.copy()
        expected.loc[[0,1],'A'] = np.nan
        assert_frame_equal(result,expected)

    def test_where_none(self):
        # GH 4667
        # setting with None changes dtype
        df = DataFrame({'series': Series(range(10))}).astype(float)
        df[df > 7] = None
        expected = DataFrame({'series': Series([0,1,2,3,4,5,6,7,np.nan,np.nan]) })
        assert_frame_equal(df, expected)

        # GH 7656
        df = DataFrame([{'A': 1, 'B': np.nan, 'C': 'Test'}, {'A': np.nan, 'B': 'Test', 'C': np.nan}])
        expected = df.where(~isnull(df), None)
        with tm.assertRaisesRegexp(TypeError, 'boolean setting on mixed-type'):
            df.where(~isnull(df), None, inplace=True)

    def test_where_align(self):

        def create():
            df = DataFrame(np.random.randn(10,3))
            df.iloc[3:5,0] = np.nan
            df.iloc[4:6,1] = np.nan
            df.iloc[5:8,2] = np.nan
            return df

        # series
        df = create()
        expected = df.fillna(df.mean())
        result = df.where(pd.notnull(df),df.mean(),axis='columns')
        assert_frame_equal(result, expected)

        df.where(pd.notnull(df),df.mean(),inplace=True,axis='columns')
        assert_frame_equal(df, expected)

        df = create().fillna(0)
        expected = df.apply(lambda x, y: x.where(x>0,y), y=df[0])
        result = df.where(df>0,df[0],axis='index')
        assert_frame_equal(result, expected)
        result = df.where(df>0,df[0],axis='rows')
        assert_frame_equal(result, expected)

        # frame
        df = create()
        expected = df.fillna(1)
        result = df.where(pd.notnull(df),DataFrame(1,index=df.index,columns=df.columns))
        assert_frame_equal(result, expected)

    def test_where_complex(self):
        # GH 6345
        expected = DataFrame([[1+1j, 2], [np.nan, 4+1j]], columns=['a', 'b'])
        df = DataFrame([[1+1j, 2], [5+1j, 4+1j]], columns=['a', 'b'])
        df[df.abs() >= 5] = np.nan
        assert_frame_equal(df,expected)

    def test_where_axis(self):
        # GH 9736
        df = DataFrame(np.random.randn(2, 2))
        mask = DataFrame([[False, False], [False, False]])
        s = Series([0, 1])

        expected = DataFrame([[0, 0], [1, 1]], dtype='float64')
        result = df.where(mask, s, axis='index')
        assert_frame_equal(result, expected)

        result = df.copy()
        result.where(mask, s, axis='index', inplace=True)
        assert_frame_equal(result, expected)

        expected = DataFrame([[0, 1], [0, 1]], dtype='float64')
        result = df.where(mask, s, axis='columns')
        assert_frame_equal(result, expected)

        result = df.copy()
        result.where(mask, s, axis='columns', inplace=True)
        assert_frame_equal(result, expected)

        # Upcast needed
        df = DataFrame([[1, 2], [3, 4]], dtype='int64')
        mask = DataFrame([[False, False], [False, False]])
        s = Series([0, np.nan])

        expected = DataFrame([[0, 0], [np.nan, np.nan]], dtype='float64')
        result = df.where(mask, s, axis='index')
        assert_frame_equal(result, expected)

        result = df.copy()
        result.where(mask, s, axis='index', inplace=True)
        assert_frame_equal(result, expected)

        expected = DataFrame([[0, np.nan], [0, np.nan]], dtype='float64')
        result = df.where(mask, s, axis='columns')
        assert_frame_equal(result, expected)

        expected = DataFrame({0 : np.array([0, 0], dtype='int64'),
                              1 : np.array([np.nan, np.nan], dtype='float64')})
        result = df.copy()
        result.where(mask, s, axis='columns', inplace=True)
        assert_frame_equal(result, expected)

        # Multiple dtypes (=> multiple Blocks)
        df = pd.concat([DataFrame(np.random.randn(10, 2)),
                     DataFrame(np.random.randint(0, 10, size=(10, 2)))],
                     ignore_index=True, axis=1)
        mask = DataFrame(False, columns=df.columns, index=df.index)
        s1 = Series(1, index=df.columns)
        s2 = Series(2, index=df.index)

        result = df.where(mask, s1, axis='columns')
        expected = DataFrame(1.0, columns=df.columns, index=df.index)
        expected[2] = expected[2].astype(int)
        expected[3] = expected[3].astype(int)
        assert_frame_equal(result, expected)

        result = df.copy()
        result.where(mask, s1, axis='columns', inplace=True)
        assert_frame_equal(result, expected)

        result = df.where(mask, s2, axis='index')
        expected = DataFrame(2.0, columns=df.columns, index=df.index)
        expected[2] = expected[2].astype(int)
        expected[3] = expected[3].astype(int)
        assert_frame_equal(result, expected)

        result = df.copy()
        result.where(mask, s2, axis='index', inplace=True)
        assert_frame_equal(result, expected)

        # DataFrame vs DataFrame
        d1 = df.copy().drop(1, axis=0)
        expected = df.copy()
        expected.loc[1, :] = np.nan

        result = df.where(mask, d1)
        assert_frame_equal(result, expected)
        result = df.where(mask, d1, axis='index')
        assert_frame_equal(result, expected)
        result = df.copy()
        result.where(mask, d1, inplace=True)
        assert_frame_equal(result, expected)
        result = df.copy()
        result.where(mask, d1, inplace=True, axis='index')
        assert_frame_equal(result, expected)

        d2 = df.copy().drop(1, axis=1)
        expected = df.copy()
        expected.loc[:, 1] = np.nan

        result = df.where(mask, d2)
        assert_frame_equal(result, expected)
        result = df.where(mask, d2, axis='columns')
        assert_frame_equal(result, expected)
        result = df.copy()
        result.where(mask, d2, inplace=True)
        assert_frame_equal(result, expected)
        result = df.copy()
        result.where(mask, d2, inplace=True, axis='columns')
        assert_frame_equal(result, expected)

    def test_mask(self):
        df = DataFrame(np.random.randn(5, 3))
        cond = df > 0

        rs = df.where(cond, np.nan)
        assert_frame_equal(rs, df.mask(df <= 0))
        assert_frame_equal(rs, df.mask(~cond))

        other = DataFrame(np.random.randn(5, 3))
        rs = df.where(cond, other)
        assert_frame_equal(rs, df.mask(df <= 0, other))
        assert_frame_equal(rs, df.mask(~cond, other))

    def test_mask_inplace(self):
        # GH8801
        df = DataFrame(np.random.randn(5, 3))
        cond = df > 0

        rdf = df.copy()

        rdf.where(cond, inplace=True)
        assert_frame_equal(rdf, df.where(cond))
        assert_frame_equal(rdf, df.mask(~cond))

        rdf = df.copy()
        rdf.where(cond, -df, inplace=True)
        assert_frame_equal(rdf, df.where(cond, -df))
        assert_frame_equal(rdf, df.mask(~cond, -df))

    def test_mask_edge_case_1xN_frame(self):
        # GH4071
        df = DataFrame([[1, 2]])
        res = df.mask(DataFrame([[True, False]]))
        expec = DataFrame([[nan, 2]])
        assert_frame_equal(res, expec)

    #----------------------------------------------------------------------
    # Transposing

    def test_transpose(self):
        frame = self.frame
        dft = frame.T
        for idx, series in compat.iteritems(dft):
            for col, value in compat.iteritems(series):
                if np.isnan(value):
                    self.assertTrue(np.isnan(frame[col][idx]))
                else:
                    self.assertEqual(value, frame[col][idx])

        # mixed type
        index, data = tm.getMixedTypeDict()
        mixed = DataFrame(data, index=index)

        mixed_T = mixed.T
        for col, s in compat.iteritems(mixed_T):
            self.assertEqual(s.dtype, np.object_)

    def test_transpose_get_view(self):
        dft = self.frame.T
        dft.values[:, 5:10] = 5

        self.assertTrue((self.frame.values[5:10] == 5).all())

    #----------------------------------------------------------------------
    # Renaming

    def test_rename(self):
        mapping = {
            'A': 'a',
            'B': 'b',
            'C': 'c',
            'D': 'd'
        }

        renamed = self.frame.rename(columns=mapping)
        renamed2 = self.frame.rename(columns=str.lower)

        assert_frame_equal(renamed, renamed2)
        assert_frame_equal(renamed2.rename(columns=str.upper),
                           self.frame, check_names=False)

        # index
        data = {
            'A': {'foo': 0, 'bar': 1}
        }

        # gets sorted alphabetical
        df = DataFrame(data)
        renamed = df.rename(index={'foo': 'bar', 'bar': 'foo'})
        self.assert_numpy_array_equal(renamed.index, ['foo', 'bar'])

        renamed = df.rename(index=str.upper)
        self.assert_numpy_array_equal(renamed.index, ['BAR', 'FOO'])

        # have to pass something
        self.assertRaises(TypeError, self.frame.rename)

        # partial columns
        renamed = self.frame.rename(columns={'C': 'foo', 'D': 'bar'})
        self.assert_numpy_array_equal(renamed.columns, ['A', 'B', 'foo', 'bar'])

        # other axis
        renamed = self.frame.T.rename(index={'C': 'foo', 'D': 'bar'})
        self.assert_numpy_array_equal(renamed.index, ['A', 'B', 'foo', 'bar'])

        # index with name
        index = Index(['foo', 'bar'], name='name')
        renamer = DataFrame(data, index=index)
        renamed = renamer.rename(index={'foo': 'bar', 'bar': 'foo'})
        self.assert_numpy_array_equal(renamed.index, ['bar', 'foo'])
        self.assertEqual(renamed.index.name, renamer.index.name)

        # MultiIndex
        tuples_index = [('foo1', 'bar1'), ('foo2', 'bar2')]
        tuples_columns = [('fizz1', 'buzz1'), ('fizz2', 'buzz2')]
        index = MultiIndex.from_tuples(tuples_index, names=['foo', 'bar'])
        columns = MultiIndex.from_tuples(tuples_columns, names=['fizz', 'buzz'])
        renamer = DataFrame([(0,0),(1,1)], index=index, columns=columns)
        renamed = renamer.rename(index={'foo1': 'foo3', 'bar2': 'bar3'},
                                 columns={'fizz1': 'fizz3', 'buzz2': 'buzz3'})
        new_index = MultiIndex.from_tuples([('foo3', 'bar1'), ('foo2', 'bar3')])
        new_columns = MultiIndex.from_tuples([('fizz3', 'buzz1'), ('fizz2', 'buzz3')])
        self.assert_numpy_array_equal(renamed.index, new_index)
        self.assert_numpy_array_equal(renamed.columns, new_columns)
        self.assertEqual(renamed.index.names, renamer.index.names)
        self.assertEqual(renamed.columns.names, renamer.columns.names)

    def test_rename_nocopy(self):
        renamed = self.frame.rename(columns={'C': 'foo'}, copy=False)
        renamed['foo'] = 1.
        self.assertTrue((self.frame['C'] == 1.).all())

    def test_rename_inplace(self):
        self.frame.rename(columns={'C': 'foo'})
        self.assertIn('C', self.frame)
        self.assertNotIn('foo', self.frame)

        c_id = id(self.frame['C'])
        frame = self.frame.copy()
        frame.rename(columns={'C': 'foo'}, inplace=True)

        self.assertNotIn('C', frame)
        self.assertIn('foo', frame)
        self.assertNotEqual(id(frame['foo']), c_id)

    def test_rename_bug(self):
        # GH 5344
        # rename set ref_locs, and set_index was not resetting
        df = DataFrame({ 0 : ['foo','bar'], 1 : ['bah','bas'], 2 : [1,2]})
        df = df.rename(columns={0 : 'a'})
        df = df.rename(columns={1 : 'b'})
        df = df.set_index(['a','b'])
        df.columns = ['2001-01-01']
        expected = DataFrame([[1],[2]],index=MultiIndex.from_tuples([('foo','bah'),('bar','bas')],
                                                                    names=['a','b']),
                             columns=['2001-01-01'])
        assert_frame_equal(df,expected)

    #----------------------------------------------------------------------
    # Time series related
    def test_diff(self):
        the_diff = self.tsframe.diff(1)

        assert_series_equal(the_diff['A'],
                            self.tsframe['A'] - self.tsframe['A'].shift(1))

        # int dtype
        a = 10000000000000000
        b = a + 1
        s = Series([a, b])

        rs = DataFrame({'s': s}).diff()
        self.assertEqual(rs.s[1], 1)

        # mixed numeric
        tf = self.tsframe.astype('float32')
        the_diff = tf.diff(1)
        assert_series_equal(the_diff['A'],
                            tf['A'] - tf['A'].shift(1))

        # issue 10907
        df = pd.DataFrame({'y': pd.Series([2]), 'z': pd.Series([3])})
        df.insert(0, 'x', 1)
        result = df.diff(axis=1)
        expected = pd.DataFrame({'x':np.nan, 'y':pd.Series(1), 'z':pd.Series(1)}).astype('float64')
        assert_frame_equal(result, expected)


    def test_diff_timedelta(self):
        # GH 4533
        df = DataFrame(dict(time=[Timestamp('20130101 9:01'),
                                  Timestamp('20130101 9:02')],
                            value=[1.0,2.0]))

        res = df.diff()
        exp = DataFrame([[pd.NaT, np.nan],
                         [Timedelta('00:01:00'), 1]],
                        columns=['time', 'value'])
        assert_frame_equal(res, exp)

    def test_diff_mixed_dtype(self):
        df = DataFrame(np.random.randn(5, 3))
        df['A'] = np.array([1, 2, 3, 4, 5], dtype=object)

        result = df.diff()
        self.assertEqual(result[0].dtype, np.float64)

    def test_diff_neg_n(self):
        rs = self.tsframe.diff(-1)
        xp = self.tsframe - self.tsframe.shift(-1)
        assert_frame_equal(rs, xp)

    def test_diff_float_n(self):
        rs = self.tsframe.diff(1.)
        xp = self.tsframe.diff(1)
        assert_frame_equal(rs, xp)

    def test_diff_axis(self):
        # GH 9727
        df = DataFrame([[1., 2.], [3., 4.]])
        assert_frame_equal(df.diff(axis=1), DataFrame([[np.nan, 1.], [np.nan, 1.]]))
        assert_frame_equal(df.diff(axis=0), DataFrame([[np.nan, np.nan], [2., 2.]]))

    def test_pct_change(self):
        rs = self.tsframe.pct_change(fill_method=None)
        assert_frame_equal(rs, self.tsframe / self.tsframe.shift(1) - 1)

        rs = self.tsframe.pct_change(2)
        filled = self.tsframe.fillna(method='pad')
        assert_frame_equal(rs, filled / filled.shift(2) - 1)

        rs = self.tsframe.pct_change(fill_method='bfill', limit=1)
        filled = self.tsframe.fillna(method='bfill', limit=1)
        assert_frame_equal(rs, filled / filled.shift(1) - 1)

        rs = self.tsframe.pct_change(freq='5D')
        filled = self.tsframe.fillna(method='pad')
        assert_frame_equal(rs, filled / filled.shift(freq='5D') - 1)

    def test_pct_change_shift_over_nas(self):
        s = Series([1., 1.5, np.nan, 2.5, 3.])

        df = DataFrame({'a': s, 'b': s})

        chg = df.pct_change()
        expected = Series([np.nan, 0.5, np.nan, 2.5 / 1.5 - 1, .2])
        edf = DataFrame({'a': expected, 'b': expected})
        assert_frame_equal(chg, edf)

    def test_shift(self):
        # naive shift
        shiftedFrame = self.tsframe.shift(5)
        self.assertTrue(shiftedFrame.index.equals(self.tsframe.index))

        shiftedSeries = self.tsframe['A'].shift(5)
        assert_series_equal(shiftedFrame['A'], shiftedSeries)

        shiftedFrame = self.tsframe.shift(-5)
        self.assertTrue(shiftedFrame.index.equals(self.tsframe.index))

        shiftedSeries = self.tsframe['A'].shift(-5)
        assert_series_equal(shiftedFrame['A'], shiftedSeries)

        # shift by 0
        unshifted = self.tsframe.shift(0)
        assert_frame_equal(unshifted, self.tsframe)

        # shift by DateOffset
        shiftedFrame = self.tsframe.shift(5, freq=datetools.BDay())
        self.assertEqual(len(shiftedFrame), len(self.tsframe))

        shiftedFrame2 = self.tsframe.shift(5, freq='B')
        assert_frame_equal(shiftedFrame, shiftedFrame2)

        d = self.tsframe.index[0]
        shifted_d = d + datetools.BDay(5)
        assert_series_equal(self.tsframe.xs(d),
                            shiftedFrame.xs(shifted_d), check_names=False)

        # shift int frame
        int_shifted = self.intframe.shift(1)

        # Shifting with PeriodIndex
        ps = tm.makePeriodFrame()
        shifted = ps.shift(1)
        unshifted = shifted.shift(-1)
        self.assertTrue(shifted.index.equals(ps.index))

        tm.assert_dict_equal(unshifted.ix[:, 0].valid(), ps.ix[:, 0],
                             compare_keys=False)

        shifted2 = ps.shift(1, 'B')
        shifted3 = ps.shift(1, datetools.bday)
        assert_frame_equal(shifted2, shifted3)
        assert_frame_equal(ps, shifted2.shift(-1, 'B'))

        assertRaisesRegexp(ValueError, 'does not match PeriodIndex freq',
                           ps.shift, freq='D')


        # shift other axis
        # GH 6371
        df = DataFrame(np.random.rand(10,5))
        expected = pd.concat([DataFrame(np.nan,index=df.index,columns=[0]),df.iloc[:,0:-1]],ignore_index=True,axis=1)
        result = df.shift(1,axis=1)
        assert_frame_equal(result,expected)

        # shift named axis
        df = DataFrame(np.random.rand(10,5))
        expected = pd.concat([DataFrame(np.nan,index=df.index,columns=[0]),df.iloc[:,0:-1]],ignore_index=True,axis=1)
        result = df.shift(1,axis='columns')
        assert_frame_equal(result,expected)

    def test_shift_bool(self):
        df = DataFrame({'high': [True, False],
                        'low': [False, False]})
        rs = df.shift(1)
        xp = DataFrame(np.array([[np.nan, np.nan],
                                 [True, False]], dtype=object),
                       columns=['high', 'low'])
        assert_frame_equal(rs, xp)

    def test_shift_categorical(self):
        # GH 9416
        s1 = pd.Series(['a', 'b', 'c'], dtype='category')
        s2 = pd.Series(['A', 'B', 'C'], dtype='category')
        df = DataFrame({'one': s1, 'two': s2})
        rs = df.shift(1)
        xp = DataFrame({'one': s1.shift(1), 'two': s2.shift(1)})
        assert_frame_equal(rs, xp)

    def test_shift_empty(self):
        # Regression test for #8019
        df = DataFrame({'foo': []})
        rs = df.shift(-1)

        assert_frame_equal(df, rs)

    def test_tshift(self):
        # PeriodIndex
        ps = tm.makePeriodFrame()
        shifted = ps.tshift(1)
        unshifted = shifted.tshift(-1)

        assert_frame_equal(unshifted, ps)

        shifted2 = ps.tshift(freq='B')
        assert_frame_equal(shifted, shifted2)

        shifted3 = ps.tshift(freq=datetools.bday)
        assert_frame_equal(shifted, shifted3)

        assertRaisesRegexp(ValueError, 'does not match', ps.tshift, freq='M')

        # DatetimeIndex
        shifted = self.tsframe.tshift(1)
        unshifted = shifted.tshift(-1)

        assert_frame_equal(self.tsframe, unshifted)

        shifted2 = self.tsframe.tshift(freq=self.tsframe.index.freq)
        assert_frame_equal(shifted, shifted2)

        inferred_ts = DataFrame(self.tsframe.values,
                                Index(np.asarray(self.tsframe.index)),
                                columns=self.tsframe.columns)
        shifted = inferred_ts.tshift(1)
        unshifted = shifted.tshift(-1)
        assert_frame_equal(shifted, self.tsframe.tshift(1))
        assert_frame_equal(unshifted, inferred_ts)

        no_freq = self.tsframe.ix[[0, 5, 7], :]
        self.assertRaises(ValueError, no_freq.tshift)

    def test_apply(self):
        # ufunc
        applied = self.frame.apply(np.sqrt)
        assert_series_equal(np.sqrt(self.frame['A']), applied['A'])

        # aggregator
        applied = self.frame.apply(np.mean)
        self.assertEqual(applied['A'], np.mean(self.frame['A']))

        d = self.frame.index[0]
        applied = self.frame.apply(np.mean, axis=1)
        self.assertEqual(applied[d], np.mean(self.frame.xs(d)))
        self.assertIs(applied.index, self.frame.index)  # want this

        # invalid axis
        df = DataFrame(
            [[1, 2, 3], [4, 5, 6], [7, 8, 9]], index=['a', 'a', 'c'])
        self.assertRaises(ValueError, df.apply, lambda x: x, 2)

        # GH9573
        df = DataFrame({'c0':['A','A','B','B'], 'c1':['C','C','D','D']})
        df = df.apply(lambda ts: ts.astype('category'))
        self.assertEqual(df.shape, (4, 2))
        self.assertTrue(isinstance(df['c0'].dtype, com.CategoricalDtype))
        self.assertTrue(isinstance(df['c1'].dtype, com.CategoricalDtype))

    def test_apply_mixed_datetimelike(self):
        # mixed datetimelike
        # GH 7778
        df = DataFrame({ 'A' : date_range('20130101',periods=3), 'B' : pd.to_timedelta(np.arange(3),unit='s') })
        result = df.apply(lambda x: x, axis=1)
        assert_frame_equal(result, df)

    def test_apply_empty(self):
        # empty
        applied = self.empty.apply(np.sqrt)
        self.assertTrue(applied.empty)

        applied = self.empty.apply(np.mean)
        self.assertTrue(applied.empty)

        no_rows = self.frame[:0]
        result = no_rows.apply(lambda x: x.mean())
        expected = Series(np.nan, index=self.frame.columns)
        assert_series_equal(result, expected)

        no_cols = self.frame.ix[:, []]
        result = no_cols.apply(lambda x: x.mean(), axis=1)
        expected = Series(np.nan, index=self.frame.index)
        assert_series_equal(result, expected)

        # 2476
        xp = DataFrame(index=['a'])
        rs = xp.apply(lambda x: x['a'], axis=1)
        assert_frame_equal(xp, rs)

        # reduce with an empty DataFrame
        x = []
        result = self.empty.apply(x.append, axis=1, reduce=False)
        assert_frame_equal(result, self.empty)
        result = self.empty.apply(x.append, axis=1, reduce=True)
        assert_series_equal(result, Series([], index=pd.Index([], dtype=object)))

        empty_with_cols = DataFrame(columns=['a', 'b', 'c'])
        result = empty_with_cols.apply(x.append, axis=1, reduce=False)
        assert_frame_equal(result, empty_with_cols)
        result = empty_with_cols.apply(x.append, axis=1, reduce=True)
        assert_series_equal(result, Series([], index=pd.Index([], dtype=object)))

        # Ensure that x.append hasn't been called
        self.assertEqual(x, [])

    def test_apply_standard_nonunique(self):
        df = DataFrame(
            [[1, 2, 3], [4, 5, 6], [7, 8, 9]], index=['a', 'a', 'c'])
        rs = df.apply(lambda s: s[0], axis=1)
        xp = Series([1, 4, 7], ['a', 'a', 'c'])
        assert_series_equal(rs, xp)

        rs = df.T.apply(lambda s: s[0], axis=0)
        assert_series_equal(rs, xp)

    def test_apply_broadcast(self):
        broadcasted = self.frame.apply(np.mean, broadcast=True)
        agged = self.frame.apply(np.mean)

        for col, ts in compat.iteritems(broadcasted):
            self.assertTrue((ts == agged[col]).all())

        broadcasted = self.frame.apply(np.mean, axis=1, broadcast=True)
        agged = self.frame.apply(np.mean, axis=1)
        for idx in broadcasted.index:
            self.assertTrue((broadcasted.xs(idx) == agged[idx]).all())

    def test_apply_raw(self):
        result0 = self.frame.apply(np.mean, raw=True)
        result1 = self.frame.apply(np.mean, axis=1, raw=True)

        expected0 = self.frame.apply(lambda x: x.values.mean())
        expected1 = self.frame.apply(lambda x: x.values.mean(), axis=1)

        assert_series_equal(result0, expected0)
        assert_series_equal(result1, expected1)

        # no reduction
        result = self.frame.apply(lambda x: x * 2, raw=True)
        expected = self.frame * 2
        assert_frame_equal(result, expected)

    def test_apply_axis1(self):
        d = self.frame.index[0]
        tapplied = self.frame.apply(np.mean, axis=1)
        self.assertEqual(tapplied[d], np.mean(self.frame.xs(d)))

    def test_apply_ignore_failures(self):
        result = self.mixed_frame._apply_standard(np.mean, 0,
                                                  ignore_failures=True)
        expected = self.mixed_frame._get_numeric_data().apply(np.mean)
        assert_series_equal(result, expected)

    def test_apply_mixed_dtype_corner(self):
        df = DataFrame({'A': ['foo'],
                        'B': [1.]})
        result = df[:0].apply(np.mean, axis=1)
        # the result here is actually kind of ambiguous, should it be a Series
        # or a DataFrame?
        expected = Series(np.nan, index=pd.Index([], dtype='int64'))
        assert_series_equal(result, expected)

        df = DataFrame({'A': ['foo'],
                        'B': [1.]})
        result = df.apply(lambda x: x['A'], axis=1)
        expected = Series(['foo'],index=[0])
        assert_series_equal(result, expected)

        result = df.apply(lambda x: x['B'], axis=1)
        expected = Series([1.],index=[0])
        assert_series_equal(result, expected)

    def test_apply_empty_infer_type(self):
        no_cols = DataFrame(index=['a', 'b', 'c'])
        no_index = DataFrame(columns=['a', 'b', 'c'])

        def _check(df, f):
            test_res = f(np.array([], dtype='f8'))
            is_reduction = not isinstance(test_res, np.ndarray)

            def _checkit(axis=0, raw=False):
                res = df.apply(f, axis=axis, raw=raw)
                if is_reduction:
                    agg_axis = df._get_agg_axis(axis)
                    tm.assertIsInstance(res, Series)
                    self.assertIs(res.index, agg_axis)
                else:
                    tm.assertIsInstance(res, DataFrame)

            _checkit()
            _checkit(axis=1)
            _checkit(raw=True)
            _checkit(axis=0, raw=True)

        _check(no_cols, lambda x: x)
        _check(no_cols, lambda x: x.mean())
        _check(no_index, lambda x: x)
        _check(no_index, lambda x: x.mean())

        result = no_cols.apply(lambda x: x.mean(), broadcast=True)
        tm.assertIsInstance(result, DataFrame)

    def test_apply_with_args_kwds(self):
        def add_some(x, howmuch=0):
            return x + howmuch

        def agg_and_add(x, howmuch=0):
            return x.mean() + howmuch

        def subtract_and_divide(x, sub, divide=1):
            return (x - sub) / divide

        result = self.frame.apply(add_some, howmuch=2)
        exp = self.frame.apply(lambda x: x + 2)
        assert_frame_equal(result, exp)

        result = self.frame.apply(agg_and_add, howmuch=2)
        exp = self.frame.apply(lambda x: x.mean() + 2)
        assert_series_equal(result, exp)

        res = self.frame.apply(subtract_and_divide, args=(2,), divide=2)
        exp = self.frame.apply(lambda x: (x - 2.) / 2.)
        assert_frame_equal(res, exp)

    def test_apply_yield_list(self):
        result = self.frame.apply(list)
        assert_frame_equal(result, self.frame)

    def test_apply_reduce_Series(self):
        self.frame.ix[::2, 'A'] = np.nan
        expected = self.frame.mean(1)
        result = self.frame.apply(np.mean, axis=1)
        assert_series_equal(result, expected)

    def test_apply_differently_indexed(self):
        df = DataFrame(np.random.randn(20, 10))

        result0 = df.apply(Series.describe, axis=0)
        expected0 = DataFrame(dict((i, v.describe())
                                   for i, v in compat.iteritems(df)),
                              columns=df.columns)
        assert_frame_equal(result0, expected0)

        result1 = df.apply(Series.describe, axis=1)
        expected1 = DataFrame(dict((i, v.describe())
                                   for i, v in compat.iteritems(df.T)),
                              columns=df.index).T
        assert_frame_equal(result1, expected1)

    def test_apply_modify_traceback(self):
        data = DataFrame({'A': ['foo', 'foo', 'foo', 'foo',
                                'bar', 'bar', 'bar', 'bar',
                                'foo', 'foo', 'foo'],
                          'B': ['one', 'one', 'one', 'two',
                                'one', 'one', 'one', 'two',
                                'two', 'two', 'one'],
                          'C': ['dull', 'dull', 'shiny', 'dull',
                                'dull', 'shiny', 'shiny', 'dull',
                                'shiny', 'shiny', 'shiny'],
                          'D': np.random.randn(11),
                          'E': np.random.randn(11),
                          'F': np.random.randn(11)})

        data.loc[4,'C'] = np.nan

        def transform(row):
            if row['C'].startswith('shin') and row['A'] == 'foo':
                row['D'] = 7
            return row

        def transform2(row):
            if (notnull(row['C']) and row['C'].startswith('shin')
                    and row['A'] == 'foo'):
                row['D'] = 7
            return row

        try:
            transformed = data.apply(transform, axis=1)
        except AttributeError as e:
            self.assertEqual(len(e.args), 2)
            self.assertEqual(e.args[1], 'occurred at index 4')
            self.assertEqual(e.args[0], "'float' object has no attribute 'startswith'")

    def test_apply_bug(self):

        # GH 6125
        import datetime
        positions = pd.DataFrame([[1, 'ABC0', 50], [1, 'YUM0', 20],
                                  [1, 'DEF0', 20], [2, 'ABC1', 50],
                                  [2, 'YUM1', 20], [2, 'DEF1', 20]],
                                 columns=['a', 'market', 'position'])
        def f(r):
            return r['market']
        expected = positions.apply(f, axis=1)

        positions = DataFrame([[datetime.datetime(2013, 1, 1), 'ABC0', 50],
                               [datetime.datetime(2013, 1, 2), 'YUM0', 20],
                               [datetime.datetime(2013, 1, 3), 'DEF0', 20],
                               [datetime.datetime(2013, 1, 4), 'ABC1', 50],
                               [datetime.datetime(2013, 1, 5), 'YUM1', 20],
                               [datetime.datetime(2013, 1, 6), 'DEF1', 20]],
                                                columns=['a', 'market', 'position'])
        result = positions.apply(f, axis=1)
        assert_series_equal(result,expected)

    def test_swapaxes(self):
        df = DataFrame(np.random.randn(10, 5))
        assert_frame_equal(df.T, df.swapaxes(0, 1))
        assert_frame_equal(df.T, df.swapaxes(1, 0))
        assert_frame_equal(df, df.swapaxes(0, 0))
        self.assertRaises(ValueError, df.swapaxes, 2, 5)

    def test_apply_convert_objects(self):
        data = DataFrame({'A': ['foo', 'foo', 'foo', 'foo',
                                'bar', 'bar', 'bar', 'bar',
                                'foo', 'foo', 'foo'],
                          'B': ['one', 'one', 'one', 'two',
                                'one', 'one', 'one', 'two',
                                'two', 'two', 'one'],
                          'C': ['dull', 'dull', 'shiny', 'dull',
                                'dull', 'shiny', 'shiny', 'dull',
                                'shiny', 'shiny', 'shiny'],
                          'D': np.random.randn(11),
                          'E': np.random.randn(11),
                          'F': np.random.randn(11)})

        result = data.apply(lambda x: x, axis=1)
        assert_frame_equal(result._convert(datetime=True), data)

    def test_apply_attach_name(self):
        result = self.frame.apply(lambda x: x.name)
        expected = Series(self.frame.columns, index=self.frame.columns)
        assert_series_equal(result, expected)

        result = self.frame.apply(lambda x: x.name, axis=1)
        expected = Series(self.frame.index, index=self.frame.index)
        assert_series_equal(result, expected)

        # non-reductions
        result = self.frame.apply(lambda x: np.repeat(x.name, len(x)))
        expected = DataFrame(np.tile(self.frame.columns,
                                     (len(self.frame.index), 1)),
                             index=self.frame.index,
                             columns=self.frame.columns)
        assert_frame_equal(result, expected)

        result = self.frame.apply(lambda x: np.repeat(x.name, len(x)),
                                  axis=1)
        expected = DataFrame(np.tile(self.frame.index,
                                     (len(self.frame.columns), 1)).T,
                             index=self.frame.index,
                             columns=self.frame.columns)
        assert_frame_equal(result, expected)

    def test_apply_multi_index(self):
        s = DataFrame([[1,2], [3,4], [5,6]])
        s.index = MultiIndex.from_arrays([['a','a','b'], ['c','d','d']])
        s.columns = ['col1','col2']
        res = s.apply(lambda x: Series({'min': min(x), 'max': max(x)}), 1)
        tm.assertIsInstance(res.index, MultiIndex)

    def test_apply_dict(self):

        # GH 8735
        A = DataFrame([['foo', 'bar'], ['spam', 'eggs']])
        A_dicts = pd.Series([dict([(0, 'foo'), (1, 'spam')]),
                             dict([(0, 'bar'), (1, 'eggs')])])
        B = DataFrame([[0, 1], [2, 3]])
        B_dicts = pd.Series([dict([(0, 0), (1, 2)]), dict([(0, 1), (1, 3)])])
        fn = lambda x: x.to_dict()

        for df, dicts in [(A, A_dicts), (B, B_dicts)]:
            reduce_true = df.apply(fn, reduce=True)
            reduce_false = df.apply(fn, reduce=False)
            reduce_none = df.apply(fn, reduce=None)

            assert_series_equal(reduce_true, dicts)
            assert_frame_equal(reduce_false, df)
            assert_series_equal(reduce_none, dicts)

    def test_applymap(self):
        applied = self.frame.applymap(lambda x: x * 2)
        assert_frame_equal(applied, self.frame * 2)
        result = self.frame.applymap(type)

        # GH #465, function returning tuples
        result = self.frame.applymap(lambda x: (x, x))
        tm.assertIsInstance(result['A'][0], tuple)

        # GH 2909, object conversion to float in constructor?
        df = DataFrame(data=[1,'a'])
        result = df.applymap(lambda x: x)
        self.assertEqual(result.dtypes[0], object)

        df = DataFrame(data=[1.,'a'])
        result = df.applymap(lambda x: x)
        self.assertEqual(result.dtypes[0], object)

        # GH2786
        df  = DataFrame(np.random.random((3,4)))
        df2 = df.copy()
        cols = ['a','a','a','a']
        df.columns = cols

        expected = df2.applymap(str)
        expected.columns = cols
        result = df.applymap(str)
        assert_frame_equal(result,expected)

        # datetime/timedelta
        df['datetime'] = Timestamp('20130101')
        df['timedelta'] = Timedelta('1 min')
        result = df.applymap(str)
        for f in ['datetime','timedelta']:
            self.assertEqual(result.loc[0,f],str(df.loc[0,f]))

    def test_filter(self):
        # items
        filtered = self.frame.filter(['A', 'B', 'E'])
        self.assertEqual(len(filtered.columns), 2)
        self.assertNotIn('E', filtered)

        filtered = self.frame.filter(['A', 'B', 'E'], axis='columns')
        self.assertEqual(len(filtered.columns), 2)
        self.assertNotIn('E', filtered)

        # other axis
        idx = self.frame.index[0:4]
        filtered = self.frame.filter(idx, axis='index')
        expected = self.frame.reindex(index=idx)
        assert_frame_equal(filtered, expected)

        # like
        fcopy = self.frame.copy()
        fcopy['AA'] = 1

        filtered = fcopy.filter(like='A')
        self.assertEqual(len(filtered.columns), 2)
        self.assertIn('AA', filtered)

        # like with ints in column names
        df = DataFrame(0., index=[0, 1, 2], columns=[0, 1, '_A', '_B'])
        filtered = df.filter(like='_')
        self.assertEqual(len(filtered.columns), 2)

        # regex with ints in column names
        # from PR #10384
        df = DataFrame(0., index=[0, 1, 2], columns=['A1', 1, 'B', 2, 'C'])
        expected = DataFrame(0., index=[0, 1, 2], columns=pd.Index([1, 2], dtype=object))
        filtered = df.filter(regex='^[0-9]+$')
        assert_frame_equal(filtered, expected)

        expected = DataFrame(0., index=[0, 1, 2], columns=[0, '0', 1, '1'])
        filtered = expected.filter(regex='^[0-9]+$')  # shouldn't remove anything
        assert_frame_equal(filtered, expected)

        # pass in None
        with assertRaisesRegexp(TypeError, 'Must pass'):
            self.frame.filter(items=None)

        # objects
        filtered = self.mixed_frame.filter(like='foo')
        self.assertIn('foo', filtered)

        # unicode columns, won't ascii-encode
        df = self.frame.rename(columns={'B': u('\u2202')})
        filtered = df.filter(like='C')
        self.assertTrue('C' in filtered)

    def test_filter_regex_search(self):
        fcopy = self.frame.copy()
        fcopy['AA'] = 1

        # regex
        filtered = fcopy.filter(regex='[A]+')
        self.assertEqual(len(filtered.columns), 2)
        self.assertIn('AA', filtered)

        # doesn't have to be at beginning
        df = DataFrame({'aBBa': [1, 2],
                        'BBaBB': [1, 2],
                        'aCCa': [1, 2],
                        'aCCaBB': [1, 2]})

        result = df.filter(regex='BB')
        exp = df[[x for x in df.columns if 'BB' in x]]
        assert_frame_equal(result, exp)

    def test_filter_corner(self):
        empty = DataFrame()

        result = empty.filter([])
        assert_frame_equal(result, empty)

        result = empty.filter(like='foo')
        assert_frame_equal(result, empty)

    def test_select(self):
        f = lambda x: x.weekday() == 2
        result = self.tsframe.select(f, axis=0)
        expected = self.tsframe.reindex(
            index=self.tsframe.index[[f(x) for x in self.tsframe.index]])
        assert_frame_equal(result, expected)

        result = self.frame.select(lambda x: x in ('B', 'D'), axis=1)
        expected = self.frame.reindex(columns=['B', 'D'])

        assert_frame_equal(result, expected, check_names=False)  # TODO should reindex check_names?

    def test_reorder_levels(self):
        index = MultiIndex(levels=[['bar'], ['one', 'two', 'three'], [0, 1]],
                           labels=[[0, 0, 0, 0, 0, 0],
                                   [0, 1, 2, 0, 1, 2],
                                   [0, 1, 0, 1, 0, 1]],
                           names=['L0', 'L1', 'L2'])
        df = DataFrame({'A': np.arange(6), 'B': np.arange(6)}, index=index)

        # no change, position
        result = df.reorder_levels([0, 1, 2])
        assert_frame_equal(df, result)

        # no change, labels
        result = df.reorder_levels(['L0', 'L1', 'L2'])
        assert_frame_equal(df, result)

        # rotate, position
        result = df.reorder_levels([1, 2, 0])
        e_idx = MultiIndex(levels=[['one', 'two', 'three'], [0, 1], ['bar']],
                           labels=[[0, 1, 2, 0, 1, 2],
                                   [0, 1, 0, 1, 0, 1],
                                   [0, 0, 0, 0, 0, 0]],
                           names=['L1', 'L2', 'L0'])
        expected = DataFrame({'A': np.arange(6), 'B': np.arange(6)},
                             index=e_idx)
        assert_frame_equal(result, expected)

        result = df.reorder_levels([0, 0, 0])
        e_idx = MultiIndex(levels=[['bar'], ['bar'], ['bar']],
                           labels=[[0, 0, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0, 0]],
                           names=['L0', 'L0', 'L0'])
        expected = DataFrame({'A': np.arange(6), 'B': np.arange(6)},
                             index=e_idx)
        assert_frame_equal(result, expected)

        result = df.reorder_levels(['L0', 'L0', 'L0'])
        assert_frame_equal(result, expected)

    def test_sort_values(self):

        # API for 9816

        # sort_index
        frame = DataFrame(np.arange(16).reshape(4, 4), index=[1, 2, 3, 4],
                          columns=['A', 'B', 'C', 'D'])

        # 9816 deprecated
        with tm.assert_produces_warning(FutureWarning):
            frame.sort(columns='A')
        with tm.assert_produces_warning(FutureWarning):
            frame.sort()

        unordered = frame.ix[[3, 2, 4, 1]]
        expected = unordered.sort_index()

        result = unordered.sort_index(axis=0)
        assert_frame_equal(result, expected)

        unordered = frame.ix[:, [2, 1, 3, 0]]
        expected = unordered.sort_index(axis=1)

        result = unordered.sort_index(axis=1)
        assert_frame_equal(result, expected)
        assert_frame_equal(result, expected)

        # sortlevel
        mi = MultiIndex.from_tuples([[1, 1, 3], [1, 1, 1]], names=list('ABC'))
        df = DataFrame([[1, 2], [3, 4]], mi)

        result = df.sort_index(level='A', sort_remaining=False)
        expected = df.sortlevel('A', sort_remaining=False)
        assert_frame_equal(result, expected)

        df = df.T
        result = df.sort_index(level='A', axis=1, sort_remaining=False)
        expected = df.sortlevel('A', axis=1, sort_remaining=False)
        assert_frame_equal(result, expected)

        # MI sort, but no by
        mi = MultiIndex.from_tuples([[1, 1, 3], [1, 1, 1]], names=list('ABC'))
        df = DataFrame([[1, 2], [3, 4]], mi)
        result = df.sort_index(sort_remaining=False)
        expected = df.sort_index()
        assert_frame_equal(result, expected)

    def test_sort_index(self):
        frame = DataFrame(np.arange(16).reshape(4, 4), index=[1, 2, 3, 4],
                          columns=['A', 'B', 'C', 'D'])

        # axis=0
        unordered = frame.ix[[3, 2, 4, 1]]
        sorted_df = unordered.sort_index(axis=0)
        expected = frame
        assert_frame_equal(sorted_df, expected)

        sorted_df = unordered.sort_index(ascending=False)
        expected = frame[::-1]
        assert_frame_equal(sorted_df, expected)

        # axis=1
        unordered = frame.ix[:, ['D', 'B', 'C', 'A']]
        sorted_df = unordered.sort_index(axis=1)
        expected = frame
        assert_frame_equal(sorted_df, expected)

        sorted_df = unordered.sort_index(axis=1, ascending=False)
        expected = frame.ix[:, ::-1]
        assert_frame_equal(sorted_df, expected)

        # by column
        sorted_df = frame.sort_values(by='A')
        indexer = frame['A'].argsort().values
        expected = frame.ix[frame.index[indexer]]
        assert_frame_equal(sorted_df, expected)

        sorted_df = frame.sort_values(by='A', ascending=False)
        indexer = indexer[::-1]
        expected = frame.ix[frame.index[indexer]]
        assert_frame_equal(sorted_df, expected)

        sorted_df = frame.sort_values(by='A', ascending=False)
        assert_frame_equal(sorted_df, expected)

        # GH4839
        sorted_df = frame.sort_values(by=['A'], ascending=[False])
        assert_frame_equal(sorted_df, expected)

        # check for now
        sorted_df = frame.sort_values(by='A')
        assert_frame_equal(sorted_df, expected[::-1])
        expected = frame.sort_values(by='A')
        assert_frame_equal(sorted_df, expected)

        expected = frame.sort_values(by=['A', 'B'], ascending=False)
        sorted_df = frame.sort_values(by=['A', 'B'])
        assert_frame_equal(sorted_df, expected[::-1])

        self.assertRaises(ValueError, lambda : frame.sort_values(by=['A','B'], axis=2, inplace=True))

        msg = 'When sorting by column, axis must be 0'
        with assertRaisesRegexp(ValueError, msg):
            frame.sort_values(by='A', axis=1)

        msg = r'Length of ascending \(5\) != length of by \(2\)'
        with assertRaisesRegexp(ValueError, msg):
            frame.sort_values(by=['A', 'B'], axis=0, ascending=[True] * 5)

    def test_sort_index_categorical_index(self):

        df = DataFrame({'A' : np.arange(6,dtype='int64'),
                        'B' : Series(list('aabbca')).astype('category',categories=list('cab')) }).set_index('B')

        result = df.sort_index()
        expected = df.iloc[[4,0,1,5,2,3]]
        assert_frame_equal(result, expected)

        result = df.sort_index(ascending=False)
        expected = df.iloc[[3,2,5,1,0,4]]
        assert_frame_equal(result, expected)

    def test_sort_nan(self):
        # GH3917
        nan = np.nan
        df = DataFrame({'A': [1, 2, nan, 1, 6, 8, 4],
                        'B': [9, nan, 5, 2, 5, 4, 5]})

        # sort one column only
        expected = DataFrame(
            {'A': [nan, 1, 1, 2, 4, 6, 8],
             'B': [5, 9, 2, nan, 5, 5, 4]},
            index=[2, 0, 3, 1, 6, 4, 5])
        sorted_df = df.sort_values(['A'], na_position='first')
        assert_frame_equal(sorted_df, expected)

        expected = DataFrame(
            {'A': [nan, 8, 6, 4, 2, 1, 1],
             'B': [5, 4, 5, 5, nan, 9, 2]},
            index=[2, 5, 4, 6, 1, 0, 3])
        sorted_df = df.sort_values(['A'], na_position='first', ascending=False)
        assert_frame_equal(sorted_df, expected)

        # na_position='last', order
        expected = DataFrame(
            {'A': [1, 1, 2, 4, 6, 8, nan],
             'B': [2, 9, nan, 5, 5, 4, 5]},
            index=[3, 0, 1, 6, 4, 5, 2])
        sorted_df = df.sort_values(['A','B'])
        assert_frame_equal(sorted_df, expected)

        # na_position='first', order
        expected = DataFrame(
            {'A': [nan, 1, 1, 2, 4, 6, 8],
             'B': [5, 2, 9, nan, 5, 5, 4]},
            index=[2, 3, 0, 1, 6, 4, 5])
        sorted_df = df.sort_values(['A','B'], na_position='first')
        assert_frame_equal(sorted_df, expected)

        # na_position='first', not order
        expected = DataFrame(
            {'A': [nan, 1, 1, 2, 4, 6, 8],
             'B': [5, 9, 2, nan, 5, 5, 4]},
            index=[2, 0, 3, 1, 6, 4, 5])
        sorted_df = df.sort_values(['A','B'], ascending=[1,0], na_position='first')
        assert_frame_equal(sorted_df, expected)

        # na_position='last', not order
        expected = DataFrame(
            {'A': [8, 6, 4, 2, 1, 1, nan],
             'B': [4, 5, 5, nan, 2, 9, 5]},
            index=[5, 4, 6, 1, 3, 0, 2])
        sorted_df = df.sort_values(['A','B'], ascending=[0,1], na_position='last')
        assert_frame_equal(sorted_df, expected)

        # Test DataFrame with nan label
        df = DataFrame({'A': [1, 2, nan, 1, 6, 8, 4],
                        'B': [9, nan, 5, 2, 5, 4, 5]},
                       index = [1, 2, 3, 4, 5, 6, nan])

        # NaN label, ascending=True, na_position='last'
        sorted_df = df.sort_index(kind='quicksort', ascending=True, na_position='last')
        expected = DataFrame({'A': [1, 2, nan, 1, 6, 8, 4],
                              'B': [9, nan, 5, 2, 5, 4, 5]},
                             index = [1, 2, 3, 4, 5, 6, nan])
        assert_frame_equal(sorted_df, expected)

        # NaN label, ascending=True, na_position='first'
        sorted_df = df.sort_index(na_position='first')
        expected = DataFrame({'A': [4, 1, 2, nan, 1, 6, 8],
                              'B': [5, 9, nan, 5, 2, 5, 4]},
                             index = [nan, 1, 2, 3, 4, 5, 6])
        assert_frame_equal(sorted_df, expected)

        # NaN label, ascending=False, na_position='last'
        sorted_df = df.sort_index(kind='quicksort', ascending=False)
        expected = DataFrame({'A': [8, 6, 1, nan, 2,   1, 4],
                              'B': [4, 5, 2, 5,   nan, 9, 5]},
                             index = [6, 5, 4, 3, 2, 1, nan])
        assert_frame_equal(sorted_df, expected)

        # NaN label, ascending=False, na_position='first'
        sorted_df = df.sort_index(kind='quicksort', ascending=False, na_position='first')
        expected = DataFrame({'A': [4, 8, 6, 1, nan, 2,   1],
                              'B': [5, 4, 5, 2, 5,   nan, 9]},
                             index = [nan, 6, 5, 4, 3, 2, 1])
        assert_frame_equal(sorted_df, expected)

    def test_stable_descending_sort(self):
        # GH #6399
        df = DataFrame([[2, 'first'], [2, 'second'], [1, 'a'], [1, 'b']],
                       columns=['sort_col', 'order'])
        sorted_df = df.sort_values(by='sort_col', kind='mergesort',
                                   ascending=False)
        assert_frame_equal(df, sorted_df)

    def test_stable_descending_multicolumn_sort(self):
        nan = np.nan
        df = DataFrame({'A': [1, 2, nan, 1, 6, 8, 4],
                        'B': [9, nan, 5, 2, 5, 4, 5]})
        # test stable mergesort
        expected = DataFrame(
            {'A': [nan, 8, 6, 4, 2, 1, 1],
             'B': [5, 4, 5, 5, nan, 2, 9]},
            index=[2, 5, 4, 6, 1, 3, 0])
        sorted_df = df.sort_values(['A','B'], ascending=[0,1], na_position='first',
                                   kind='mergesort')
        assert_frame_equal(sorted_df, expected)

        expected = DataFrame(
            {'A': [nan, 8, 6, 4, 2, 1, 1],
             'B': [5, 4, 5, 5, nan, 9, 2]},
            index=[2, 5, 4, 6, 1, 0, 3])
        sorted_df = df.sort_values(['A','B'], ascending=[0,0], na_position='first',
                                   kind='mergesort')
        assert_frame_equal(sorted_df, expected)

    def test_sort_index_multicolumn(self):
        import random
        A = np.arange(5).repeat(20)
        B = np.tile(np.arange(5), 20)
        random.shuffle(A)
        random.shuffle(B)
        frame = DataFrame({'A': A, 'B': B,
                           'C': np.random.randn(100)})

        # use .sort_values #9816
        with tm.assert_produces_warning(FutureWarning):
            frame.sort_index(by=['A', 'B'])
        result = frame.sort_values(by=['A', 'B'])
        indexer = np.lexsort((frame['B'], frame['A']))
        expected = frame.take(indexer)
        assert_frame_equal(result, expected)

        # use .sort_values #9816
        with tm.assert_produces_warning(FutureWarning):
            frame.sort_index(by=['A', 'B'], ascending=False)
        result = frame.sort_values(by=['A', 'B'], ascending=False)
        indexer = np.lexsort((frame['B'].rank(ascending=False),
                              frame['A'].rank(ascending=False)))
        expected = frame.take(indexer)
        assert_frame_equal(result, expected)

        # use .sort_values #9816
        with tm.assert_produces_warning(FutureWarning):
            frame.sort_index(by=['B', 'A'])
        result = frame.sort_values(by=['B', 'A'])
        indexer = np.lexsort((frame['A'], frame['B']))
        expected = frame.take(indexer)
        assert_frame_equal(result, expected)

    def test_sort_index_inplace(self):
        frame = DataFrame(np.random.randn(4, 4), index=[1, 2, 3, 4],
                          columns=['A', 'B', 'C', 'D'])

        # axis=0
        unordered = frame.ix[[3, 2, 4, 1]]
        a_id = id(unordered['A'])
        df = unordered.copy()
        df.sort_index(inplace=True)
        expected = frame
        assert_frame_equal(df, expected)
        self.assertNotEqual(a_id, id(df['A']))

        df = unordered.copy()
        df.sort_index(ascending=False, inplace=True)
        expected = frame[::-1]
        assert_frame_equal(df, expected)

        # axis=1
        unordered = frame.ix[:, ['D', 'B', 'C', 'A']]
        df = unordered.copy()
        df.sort_index(axis=1, inplace=True)
        expected = frame
        assert_frame_equal(df, expected)

        df = unordered.copy()
        df.sort_index(axis=1, ascending=False, inplace=True)
        expected = frame.ix[:, ::-1]
        assert_frame_equal(df, expected)

    def test_sort_index_different_sortorder(self):
        A = np.arange(20).repeat(5)
        B = np.tile(np.arange(5), 20)

        indexer = np.random.permutation(100)
        A = A.take(indexer)
        B = B.take(indexer)

        df = DataFrame({'A': A, 'B': B,
                        'C': np.random.randn(100)})

        # use .sort_values #9816
        with tm.assert_produces_warning(FutureWarning):
            df.sort_index(by=['A', 'B'], ascending=[1, 0])
        result = df.sort_values(by=['A', 'B'], ascending=[1, 0])

        ex_indexer = np.lexsort((df.B.max() - df.B, df.A))
        expected = df.take(ex_indexer)
        assert_frame_equal(result, expected)

        # test with multiindex, too
        idf = df.set_index(['A', 'B'])

        result = idf.sort_index(ascending=[1, 0])
        expected = idf.take(ex_indexer)
        assert_frame_equal(result, expected)

        # also, Series!
        result = idf['C'].sort_index(ascending=[1, 0])
        assert_series_equal(result, expected['C'])

    def test_sort_inplace(self):
        frame = DataFrame(np.random.randn(4, 4), index=[1, 2, 3, 4],
                          columns=['A', 'B', 'C', 'D'])

        sorted_df = frame.copy()
        sorted_df.sort_values(by='A', inplace=True)
        expected = frame.sort_values(by='A')
        assert_frame_equal(sorted_df, expected)

        sorted_df = frame.copy()
        sorted_df.sort_values(by='A', ascending=False, inplace=True)
        expected = frame.sort_values(by='A', ascending=False)
        assert_frame_equal(sorted_df, expected)

        sorted_df = frame.copy()
        sorted_df.sort_values(by=['A', 'B'], ascending=False, inplace=True)
        expected = frame.sort_values(by=['A', 'B'], ascending=False)
        assert_frame_equal(sorted_df, expected)

    def test_sort_index_duplicates(self):

        ### with 9816, these are all translated to .sort_values

        df = DataFrame([lrange(5,9), lrange(4)],
                       columns=['a', 'a', 'b', 'b'])

        with assertRaisesRegexp(ValueError, 'duplicate'):
            # use .sort_values #9816
            with tm.assert_produces_warning(FutureWarning):
                df.sort_index(by='a')
        with assertRaisesRegexp(ValueError, 'duplicate'):
                df.sort_values(by='a')

        with assertRaisesRegexp(ValueError, 'duplicate'):
            # use .sort_values #9816
            with tm.assert_produces_warning(FutureWarning):
                df.sort_index(by=['a'])
        with assertRaisesRegexp(ValueError, 'duplicate'):
            df.sort_values(by=['a'])

        with assertRaisesRegexp(ValueError, 'duplicate'):
            # use .sort_values #9816
            with tm.assert_produces_warning(FutureWarning):
                # multi-column 'by' is separate codepath
                df.sort_index(by=['a', 'b'])
        with assertRaisesRegexp(ValueError, 'duplicate'):
            # multi-column 'by' is separate codepath
            df.sort_values(by=['a', 'b'])

        # with multi-index
        # GH4370
        df = DataFrame(np.random.randn(4,2),columns=MultiIndex.from_tuples([('a',0),('a',1)]))
        with assertRaisesRegexp(ValueError, 'levels'):
            # use .sort_values #9816
            with tm.assert_produces_warning(FutureWarning):
                df.sort_index(by='a')
        with assertRaisesRegexp(ValueError, 'levels'):
            df.sort_values(by='a')

        # convert tuples to a list of tuples
        # use .sort_values #9816
        with tm.assert_produces_warning(FutureWarning):
            df.sort_index(by=[('a',1)])
        expected = df.sort_values(by=[('a',1)])

        # use .sort_values #9816
        with tm.assert_produces_warning(FutureWarning):
            df.sort_index(by=('a',1))
        result = df.sort_values(by=('a',1))
        assert_frame_equal(result, expected)

    def test_sortlevel(self):
        mi = MultiIndex.from_tuples([[1, 1, 3], [1, 1, 1]], names=list('ABC'))
        df = DataFrame([[1, 2], [3, 4]], mi)
        res = df.sortlevel('A', sort_remaining=False)
        assert_frame_equal(df, res)

        res = df.sortlevel(['A', 'B'], sort_remaining=False)
        assert_frame_equal(df, res)

    def test_sort_datetimes(self):

        # GH 3461, argsort / lexsort differences for a datetime column
        df = DataFrame(['a','a','a','b','c','d','e','f','g'],
                       columns=['A'],
                       index=date_range('20130101',periods=9))
        dts = [Timestamp(x)
               for x in  ['2004-02-11','2004-01-21','2004-01-26',
                          '2005-09-20','2010-10-04','2009-05-12',
                          '2008-11-12','2010-09-28','2010-09-28']]
        df['B'] = dts[::2] + dts[1::2]
        df['C'] = 2.
        df['A1'] = 3.

        df1 = df.sort_values(by='A')
        df2 = df.sort_values(by=['A'])
        assert_frame_equal(df1,df2)

        df1 = df.sort_values(by='B')
        df2 = df.sort_values(by=['B'])
        assert_frame_equal(df1,df2)

    def test_frame_column_inplace_sort_exception(self):
        s = self.frame['A']
        with assertRaisesRegexp(ValueError, "This Series is a view"):
            s.sort_values(inplace=True)

        cp = s.copy()
        cp.sort_values() # it works!

    def test_combine_first(self):
        # disjoint
        head, tail = self.frame[:5], self.frame[5:]

        combined = head.combine_first(tail)
        reordered_frame = self.frame.reindex(combined.index)
        assert_frame_equal(combined, reordered_frame)
        self.assertTrue(tm.equalContents(combined.columns, self.frame.columns))
        assert_series_equal(combined['A'], reordered_frame['A'])

        # same index
        fcopy = self.frame.copy()
        fcopy['A'] = 1
        del fcopy['C']

        fcopy2 = self.frame.copy()
        fcopy2['B'] = 0
        del fcopy2['D']

        combined = fcopy.combine_first(fcopy2)

        self.assertTrue((combined['A'] == 1).all())
        assert_series_equal(combined['B'], fcopy['B'])
        assert_series_equal(combined['C'], fcopy2['C'])
        assert_series_equal(combined['D'], fcopy['D'])

        # overlap
        head, tail = reordered_frame[:10].copy(), reordered_frame
        head['A'] = 1

        combined = head.combine_first(tail)
        self.assertTrue((combined['A'][:10] == 1).all())

        # reverse overlap
        tail['A'][:10] = 0
        combined = tail.combine_first(head)
        self.assertTrue((combined['A'][:10] == 0).all())

        # no overlap
        f = self.frame[:10]
        g = self.frame[10:]
        combined = f.combine_first(g)
        assert_series_equal(combined['A'].reindex(f.index), f['A'])
        assert_series_equal(combined['A'].reindex(g.index), g['A'])

        # corner cases
        comb = self.frame.combine_first(self.empty)
        assert_frame_equal(comb, self.frame)

        comb = self.empty.combine_first(self.frame)
        assert_frame_equal(comb, self.frame)

        comb = self.frame.combine_first(DataFrame(index=["faz", "boo"]))
        self.assertTrue("faz" in comb.index)

        # #2525
        df = DataFrame({'a': [1]}, index=[datetime(2012, 1, 1)])
        df2 = DataFrame({}, columns=['b'])
        result = df.combine_first(df2)
        self.assertTrue('b' in result)

    def test_combine_first_mixed_bug(self):
        idx = Index(['a', 'b', 'c', 'e'])
        ser1 = Series([5.0, -9.0, 4.0, 100.], index=idx)
        ser2 = Series(['a', 'b', 'c', 'e'], index=idx)
        ser3 = Series([12, 4, 5, 97], index=idx)

        frame1 = DataFrame({"col0": ser1,
                            "col2": ser2,
                            "col3": ser3})

        idx = Index(['a', 'b', 'c', 'f'])
        ser1 = Series([5.0, -9.0, 4.0, 100.], index=idx)
        ser2 = Series(['a', 'b', 'c', 'f'], index=idx)
        ser3 = Series([12, 4, 5, 97], index=idx)

        frame2 = DataFrame({"col1": ser1,
                            "col2": ser2,
                            "col5": ser3})

        combined = frame1.combine_first(frame2)
        self.assertEqual(len(combined.columns), 5)

        # gh 3016 (same as in update)
        df = DataFrame([[1.,2.,False, True],[4.,5.,True,False]],
                       columns=['A','B','bool1','bool2'])

        other = DataFrame([[45,45]],index=[0],columns=['A','B'])
        result = df.combine_first(other)
        assert_frame_equal(result, df)

        df.ix[0,'A'] = np.nan
        result = df.combine_first(other)
        df.ix[0,'A'] = 45
        assert_frame_equal(result, df)

        # doc example
        df1 = DataFrame({'A' : [1., np.nan, 3., 5., np.nan],
                         'B' : [np.nan, 2., 3., np.nan, 6.]})

        df2 = DataFrame({'A' : [5., 2., 4., np.nan, 3., 7.],
                         'B' : [np.nan, np.nan, 3., 4., 6., 8.]})

        result = df1.combine_first(df2)
        expected = DataFrame({ 'A' : [1,2,3,5,3,7.], 'B' : [np.nan,2,3,4,6,8] })
        assert_frame_equal(result,expected)

        # GH3552, return object dtype with bools
        df1 = DataFrame([[np.nan, 3.,True], [-4.6, np.nan, True], [np.nan, 7., False]])
        df2 = DataFrame([[-42.6, np.nan, True], [-5., 1.6, False]], index=[1, 2])

        result = df1.combine_first(df2)[2]
        expected = Series([True, True, False], name=2)
        assert_series_equal(result, expected)

        # GH 3593, converting datetime64[ns] incorrecly
        df0 = DataFrame({"a":[datetime(2000, 1, 1), datetime(2000, 1, 2), datetime(2000, 1, 3)]})
        df1 = DataFrame({"a":[None, None, None]})
        df2 = df1.combine_first(df0)
        assert_frame_equal(df2, df0)

        df2 = df0.combine_first(df1)
        assert_frame_equal(df2, df0)

        df0 = DataFrame({"a":[datetime(2000, 1, 1), datetime(2000, 1, 2), datetime(2000, 1, 3)]})
        df1 = DataFrame({"a":[datetime(2000, 1, 2), None, None]})
        df2 = df1.combine_first(df0)
        result = df0.copy()
        result.iloc[0,:] = df1.iloc[0,:]
        assert_frame_equal(df2, result)

        df2 = df0.combine_first(df1)
        assert_frame_equal(df2, df0)

    def test_update(self):
        df = DataFrame([[1.5, nan, 3.],
                        [1.5, nan, 3.],
                        [1.5, nan, 3],
                        [1.5, nan, 3]])

        other = DataFrame([[3.6, 2., np.nan],
                           [np.nan, np.nan, 7]], index=[1, 3])

        df.update(other)

        expected = DataFrame([[1.5, nan, 3],
                              [3.6, 2, 3],
                              [1.5, nan, 3],
                              [1.5, nan, 7.]])
        assert_frame_equal(df, expected)

    def test_update_dtypes(self):

        # gh 3016
        df = DataFrame([[1.,2.,False, True],[4.,5.,True,False]],
                       columns=['A','B','bool1','bool2'])

        other = DataFrame([[45,45]],index=[0],columns=['A','B'])
        df.update(other)

        expected = DataFrame([[45.,45.,False, True],[4.,5.,True,False]],
                             columns=['A','B','bool1','bool2'])
        assert_frame_equal(df, expected)

    def test_update_nooverwrite(self):
        df = DataFrame([[1.5, nan, 3.],
                        [1.5, nan, 3.],
                        [1.5, nan, 3],
                        [1.5, nan, 3]])

        other = DataFrame([[3.6, 2., np.nan],
                           [np.nan, np.nan, 7]], index=[1, 3])

        df.update(other, overwrite=False)

        expected = DataFrame([[1.5, nan, 3],
                              [1.5, 2, 3],
                              [1.5, nan, 3],
                              [1.5, nan, 3.]])
        assert_frame_equal(df, expected)

    def test_update_filtered(self):
        df = DataFrame([[1.5, nan, 3.],
                        [1.5, nan, 3.],
                        [1.5, nan, 3],
                        [1.5, nan, 3]])

        other = DataFrame([[3.6, 2., np.nan],
                           [np.nan, np.nan, 7]], index=[1, 3])

        df.update(other, filter_func=lambda x: x > 2)

        expected = DataFrame([[1.5, nan, 3],
                              [1.5, nan, 3],
                              [1.5, nan, 3],
                              [1.5, nan, 7.]])
        assert_frame_equal(df, expected)

    def test_update_raise(self):
        df = DataFrame([[1.5, 1, 3.],
                        [1.5, nan, 3.],
                        [1.5, nan, 3],
                        [1.5, nan, 3]])

        other = DataFrame([[2., nan],
                           [nan, 7]], index=[1, 3], columns=[1, 2])
        with assertRaisesRegexp(ValueError, "Data overlaps"):
            df.update(other, raise_conflict=True)

    def test_update_from_non_df(self):
        d = {'a': Series([1, 2, 3, 4]), 'b': Series([5, 6, 7, 8])}
        df = DataFrame(d)

        d['a'] = Series([5, 6, 7, 8])
        df.update(d)

        expected = DataFrame(d)

        assert_frame_equal(df, expected)

        d = {'a': [1, 2, 3, 4], 'b': [5, 6, 7, 8]}
        df = DataFrame(d)

        d['a'] = [5, 6, 7, 8]
        df.update(d)

        expected = DataFrame(d)

        assert_frame_equal(df, expected)

    def test_clip(self):
        median = self.frame.median().median()

        capped = self.frame.clip_upper(median)
        self.assertFalse((capped.values > median).any())

        floored = self.frame.clip_lower(median)
        self.assertFalse((floored.values < median).any())

        double = self.frame.clip(upper=median, lower=median)
        self.assertFalse((double.values != median).any())

    def test_dataframe_clip(self):
        # GH #2747
        df = DataFrame(np.random.randn(1000,2))

        for lb, ub in [(-1,1),(1,-1)]:
            clipped_df = df.clip(lb, ub)

            lb, ub = min(lb,ub), max(ub,lb)
            lb_mask = df.values <= lb
            ub_mask = df.values >= ub
            mask = ~lb_mask & ~ub_mask
            self.assertTrue((clipped_df.values[lb_mask] == lb).all() == True)
            self.assertTrue((clipped_df.values[ub_mask] == ub).all() == True)
            self.assertTrue((clipped_df.values[mask] ==
                             df.values[mask]).all() == True)

    def test_clip_against_series(self):
        # GH #6966

        df = DataFrame(np.random.randn(1000, 2))
        lb = Series(np.random.randn(1000))
        ub = lb + 1

        clipped_df = df.clip(lb, ub, axis=0)

        for i in range(2):
            lb_mask = df.iloc[:, i] <= lb
            ub_mask = df.iloc[:, i] >= ub
            mask = ~lb_mask & ~ub_mask

            result = clipped_df.loc[lb_mask, i]
            assert_series_equal(result, lb[lb_mask], check_names=False)
            self.assertEqual(result.name, i)

            result = clipped_df.loc[ub_mask, i]
            assert_series_equal(result, ub[ub_mask], check_names=False)
            self.assertEqual(result.name, i)

            assert_series_equal(clipped_df.loc[mask, i], df.loc[mask, i])

    def test_clip_against_frame(self):
        df = DataFrame(np.random.randn(1000, 2))
        lb = DataFrame(np.random.randn(1000, 2))
        ub = lb + 1

        clipped_df = df.clip(lb, ub)

        lb_mask = df <= lb
        ub_mask = df >= ub
        mask = ~lb_mask & ~ub_mask

        assert_frame_equal(clipped_df[lb_mask], lb[lb_mask])
        assert_frame_equal(clipped_df[ub_mask], ub[ub_mask])
        assert_frame_equal(clipped_df[mask], df[mask])

    def test_get_X_columns(self):
        # numeric and object columns

        df = DataFrame({'a': [1, 2, 3],
                        'b' : [True, False, True],
                        'c': ['foo', 'bar', 'baz'],
                        'd': [None, None, None],
                        'e': [3.14, 0.577, 2.773]})

        self.assert_numpy_array_equal(df._get_numeric_data().columns,
                                      ['a', 'b', 'e'])

    def test_is_mixed_type(self):
        self.assertFalse(self.frame._is_mixed_type)
        self.assertTrue(self.mixed_frame._is_mixed_type)

    def test_get_numeric_data(self):
        intname = np.dtype(np.int_).name
        floatname = np.dtype(np.float_).name
        datetime64name = np.dtype('M8[ns]').name
        objectname = np.dtype(np.object_).name

        df = DataFrame({'a': 1., 'b': 2, 'c': 'foo', 'f' : Timestamp('20010102')},
                       index=np.arange(10))
        result = df.get_dtype_counts()
        expected = Series({'int64': 1, 'float64' : 1, datetime64name: 1, objectname : 1})
        result.sort_index()
        expected.sort_index()
        assert_series_equal(result, expected)

        df = DataFrame({'a': 1., 'b': 2, 'c': 'foo',
                        'd' : np.array([1.]*10,dtype='float32'),
                        'e' : np.array([1]*10,dtype='int32'),
                        'f' : np.array([1]*10,dtype='int16'),
                        'g' : Timestamp('20010102')},
                       index=np.arange(10))

        result = df._get_numeric_data()
        expected = df.ix[:, ['a', 'b','d','e','f']]
        assert_frame_equal(result, expected)

        only_obj = df.ix[:, ['c','g']]
        result = only_obj._get_numeric_data()
        expected = df.ix[:, []]
        assert_frame_equal(result, expected)

        df = DataFrame.from_dict({'a':[1,2], 'b':['foo','bar'],'c':[np.pi,np.e]})
        result = df._get_numeric_data()
        expected = DataFrame.from_dict({'a':[1,2], 'c':[np.pi,np.e]})
        assert_frame_equal(result, expected)

        df = result.copy()
        result = df._get_numeric_data()
        expected = df
        assert_frame_equal(result, expected)

    def test_round(self):

        # GH 2665

        # Test that rounding an empty DataFrame does nothing
        df = DataFrame()
        assert_frame_equal(df, df.round())

        # Here's the test frame we'll be working with
        df = DataFrame(
            {'col1': [1.123, 2.123, 3.123], 'col2': [1.234, 2.234, 3.234]})

        # Default round to integer (i.e. decimals=0)
        expected_rounded = DataFrame(
            {'col1': [1., 2., 3.], 'col2': [1., 2., 3.]})
        assert_frame_equal(df.round(), expected_rounded)

        # Round with an integer
        decimals = 2
        expected_rounded = DataFrame(
            {'col1': [1.12, 2.12, 3.12], 'col2': [1.23, 2.23, 3.23]})
        assert_frame_equal(df.round(decimals), expected_rounded)

        # This should also work with np.round (since np.round dispatches to
        # df.round)
        assert_frame_equal(np.round(df, decimals), expected_rounded)

        # Round with a list
        round_list = [1, 2]
        with self.assertRaises(TypeError):
            df.round(round_list)

        # Round with a dictionary
        expected_rounded = DataFrame(
            {'col1': [1.1, 2.1, 3.1], 'col2': [1.23, 2.23, 3.23]})
        round_dict = {'col1': 1, 'col2': 2}
        assert_frame_equal(df.round(round_dict), expected_rounded)

        # Incomplete dict
        expected_partially_rounded = DataFrame(
            {'col1': [1.123, 2.123, 3.123], 'col2': [1.2, 2.2, 3.2]})
        partial_round_dict = {'col2': 1}
        assert_frame_equal(
            df.round(partial_round_dict), expected_partially_rounded)

        # Dict with unknown elements
        wrong_round_dict = {'col3': 2, 'col2': 1}
        assert_frame_equal(
            df.round(wrong_round_dict), expected_partially_rounded)

        # float input to `decimals`
        non_int_round_dict = {'col1': 1, 'col2': 0.5}
        with self.assertRaises(TypeError):
            df.round(non_int_round_dict)

        # String input
        non_int_round_dict = {'col1': 1, 'col2': 'foo'}
        with self.assertRaises(TypeError):
            df.round(non_int_round_dict)

        non_int_round_Series = Series(non_int_round_dict)
        with self.assertRaises(TypeError):
            df.round(non_int_round_Series)

        # List input
        non_int_round_dict = {'col1': 1, 'col2': [1, 2]}
        with self.assertRaises(TypeError):
            df.round(non_int_round_dict)

        non_int_round_Series = Series(non_int_round_dict)
        with self.assertRaises(TypeError):
            df.round(non_int_round_Series)

        # Non integer Series inputs
        non_int_round_Series = Series(non_int_round_dict)
        with self.assertRaises(TypeError):
            df.round(non_int_round_Series)

        non_int_round_Series = Series(non_int_round_dict)
        with self.assertRaises(TypeError):
            df.round(non_int_round_Series)

        # Negative numbers
        negative_round_dict = {'col1': -1, 'col2': -2}
        big_df = df * 100
        expected_neg_rounded = DataFrame(
            {'col1': [110., 210, 310], 'col2': [100., 200, 300]})
        assert_frame_equal(
            big_df.round(negative_round_dict), expected_neg_rounded)

        # nan in Series round
        nan_round_Series = Series({'col1': nan, 'col2':1})
        expected_nan_round = DataFrame(
                {'col1': [1.123, 2.123, 3.123], 'col2': [1.2, 2.2, 3.2]})
        if sys.version < LooseVersion('2.7'):
            # Rounding with decimal is a ValueError in Python < 2.7
            with self.assertRaises(ValueError):
                df.round(nan_round_Series)
        else:
            with self.assertRaises(TypeError):
                df.round(nan_round_Series)

        # Make sure this doesn't break existing Series.round
        assert_series_equal(df['col1'].round(1), expected_rounded['col1'])

        # named columns
        # GH 11986
        decimals = 2
        expected_rounded = DataFrame(
            {'col1': [1.12, 2.12, 3.12], 'col2': [1.23, 2.23, 3.23]})
        df.columns.name = "cols"
        expected_rounded.columns.name = "cols"
        assert_frame_equal(df.round(decimals), expected_rounded)

        # interaction of named columns & series
        assert_series_equal(df['col1'].round(decimals),
                            expected_rounded['col1'])
        assert_series_equal(df.round(decimals)['col1'],
                            expected_rounded['col1'])

    def test_round_mixed_type(self):
        # GH11885
        df = DataFrame({'col1': [1.1, 2.2, 3.3, 4.4],
                        'col2': ['1', 'a', 'c', 'f'],
                        'col3': date_range('20111111', periods=4)})
        round_0 = DataFrame({'col1': [1., 2., 3., 4.],
                             'col2': ['1', 'a', 'c', 'f'],
                             'col3': date_range('20111111', periods=4)})
        assert_frame_equal(df.round(), round_0)
        assert_frame_equal(df.round(1), df)
        assert_frame_equal(df.round({'col1': 1}), df)
        assert_frame_equal(df.round({'col1': 0}), round_0)
        assert_frame_equal(df.round({'col1': 0, 'col2': 1}), round_0)
        assert_frame_equal(df.round({'col3': 1}), df)

    def test_round_issue(self):
        # GH11611

        df = pd.DataFrame(np.random.random([3, 3]), columns=['A', 'B', 'C'],
                          index=['first', 'second', 'third'])

        dfs = pd.concat((df, df), axis=1)
        rounded = dfs.round()
        self.assertTrue(rounded.index.equals(dfs.index))

        decimals = pd.Series([1, 0, 2], index=['A', 'B', 'A'])
        self.assertRaises(ValueError, df.round, decimals)

    def test_built_in_round(self):
        if not compat.PY3:
            raise nose.SkipTest("build in round cannot be overriden "
                                "prior to Python 3")

        # GH11763
        # Here's the test frame we'll be working with
        df = DataFrame(
            {'col1': [1.123, 2.123, 3.123], 'col2': [1.234, 2.234, 3.234]})

        # Default round to integer (i.e. decimals=0)
        expected_rounded = DataFrame(
            {'col1': [1., 2., 3.], 'col2': [1., 2., 3.]})
        assert_frame_equal(round(df), expected_rounded)

    def test_axis_aliases(self):

        f = self.frame

        # reg name
        expected = f.sum(axis=0)
        result = f.sum(axis='index')
        assert_series_equal(result, expected)

        expected = f.sum(axis=1)
        result = f.sum(axis='columns')
        assert_series_equal(result, expected)

    def test_combine_first_mixed(self):
        a = Series(['a', 'b'], index=lrange(2))
        b = Series(lrange(2), index=lrange(2))
        f = DataFrame({'A': a, 'B': b})

        a = Series(['a', 'b'], index=lrange(5, 7))
        b = Series(lrange(2), index=lrange(5, 7))
        g = DataFrame({'A': a, 'B': b})

        combined = f.combine_first(g)

    def test_more_asMatrix(self):
        values = self.mixed_frame.as_matrix()
        self.assertEqual(values.shape[1], len(self.mixed_frame.columns))

    def test_reindex_boolean(self):
        frame = DataFrame(np.ones((10, 2), dtype=bool),
                          index=np.arange(0, 20, 2),
                          columns=[0, 2])

        reindexed = frame.reindex(np.arange(10))
        self.assertEqual(reindexed.values.dtype, np.object_)
        self.assertTrue(isnull(reindexed[0][1]))

        reindexed = frame.reindex(columns=lrange(3))
        self.assertEqual(reindexed.values.dtype, np.object_)
        self.assertTrue(isnull(reindexed[1]).all())

    def test_reindex_objects(self):
        reindexed = self.mixed_frame.reindex(columns=['foo', 'A', 'B'])
        self.assertIn('foo', reindexed)

        reindexed = self.mixed_frame.reindex(columns=['A', 'B'])
        self.assertNotIn('foo', reindexed)

    def test_reindex_corner(self):
        index = Index(['a', 'b', 'c'])
        dm = self.empty.reindex(index=[1, 2, 3])
        reindexed = dm.reindex(columns=index)
        self.assertTrue(reindexed.columns.equals(index))

        # ints are weird

        smaller = self.intframe.reindex(columns=['A', 'B', 'E'])
        self.assertEqual(smaller['E'].dtype, np.float64)

    def test_reindex_axis(self):
        cols = ['A', 'B', 'E']
        reindexed1 = self.intframe.reindex_axis(cols, axis=1)
        reindexed2 = self.intframe.reindex(columns=cols)
        assert_frame_equal(reindexed1, reindexed2)

        rows = self.intframe.index[0:5]
        reindexed1 = self.intframe.reindex_axis(rows, axis=0)
        reindexed2 = self.intframe.reindex(index=rows)
        assert_frame_equal(reindexed1, reindexed2)

        self.assertRaises(ValueError, self.intframe.reindex_axis, rows, axis=2)

        # no-op case
        cols = self.frame.columns.copy()
        newFrame = self.frame.reindex_axis(cols, axis=1)
        assert_frame_equal(newFrame, self.frame)

    def test_reindex_with_nans(self):
        df = DataFrame([[1, 2], [3, 4], [np.nan, np.nan], [7, 8], [9, 10]],
                       columns=['a', 'b'],
                       index=[100.0, 101.0, np.nan, 102.0, 103.0])

        result = df.reindex(index=[101.0, 102.0, 103.0])
        expected = df.iloc[[1, 3, 4]]
        assert_frame_equal(result, expected)

        result = df.reindex(index=[103.0])
        expected = df.iloc[[4]]
        assert_frame_equal(result, expected)

        result = df.reindex(index=[101.0])
        expected = df.iloc[[1]]
        assert_frame_equal(result, expected)

    def test_reindex_multi(self):
        df = DataFrame(np.random.randn(3, 3))

        result = df.reindex(lrange(4), lrange(4))
        expected = df.reindex(lrange(4)).reindex(columns=lrange(4))

        assert_frame_equal(result, expected)

        df = DataFrame(np.random.randint(0, 10, (3, 3)))

        result = df.reindex(lrange(4), lrange(4))
        expected = df.reindex(lrange(4)).reindex(columns=lrange(4))

        assert_frame_equal(result, expected)

        df = DataFrame(np.random.randint(0, 10, (3, 3)))

        result = df.reindex(lrange(2), lrange(2))
        expected = df.reindex(lrange(2)).reindex(columns=lrange(2))

        assert_frame_equal(result, expected)

        df = DataFrame(np.random.randn(5, 3) + 1j, columns=['a', 'b', 'c'])

        result = df.reindex(index=[0, 1], columns=['a', 'b'])
        expected = df.reindex([0, 1]).reindex(columns=['a', 'b'])

        assert_frame_equal(result, expected)

    def test_rename_objects(self):
        renamed = self.mixed_frame.rename(columns=str.upper)
        self.assertIn('FOO', renamed)
        self.assertNotIn('foo', renamed)

    def test_fill_corner(self):
        self.mixed_frame.ix[5:20,'foo'] = nan
        self.mixed_frame.ix[-10:,'A'] = nan

        filled = self.mixed_frame.fillna(value=0)
        self.assertTrue((filled.ix[5:20,'foo'] == 0).all())
        del self.mixed_frame['foo']

        empty_float = self.frame.reindex(columns=[])
        result = empty_float.fillna(value=0)

    def test_repr_with_mi_nat(self):
        df = DataFrame({'X': [1, 2]},
                       index=[[pd.NaT, pd.Timestamp('20130101')], ['a', 'b']])
        res = repr(df)
        exp = '              X\nNaT        a  1\n2013-01-01 b  2'
        nose.tools.assert_equal(res, exp)

    def test_reset_index(self):
        stacked = self.frame.stack()[::2]
        stacked = DataFrame({'foo': stacked, 'bar': stacked})

        names = ['first', 'second']
        stacked.index.names = names
        deleveled = stacked.reset_index()
        for i, (lev, lab) in enumerate(zip(stacked.index.levels,
                                           stacked.index.labels)):
            values = lev.take(lab)
            name = names[i]
            assert_almost_equal(values, deleveled[name])

        stacked.index.names = [None, None]
        deleveled2 = stacked.reset_index()
        self.assert_numpy_array_equal(deleveled['first'],
                                      deleveled2['level_0'])
        self.assert_numpy_array_equal(deleveled['second'],
                                      deleveled2['level_1'])

        # default name assigned
        rdf = self.frame.reset_index()
        self.assert_numpy_array_equal(rdf['index'], self.frame.index.values)

        # default name assigned, corner case
        df = self.frame.copy()
        df['index'] = 'foo'
        rdf = df.reset_index()
        self.assert_numpy_array_equal(rdf['level_0'], self.frame.index.values)

        # but this is ok
        self.frame.index.name = 'index'
        deleveled = self.frame.reset_index()
        self.assert_numpy_array_equal(deleveled['index'],
                                      self.frame.index.values)
        self.assert_numpy_array_equal(deleveled.index,
                                      np.arange(len(deleveled)))

        # preserve column names
        self.frame.columns.name = 'columns'
        resetted = self.frame.reset_index()
        self.assertEqual(resetted.columns.name, 'columns')

        # only remove certain columns
        frame = self.frame.reset_index().set_index(['index', 'A', 'B'])
        rs = frame.reset_index(['A', 'B'])

        assert_frame_equal(rs, self.frame, check_names=False)  # TODO should reset_index check_names ?

        rs = frame.reset_index(['index', 'A', 'B'])
        assert_frame_equal(rs, self.frame.reset_index(), check_names=False)

        rs = frame.reset_index(['index', 'A', 'B'])
        assert_frame_equal(rs, self.frame.reset_index(), check_names=False)

        rs = frame.reset_index('A')
        xp = self.frame.reset_index().set_index(['index', 'B'])
        assert_frame_equal(rs, xp, check_names=False)

        # test resetting in place
        df = self.frame.copy()
        resetted = self.frame.reset_index()
        df.reset_index(inplace=True)
        assert_frame_equal(df, resetted, check_names=False)

        frame = self.frame.reset_index().set_index(['index', 'A', 'B'])
        rs = frame.reset_index('A', drop=True)
        xp = self.frame.copy()
        del xp['A']
        xp = xp.set_index(['B'], append=True)
        assert_frame_equal(rs, xp, check_names=False)

    def test_reset_index_right_dtype(self):
        time = np.arange(0.0, 10, np.sqrt(2) / 2)
        s1 = Series((9.81 * time ** 2) / 2,
                    index=Index(time, name='time'),
                    name='speed')
        df = DataFrame(s1)

        resetted = s1.reset_index()
        self.assertEqual(resetted['time'].dtype, np.float64)

        resetted = df.reset_index()
        self.assertEqual(resetted['time'].dtype, np.float64)

    def test_reset_index_multiindex_col(self):
        vals = np.random.randn(3, 3).astype(object)
        idx = ['x', 'y', 'z']
        full = np.hstack(([[x] for x in idx], vals))
        df = DataFrame(vals, Index(idx, name='a'),
                       columns=[['b', 'b', 'c'], ['mean', 'median', 'mean']])
        rs = df.reset_index()
        xp = DataFrame(full, columns=[['a', 'b', 'b', 'c'],
                                      ['', 'mean', 'median', 'mean']])
        assert_frame_equal(rs, xp)

        rs = df.reset_index(col_fill=None)
        xp = DataFrame(full, columns=[['a', 'b', 'b', 'c'],
                                      ['a', 'mean', 'median', 'mean']])
        assert_frame_equal(rs, xp)

        rs = df.reset_index(col_level=1, col_fill='blah')
        xp = DataFrame(full, columns=[['blah', 'b', 'b', 'c'],
                                      ['a', 'mean', 'median', 'mean']])
        assert_frame_equal(rs, xp)

        df = DataFrame(vals,
                       MultiIndex.from_arrays([[0, 1, 2], ['x', 'y', 'z']],
                                              names=['d', 'a']),
                       columns=[['b', 'b', 'c'], ['mean', 'median', 'mean']])
        rs = df.reset_index('a', )
        xp = DataFrame(full, Index([0, 1, 2], name='d'),
                       columns=[['a', 'b', 'b', 'c'],
                                ['', 'mean', 'median', 'mean']])
        assert_frame_equal(rs, xp)

        rs = df.reset_index('a', col_fill=None)
        xp = DataFrame(full, Index(lrange(3), name='d'),
                       columns=[['a', 'b', 'b', 'c'],
                                ['a', 'mean', 'median', 'mean']])
        assert_frame_equal(rs, xp)

        rs = df.reset_index('a', col_fill='blah', col_level=1)
        xp = DataFrame(full, Index(lrange(3), name='d'),
                       columns=[['blah', 'b', 'b', 'c'],
                                ['a', 'mean', 'median', 'mean']])
        assert_frame_equal(rs, xp)

    def test_reset_index_with_datetimeindex_cols(self):
        # GH5818
        #
        df = pd.DataFrame([[1, 2], [3, 4]],
                          columns=pd.date_range('1/1/2013', '1/2/2013'),
                          index=['A', 'B'])

        result = df.reset_index()
        expected = pd.DataFrame([['A', 1, 2], ['B', 3, 4]],
                          columns=['index', datetime(2013, 1, 1),
                                   datetime(2013, 1, 2)])
        assert_frame_equal(result, expected)

    #----------------------------------------------------------------------
    # Tests to cope with refactored internals
    def test_as_matrix_numeric_cols(self):
        self.frame['foo'] = 'bar'

        values = self.frame.as_matrix(['A', 'B', 'C', 'D'])
        self.assertEqual(values.dtype, np.float64)

    def test_as_matrix_lcd(self):

        # mixed lcd
        values = self.mixed_float.as_matrix(['A', 'B', 'C', 'D'])
        self.assertEqual(values.dtype, np.float64)

        values = self.mixed_float.as_matrix(['A', 'B', 'C' ])
        self.assertEqual(values.dtype, np.float32)

        values = self.mixed_float.as_matrix(['C'])
        self.assertEqual(values.dtype, np.float16)

        values = self.mixed_int.as_matrix(['A','B','C','D'])
        self.assertEqual(values.dtype, np.int64)

        values = self.mixed_int.as_matrix(['A','D'])
        self.assertEqual(values.dtype, np.int64)

        # guess all ints are cast to uints....
        values = self.mixed_int.as_matrix(['A','B','C'])
        self.assertEqual(values.dtype, np.int64)

        values = self.mixed_int.as_matrix(['A','C'])
        self.assertEqual(values.dtype, np.int32)

        values = self.mixed_int.as_matrix(['C','D'])
        self.assertEqual(values.dtype, np.int64)

        values = self.mixed_int.as_matrix(['A'])
        self.assertEqual(values.dtype, np.int32)

        values = self.mixed_int.as_matrix(['C'])
        self.assertEqual(values.dtype, np.uint8)

    def test_constructor_with_convert(self):
        # this is actually mostly a test of lib.maybe_convert_objects
        # #2845
        df = DataFrame({'A' : [2**63-1] })
        result = df['A']
        expected = Series(np.asarray([2**63-1], np.int64), name='A')
        assert_series_equal(result, expected)

        df = DataFrame({'A' : [2**63] })
        result = df['A']
        expected = Series(np.asarray([2**63], np.object_), name='A')
        assert_series_equal(result, expected)

        df = DataFrame({'A' : [datetime(2005, 1, 1), True] })
        result = df['A']
        expected = Series(np.asarray([datetime(2005, 1, 1), True], np.object_),
                          name='A')
        assert_series_equal(result, expected)

        df = DataFrame({'A' : [None, 1] })
        result = df['A']
        expected = Series(np.asarray([np.nan, 1], np.float_), name='A')
        assert_series_equal(result, expected)

        df = DataFrame({'A' : [1.0, 2] })
        result = df['A']
        expected = Series(np.asarray([1.0, 2], np.float_), name='A')
        assert_series_equal(result, expected)

        df = DataFrame({'A' : [1.0+2.0j, 3] })
        result = df['A']
        expected = Series(np.asarray([1.0+2.0j, 3], np.complex_), name='A')
        assert_series_equal(result, expected)

        df = DataFrame({'A' : [1.0+2.0j, 3.0] })
        result = df['A']
        expected = Series(np.asarray([1.0+2.0j, 3.0], np.complex_), name='A')
        assert_series_equal(result, expected)

        df = DataFrame({'A' : [1.0+2.0j, True] })
        result = df['A']
        expected = Series(np.asarray([1.0+2.0j, True], np.object_), name='A')
        assert_series_equal(result, expected)

        df = DataFrame({'A' : [1.0, None] })
        result = df['A']
        expected = Series(np.asarray([1.0, np.nan], np.float_), name='A')
        assert_series_equal(result, expected)

        df = DataFrame({'A' : [1.0+2.0j, None] })
        result = df['A']
        expected = Series(np.asarray([1.0+2.0j, np.nan], np.complex_), name='A')
        assert_series_equal(result, expected)

        df = DataFrame({'A' : [2.0, 1, True, None] })
        result = df['A']
        expected = Series(np.asarray([2.0, 1, True, None], np.object_), name='A')
        assert_series_equal(result, expected)

        df = DataFrame({'A' : [2.0, 1, datetime(2006, 1, 1), None] })
        result = df['A']
        expected = Series(np.asarray([2.0, 1, datetime(2006, 1, 1),
                                      None], np.object_), name='A')
        assert_series_equal(result, expected)

    def test_construction_with_mixed(self):
        # test construction edge cases with mixed types

        # f7u12, this does not work without extensive workaround
        data = [[datetime(2001, 1, 5), nan, datetime(2001, 1, 2)],
                [datetime(2000, 1, 2), datetime(2000, 1, 3),
                 datetime(2000, 1, 1)]]
        df = DataFrame(data)

        # check dtypes
        result = df.get_dtype_counts().sort_values()
        expected = Series({ 'datetime64[ns]' : 3 })

        # mixed-type frames
        self.mixed_frame['datetime'] = datetime.now()
        self.mixed_frame['timedelta'] = timedelta(days=1,seconds=1)
        self.assertEqual(self.mixed_frame['datetime'].dtype, 'M8[ns]')
        self.assertEqual(self.mixed_frame['timedelta'].dtype, 'm8[ns]')
        result = self.mixed_frame.get_dtype_counts().sort_values()
        expected = Series({ 'float64' : 4,
                            'object' : 1,
                            'datetime64[ns]' : 1,
                            'timedelta64[ns]' : 1}).sort_values()
        assert_series_equal(result,expected)

    def test_construction_with_conversions(self):

        # convert from a numpy array of non-ns timedelta64
        arr = np.array([1,2,3],dtype='timedelta64[s]')
        s = Series(arr)
        expected = Series(timedelta_range('00:00:01',periods=3,freq='s'))
        assert_series_equal(s,expected)

        df = DataFrame(index=range(3))
        df['A'] = arr
        expected = DataFrame({'A' : timedelta_range('00:00:01',periods=3,freq='s')},
                             index=range(3))
        assert_frame_equal(df,expected)

        # convert from a numpy array of non-ns datetime64
        #### note that creating a numpy datetime64 is in LOCAL time!!!!
        #### seems to work for M8[D], but not for M8[s]

        s = Series(np.array(['2013-01-01','2013-01-02','2013-01-03'],dtype='datetime64[D]'))
        assert_series_equal(s,Series(date_range('20130101',periods=3,freq='D')))
        #s = Series(np.array(['2013-01-01 00:00:01','2013-01-01 00:00:02','2013-01-01 00:00:03'],dtype='datetime64[s]'))
        #assert_series_equal(s,date_range('20130101 00:00:01',period=3,freq='s'))

        expected = DataFrame({
            'dt1' : Timestamp('20130101'),
            'dt2' : date_range('20130101',periods=3),
            #'dt3' : date_range('20130101 00:00:01',periods=3,freq='s'),
            },index=range(3))


        df = DataFrame(index=range(3))
        df['dt1'] = np.datetime64('2013-01-01')
        df['dt2'] = np.array(['2013-01-01','2013-01-02','2013-01-03'],dtype='datetime64[D]')
        #df['dt3'] = np.array(['2013-01-01 00:00:01','2013-01-01 00:00:02','2013-01-01 00:00:03'],dtype='datetime64[s]')
        assert_frame_equal(df, expected)

    def test_constructor_frame_copy(self):
        cop = DataFrame(self.frame, copy=True)
        cop['A'] = 5
        self.assertTrue((cop['A'] == 5).all())
        self.assertFalse((self.frame['A'] == 5).all())

    def test_constructor_ndarray_copy(self):
        df = DataFrame(self.frame.values)

        self.frame.values[5] = 5
        self.assertTrue((df.values[5] == 5).all())

        df = DataFrame(self.frame.values, copy=True)
        self.frame.values[6] = 6
        self.assertFalse((df.values[6] == 6).all())

    def test_constructor_series_copy(self):
        series = self.frame._series

        df = DataFrame({'A': series['A']})
        df['A'][:] = 5

        self.assertFalse((series['A'] == 5).all())

    def test_constructor_compound_dtypes(self):
        # GH 5191
        # compound dtypes should raise not-implementederror

        def f(dtype):
            return DataFrame(data = list(itertools.repeat((datetime(2001, 1, 1), "aa", 20), 9)),
                             columns=["A", "B", "C"], dtype=dtype)

        self.assertRaises(NotImplementedError, f, [("A","datetime64[h]"), ("B","str"), ("C","int32")])

        # these work (though results may be unexpected)
        f('int64')
        f('float64')

        # 10822
        # invalid error message on dt inference
        if not is_platform_windows():
            f('M8[ns]')

    def test_assign_columns(self):
        self.frame['hi'] = 'there'

        frame = self.frame.copy()
        frame.columns = ['foo', 'bar', 'baz', 'quux', 'foo2']
        assert_series_equal(self.frame['C'], frame['baz'], check_names=False)
        assert_series_equal(self.frame['hi'], frame['foo2'], check_names=False)

    def test_columns_with_dups(self):

        # GH 3468 related

        # basic
        df = DataFrame([[1,2]], columns=['a','a'])
        df.columns = ['a','a.1']
        str(df)
        expected = DataFrame([[1,2]], columns=['a','a.1'])
        assert_frame_equal(df, expected)

        df = DataFrame([[1,2,3]], columns=['b','a','a'])
        df.columns = ['b','a','a.1']
        str(df)
        expected = DataFrame([[1,2,3]], columns=['b','a','a.1'])
        assert_frame_equal(df, expected)

        # with a dup index
        df = DataFrame([[1,2]], columns=['a','a'])
        df.columns = ['b','b']
        str(df)
        expected = DataFrame([[1,2]], columns=['b','b'])
        assert_frame_equal(df, expected)

        # multi-dtype
        df = DataFrame([[1,2,1.,2.,3.,'foo','bar']], columns=['a','a','b','b','d','c','c'])
        df.columns = list('ABCDEFG')
        str(df)
        expected = DataFrame([[1,2,1.,2.,3.,'foo','bar']], columns=list('ABCDEFG'))
        assert_frame_equal(df, expected)

        # this is an error because we cannot disambiguate the dup columns
        self.assertRaises(Exception, lambda x: DataFrame([[1,2,'foo','bar']], columns=['a','a','a','a']))

        # dups across blocks
        df_float  = DataFrame(np.random.randn(10, 3),dtype='float64')
        df_int    = DataFrame(np.random.randn(10, 3),dtype='int64')
        df_bool   = DataFrame(True,index=df_float.index,columns=df_float.columns)
        df_object = DataFrame('foo',index=df_float.index,columns=df_float.columns)
        df_dt     = DataFrame(Timestamp('20010101'),index=df_float.index,columns=df_float.columns)
        df        = pd.concat([ df_float, df_int, df_bool, df_object, df_dt ], axis=1)

        self.assertEqual(len(df._data._blknos), len(df.columns))
        self.assertEqual(len(df._data._blklocs), len(df.columns))

        # testing iget
        for i in range(len(df.columns)):
            df.iloc[:,i]

        # dup columns across dtype GH 2079/2194
        vals = [[1, -1, 2.], [2, -2, 3.]]
        rs = DataFrame(vals, columns=['A', 'A', 'B'])
        xp = DataFrame(vals)
        xp.columns = ['A', 'A', 'B']
        assert_frame_equal(rs, xp)

    def test_insert_column_bug_4032(self):

        # GH4032, inserting a column and renaming causing errors
        df = DataFrame({'b': [1.1, 2.2]})
        df = df.rename(columns={})
        df.insert(0, 'a', [1, 2])

        result = df.rename(columns={})
        str(result)
        expected = DataFrame([[1,1.1],[2, 2.2]],columns=['a','b'])
        assert_frame_equal(result,expected)
        df.insert(0, 'c', [1.3, 2.3])

        result = df.rename(columns={})
        str(result)

        expected = DataFrame([[1.3,1,1.1],[2.3,2, 2.2]],columns=['c','a','b'])
        assert_frame_equal(result,expected)

    def test_cast_internals(self):
        casted = DataFrame(self.frame._data, dtype=int)
        expected = DataFrame(self.frame._series, dtype=int)
        assert_frame_equal(casted, expected)

        casted = DataFrame(self.frame._data, dtype=np.int32)
        expected = DataFrame(self.frame._series, dtype=np.int32)
        assert_frame_equal(casted, expected)

    def test_consolidate(self):
        self.frame['E'] = 7.
        consolidated = self.frame.consolidate()
        self.assertEqual(len(consolidated._data.blocks), 1)

        # Ensure copy, do I want this?
        recons = consolidated.consolidate()
        self.assertIsNot(recons, consolidated)
        assert_frame_equal(recons, consolidated)

        self.frame['F'] = 8.
        self.assertEqual(len(self.frame._data.blocks), 3)
        self.frame.consolidate(inplace=True)
        self.assertEqual(len(self.frame._data.blocks), 1)

    def test_consolidate_inplace(self):
        frame = self.frame.copy()

        # triggers in-place consolidation
        for letter in range(ord('A'), ord('Z')):
            self.frame[chr(letter)] = chr(letter)

    def test_as_matrix_consolidate(self):
        self.frame['E'] = 7.
        self.assertFalse(self.frame._data.is_consolidated())
        _ = self.frame.as_matrix()
        self.assertTrue(self.frame._data.is_consolidated())

    def test_modify_values(self):
        self.frame.values[5] = 5
        self.assertTrue((self.frame.values[5] == 5).all())

        # unconsolidated
        self.frame['E'] = 7.
        self.frame.values[6] = 6
        self.assertTrue((self.frame.values[6] == 6).all())

    def test_boolean_set_uncons(self):
        self.frame['E'] = 7.

        expected = self.frame.values.copy()
        expected[expected > 1] = 2

        self.frame[self.frame > 1] = 2
        assert_almost_equal(expected, self.frame.values)

    def test_xs_view(self):
        # in 0.14 this will return a view if possible a copy otherwise, but
        # this is numpy dependent

        dm = DataFrame(np.arange(20.).reshape(4, 5),
                       index=lrange(4), columns=lrange(5))

        dm.xs(2)[:] = 10
        self.assertTrue((dm.xs(2) == 10).all())

    def test_boolean_indexing(self):
        idx = lrange(3)
        cols = ['A','B','C']
        df1 = DataFrame(index=idx, columns=cols,
                        data=np.array([[0.0, 0.5, 1.0],
                                       [1.5, 2.0, 2.5],
                                       [3.0, 3.5, 4.0]],
                                      dtype=float))
        df2 = DataFrame(index=idx, columns=cols,
                        data=np.ones((len(idx), len(cols))))

        expected = DataFrame(index=idx, columns=cols,
                             data=np.array([[0.0, 0.5, 1.0],
                                            [1.5, 2.0, -1],
                                            [-1, -1, -1]], dtype=float))

        df1[df1 > 2.0 * df2] = -1
        assert_frame_equal(df1, expected)
        with assertRaisesRegexp(ValueError, 'Item wrong length'):
            df1[df1.index[:-1] > 2] = -1

    def test_boolean_indexing_mixed(self):
        df = DataFrame(
            {long(0): {35: np.nan, 40: np.nan, 43: np.nan, 49: np.nan, 50: np.nan},
             long(1): {35: np.nan,
                  40: 0.32632316859446198,
                  43: np.nan,
                  49: 0.32632316859446198,
                  50: 0.39114724480578139},
             long(2): {35: np.nan, 40: np.nan, 43: 0.29012581014105987, 49: np.nan, 50: np.nan},
             long(3): {35: np.nan, 40: np.nan, 43: np.nan, 49: np.nan, 50: np.nan},
             long(4): {35: 0.34215328467153283, 40: np.nan, 43: np.nan, 49: np.nan, 50: np.nan},
             'y': {35: 0, 40: 0, 43: 0, 49: 0, 50: 1}})

        # mixed int/float ok
        df2 = df.copy()
        df2[df2>0.3] = 1
        expected = df.copy()
        expected.loc[40,1] = 1
        expected.loc[49,1] = 1
        expected.loc[50,1] = 1
        expected.loc[35,4] = 1
        assert_frame_equal(df2,expected)

        df['foo'] = 'test'
        with tm.assertRaisesRegexp(TypeError, 'boolean setting on mixed-type'):
            df[df > 0.3] = 1

    def test_fillna_col_reordering(self):
        idx = lrange(20)
        cols = ["COL." + str(i) for i in range(5, 0, -1)]
        data = np.random.rand(20, 5)
        df = DataFrame(index=lrange(20), columns=cols, data=data)
        filled = df.fillna(method='ffill')
        self.assertEqual(df.columns.tolist(), filled.columns.tolist())

    def test_take(self):

        # homogeneous
        #----------------------------------------
        order  = [3, 1, 2, 0]
        for df in [self.frame]:

            result = df.take(order, axis=0)
            expected = df.reindex(df.index.take(order))
            assert_frame_equal(result, expected)

            # axis = 1
            result = df.take(order, axis=1)
            expected = df.ix[:, ['D', 'B', 'C', 'A']]
            assert_frame_equal(result, expected, check_names=False)

        # neg indicies
        order = [2,1,-1]
        for df in [self.frame]:

            result = df.take(order, axis=0)
            expected = df.reindex(df.index.take(order))
            assert_frame_equal(result, expected)

            # axis = 1
            result = df.take(order, axis=1)
            expected = df.ix[:, ['C', 'B', 'D']]
            assert_frame_equal(result, expected, check_names=False)

        # illegal indices
        self.assertRaises(IndexError, df.take, [3,1,2,30], axis=0)
        self.assertRaises(IndexError, df.take, [3,1,2,-31], axis=0)
        self.assertRaises(IndexError, df.take, [3,1,2,5], axis=1)
        self.assertRaises(IndexError, df.take, [3,1,2,-5], axis=1)

        # mixed-dtype
        #----------------------------------------
        order  = [4, 1, 2, 0, 3]
        for df in [self.mixed_frame]:

            result = df.take(order, axis=0)
            expected = df.reindex(df.index.take(order))
            assert_frame_equal(result, expected)

            # axis = 1
            result = df.take(order, axis=1)
            expected = df.ix[:, ['foo', 'B', 'C', 'A', 'D']]
            assert_frame_equal(result, expected)

        # neg indicies
        order = [4,1,-2]
        for df in [self.mixed_frame]:

            result = df.take(order, axis=0)
            expected = df.reindex(df.index.take(order))
            assert_frame_equal(result, expected)

            # axis = 1
            result = df.take(order, axis=1)
            expected = df.ix[:, ['foo', 'B', 'D']]
            assert_frame_equal(result, expected)

        # by dtype
        order = [1, 2, 0, 3]
        for df in [self.mixed_float,self.mixed_int]:

            result = df.take(order, axis=0)
            expected = df.reindex(df.index.take(order))
            assert_frame_equal(result, expected)

            # axis = 1
            result = df.take(order, axis=1)
            expected = df.ix[:, ['B', 'C', 'A', 'D']]
            assert_frame_equal(result, expected)

    def test_iterkv_deprecation(self):
        with tm.assert_produces_warning(FutureWarning):
            self.mixed_float.iterkv()

    def test_iterkv_names(self):
        for k, v in compat.iteritems(self.mixed_frame):
            self.assertEqual(v.name, k)

    def test_series_put_names(self):
        series = self.mixed_frame._series
        for k, v in compat.iteritems(series):
            self.assertEqual(v.name, k)

    def test_dot(self):
        a = DataFrame(np.random.randn(3, 4), index=['a', 'b', 'c'],
                      columns=['p', 'q', 'r', 's'])
        b = DataFrame(np.random.randn(4, 2), index=['p', 'q', 'r', 's'],
                      columns=['one', 'two'])

        result = a.dot(b)
        expected = DataFrame(np.dot(a.values, b.values),
                             index=['a', 'b', 'c'],
                             columns=['one', 'two'])
        # Check alignment
        b1 = b.reindex(index=reversed(b.index))
        result = a.dot(b)
        assert_frame_equal(result, expected)

        # Check series argument
        result = a.dot(b['one'])
        assert_series_equal(result, expected['one'], check_names=False)
        self.assertTrue(result.name is None)

        result = a.dot(b1['one'])
        assert_series_equal(result, expected['one'], check_names=False)
        self.assertTrue(result.name is None)

        # can pass correct-length arrays
        row = a.ix[0].values

        result = a.dot(row)
        exp = a.dot(a.ix[0])
        assert_series_equal(result, exp)

        with assertRaisesRegexp(ValueError, 'Dot product shape mismatch'):
            a.dot(row[:-1])

        a = np.random.rand(1, 5)
        b = np.random.rand(5, 1)
        A = DataFrame(a)
        B = DataFrame(b)

        # it works
        result = A.dot(b)

        # unaligned
        df = DataFrame(randn(3, 4), index=[1, 2, 3], columns=lrange(4))
        df2 = DataFrame(randn(5, 3), index=lrange(5), columns=[1, 2, 3])

        assertRaisesRegexp(ValueError, 'aligned', df.dot, df2)

    def test_stale_cached_series_bug_473(self):

        # this is chained, but ok
        with option_context('chained_assignment',None):
            Y = DataFrame(np.random.random((4, 4)), index=('a', 'b', 'c', 'd'),
                          columns=('e', 'f', 'g', 'h'))
            repr(Y)
            Y['e'] = Y['e'].astype('object')
            Y['g']['c'] = np.NaN
            repr(Y)
            result = Y.sum()
            exp = Y['g'].sum()
            self.assertTrue(isnull(Y['g']['c']))

    def test_index_namedtuple(self):
        from collections import namedtuple
        IndexType = namedtuple("IndexType", ["a", "b"])
        idx1 = IndexType("foo", "bar")
        idx2 = IndexType("baz", "bof")
        index = Index([idx1, idx2],
                      name="composite_index", tupleize_cols=False)
        df = DataFrame([(1, 2), (3, 4)], index=index, columns=["A", "B"])
        result = df.ix[IndexType("foo", "bar")]["A"]
        self.assertEqual(result, 1)

    def test_empty_nonzero(self):
        df = DataFrame([1, 2, 3])
        self.assertFalse(df.empty)
        df = DataFrame(index=['a', 'b'], columns=['c', 'd']).dropna()
        self.assertTrue(df.empty)
        self.assertTrue(df.T.empty)

    def test_consolidate_datetime64(self):
        # numpy vstack bug

        data = """\
starting,ending,measure
2012-06-21 00:00,2012-06-23 07:00,77
2012-06-23 07:00,2012-06-23 16:30,65
2012-06-23 16:30,2012-06-25 08:00,77
2012-06-25 08:00,2012-06-26 12:00,0
2012-06-26 12:00,2012-06-27 08:00,77
"""
        df = read_csv(StringIO(data), parse_dates=[0, 1])

        ser_starting = df.starting
        ser_starting.index = ser_starting.values
        ser_starting = ser_starting.tz_localize('US/Eastern')
        ser_starting = ser_starting.tz_convert('UTC')

        ser_ending = df.ending
        ser_ending.index = ser_ending.values
        ser_ending = ser_ending.tz_localize('US/Eastern')
        ser_ending = ser_ending.tz_convert('UTC')

        df.starting = ser_starting.index
        df.ending = ser_ending.index

        tm.assert_index_equal(pd.DatetimeIndex(df.starting), ser_starting.index)
        tm.assert_index_equal(pd.DatetimeIndex(df.ending), ser_ending.index)

    def test_strange_column_corruption_issue(self):

        df = DataFrame(index=[0, 1])
        df[0] = nan
        wasCol = {}
        # uncommenting these makes the results match
        # for col in xrange(100, 200):
        #    wasCol[col] = 1
        #    df[col] = nan

        for i, dt in enumerate(df.index):
            for col in range(100, 200):
                if not col in wasCol:
                    wasCol[col] = 1
                    df[col] = nan
                df[col][dt] = i

        myid = 100

        first = len(df.ix[isnull(df[myid]), [myid]])
        second = len(df.ix[isnull(df[myid]), [myid]])
        self.assertTrue(first == second == 0)

    def test_inplace_return_self(self):
        # re #1893

        data = DataFrame({'a': ['foo', 'bar', 'baz', 'qux'],
                          'b': [0, 0, 1, 1],
                          'c': [1, 2, 3, 4]})

        def _check_f(base, f):
            result = f(base)
            self.assertTrue(result is None)

        # -----DataFrame-----

        # set_index
        f = lambda x: x.set_index('a', inplace=True)
        _check_f(data.copy(), f)

        # reset_index
        f = lambda x: x.reset_index(inplace=True)
        _check_f(data.set_index('a'), f)

        # drop_duplicates
        f = lambda x: x.drop_duplicates(inplace=True)
        _check_f(data.copy(), f)

        # sort
        f = lambda x: x.sort_values('b', inplace=True)
        _check_f(data.copy(), f)

        # sort_index
        f = lambda x: x.sort_index(inplace=True)
        _check_f(data.copy(), f)

        # sortlevel
        f = lambda x: x.sortlevel(0, inplace=True)
        _check_f(data.set_index(['a', 'b']), f)

        # fillna
        f = lambda x: x.fillna(0, inplace=True)
        _check_f(data.copy(), f)

        # replace
        f = lambda x: x.replace(1, 0, inplace=True)
        _check_f(data.copy(), f)

        # rename
        f = lambda x: x.rename({1: 'foo'}, inplace=True)
        _check_f(data.copy(), f)

        # -----Series-----
        d = data.copy()['c']

        # reset_index
        f = lambda x: x.reset_index(inplace=True, drop=True)
        _check_f(data.set_index('a')['c'], f)

        # fillna
        f = lambda x: x.fillna(0, inplace=True)
        _check_f(d.copy(), f)

        # replace
        f = lambda x: x.replace(1, 0, inplace=True)
        _check_f(d.copy(), f)

        # rename
        f = lambda x: x.rename({1: 'foo'}, inplace=True)
        _check_f(d.copy(), f)

    def test_isin(self):
        # GH #4211
        df = DataFrame({'vals': [1, 2, 3, 4], 'ids': ['a', 'b', 'f', 'n'],
                        'ids2': ['a', 'n', 'c', 'n']},
                        index=['foo', 'bar', 'baz', 'qux'])
        other = ['a', 'b', 'c']

        result = df.isin(other)
        expected = DataFrame([df.loc[s].isin(other) for s in df.index])
        assert_frame_equal(result, expected)

    def test_isin_empty(self):
        df = DataFrame({'A': ['a', 'b', 'c'], 'B': ['a', 'e', 'f']})
        result = df.isin([])
        expected = pd.DataFrame(False, df.index, df.columns)
        assert_frame_equal(result, expected)

    def test_isin_dict(self):
        df = DataFrame({'A': ['a', 'b', 'c'], 'B': ['a', 'e', 'f']})
        d = {'A': ['a']}

        expected = DataFrame(False, df.index, df.columns)
        expected.loc[0, 'A'] = True

        result = df.isin(d)
        assert_frame_equal(result, expected)

        # non unique columns
        df = DataFrame({'A': ['a', 'b', 'c'], 'B': ['a', 'e', 'f']})
        df.columns = ['A', 'A']
        expected = DataFrame(False, df.index, df.columns)
        expected.loc[0, 'A'] = True
        result = df.isin(d)
        assert_frame_equal(result, expected)

    def test_isin_with_string_scalar(self):
        #GH4763
        df = DataFrame({'vals': [1, 2, 3, 4], 'ids': ['a', 'b', 'f', 'n'],
                        'ids2': ['a', 'n', 'c', 'n']},
                        index=['foo', 'bar', 'baz', 'qux'])
        with tm.assertRaises(TypeError):
            df.isin('a')

        with tm.assertRaises(TypeError):
            df.isin('aaa')

    def test_isin_df(self):
        df1 = DataFrame({'A': [1, 2, 3, 4], 'B': [2, np.nan, 4, 4]})
        df2 = DataFrame({'A': [0, 2, 12, 4], 'B': [2, np.nan, 4, 5]})
        expected = DataFrame(False, df1.index, df1.columns)
        result = df1.isin(df2)
        expected['A'].loc[[1, 3]] = True
        expected['B'].loc[[0, 2]] = True
        assert_frame_equal(result, expected)

        # partial overlapping columns
        df2.columns = ['A', 'C']
        result = df1.isin(df2)
        expected['B'] = False
        assert_frame_equal(result, expected)

    def test_isin_df_dupe_values(self):
        df1 = DataFrame({'A': [1, 2, 3, 4], 'B': [2, np.nan, 4, 4]})
        # just cols duped
        df2 = DataFrame([[0, 2], [12, 4], [2, np.nan], [4, 5]],
                        columns=['B', 'B'])
        with tm.assertRaises(ValueError):
            df1.isin(df2)

        # just index duped
        df2 = DataFrame([[0, 2], [12, 4], [2, np.nan], [4, 5]],
                        columns=['A', 'B'], index=[0, 0, 1, 1])
        with tm.assertRaises(ValueError):
            df1.isin(df2)

        # cols and index:
        df2.columns = ['B', 'B']
        with tm.assertRaises(ValueError):
            df1.isin(df2)

    def test_isin_dupe_self(self):
        other = DataFrame({'A': [1, 0, 1, 0], 'B': [1, 1, 0, 0]})
        df = DataFrame([[1, 1], [1, 0], [0, 0]], columns=['A','A'])
        result = df.isin(other)
        expected = DataFrame(False, index=df.index, columns=df.columns)
        expected.loc[0] = True
        expected.iloc[1, 1] = True
        assert_frame_equal(result, expected)

    def test_isin_against_series(self):
        df = pd.DataFrame({'A': [1, 2, 3, 4], 'B': [2, np.nan, 4, 4]},
                          index=['a', 'b', 'c', 'd'])
        s = pd.Series([1, 3, 11, 4], index=['a', 'b', 'c', 'd'])
        expected = DataFrame(False, index=df.index, columns=df.columns)
        expected['A'].loc['a'] = True
        expected.loc['d'] = True
        result = df.isin(s)
        assert_frame_equal(result, expected)

    def test_isin_multiIndex(self):
        idx = MultiIndex.from_tuples([(0, 'a', 'foo'), (0, 'a', 'bar'),
                                      (0, 'b', 'bar'), (0, 'b', 'baz'),
                                      (2, 'a', 'foo'), (2, 'a', 'bar'),
                                      (2, 'c', 'bar'), (2, 'c', 'baz'),
                                      (1, 'b', 'foo'), (1, 'b', 'bar'),
                                      (1, 'c', 'bar'), (1, 'c', 'baz')])
        df1 = DataFrame({'A': np.ones(12),
                         'B': np.zeros(12)}, index=idx)
        df2 = DataFrame({'A': [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
                         'B': [1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1]})
        # against regular index
        expected = DataFrame(False, index=df1.index, columns=df1.columns)
        result = df1.isin(df2)
        assert_frame_equal(result, expected)

        df2.index = idx
        expected = df2.values.astype(np.bool)
        expected[:, 1] = ~expected[:, 1]
        expected = DataFrame(expected, columns=['A', 'B'], index=idx)

        result = df1.isin(df2)
        assert_frame_equal(result, expected)

    def test_concat_empty_dataframe_dtypes(self):
        df = DataFrame(columns=list("abc"))
        df['a'] = df['a'].astype(np.bool_)
        df['b'] = df['b'].astype(np.int32)
        df['c'] = df['c'].astype(np.float64)

        result = pd.concat([df, df])
        self.assertEqual(result['a'].dtype, np.bool_)
        self.assertEqual(result['b'].dtype, np.int32)
        self.assertEqual(result['c'].dtype, np.float64)

        result = pd.concat([df, df.astype(np.float64)])
        self.assertEqual(result['a'].dtype, np.object_)
        self.assertEqual(result['b'].dtype, np.float64)
        self.assertEqual(result['c'].dtype, np.float64)

    def test_empty_frame_dtypes_ftypes(self):
        empty_df = pd.DataFrame()
        assert_series_equal(empty_df.dtypes, pd.Series(dtype=np.object))
        assert_series_equal(empty_df.ftypes, pd.Series(dtype=np.object))

        nocols_df = pd.DataFrame(index=[1,2,3])
        assert_series_equal(nocols_df.dtypes, pd.Series(dtype=np.object))
        assert_series_equal(nocols_df.ftypes, pd.Series(dtype=np.object))

        norows_df = pd.DataFrame(columns=list("abc"))
        assert_series_equal(norows_df.dtypes, pd.Series(np.object, index=list("abc")))
        assert_series_equal(norows_df.ftypes, pd.Series('object:dense', index=list("abc")))

        norows_int_df = pd.DataFrame(columns=list("abc")).astype(np.int32)
        assert_series_equal(norows_int_df.dtypes, pd.Series(np.dtype('int32'), index=list("abc")))
        assert_series_equal(norows_int_df.ftypes, pd.Series('int32:dense', index=list("abc")))

        odict = OrderedDict
        df = pd.DataFrame(odict([('a', 1), ('b', True), ('c', 1.0)]), index=[1, 2, 3])
        assert_series_equal(df.dtypes, pd.Series(odict([('a', np.int64),
                                                        ('b', np.bool),
                                                        ('c', np.float64)])))
        assert_series_equal(df.ftypes, pd.Series(odict([('a', 'int64:dense'),
                                                        ('b', 'bool:dense'),
                                                        ('c', 'float64:dense')])))

        # same but for empty slice of df
        assert_series_equal(df[:0].dtypes, pd.Series(odict([('a', np.int64),
                                                            ('b', np.bool),
                                                            ('c', np.float64)])))
        assert_series_equal(df[:0].ftypes, pd.Series(odict([('a', 'int64:dense'),
                                                            ('b', 'bool:dense'),
                                                            ('c', 'float64:dense')])))

    def test_dtypes_are_correct_after_column_slice(self):
        # GH6525
        df = pd.DataFrame(index=range(5), columns=list("abc"), dtype=np.float_)
        odict = OrderedDict
        assert_series_equal(df.dtypes,
                            pd.Series(odict([('a', np.float_), ('b', np.float_),
                                             ('c', np.float_),])))
        assert_series_equal(df.iloc[:,2:].dtypes,
                            pd.Series(odict([('c', np.float_)])))
        assert_series_equal(df.dtypes,
                            pd.Series(odict([('a', np.float_), ('b', np.float_),
                                             ('c', np.float_),])))

    def test_set_index_names(self):
        df = pd.util.testing.makeDataFrame()
        df.index.name = 'name'

        self.assertEqual(df.set_index(df.index).index.names, ['name'])

        mi = MultiIndex.from_arrays(df[['A', 'B']].T.values, names=['A', 'B'])
        mi2 = MultiIndex.from_arrays(df[['A', 'B', 'A', 'B']].T.values,
                                     names=['A', 'B', 'A', 'B'])

        df = df.set_index(['A', 'B'])

        self.assertEqual(df.set_index(df.index).index.names, ['A', 'B'])

        # Check that set_index isn't converting a MultiIndex into an Index
        self.assertTrue(isinstance(df.set_index(df.index).index, MultiIndex))

        # Check actual equality
        tm.assert_index_equal(df.set_index(df.index).index, mi)

        # Check that [MultiIndex, MultiIndex] yields a MultiIndex rather
        # than a pair of tuples
        self.assertTrue(isinstance(df.set_index([df.index, df.index]).index, MultiIndex))

        # Check equality
        tm.assert_index_equal(df.set_index([df.index, df.index]).index, mi2)

    def test_select_dtypes_include(self):
        df = DataFrame({'a': list('abc'),
                        'b': list(range(1, 4)),
                        'c': np.arange(3, 6).astype('u1'),
                        'd': np.arange(4.0, 7.0, dtype='float64'),
                        'e': [True, False, True],
                        'f': pd.Categorical(list('abc'))})
        ri = df.select_dtypes(include=[np.number])
        ei = df[['b', 'c', 'd']]
        assert_frame_equal(ri, ei)

        ri = df.select_dtypes(include=[np.number, 'category'])
        ei = df[['b', 'c', 'd', 'f']]
        assert_frame_equal(ri, ei)

    def test_select_dtypes_exclude(self):
        df = DataFrame({'a': list('abc'),
                        'b': list(range(1, 4)),
                        'c': np.arange(3, 6).astype('u1'),
                        'd': np.arange(4.0, 7.0, dtype='float64'),
                        'e': [True, False, True]})
        re = df.select_dtypes(exclude=[np.number])
        ee = df[['a', 'e']]
        assert_frame_equal(re, ee)

    def test_select_dtypes_exclude_include(self):
        df = DataFrame({'a': list('abc'),
                        'b': list(range(1, 4)),
                        'c': np.arange(3, 6).astype('u1'),
                        'd': np.arange(4.0, 7.0, dtype='float64'),
                        'e': [True, False, True],
                        'f': pd.date_range('now', periods=3).values})
        exclude = np.datetime64,
        include = np.bool_, 'integer'
        r = df.select_dtypes(include=include, exclude=exclude)
        e = df[['b', 'c', 'e']]
        assert_frame_equal(r, e)

        exclude = 'datetime',
        include = 'bool', 'int64', 'int32'
        r = df.select_dtypes(include=include, exclude=exclude)
        e = df[['b', 'e']]
        assert_frame_equal(r, e)

    def test_select_dtypes_not_an_attr_but_still_valid_dtype(self):
        df = DataFrame({'a': list('abc'),
                        'b': list(range(1, 4)),
                        'c': np.arange(3, 6).astype('u1'),
                        'd': np.arange(4.0, 7.0, dtype='float64'),
                        'e': [True, False, True],
                        'f': pd.date_range('now', periods=3).values})
        df['g'] = df.f.diff()
        assert not hasattr(np, 'u8')
        r = df.select_dtypes(include=['i8', 'O'], exclude=['timedelta'])
        e = df[['a', 'b']]
        assert_frame_equal(r, e)

        r = df.select_dtypes(include=['i8', 'O', 'timedelta64[ns]'])
        e = df[['a', 'b', 'g']]
        assert_frame_equal(r, e)

    def test_select_dtypes_empty(self):
        df = DataFrame({'a': list('abc'), 'b': list(range(1, 4))})
        with tm.assertRaisesRegexp(ValueError, 'at least one of include or '
                                   'exclude must be nonempty'):
            df.select_dtypes()

    def test_select_dtypes_raises_on_string(self):
        df = DataFrame({'a': list('abc'), 'b': list(range(1, 4))})
        with tm.assertRaisesRegexp(TypeError, 'include and exclude .+ non-'):
            df.select_dtypes(include='object')
        with tm.assertRaisesRegexp(TypeError, 'include and exclude .+ non-'):
            df.select_dtypes(exclude='object')
        with tm.assertRaisesRegexp(TypeError, 'include and exclude .+ non-'):
            df.select_dtypes(include=int, exclude='object')

    def test_select_dtypes_bad_datetime64(self):
        df = DataFrame({'a': list('abc'),
                        'b': list(range(1, 4)),
                        'c': np.arange(3, 6).astype('u1'),
                        'd': np.arange(4.0, 7.0, dtype='float64'),
                        'e': [True, False, True],
                        'f': pd.date_range('now', periods=3).values})
        with tm.assertRaisesRegexp(ValueError, '.+ is too specific'):
            df.select_dtypes(include=['datetime64[D]'])

        with tm.assertRaisesRegexp(ValueError, '.+ is too specific'):
            df.select_dtypes(exclude=['datetime64[as]'])

    def test_select_dtypes_str_raises(self):
        df = DataFrame({'a': list('abc'),
                        'g': list(u('abc')),
                        'b': list(range(1, 4)),
                        'c': np.arange(3, 6).astype('u1'),
                        'd': np.arange(4.0, 7.0, dtype='float64'),
                        'e': [True, False, True],
                        'f': pd.date_range('now', periods=3).values})
        string_dtypes = set((str, 'str', np.string_, 'S1',
                             'unicode', np.unicode_, 'U1'))
        try:
            string_dtypes.add(unicode)
        except NameError:
            pass
        for dt in string_dtypes:
            with tm.assertRaisesRegexp(TypeError,
                                       'string dtypes are not allowed'):
                df.select_dtypes(include=[dt])
            with tm.assertRaisesRegexp(TypeError,
                                       'string dtypes are not allowed'):
                df.select_dtypes(exclude=[dt])

    def test_select_dtypes_bad_arg_raises(self):
        df = DataFrame({'a': list('abc'),
                        'g': list(u('abc')),
                        'b': list(range(1, 4)),
                        'c': np.arange(3, 6).astype('u1'),
                        'd': np.arange(4.0, 7.0, dtype='float64'),
                        'e': [True, False, True],
                        'f': pd.date_range('now', periods=3).values})
        with tm.assertRaisesRegexp(TypeError, 'data type.*not understood'):
            df.select_dtypes(['blargy, blarg, blarg'])

    def test_select_dtypes_typecodes(self):
        # GH 11990
        df = mkdf(30, 3, data_gen_f=lambda x, y: np.random.random())
        expected = df
        FLOAT_TYPES = list(np.typecodes['AllFloat'])
        assert_frame_equal(df.select_dtypes(FLOAT_TYPES), expected)

    def test_assign(self):
        df = DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})
        original = df.copy()
        result = df.assign(C=df.B / df.A)
        expected = df.copy()
        expected['C'] = [4, 2.5, 2]
        assert_frame_equal(result, expected)

        # lambda syntax
        result = df.assign(C=lambda x: x.B / x.A)
        assert_frame_equal(result, expected)

        # original is unmodified
        assert_frame_equal(df, original)

        # Non-Series array-like
        result = df.assign(C=[4, 2.5, 2])
        assert_frame_equal(result, expected)
        # original is unmodified
        assert_frame_equal(df, original)

        result = df.assign(B=df.B / df.A)
        expected = expected.drop('B', axis=1).rename(columns={'C': 'B'})
        assert_frame_equal(result, expected)

        # overwrite
        result = df.assign(A=df.A + df.B)
        expected = df.copy()
        expected['A'] = [5, 7, 9]
        assert_frame_equal(result, expected)

        # lambda
        result = df.assign(A=lambda x: x.A + x.B)
        assert_frame_equal(result, expected)

    def test_assign_multiple(self):
        df = DataFrame([[1, 4], [2, 5], [3, 6]], columns=['A', 'B'])
        result = df.assign(C=[7, 8, 9], D=df.A, E=lambda x: x.B)
        expected = DataFrame([[1, 4, 7, 1, 4], [2, 5, 8, 2, 5],
                              [3, 6, 9, 3, 6]], columns=list('ABCDE'))
        assert_frame_equal(result, expected)

    def test_assign_alphabetical(self):
        # GH 9818
        df = DataFrame([[1, 2], [3, 4]], columns=['A', 'B'])
        result = df.assign(D=df.A + df.B, C=df.A - df.B)
        expected = DataFrame([[1, 2, -1, 3], [3, 4, -1, 7]],
                             columns=list('ABCD'))
        assert_frame_equal(result, expected)
        result = df.assign(C=df.A - df.B, D=df.A + df.B)
        assert_frame_equal(result, expected)

    def test_assign_bad(self):
        df = DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})
        # non-keyword argument
        with tm.assertRaises(TypeError):
            df.assign(lambda x: x.A)
        with tm.assertRaises(AttributeError):
            df.assign(C=df.A, D=df.A + df.C)
        with tm.assertRaises(KeyError):
            df.assign(C=lambda df: df.A, D=lambda df: df['A'] + df['C'])
        with tm.assertRaises(KeyError):
            df.assign(C=df.A, D=lambda x: x['A'] + x['C'])

    def test_dataframe_metadata(self):

        df = SubclassedDataFrame({'X': [1, 2, 3], 'Y': [1, 2, 3]},
                                 index=['a', 'b', 'c'])
        df.testattr = 'XXX'

        self.assertEqual(df.testattr, 'XXX')
        self.assertEqual(df[['X']].testattr, 'XXX')
        self.assertEqual(df.loc[['a', 'b'], :].testattr, 'XXX')
        self.assertEqual(df.iloc[[0, 1], :].testattr, 'XXX')

        # GH9776
        self.assertEqual(df.iloc[0:1, :].testattr, 'XXX')

        # GH10553
        unpickled = self.round_trip_pickle(df)
        assert_frame_equal(df, unpickled)
        self.assertEqual(df._metadata, unpickled._metadata)
        self.assertEqual(df.testattr, unpickled.testattr)

    def test_to_panel_expanddim(self):
        # GH 9762

        class SubclassedFrame(DataFrame):
            @property
            def _constructor_expanddim(self):
                return SubclassedPanel

        class SubclassedPanel(Panel):
            pass

        index = MultiIndex.from_tuples([(0, 0), (0, 1), (0, 2)])
        df = SubclassedFrame({'X':[1, 2, 3], 'Y': [4, 5, 6]}, index=index)
        result = df.to_panel()
        self.assertTrue(isinstance(result, SubclassedPanel))
        expected = SubclassedPanel([[[1, 2, 3]], [[4, 5, 6]]],
                                   items=['X', 'Y'], major_axis=[0],
                                   minor_axis=[0, 1, 2],
                                   dtype='int64')
        tm.assert_panel_equal(result, expected)

    def test_subclass_attr_err_propagation(self):
        # GH 11808
        class A(DataFrame):

            @property
            def bar(self):
                return self.i_dont_exist
        with tm.assertRaisesRegexp(AttributeError, '.*i_dont_exist.*'):
            A().bar


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
