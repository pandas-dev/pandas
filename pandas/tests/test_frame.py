# -*- coding: utf-8 -*-

from __future__ import print_function
# pylint: disable-msg=W0612,E1101
from copy import deepcopy
from datetime import datetime, timedelta, time, date
import sys
import operator
import re
import csv
import nose
import functools
import itertools
from itertools import product, permutations
from distutils.version import LooseVersion

from pandas.compat import(
    map, zip, range, long, lrange, lmap, lzip,
    OrderedDict, u, StringIO, string_types,
    is_platform_windows
)
from pandas import compat

from numpy import random, nan, inf
from numpy.random import randn
import numpy as np
import numpy.ma as ma
import numpy.ma.mrecords as mrecords

import pandas.core.nanops as nanops
import pandas.core.common as com
import pandas.core.format as fmt
import pandas.core.datetools as datetools
from pandas import (DataFrame, Index, Series, Panel, notnull, isnull,
                    MultiIndex, DatetimeIndex, Timestamp, date_range,
                    read_csv, timedelta_range, Timedelta, CategoricalIndex,
                    option_context, period_range)
from pandas.core.dtypes import DatetimeTZDtype
import pandas as pd
from pandas.parser import CParserError
from pandas.util.misc import is_little_endian

from pandas.util.testing import (assert_almost_equal,
                                 assert_numpy_array_equal,
                                 assert_series_equal,
                                 assert_frame_equal,
                                 assertRaisesRegexp,
                                 assertRaises,
                                 makeCustomDataframe as mkdf,
                                 ensure_clean,
                                 SubclassedDataFrame)
from pandas.core.indexing import IndexingError
from pandas.core.common import PandasError

import pandas.util.testing as tm
import pandas.lib as lib

from numpy.testing.decorators import slow

#---------------------------------------------------------------------
# DataFrame test cases

JOIN_TYPES = ['inner', 'outer', 'left', 'right']
MIXED_FLOAT_DTYPES = ['float16','float32','float64']
MIXED_INT_DTYPES   = ['uint8','uint16','uint32','uint64','int8','int16',
                      'int32','int64']

def _check_mixed_float(df, dtype = None):

    # float16 are most likely to be upcasted to float32
    dtypes = dict(A = 'float32', B = 'float32', C = 'float16', D = 'float64')
    if isinstance(dtype, compat.string_types):
        dtypes = dict([ (k,dtype) for k, v in dtypes.items() ])
    elif isinstance(dtype, dict):
        dtypes.update(dtype)
    if dtypes.get('A'):
        assert(df.dtypes['A'] == dtypes['A'])
    if dtypes.get('B'):
        assert(df.dtypes['B'] == dtypes['B'])
    if dtypes.get('C'):
        assert(df.dtypes['C'] == dtypes['C'])
    if dtypes.get('D'):
        assert(df.dtypes['D'] == dtypes['D'])


def _check_mixed_int(df, dtype = None):
    dtypes = dict(A = 'int32', B = 'uint64', C = 'uint8', D = 'int64')
    if isinstance(dtype, compat.string_types):
        dtypes = dict([ (k,dtype) for k, v in dtypes.items() ])
    elif isinstance(dtype, dict):
        dtypes.update(dtype)
    if dtypes.get('A'):
        assert(df.dtypes['A'] == dtypes['A'])
    if dtypes.get('B'):
        assert(df.dtypes['B'] == dtypes['B'])
    if dtypes.get('C'):
        assert(df.dtypes['C'] == dtypes['C'])
    if dtypes.get('D'):
        assert(df.dtypes['D'] == dtypes['D'])


class CheckIndexing(object):

    _multiprocess_can_split_ = True

    def test_getitem(self):
        # slicing
        sl = self.frame[:20]
        self.assertEqual(20, len(sl.index))

        # column access

        for _, series in compat.iteritems(sl):
            self.assertEqual(20, len(series.index))
            self.assertTrue(tm.equalContents(series.index, sl.index))

        for key, _ in compat.iteritems(self.frame._series):
            self.assertIsNotNone(self.frame[key])

        self.assertNotIn('random', self.frame)
        with assertRaisesRegexp(KeyError, 'random'):
            self.frame['random']

        df = self.frame.copy()
        df['$10'] = randn(len(df))
        ad = randn(len(df))
        df['@awesome_domain'] = ad
        self.assertRaises(KeyError, df.__getitem__, 'df["$10"]')
        res = df['@awesome_domain']
        assert_numpy_array_equal(ad, res.values)

    def test_getitem_dupe_cols(self):
        df = DataFrame([[1, 2, 3], [4, 5, 6]], columns=['a', 'a', 'b'])
        try:
            df[['baf']]
        except KeyError:
            pass
        else:
            self.fail("Dataframe failed to raise KeyError")

    def test_get(self):
        b = self.frame.get('B')
        assert_series_equal(b, self.frame['B'])

        self.assertIsNone(self.frame.get('foo'))
        assert_series_equal(self.frame.get('foo', self.frame['B']),
                            self.frame['B'])
        # None
        # GH 5652
        for df in [DataFrame(), DataFrame(columns=list('AB')), DataFrame(columns=list('AB'),index=range(3)) ]:
            result = df.get(None)
            self.assertIsNone(result)

    def test_getitem_iterator(self):
        idx = iter(['A', 'B', 'C'])
        result = self.frame.ix[:, idx]
        expected = self.frame.ix[:, ['A', 'B', 'C']]
        assert_frame_equal(result, expected)

    def test_getitem_list(self):
        self.frame.columns.name = 'foo'

        result = self.frame[['B', 'A']]
        result2 = self.frame[Index(['B', 'A'])]

        expected = self.frame.ix[:, ['B', 'A']]
        expected.columns.name = 'foo'

        assert_frame_equal(result, expected)
        assert_frame_equal(result2, expected)

        self.assertEqual(result.columns.name, 'foo')

        with assertRaisesRegexp(KeyError, 'not in index'):
            self.frame[['B', 'A', 'food']]
        with assertRaisesRegexp(KeyError, 'not in index'):
            self.frame[Index(['B', 'A', 'foo'])]

        # tuples
        df = DataFrame(randn(8, 3),
                       columns=Index([('foo', 'bar'), ('baz', 'qux'),
                                      ('peek', 'aboo')], name=['sth', 'sth2']))

        result = df[[('foo', 'bar'), ('baz', 'qux')]]
        expected = df.ix[:, :2]
        assert_frame_equal(result, expected)
        self.assertEqual(result.columns.names, ['sth', 'sth2'])

    def test_setitem_list(self):

        self.frame['E'] = 'foo'
        data = self.frame[['A', 'B']]
        self.frame[['B', 'A']] = data

        assert_series_equal(self.frame['B'], data['A'], check_names=False)
        assert_series_equal(self.frame['A'], data['B'], check_names=False)

        with assertRaisesRegexp(ValueError, 'Columns must be same length as key'):
            data[['A']] = self.frame[['A', 'B']]
        with assertRaisesRegexp(ValueError, 'Length of values does not match '
                                'length of index'):
            data['A'] = range(len(data.index) - 1)

        df = DataFrame(0, lrange(3), ['tt1', 'tt2'], dtype=np.int_)
        df.ix[1, ['tt1', 'tt2']] = [1, 2]

        result = df.ix[1, ['tt1', 'tt2']]
        expected = Series([1, 2], df.columns, dtype=np.int_, name=1)
        assert_series_equal(result, expected)

        df['tt1'] = df['tt2'] = '0'
        df.ix[1, ['tt1', 'tt2']] = ['1', '2']
        result = df.ix[1, ['tt1', 'tt2']]
        expected = Series(['1', '2'], df.columns, name=1)
        assert_series_equal(result, expected)

    def test_setitem_list_not_dataframe(self):
        data = np.random.randn(len(self.frame), 2)
        self.frame[['A', 'B']] = data
        assert_almost_equal(self.frame[['A', 'B']].values, data)

    def test_setitem_list_of_tuples(self):
        tuples = lzip(self.frame['A'], self.frame['B'])
        self.frame['tuples'] = tuples

        result = self.frame['tuples']
        expected = Series(tuples, index=self.frame.index, name='tuples')
        assert_series_equal(result, expected)

    def test_setitem_mulit_index(self):
        # GH7655, test that assigning to a sub-frame of a frame
        # with multi-index columns aligns both rows and columns
        it = ['jim', 'joe', 'jolie'], ['first', 'last'], \
             ['left', 'center', 'right']

        cols = MultiIndex.from_product(it)
        index = pd.date_range('20141006',periods=20)
        vals = np.random.randint(1, 1000, (len(index), len(cols)))
        df = pd.DataFrame(vals, columns=cols, index=index)

        i, j = df.index.values.copy(), it[-1][:]

        np.random.shuffle(i)
        df['jim'] = df['jolie'].loc[i, ::-1]
        assert_frame_equal(df['jim'], df['jolie'])

        np.random.shuffle(j)
        df[('joe', 'first')] = df[('jolie', 'last')].loc[i, j]
        assert_frame_equal(df[('joe', 'first')], df[('jolie', 'last')])

        np.random.shuffle(j)
        df[('joe', 'last')] = df[('jolie', 'first')].loc[i, j]
        assert_frame_equal(df[('joe', 'last')], df[('jolie', 'first')])

    def test_inplace_ops_alignment(self):

        # inplace ops / ops alignment
        # GH 8511

        columns = list('abcdefg')
        X_orig = DataFrame(np.arange(10*len(columns)).reshape(-1,len(columns)), columns=columns, index=range(10))
        Z = 100*X_orig.iloc[:,1:-1].copy()
        block1 = list('bedcf')
        subs = list('bcdef')

        # add
        X = X_orig.copy()
        result1 = (X[block1] + Z).reindex(columns=subs)

        X[block1] += Z
        result2 = X.reindex(columns=subs)

        X = X_orig.copy()
        result3 = (X[block1] + Z[block1]).reindex(columns=subs)

        X[block1] += Z[block1]
        result4 = X.reindex(columns=subs)

        assert_frame_equal(result1, result2)
        assert_frame_equal(result1, result3)
        assert_frame_equal(result1, result4)

        # sub
        X = X_orig.copy()
        result1 = (X[block1] - Z).reindex(columns=subs)

        X[block1] -= Z
        result2 = X.reindex(columns=subs)

        X = X_orig.copy()
        result3 = (X[block1] - Z[block1]).reindex(columns=subs)

        X[block1] -= Z[block1]
        result4 = X.reindex(columns=subs)

        assert_frame_equal(result1, result2)
        assert_frame_equal(result1, result3)
        assert_frame_equal(result1, result4)

    def test_inplace_ops_identity(self):

        # GH 5104
        # make sure that we are actually changing the object
        s_orig = Series([1, 2, 3])
        df_orig = DataFrame(np.random.randint(0,5,size=10).reshape(-1,5))

        # no dtype change
        s = s_orig.copy()
        s2 = s
        s += 1
        assert_series_equal(s,s2)
        assert_series_equal(s_orig+1,s)
        self.assertIs(s,s2)
        self.assertIs(s._data,s2._data)

        df = df_orig.copy()
        df2 = df
        df += 1
        assert_frame_equal(df,df2)
        assert_frame_equal(df_orig+1,df)
        self.assertIs(df,df2)
        self.assertIs(df._data,df2._data)

        # dtype change
        s = s_orig.copy()
        s2 = s
        s += 1.5
        assert_series_equal(s,s2)
        assert_series_equal(s_orig+1.5,s)

        df = df_orig.copy()
        df2 = df
        df += 1.5
        assert_frame_equal(df,df2)
        assert_frame_equal(df_orig+1.5,df)
        self.assertIs(df,df2)
        self.assertIs(df._data,df2._data)

        # mixed dtype
        arr = np.random.randint(0,10,size=5)
        df_orig = DataFrame({'A' : arr.copy(), 'B' : 'foo'})
        df = df_orig.copy()
        df2 = df
        df['A'] += 1
        expected = DataFrame({'A' : arr.copy()+1, 'B' : 'foo'})
        assert_frame_equal(df,expected)
        assert_frame_equal(df2,expected)
        self.assertIs(df._data,df2._data)

        df = df_orig.copy()
        df2 = df
        df['A'] += 1.5
        expected = DataFrame({'A' : arr.copy()+1.5, 'B' : 'foo'})
        assert_frame_equal(df,expected)
        assert_frame_equal(df2,expected)
        self.assertIs(df._data,df2._data)

    def test_getitem_boolean(self):
        # boolean indexing
        d = self.tsframe.index[10]
        indexer = self.tsframe.index > d
        indexer_obj = indexer.astype(object)

        subindex = self.tsframe.index[indexer]
        subframe = self.tsframe[indexer]

        self.assert_numpy_array_equal(subindex, subframe.index)
        with assertRaisesRegexp(ValueError, 'Item wrong length'):
            self.tsframe[indexer[:-1]]

        subframe_obj = self.tsframe[indexer_obj]
        assert_frame_equal(subframe_obj, subframe)

        with tm.assertRaisesRegexp(ValueError, 'boolean values only'):
            self.tsframe[self.tsframe]

        # test that Series work
        indexer_obj = Series(indexer_obj, self.tsframe.index)

        subframe_obj = self.tsframe[indexer_obj]
        assert_frame_equal(subframe_obj, subframe)

        # test that Series indexers reindex
        import warnings
        warnings.filterwarnings(action='ignore', category=UserWarning)

        indexer_obj = indexer_obj.reindex(self.tsframe.index[::-1])

        subframe_obj = self.tsframe[indexer_obj]
        assert_frame_equal(subframe_obj, subframe)

        warnings.filterwarnings(action='default', category=UserWarning)

        # test df[df > 0]
        for df in [ self.tsframe, self.mixed_frame, self.mixed_float, self.mixed_int ]:

            data = df._get_numeric_data()
            bif = df[df > 0]
            bifw = DataFrame(dict([ (c,np.where(data[c] > 0, data[c], np.nan)) for c in data.columns ]),
                             index=data.index, columns=data.columns)

            # add back other columns to compare
            for c in df.columns:
                if c not in bifw:
                    bifw[c] = df[c]
            bifw = bifw.reindex(columns = df.columns)

            assert_frame_equal(bif, bifw, check_dtype=False)
            for c in df.columns:
                if bif[c].dtype != bifw[c].dtype:
                    self.assertEqual(bif[c].dtype, df[c].dtype)

    def test_getitem_boolean_casting(self):

        # don't upcast if we don't need to
        df = self.tsframe.copy()
        df['E'] = 1
        df['E'] = df['E'].astype('int32')
        df['E1'] = df['E'].copy()
        df['F'] = 1
        df['F'] = df['F'].astype('int64')
        df['F1'] = df['F'].copy()

        casted = df[df>0]
        result = casted.get_dtype_counts()
        expected = Series({'float64': 4, 'int32' : 2, 'int64' : 2})
        assert_series_equal(result, expected)

        # int block splitting
        df.ix[1:3,['E1','F1']] = 0
        casted = df[df>0]
        result = casted.get_dtype_counts()
        expected = Series({'float64': 6, 'int32' : 1, 'int64' : 1})
        assert_series_equal(result, expected)

        # where dtype conversions
        # GH 3733
        df = DataFrame(data = np.random.randn(100, 50))
        df = df.where(df > 0) # create nans
        bools = df > 0
        mask = isnull(df)
        expected = bools.astype(float).mask(mask)
        result = bools.mask(mask)
        assert_frame_equal(result,expected)

    def test_getitem_boolean_list(self):
        df = DataFrame(np.arange(12).reshape(3, 4))

        def _checkit(lst):
            result = df[lst]
            expected = df.ix[df.index[lst]]
            assert_frame_equal(result, expected)

        _checkit([True, False, True])
        _checkit([True, True, True])
        _checkit([False, False, False])

    def test_getitem_boolean_iadd(self):
        arr = randn(5, 5)

        df = DataFrame(arr.copy(), columns = ['A','B','C','D','E'])

        df[df < 0] += 1
        arr[arr < 0] += 1

        assert_almost_equal(df.values, arr)

    def test_boolean_index_empty_corner(self):
        # #2096
        blah = DataFrame(np.empty([0, 1]), columns=['A'],
                         index=DatetimeIndex([]))

        # both of these should succeed trivially
        k = np.array([], bool)

        blah[k]
        blah[k] = 0

    def test_getitem_ix_mixed_integer(self):
        df = DataFrame(np.random.randn(4, 3),
                       index=[1, 10, 'C', 'E'], columns=[1, 2, 3])

        result = df.ix[:-1]
        expected = df.ix[df.index[:-1]]
        assert_frame_equal(result, expected)

        result = df.ix[[1, 10]]
        expected = df.ix[Index([1, 10], dtype=object)]
        assert_frame_equal(result, expected)

    def test_getitem_setitem_ix_negative_integers(self):
        result = self.frame.ix[:, -1]
        assert_series_equal(result, self.frame['D'])

        result = self.frame.ix[:, [-1]]
        assert_frame_equal(result, self.frame[['D']])

        result = self.frame.ix[:, [-1, -2]]
        assert_frame_equal(result, self.frame[['D', 'C']])

        self.frame.ix[:, [-1]] = 0
        self.assertTrue((self.frame['D'] == 0).all())

        df = DataFrame(np.random.randn(8, 4))
        self.assertTrue(isnull(df.ix[:, [-1]].values).all())

        # #1942
        a = DataFrame(randn(20, 2), index=[chr(x + 65) for x in range(20)])
        a.ix[-1] = a.ix[-2]

        assert_series_equal(a.ix[-1], a.ix[-2], check_names=False)
        self.assertEqual(a.ix[-1].name, 'T')
        self.assertEqual(a.ix[-2].name, 'S')

    def test_getattr(self):
        tm.assert_series_equal(self.frame.A, self.frame['A'])
        self.assertRaises(AttributeError, getattr, self.frame,
                          'NONEXISTENT_NAME')

    def test_setattr_column(self):
        df = DataFrame({'foobar': 1}, index=lrange(10))

        df.foobar = 5
        self.assertTrue((df.foobar == 5).all())

    def test_setitem(self):
        # not sure what else to do here
        series = self.frame['A'][::2]
        self.frame['col5'] = series
        self.assertIn('col5', self.frame)
        tm.assert_dict_equal(series, self.frame['col5'],
                             compare_keys=False)

        series = self.frame['A']
        self.frame['col6'] = series
        tm.assert_dict_equal(series, self.frame['col6'],
                             compare_keys=False)

        with tm.assertRaises(KeyError):
            self.frame[randn(len(self.frame) + 1)] = 1

        # set ndarray
        arr = randn(len(self.frame))
        self.frame['col9'] = arr
        self.assertTrue((self.frame['col9'] == arr).all())

        self.frame['col7'] = 5
        assert((self.frame['col7'] == 5).all())

        self.frame['col0'] = 3.14
        assert((self.frame['col0'] == 3.14).all())

        self.frame['col8'] = 'foo'
        assert((self.frame['col8'] == 'foo').all())

        # this is partially a view (e.g. some blocks are view)
        # so raise/warn
        smaller = self.frame[:2]
        def f():
            smaller['col10'] = ['1', '2']
        self.assertRaises(com.SettingWithCopyError, f)
        self.assertEqual(smaller['col10'].dtype, np.object_)
        self.assertTrue((smaller['col10'] == ['1', '2']).all())

        # with a dtype
        for dtype in ['int32','int64','float32','float64']:
            self.frame[dtype] = np.array(arr,dtype=dtype)
            self.assertEqual(self.frame[dtype].dtype.name, dtype)

        # dtype changing GH4204
        df = DataFrame([[0,0]])
        df.iloc[0] = np.nan
        expected = DataFrame([[np.nan,np.nan]])
        assert_frame_equal(df,expected)

        df = DataFrame([[0,0]])
        df.loc[0] = np.nan
        assert_frame_equal(df,expected)

    def test_setitem_tuple(self):
        self.frame['A', 'B'] = self.frame['A']
        assert_series_equal(self.frame['A', 'B'], self.frame['A'], check_names=False)

    def test_setitem_always_copy(self):
        s = self.frame['A'].copy()
        self.frame['E'] = s

        self.frame['E'][5:10] = nan
        self.assertTrue(notnull(s[5:10]).all())

    def test_setitem_boolean(self):
        df = self.frame.copy()
        values = self.frame.values

        df[df['A'] > 0] = 4
        values[values[:, 0] > 0] = 4
        assert_almost_equal(df.values, values)

        # test that column reindexing works
        series = df['A'] == 4
        series = series.reindex(df.index[::-1])
        df[series] = 1
        values[values[:, 0] == 4] = 1
        assert_almost_equal(df.values, values)

        df[df > 0] = 5
        values[values > 0] = 5
        assert_almost_equal(df.values, values)

        df[df == 5] = 0
        values[values == 5] = 0
        assert_almost_equal(df.values, values)

        # a df that needs alignment first
        df[df[:-1] < 0] = 2
        np.putmask(values[:-1], values[:-1] < 0, 2)
        assert_almost_equal(df.values, values)

        # indexed with same shape but rows-reversed df
        df[df[::-1] == 2] = 3
        values[values == 2] = 3
        assert_almost_equal(df.values, values)

        with assertRaisesRegexp(TypeError, 'Must pass DataFrame with boolean '
                                'values only'):
            df[df * 0] = 2

        # index with DataFrame
        mask = df > np.abs(df)
        expected = df.copy()
        df[df > np.abs(df)] = nan
        expected.values[mask.values] = nan
        assert_frame_equal(df, expected)

        # set from DataFrame
        expected = df.copy()
        df[df > np.abs(df)] = df * 2
        np.putmask(expected.values, mask.values, df.values * 2)
        assert_frame_equal(df, expected)

    def test_setitem_cast(self):
        self.frame['D'] = self.frame['D'].astype('i8')
        self.assertEqual(self.frame['D'].dtype, np.int64)

        # #669, should not cast?
        # this is now set to int64, which means a replacement of the column to
        # the value dtype (and nothing to do with the existing dtype)
        self.frame['B'] = 0
        self.assertEqual(self.frame['B'].dtype, np.int64)

        # cast if pass array of course
        self.frame['B'] = np.arange(len(self.frame))
        self.assertTrue(issubclass(self.frame['B'].dtype.type, np.integer))

        self.frame['foo'] = 'bar'
        self.frame['foo'] = 0
        self.assertEqual(self.frame['foo'].dtype, np.int64)

        self.frame['foo'] = 'bar'
        self.frame['foo'] = 2.5
        self.assertEqual(self.frame['foo'].dtype, np.float64)

        self.frame['something'] = 0
        self.assertEqual(self.frame['something'].dtype, np.int64)
        self.frame['something'] = 2
        self.assertEqual(self.frame['something'].dtype, np.int64)
        self.frame['something'] = 2.5
        self.assertEqual(self.frame['something'].dtype, np.float64)

        # GH 7704
        # dtype conversion on setting
        df = DataFrame(np.random.rand(30, 3), columns=tuple('ABC'))
        df['event'] = np.nan
        df.loc[10,'event'] = 'foo'
        result = df.get_dtype_counts().sort_values()
        expected = Series({'float64' : 3, 'object' : 1 }).sort_values()
        assert_series_equal(result, expected)

    def test_setitem_boolean_column(self):
        expected = self.frame.copy()
        mask = self.frame['A'] > 0

        self.frame.ix[mask, 'B'] = 0
        expected.values[mask.values, 1] = 0

        assert_frame_equal(self.frame, expected)

    def test_setitem_corner(self):
        # corner case
        df = DataFrame({'B': [1., 2., 3.],
                        'C': ['a', 'b', 'c']},
                       index=np.arange(3))
        del df['B']
        df['B'] = [1., 2., 3.]
        self.assertIn('B', df)
        self.assertEqual(len(df.columns), 2)

        df['A'] = 'beginning'
        df['E'] = 'foo'
        df['D'] = 'bar'
        df[datetime.now()] = 'date'
        df[datetime.now()] = 5.

        # what to do when empty frame with index
        dm = DataFrame(index=self.frame.index)
        dm['A'] = 'foo'
        dm['B'] = 'bar'
        self.assertEqual(len(dm.columns), 2)
        self.assertEqual(dm.values.dtype, np.object_)

        # upcast
        dm['C'] = 1
        self.assertEqual(dm['C'].dtype, np.int64)

        dm['E'] = 1.
        self.assertEqual(dm['E'].dtype, np.float64)

        # set existing column
        dm['A'] = 'bar'
        self.assertEqual('bar', dm['A'][0])

        dm = DataFrame(index=np.arange(3))
        dm['A'] = 1
        dm['foo'] = 'bar'
        del dm['foo']
        dm['foo'] = 'bar'
        self.assertEqual(dm['foo'].dtype, np.object_)

        dm['coercable'] = ['1', '2', '3']
        self.assertEqual(dm['coercable'].dtype, np.object_)

    def test_setitem_corner2(self):
        data = {"title": ['foobar', 'bar', 'foobar'] + ['foobar'] * 17,
                "cruft": np.random.random(20)}

        df = DataFrame(data)
        ix = df[df['title'] == 'bar'].index

        df.ix[ix, ['title']] = 'foobar'
        df.ix[ix, ['cruft']] = 0

        assert(df.ix[1, 'title'] == 'foobar')
        assert(df.ix[1, 'cruft'] == 0)

    def test_setitem_ambig(self):
        # difficulties with mixed-type data
        from decimal import Decimal

        # created as float type
        dm = DataFrame(index=lrange(3), columns=lrange(3))

        coercable_series = Series([Decimal(1) for _ in range(3)],
                                  index=lrange(3))
        uncoercable_series = Series(['foo', 'bzr', 'baz'], index=lrange(3))

        dm[0] = np.ones(3)
        self.assertEqual(len(dm.columns), 3)
        # self.assertIsNone(dm.objects)

        dm[1] = coercable_series
        self.assertEqual(len(dm.columns), 3)
        # self.assertIsNone(dm.objects)

        dm[2] = uncoercable_series
        self.assertEqual(len(dm.columns), 3)
        # self.assertIsNotNone(dm.objects)
        self.assertEqual(dm[2].dtype, np.object_)

    def test_setitem_clear_caches(self):
        # GH #304
        df = DataFrame({'x': [1.1, 2.1, 3.1, 4.1], 'y': [5.1, 6.1, 7.1, 8.1]},
                       index=[0, 1, 2, 3])
        df.insert(2, 'z', np.nan)

        # cache it
        foo = df['z']

        df.ix[2:, 'z'] = 42

        expected = Series([np.nan, np.nan, 42, 42], index=df.index, name='z')
        self.assertIsNot(df['z'], foo)
        assert_series_equal(df['z'], expected)

    def test_setitem_None(self):
        # GH #766
        self.frame[None] = self.frame['A']
        assert_series_equal(self.frame.iloc[:,-1], self.frame['A'], check_names=False)
        assert_series_equal(self.frame.loc[:,None], self.frame['A'], check_names=False)
        assert_series_equal(self.frame[None], self.frame['A'], check_names=False)
        repr(self.frame)

    def test_setitem_empty(self):
        # GH 9596
        df = pd.DataFrame({'a': ['1', '2', '3'],
                           'b': ['11', '22', '33'],
                           'c': ['111', '222', '333']})

        result = df.copy()
        result.loc[result.b.isnull(), 'a'] = result.a
        assert_frame_equal(result, df)

    def test_setitem_empty_frame_with_boolean(self):
        # Test for issue #10126

        for dtype in ('float', 'int64'):
            for df in [
                    pd.DataFrame(dtype=dtype),
                    pd.DataFrame(dtype=dtype, index=[1]),
                    pd.DataFrame(dtype=dtype, columns=['A']),
            ]:
                df2 = df.copy()
                df[df > df2] = 47
                assert_frame_equal(df, df2)

    def test_delitem_corner(self):
        f = self.frame.copy()
        del f['D']
        self.assertEqual(len(f.columns), 3)
        self.assertRaises(KeyError, f.__delitem__, 'D')
        del f['B']
        self.assertEqual(len(f.columns), 2)

    def test_getitem_fancy_2d(self):
        f = self.frame
        ix = f.ix

        assert_frame_equal(ix[:, ['B', 'A']], f.reindex(columns=['B', 'A']))

        subidx = self.frame.index[[5, 4, 1]]
        assert_frame_equal(ix[subidx, ['B', 'A']],
                           f.reindex(index=subidx, columns=['B', 'A']))

        # slicing rows, etc.
        assert_frame_equal(ix[5:10], f[5:10])
        assert_frame_equal(ix[5:10, :], f[5:10])
        assert_frame_equal(ix[:5, ['A', 'B']],
                           f.reindex(index=f.index[:5], columns=['A', 'B']))

        # slice rows with labels, inclusive!
        expected = ix[5:11]
        result = ix[f.index[5]:f.index[10]]
        assert_frame_equal(expected, result)

        # slice columns
        assert_frame_equal(ix[:, :2], f.reindex(columns=['A', 'B']))

        # get view
        exp = f.copy()
        ix[5:10].values[:] = 5
        exp.values[5:10] = 5
        assert_frame_equal(f, exp)

        self.assertRaises(ValueError, ix.__getitem__, f > 0.5)

    def test_slice_floats(self):
        index = [52195.504153, 52196.303147, 52198.369883]
        df = DataFrame(np.random.rand(3, 2), index=index)

        s1 = df.ix[52195.1:52196.5]
        self.assertEqual(len(s1), 2)

        s1 = df.ix[52195.1:52196.6]
        self.assertEqual(len(s1), 2)

        s1 = df.ix[52195.1:52198.9]
        self.assertEqual(len(s1), 3)

    def test_getitem_fancy_slice_integers_step(self):
        df = DataFrame(np.random.randn(10, 5))

        # this is OK
        result = df.ix[:8:2]
        df.ix[:8:2] = np.nan
        self.assertTrue(isnull(df.ix[:8:2]).values.all())

    def test_getitem_setitem_integer_slice_keyerrors(self):
        df = DataFrame(np.random.randn(10, 5), index=lrange(0, 20, 2))

        # this is OK
        cp = df.copy()
        cp.ix[4:10] = 0
        self.assertTrue((cp.ix[4:10] == 0).values.all())

        # so is this
        cp = df.copy()
        cp.ix[3:11] = 0
        self.assertTrue((cp.ix[3:11] == 0).values.all())

        result = df.ix[4:10]
        result2 = df.ix[3:11]
        expected = df.reindex([4, 6, 8, 10])

        assert_frame_equal(result, expected)
        assert_frame_equal(result2, expected)

        # non-monotonic, raise KeyError
        df2 = df.iloc[lrange(5) + lrange(5, 10)[::-1]]
        self.assertRaises(KeyError, df2.ix.__getitem__, slice(3, 11))
        self.assertRaises(KeyError, df2.ix.__setitem__, slice(3, 11), 0)

    def test_setitem_fancy_2d(self):
        f = self.frame
        ix = f.ix

        # case 1
        frame = self.frame.copy()
        expected = frame.copy()
        frame.ix[:, ['B', 'A']] = 1
        expected['B'] = 1.
        expected['A'] = 1.
        assert_frame_equal(frame, expected)

        # case 2
        frame = self.frame.copy()
        frame2 = self.frame.copy()

        expected = frame.copy()

        subidx = self.frame.index[[5, 4, 1]]
        values = randn(3, 2)

        frame.ix[subidx, ['B', 'A']] = values
        frame2.ix[[5, 4, 1], ['B', 'A']] = values

        expected['B'].ix[subidx] = values[:, 0]
        expected['A'].ix[subidx] = values[:, 1]

        assert_frame_equal(frame, expected)
        assert_frame_equal(frame2, expected)

        # case 3: slicing rows, etc.
        frame = self.frame.copy()

        expected1 = self.frame.copy()
        frame.ix[5:10] = 1.
        expected1.values[5:10] = 1.
        assert_frame_equal(frame, expected1)

        expected2 = self.frame.copy()
        arr = randn(5, len(frame.columns))
        frame.ix[5:10] = arr
        expected2.values[5:10] = arr
        assert_frame_equal(frame, expected2)

        # case 4
        frame = self.frame.copy()
        frame.ix[5:10, :] = 1.
        assert_frame_equal(frame, expected1)
        frame.ix[5:10, :] = arr
        assert_frame_equal(frame, expected2)

        # case 5
        frame = self.frame.copy()
        frame2 = self.frame.copy()

        expected = self.frame.copy()
        values = randn(5, 2)

        frame.ix[:5, ['A', 'B']] = values
        expected['A'][:5] = values[:, 0]
        expected['B'][:5] = values[:, 1]
        assert_frame_equal(frame, expected)

        frame2.ix[:5, [0, 1]] = values
        assert_frame_equal(frame2, expected)

        # case 6: slice rows with labels, inclusive!
        frame = self.frame.copy()
        expected = self.frame.copy()

        frame.ix[frame.index[5]:frame.index[10]] = 5.
        expected.values[5:11] = 5
        assert_frame_equal(frame, expected)

        # case 7: slice columns
        frame = self.frame.copy()
        frame2 = self.frame.copy()
        expected = self.frame.copy()

        # slice indices
        frame.ix[:, 1:3] = 4.
        expected.values[:, 1:3] = 4.
        assert_frame_equal(frame, expected)

        # slice with labels
        frame.ix[:, 'B':'C'] = 4.
        assert_frame_equal(frame, expected)

        # new corner case of boolean slicing / setting
        frame = DataFrame(lzip([2, 3, 9, 6, 7], [np.nan] * 5),
                          columns=['a', 'b'])
        lst = [100]
        lst.extend([np.nan] * 4)
        expected = DataFrame(lzip([100, 3, 9, 6, 7], lst),
                             columns=['a', 'b'])
        frame[frame['a'] == 2] = 100
        assert_frame_equal(frame, expected)

    def test_fancy_getitem_slice_mixed(self):
        sliced = self.mixed_frame.ix[:, -3:]
        self.assertEqual(sliced['D'].dtype, np.float64)

        # get view with single block
        # setting it triggers setting with copy
        sliced = self.frame.ix[:, -3:]
        def f():
            sliced['C'] = 4.
        self.assertRaises(com.SettingWithCopyError, f)
        self.assertTrue((self.frame['C'] == 4).all())

    def test_fancy_setitem_int_labels(self):
        # integer index defers to label-based indexing

        df = DataFrame(np.random.randn(10, 5), index=np.arange(0, 20, 2))

        tmp = df.copy()
        exp = df.copy()
        tmp.ix[[0, 2, 4]] = 5
        exp.values[:3] = 5
        assert_frame_equal(tmp, exp)

        tmp = df.copy()
        exp = df.copy()
        tmp.ix[6] = 5
        exp.values[3] = 5
        assert_frame_equal(tmp, exp)

        tmp = df.copy()
        exp = df.copy()
        tmp.ix[:, 2] = 5

        # tmp correctly sets the dtype
        # so match the exp way
        exp[2] = 5
        assert_frame_equal(tmp, exp)

    def test_fancy_getitem_int_labels(self):
        df = DataFrame(np.random.randn(10, 5), index=np.arange(0, 20, 2))

        result = df.ix[[4, 2, 0], [2, 0]]
        expected = df.reindex(index=[4, 2, 0], columns=[2, 0])
        assert_frame_equal(result, expected)

        result = df.ix[[4, 2, 0]]
        expected = df.reindex(index=[4, 2, 0])
        assert_frame_equal(result, expected)

        result = df.ix[4]
        expected = df.xs(4)
        assert_series_equal(result, expected)

        result = df.ix[:, 3]
        expected = df[3]
        assert_series_equal(result, expected)

    def test_fancy_index_int_labels_exceptions(self):
        df = DataFrame(np.random.randn(10, 5), index=np.arange(0, 20, 2))

        # labels that aren't contained
        self.assertRaises(KeyError, df.ix.__setitem__,
                          ([0, 1, 2], [2, 3, 4]), 5)

        # try to set indices not contained in frame
        self.assertRaises(KeyError,
                          self.frame.ix.__setitem__,
                          ['foo', 'bar', 'baz'], 1)
        self.assertRaises(KeyError,
                          self.frame.ix.__setitem__,
                          (slice(None, None), ['E']), 1)

        # partial setting now allows this GH2578
        #self.assertRaises(KeyError,
        #                  self.frame.ix.__setitem__,
        #                  (slice(None, None), 'E'), 1)

    def test_setitem_fancy_mixed_2d(self):
        self.mixed_frame.ix[:5, ['C', 'B', 'A']] = 5
        result = self.mixed_frame.ix[:5, ['C', 'B', 'A']]
        self.assertTrue((result.values == 5).all())

        self.mixed_frame.ix[5] = np.nan
        self.assertTrue(isnull(self.mixed_frame.ix[5]).all())

        self.mixed_frame.ix[5] = self.mixed_frame.ix[6]
        assert_series_equal(self.mixed_frame.ix[5], self.mixed_frame.ix[6],
                            check_names=False)

        # #1432
        df = DataFrame({1: [1., 2., 3.],
                        2: [3, 4, 5]})
        self.assertTrue(df._is_mixed_type)

        df.ix[1] = [5, 10]

        expected = DataFrame({1: [1., 5., 3.],
                              2: [3, 10, 5]})

        assert_frame_equal(df, expected)

    def test_ix_align(self):
        b = Series(randn(10), name=0).sort_values()
        df_orig = DataFrame(randn(10, 4))
        df = df_orig.copy()

        df.ix[:, 0] = b
        assert_series_equal(df.ix[:, 0].reindex(b.index), b)

        dft = df_orig.T
        dft.ix[0, :] = b
        assert_series_equal(dft.ix[0, :].reindex(b.index), b)

        df = df_orig.copy()
        df.ix[:5, 0] = b
        s = df.ix[:5, 0]
        assert_series_equal(s, b.reindex(s.index))

        dft = df_orig.T
        dft.ix[0, :5] = b
        s = dft.ix[0, :5]
        assert_series_equal(s, b.reindex(s.index))

        df = df_orig.copy()
        idx = [0, 1, 3, 5]
        df.ix[idx, 0] = b
        s = df.ix[idx, 0]
        assert_series_equal(s, b.reindex(s.index))

        dft = df_orig.T
        dft.ix[0, idx] = b
        s = dft.ix[0, idx]
        assert_series_equal(s, b.reindex(s.index))

    def test_ix_frame_align(self):
        b = DataFrame(np.random.randn(3, 4))
        df_orig = DataFrame(randn(10, 4))
        df = df_orig.copy()

        df.ix[:3] = b
        out = b.ix[:3]
        assert_frame_equal(out, b)

        b.sort_index(inplace=True)

        df = df_orig.copy()
        df.ix[[0, 1, 2]] = b
        out = df.ix[[0, 1, 2]].reindex(b.index)
        assert_frame_equal(out, b)

        df = df_orig.copy()
        df.ix[:3] = b
        out = df.ix[:3]
        assert_frame_equal(out, b.reindex(out.index))

    def test_getitem_setitem_non_ix_labels(self):
        df = tm.makeTimeDataFrame()

        start, end = df.index[[5, 10]]

        result = df.ix[start:end]
        result2 = df[start:end]
        expected = df[5:11]
        assert_frame_equal(result, expected)
        assert_frame_equal(result2, expected)

        result = df.copy()
        result.ix[start:end] = 0
        result2 = df.copy()
        result2[start:end] = 0
        expected = df.copy()
        expected[5:11] = 0
        assert_frame_equal(result, expected)
        assert_frame_equal(result2, expected)

    def test_ix_multi_take(self):
        df = DataFrame(np.random.randn(3, 2))
        rs = df.ix[df.index == 0, :]
        xp = df.reindex([0])
        assert_frame_equal(rs, xp)

        """ #1321
        df = DataFrame(np.random.randn(3, 2))
        rs = df.ix[df.index==0, df.columns==1]
        xp = df.reindex([0], [1])
        assert_frame_equal(rs, xp)
        """

    def test_ix_multi_take_nonint_index(self):
        df = DataFrame(np.random.randn(3, 2), index=['x', 'y', 'z'],
                       columns=['a', 'b'])
        rs = df.ix[[0], [0]]
        xp = df.reindex(['x'], columns=['a'])
        assert_frame_equal(rs, xp)

    def test_ix_multi_take_multiindex(self):
        df = DataFrame(np.random.randn(3, 2), index=['x', 'y', 'z'],
                       columns=[['a', 'b'], ['1', '2']])
        rs = df.ix[[0], [0]]
        xp = df.reindex(['x'], columns=[('a', '1')])
        assert_frame_equal(rs, xp)

    def test_ix_dup(self):
        idx = Index(['a', 'a', 'b', 'c', 'd', 'd'])
        df = DataFrame(np.random.randn(len(idx), 3), idx)

        sub = df.ix[:'d']
        assert_frame_equal(sub, df)

        sub = df.ix['a':'c']
        assert_frame_equal(sub, df.ix[0:4])

        sub = df.ix['b':'d']
        assert_frame_equal(sub, df.ix[2:])

    def test_getitem_fancy_1d(self):
        f = self.frame
        ix = f.ix

        # return self if no slicing...for now
        self.assertIs(ix[:, :], f)

        # low dimensional slice
        xs1 = ix[2, ['C', 'B', 'A']]
        xs2 = f.xs(f.index[2]).reindex(['C', 'B', 'A'])
        assert_series_equal(xs1, xs2)

        ts1 = ix[5:10, 2]
        ts2 = f[f.columns[2]][5:10]
        assert_series_equal(ts1, ts2)

        # positional xs
        xs1 = ix[0]
        xs2 = f.xs(f.index[0])
        assert_series_equal(xs1, xs2)

        xs1 = ix[f.index[5]]
        xs2 = f.xs(f.index[5])
        assert_series_equal(xs1, xs2)

        # single column
        assert_series_equal(ix[:, 'A'], f['A'])

        # return view
        exp = f.copy()
        exp.values[5] = 4
        ix[5][:] = 4
        assert_frame_equal(exp, f)

        exp.values[:, 1] = 6
        ix[:, 1][:] = 6
        assert_frame_equal(exp, f)

        # slice of mixed-frame
        xs = self.mixed_frame.ix[5]
        exp = self.mixed_frame.xs(self.mixed_frame.index[5])
        assert_series_equal(xs, exp)

    def test_setitem_fancy_1d(self):

        # case 1: set cross-section for indices
        frame = self.frame.copy()
        expected = self.frame.copy()

        frame.ix[2, ['C', 'B', 'A']] = [1., 2., 3.]
        expected['C'][2] = 1.
        expected['B'][2] = 2.
        expected['A'][2] = 3.
        assert_frame_equal(frame, expected)

        frame2 = self.frame.copy()
        frame2.ix[2, [3, 2, 1]] = [1., 2., 3.]
        assert_frame_equal(frame, expected)

        # case 2, set a section of a column
        frame = self.frame.copy()
        expected = self.frame.copy()

        vals = randn(5)
        expected.values[5:10, 2] = vals
        frame.ix[5:10, 2] = vals
        assert_frame_equal(frame, expected)

        frame2 = self.frame.copy()
        frame2.ix[5:10, 'B'] = vals
        assert_frame_equal(frame, expected)

        # case 3: full xs
        frame = self.frame.copy()
        expected = self.frame.copy()

        frame.ix[4] = 5.
        expected.values[4] = 5.
        assert_frame_equal(frame, expected)

        frame.ix[frame.index[4]] = 6.
        expected.values[4] = 6.
        assert_frame_equal(frame, expected)

        # single column
        frame = self.frame.copy()
        expected = self.frame.copy()

        frame.ix[:, 'A'] = 7.
        expected['A'] = 7.
        assert_frame_equal(frame, expected)

    def test_getitem_fancy_scalar(self):
        f = self.frame
        ix = f.ix
        # individual value
        for col in f.columns:
            ts = f[col]
            for idx in f.index[::5]:
                assert_almost_equal(ix[idx, col], ts[idx])

    def test_setitem_fancy_scalar(self):
        f = self.frame
        expected = self.frame.copy()
        ix = f.ix
        # individual value
        for j, col in enumerate(f.columns):
            ts = f[col]
            for idx in f.index[::5]:
                i = f.index.get_loc(idx)
                val = randn()
                expected.values[i, j] = val
                ix[idx, col] = val
                assert_frame_equal(f, expected)

    def test_getitem_fancy_boolean(self):
        f = self.frame
        ix = f.ix

        expected = f.reindex(columns=['B', 'D'])
        result = ix[:, [False, True, False, True]]
        assert_frame_equal(result, expected)

        expected = f.reindex(index=f.index[5:10], columns=['B', 'D'])
        result = ix[5:10, [False, True, False, True]]
        assert_frame_equal(result, expected)

        boolvec = f.index > f.index[7]
        expected = f.reindex(index=f.index[boolvec])
        result = ix[boolvec]
        assert_frame_equal(result, expected)
        result = ix[boolvec, :]
        assert_frame_equal(result, expected)

        result = ix[boolvec, 2:]
        expected = f.reindex(index=f.index[boolvec],
                             columns=['C', 'D'])
        assert_frame_equal(result, expected)

    def test_setitem_fancy_boolean(self):
        # from 2d, set with booleans
        frame = self.frame.copy()
        expected = self.frame.copy()

        mask = frame['A'] > 0
        frame.ix[mask] = 0.
        expected.values[mask.values] = 0.
        assert_frame_equal(frame, expected)

        frame = self.frame.copy()
        expected = self.frame.copy()
        frame.ix[mask, ['A', 'B']] = 0.
        expected.values[mask.values, :2] = 0.
        assert_frame_equal(frame, expected)

    def test_getitem_fancy_ints(self):
        result = self.frame.ix[[1, 4, 7]]
        expected = self.frame.ix[self.frame.index[[1, 4, 7]]]
        assert_frame_equal(result, expected)

        result = self.frame.ix[:, [2, 0, 1]]
        expected = self.frame.ix[:, self.frame.columns[[2, 0, 1]]]
        assert_frame_equal(result, expected)

    def test_getitem_setitem_fancy_exceptions(self):
        ix = self.frame.ix
        with assertRaisesRegexp(IndexingError, 'Too many indexers'):
            ix[:, :, :]

        with assertRaises(IndexingError):
            ix[:, :, :] = 1

    def test_getitem_setitem_boolean_misaligned(self):
        # boolean index misaligned labels
        mask = self.frame['A'][::-1] > 1

        result = self.frame.ix[mask]
        expected = self.frame.ix[mask[::-1]]
        assert_frame_equal(result, expected)

        cp = self.frame.copy()
        expected = self.frame.copy()
        cp.ix[mask] = 0
        expected.ix[mask] = 0
        assert_frame_equal(cp, expected)

    def test_getitem_setitem_boolean_multi(self):
        df = DataFrame(np.random.randn(3, 2))

        # get
        k1 = np.array([True, False, True])
        k2 = np.array([False, True])
        result = df.ix[k1, k2]
        expected = df.ix[[0, 2], [1]]
        assert_frame_equal(result, expected)

        expected = df.copy()
        df.ix[np.array([True, False, True]),
              np.array([False, True])] = 5
        expected.ix[[0, 2], [1]] = 5
        assert_frame_equal(df, expected)

    def test_getitem_setitem_float_labels(self):
        index = Index([1.5, 2, 3, 4, 5])
        df = DataFrame(np.random.randn(5, 5), index=index)

        result = df.ix[1.5:4]
        expected = df.reindex([1.5, 2, 3, 4])
        assert_frame_equal(result, expected)
        self.assertEqual(len(result), 4)

        result = df.ix[4:5]
        expected = df.reindex([4, 5])
        assert_frame_equal(result, expected)
        self.assertEqual(len(result), 2)

        # loc_float changes this to work properly
        result = df.ix[1:2]
        expected = df.iloc[0:2]
        assert_frame_equal(result, expected)

        df.ix[1:2] = 0
        result = df[1:2]
        self.assertTrue((result==0).all().all())

        # #2727
        index = Index([1.0, 2.5, 3.5, 4.5, 5.0])
        df = DataFrame(np.random.randn(5, 5), index=index)

        # positional slicing only via iloc!
        # stacklevel=False -> needed stacklevel depends on index type
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            result = df.iloc[1.0:5]

        expected = df.reindex([2.5, 3.5, 4.5, 5.0])
        assert_frame_equal(result, expected)
        self.assertEqual(len(result), 4)

        result = df.iloc[4:5]
        expected = df.reindex([5.0])
        assert_frame_equal(result, expected)
        self.assertEqual(len(result), 1)

        # GH 4892, float indexers in iloc are deprecated
        import warnings
        warnings.filterwarnings(action='error', category=FutureWarning)

        cp = df.copy()
        def f():
            cp.iloc[1.0:5] = 0
        self.assertRaises(FutureWarning, f)
        def f():
            result = cp.iloc[1.0:5] == 0
        self.assertRaises(FutureWarning, f)
        self.assertTrue(result.values.all())
        self.assertTrue((cp.iloc[0:1] == df.iloc[0:1]).values.all())

        warnings.filterwarnings(action='default', category=FutureWarning)

        cp = df.copy()
        cp.iloc[4:5] = 0
        self.assertTrue((cp.iloc[4:5] == 0).values.all())
        self.assertTrue((cp.iloc[0:4] == df.iloc[0:4]).values.all())

        # float slicing
        result = df.ix[1.0:5]
        expected = df
        assert_frame_equal(result, expected)
        self.assertEqual(len(result), 5)

        result = df.ix[1.1:5]
        expected = df.reindex([2.5, 3.5, 4.5, 5.0])
        assert_frame_equal(result, expected)
        self.assertEqual(len(result), 4)

        result = df.ix[4.51:5]
        expected = df.reindex([5.0])
        assert_frame_equal(result, expected)
        self.assertEqual(len(result), 1)

        result = df.ix[1.0:5.0]
        expected = df.reindex([1.0, 2.5, 3.5, 4.5, 5.0])
        assert_frame_equal(result, expected)
        self.assertEqual(len(result), 5)

        cp = df.copy()
        cp.ix[1.0:5.0] = 0
        result = cp.ix[1.0:5.0]
        self.assertTrue((result == 0).values.all())

    def test_setitem_single_column_mixed(self):
        df = DataFrame(randn(5, 3), index=['a', 'b', 'c', 'd', 'e'],
                       columns=['foo', 'bar', 'baz'])
        df['str'] = 'qux'
        df.ix[::2, 'str'] = nan
        expected = [nan, 'qux', nan, 'qux', nan]
        assert_almost_equal(df['str'].values, expected)

    def test_setitem_single_column_mixed_datetime(self):
        df = DataFrame(randn(5, 3), index=['a', 'b', 'c', 'd', 'e'],
                       columns=['foo', 'bar', 'baz'])

        df['timestamp'] = Timestamp('20010102')

        # check our dtypes
        result = df.get_dtype_counts()
        expected = Series({'float64': 3, 'datetime64[ns]': 1})
        assert_series_equal(result, expected)

        # set an allowable datetime64 type
        from pandas import tslib
        df.ix['b', 'timestamp'] = tslib.iNaT
        self.assertTrue(com.isnull(df.ix['b', 'timestamp']))

        # allow this syntax
        df.ix['c', 'timestamp'] = nan
        self.assertTrue(com.isnull(df.ix['c', 'timestamp']))

        # allow this syntax
        df.ix['d', :] = nan
        self.assertTrue(com.isnull(df.ix['c', :]).all() == False)

        # as of GH 3216 this will now work!
        # try to set with a list like item
        #self.assertRaises(
        #    Exception, df.ix.__setitem__, ('d', 'timestamp'), [nan])

    def test_setitem_frame(self):
        piece = self.frame.ix[:2, ['A', 'B']]
        self.frame.ix[-2:, ['A', 'B']] = piece.values
        assert_almost_equal(self.frame.ix[-2:, ['A', 'B']].values,
                            piece.values)

        # GH 3216

        # already aligned
        f = self.mixed_frame.copy()
        piece = DataFrame([[ 1, 2], [3, 4]], index=f.index[0:2],columns=['A', 'B'])
        key = (slice(None,2), ['A', 'B'])
        f.ix[key] = piece
        assert_almost_equal(f.ix[0:2, ['A', 'B']].values,
                            piece.values)

        # rows unaligned
        f = self.mixed_frame.copy()
        piece = DataFrame([[ 1, 2 ], [3, 4], [5, 6], [7, 8]], index=list(f.index[0:2]) + ['foo','bar'],columns=['A', 'B'])
        key = (slice(None,2), ['A', 'B'])
        f.ix[key] = piece
        assert_almost_equal(f.ix[0:2:, ['A', 'B']].values,
                            piece.values[0:2])

        # key is unaligned with values
        f = self.mixed_frame.copy()
        piece = f.ix[:2, ['A']]
        piece.index = f.index[-2:]
        key = (slice(-2, None), ['A', 'B'])
        f.ix[key] = piece
        piece['B'] = np.nan
        assert_almost_equal(f.ix[-2:, ['A', 'B']].values,
                            piece.values)

        # ndarray
        f = self.mixed_frame.copy()
        piece = self.mixed_frame.ix[:2, ['A', 'B']]
        key = (slice(-2, None), ['A', 'B'])
        f.ix[key] = piece.values
        assert_almost_equal(f.ix[-2:, ['A', 'B']].values,
                            piece.values)


        # needs upcasting
        df = DataFrame([[1,2,'foo'],[3,4,'bar']],columns=['A','B','C'])
        df2 = df.copy()
        df2.ix[:,['A','B']] = df.ix[:,['A','B']]+0.5
        expected = df.reindex(columns=['A','B'])
        expected += 0.5
        expected['C'] = df['C']
        assert_frame_equal(df2, expected)

    def test_setitem_frame_align(self):
        piece = self.frame.ix[:2, ['A', 'B']]
        piece.index = self.frame.index[-2:]
        piece.columns = ['A', 'B']
        self.frame.ix[-2:, ['A', 'B']] = piece
        assert_almost_equal(self.frame.ix[-2:, ['A', 'B']].values,
                            piece.values)

    def test_setitem_fancy_exceptions(self):
        pass

    def test_getitem_boolean_missing(self):
        pass

    def test_setitem_boolean_missing(self):
        pass

    def test_getitem_setitem_ix_duplicates(self):
        # #1201
        df = DataFrame(np.random.randn(5, 3),
                       index=['foo', 'foo', 'bar', 'baz', 'bar'])

        result = df.ix['foo']
        expected = df[:2]
        assert_frame_equal(result, expected)

        result = df.ix['bar']
        expected = df.ix[[2, 4]]
        assert_frame_equal(result, expected)

        result = df.ix['baz']
        expected = df.ix[3]
        assert_series_equal(result, expected)

    def test_getitem_ix_boolean_duplicates_multiple(self):
        # #1201
        df = DataFrame(np.random.randn(5, 3),
                       index=['foo', 'foo', 'bar', 'baz', 'bar'])

        result = df.ix[['bar']]
        exp = df.ix[[2, 4]]
        assert_frame_equal(result, exp)

        result = df.ix[df[1] > 0]
        exp = df[df[1] > 0]
        assert_frame_equal(result, exp)

        result = df.ix[df[0] > 0]
        exp = df[df[0] > 0]
        assert_frame_equal(result, exp)

    def test_getitem_setitem_ix_bool_keyerror(self):
        # #2199
        df = DataFrame({'a': [1, 2, 3]})

        self.assertRaises(KeyError, df.ix.__getitem__, False)
        self.assertRaises(KeyError, df.ix.__getitem__, True)

        self.assertRaises(KeyError, df.ix.__setitem__, False, 0)
        self.assertRaises(KeyError, df.ix.__setitem__, True, 0)

    def test_getitem_list_duplicates(self):
        # #1943
        df = DataFrame(np.random.randn(4, 4), columns=list('AABC'))
        df.columns.name = 'foo'

        result = df[['B', 'C']]
        self.assertEqual(result.columns.name, 'foo')

        expected = df.ix[:, 2:]
        assert_frame_equal(result, expected)

    def test_get_value(self):
        for idx in self.frame.index:
            for col in self.frame.columns:
                result = self.frame.get_value(idx, col)
                expected = self.frame[col][idx]
                assert_almost_equal(result, expected)

    def test_iteritems(self):
        df = DataFrame([[1, 2, 3], [4, 5, 6]], columns=['a', 'a', 'b'])
        for k, v in compat.iteritems(df):
            self.assertEqual(type(v), Series)

    def test_lookup(self):
        def alt(df, rows, cols):
            result = []
            for r, c in zip(rows, cols):
                result.append(df.get_value(r, c))
            return result

        def testit(df):
            rows = list(df.index) * len(df.columns)
            cols = list(df.columns) * len(df.index)
            result = df.lookup(rows, cols)
            expected = alt(df, rows, cols)
            assert_almost_equal(result, expected)

        testit(self.mixed_frame)
        testit(self.frame)

        df = DataFrame({'label': ['a', 'b', 'a', 'c'],
                        'mask_a': [True, True, False, True],
                        'mask_b': [True, False, False, False],
                        'mask_c': [False, True, False, True]})
        df['mask'] = df.lookup(df.index, 'mask_' + df['label'])
        exp_mask = alt(df, df.index, 'mask_' + df['label'])
        assert_almost_equal(df['mask'], exp_mask)
        self.assertEqual(df['mask'].dtype, np.bool_)

        with tm.assertRaises(KeyError):
            self.frame.lookup(['xyz'], ['A'])

        with tm.assertRaises(KeyError):
            self.frame.lookup([self.frame.index[0]], ['xyz'])

        with tm.assertRaisesRegexp(ValueError, 'same size'):
            self.frame.lookup(['a', 'b', 'c'], ['a'])

    def test_set_value(self):
        for idx in self.frame.index:
            for col in self.frame.columns:
                self.frame.set_value(idx, col, 1)
                assert_almost_equal(self.frame[col][idx], 1)

    def test_set_value_resize(self):

        res = self.frame.set_value('foobar', 'B', 0)
        self.assertIs(res, self.frame)
        self.assertEqual(res.index[-1], 'foobar')
        self.assertEqual(res.get_value('foobar', 'B'), 0)

        self.frame.loc['foobar','qux'] = 0
        self.assertEqual(self.frame.get_value('foobar', 'qux'), 0)

        res = self.frame.copy()
        res3 = res.set_value('foobar', 'baz', 'sam')
        self.assertEqual(res3['baz'].dtype, np.object_)

        res = self.frame.copy()
        res3 = res.set_value('foobar', 'baz', True)
        self.assertEqual(res3['baz'].dtype, np.object_)

        res = self.frame.copy()
        res3 = res.set_value('foobar', 'baz', 5)
        self.assertTrue(com.is_float_dtype(res3['baz']))
        self.assertTrue(isnull(res3['baz'].drop(['foobar'])).all())
        self.assertRaises(ValueError, res3.set_value, 'foobar', 'baz', 'sam')

    def test_set_value_with_index_dtype_change(self):
        df_orig = DataFrame(randn(3, 3), index=lrange(3), columns=list('ABC'))

        # this is actually ambiguous as the 2 is interpreted as a positional
        # so column is not created
        df = df_orig.copy()
        df.set_value('C', 2, 1.0)
        self.assertEqual(list(df.index), list(df_orig.index) + ['C'])
        #self.assertEqual(list(df.columns), list(df_orig.columns) + [2])

        df = df_orig.copy()
        df.loc['C', 2] = 1.0
        self.assertEqual(list(df.index), list(df_orig.index) + ['C'])
        #self.assertEqual(list(df.columns), list(df_orig.columns) + [2])

        # create both new
        df = df_orig.copy()
        df.set_value('C', 'D', 1.0)
        self.assertEqual(list(df.index), list(df_orig.index) + ['C'])
        self.assertEqual(list(df.columns), list(df_orig.columns) + ['D'])

        df = df_orig.copy()
        df.loc['C', 'D'] = 1.0
        self.assertEqual(list(df.index), list(df_orig.index) + ['C'])
        self.assertEqual(list(df.columns), list(df_orig.columns) + ['D'])

    def test_get_set_value_no_partial_indexing(self):
        # partial w/ MultiIndex raise exception
        index = MultiIndex.from_tuples([(0, 1), (0, 2), (1, 1), (1, 2)])
        df = DataFrame(index=index, columns=lrange(4))
        self.assertRaises(KeyError, df.get_value, 0, 1)
        # self.assertRaises(KeyError, df.set_value, 0, 1, 0)

    def test_single_element_ix_dont_upcast(self):
        self.frame['E'] = 1
        self.assertTrue(issubclass(self.frame['E'].dtype.type,
                                   (int, np.integer)))

        result = self.frame.ix[self.frame.index[5], 'E']
        self.assertTrue(com.is_integer(result))

    def test_irow(self):
        df = DataFrame(np.random.randn(10, 4), index=lrange(0, 20, 2))

        # 10711, deprecated
        with tm.assert_produces_warning(FutureWarning):
            df.irow(1)

        result = df.iloc[1]
        exp = df.ix[2]
        assert_series_equal(result, exp)

        result = df.iloc[2]
        exp = df.ix[4]
        assert_series_equal(result, exp)

        # slice
        result = df.iloc[slice(4, 8)]
        expected = df.ix[8:14]
        assert_frame_equal(result, expected)

        # verify slice is view
        # setting it makes it raise/warn
        def f():
            result[2] = 0.
        self.assertRaises(com.SettingWithCopyError, f)
        exp_col = df[2].copy()
        exp_col[4:8] = 0.
        assert_series_equal(df[2], exp_col)

        # list of integers
        result = df.iloc[[1, 2, 4, 6]]
        expected = df.reindex(df.index[[1, 2, 4, 6]])
        assert_frame_equal(result, expected)

    def test_icol(self):

        df = DataFrame(np.random.randn(4, 10), columns=lrange(0, 20, 2))

        # 10711, deprecated
        with tm.assert_produces_warning(FutureWarning):
            df.icol(1)

        result = df.iloc[:, 1]
        exp = df.ix[:, 2]
        assert_series_equal(result, exp)

        result = df.iloc[:, 2]
        exp = df.ix[:, 4]
        assert_series_equal(result, exp)

        # slice
        result = df.iloc[:, slice(4, 8)]
        expected = df.ix[:, 8:14]
        assert_frame_equal(result, expected)

        # verify slice is view
        # and that we are setting a copy
        def f():
            result[8] = 0.
        self.assertRaises(com.SettingWithCopyError, f)
        self.assertTrue((df[8] == 0).all())

        # list of integers
        result = df.iloc[:, [1, 2, 4, 6]]
        expected = df.reindex(columns=df.columns[[1, 2, 4, 6]])
        assert_frame_equal(result, expected)

    def test_irow_icol_duplicates(self):
        # 10711, deprecated

        df = DataFrame(np.random.rand(3, 3), columns=list('ABC'),
                       index=list('aab'))

        result = df.iloc[0]
        result2 = df.ix[0]
        tm.assertIsInstance(result, Series)
        assert_almost_equal(result.values, df.values[0])
        assert_series_equal(result, result2)

        result = df.T.iloc[:, 0]
        result2 = df.T.ix[:, 0]
        tm.assertIsInstance(result, Series)
        assert_almost_equal(result.values, df.values[0])
        assert_series_equal(result, result2)

        # multiindex
        df = DataFrame(np.random.randn(3, 3), columns=[['i', 'i', 'j'],
                                                       ['A', 'A', 'B']],
                       index=[['i', 'i', 'j'], ['X', 'X', 'Y']])
        rs = df.iloc[0]
        xp = df.ix[0]
        assert_series_equal(rs, xp)

        rs = df.iloc[:, 0]
        xp = df.T.ix[0]
        assert_series_equal(rs, xp)

        rs = df.iloc[:, [0]]
        xp = df.ix[:, [0]]
        assert_frame_equal(rs, xp)

        # #2259
        df = DataFrame([[1, 2, 3], [4, 5, 6]], columns=[1, 1, 2])
        result = df.iloc[:, [0]]
        expected = df.take([0], axis=1)
        assert_frame_equal(result, expected)

    def test_icol_sparse_propegate_fill_value(self):
        from pandas.sparse.api import SparseDataFrame
        df = SparseDataFrame({'A': [999, 1]}, default_fill_value=999)
        self.assertTrue(len(df['A'].sp_values) == len(df.iloc[:, 0].sp_values))

    def test_iget_value(self):
        # 10711 deprecated

        with tm.assert_produces_warning(FutureWarning):
            self.frame.iget_value(0,0)

        for i, row in enumerate(self.frame.index):
            for j, col in enumerate(self.frame.columns):
                result = self.frame.iat[i,j]
                expected = self.frame.at[row, col]
                assert_almost_equal(result, expected)

    def test_nested_exception(self):
        # Ignore the strange way of triggering the problem
        # (which may get fixed), it's just a way to trigger
        # the issue or reraising an outer exception without
        # a named argument
        df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [7, 8,
                       9]}).set_index(["a", "b"])
        l = list(df.index)
        l[0] = ["a", "b"]
        df.index = l

        try:
            repr(df)
        except Exception as e:
            self.assertNotEqual(type(e), UnboundLocalError)

    def test_reindex_methods(self):
        df = pd.DataFrame({'x': list(range(5))})
        target = np.array([-0.1, 0.9, 1.1, 1.5])

        for method, expected_values in [('nearest', [0, 1, 1, 2]),
                                        ('pad', [np.nan, 0, 1, 1]),
                                        ('backfill', [0, 1, 2, 2])]:
            expected = pd.DataFrame({'x': expected_values}, index=target)
            actual = df.reindex(target, method=method)
            assert_frame_equal(expected, actual)

            actual = df.reindex_like(df, method=method, tolerance=0)
            assert_frame_equal(df, actual)

            actual = df.reindex(target, method=method, tolerance=1)
            assert_frame_equal(expected, actual)

            e2 = expected[::-1]
            actual = df.reindex(target[::-1], method=method)
            assert_frame_equal(e2, actual)

            new_order = [3, 0, 2, 1]
            e2 = expected.iloc[new_order]
            actual = df.reindex(target[new_order], method=method)
            assert_frame_equal(e2, actual)

            switched_method = ('pad' if method == 'backfill'
                               else 'backfill' if method == 'pad'
                               else method)
            actual = df[::-1].reindex(target, method=switched_method)
            assert_frame_equal(expected, actual)

        expected = pd.DataFrame({'x': [0, 1, 1, np.nan]}, index=target)
        actual = df.reindex(target, method='nearest', tolerance=0.2)
        assert_frame_equal(expected, actual)

    def test_non_monotonic_reindex_methods(self):
        dr = pd.date_range('2013-08-01', periods=6, freq='B')
        data = np.random.randn(6,1)
        df = pd.DataFrame(data, index=dr, columns=list('A'))
        df_rev = pd.DataFrame(data, index=dr[[3, 4, 5] + [0, 1, 2]],
                              columns=list('A'))
        # index is not monotonic increasing or decreasing
        self.assertRaises(ValueError, df_rev.reindex, df.index, method='pad')
        self.assertRaises(ValueError, df_rev.reindex, df.index, method='ffill')
        self.assertRaises(ValueError, df_rev.reindex, df.index, method='bfill')
        self.assertRaises(ValueError, df_rev.reindex, df.index, method='nearest')

    def test_reindex_level(self):
        from itertools import permutations
        icol = ['jim', 'joe', 'jolie']

        def verify_first_level(df, level, idx):
            f = lambda val: np.nonzero(df[level] == val)[0]
            i = np.concatenate(list(map(f, idx)))
            left = df.set_index(icol).reindex(idx, level=level)
            right = df.iloc[i].set_index(icol)
            assert_frame_equal(left, right)

        def verify(df, level, idx, indexer):
            left = df.set_index(icol).reindex(idx, level=level)
            right = df.iloc[indexer].set_index(icol)
            assert_frame_equal(left, right)

        df = pd.DataFrame({'jim':list('B' * 4 + 'A' * 2 + 'C' * 3),
                           'joe':list('abcdeabcd')[::-1],
                           'jolie':[10, 20, 30] * 3,
                           'joline': np.random.randint(0, 1000, 9)})

        target = [['C', 'B', 'A'], ['F', 'C', 'A', 'D'], ['A'], ['D', 'F'],
                  ['A', 'B', 'C'], ['C', 'A', 'B'], ['C', 'B'], ['C', 'A'],
                  ['A', 'B'], ['B', 'A', 'C'], ['A', 'C', 'B']]

        for idx in target:
            verify_first_level(df, 'jim', idx)

        verify(df, 'joe', list('abcde'), [3, 2, 1, 0, 5, 4, 8, 7, 6])
        verify(df, 'joe', list('abcd'),  [3, 2, 1, 0, 5, 8, 7, 6])
        verify(df, 'joe', list('abc'),   [3, 2, 1, 8, 7, 6])
        verify(df, 'joe', list('eca'),   [1, 3, 4, 6, 8])
        verify(df, 'joe', list('edc'),   [0, 1, 4, 5, 6])
        verify(df, 'joe', list('eadbc'), [3, 0, 2, 1, 4, 5, 8, 7, 6])
        verify(df, 'joe', list('edwq'),  [0, 4, 5])
        verify(df, 'joe', list('wq'),    [])

        df = DataFrame({'jim':['mid'] * 5 + ['btm'] * 8 + ['top'] * 7,
                        'joe':['3rd'] * 2 + ['1st'] * 3 + ['2nd'] * 3 +
                              ['1st'] * 2 + ['3rd'] * 3 + ['1st'] * 2 +
                              ['3rd'] * 3 + ['2nd'] * 2,
                        # this needs to be jointly unique with jim and joe or
                        # reindexing will fail ~1.5% of the time, this works
                        # out to needing unique groups of same size as joe
                        'jolie': np.concatenate([np.random.choice(1000, x, replace=False)
                                                 for x in [2, 3, 3, 2, 3, 2, 3, 2]]),
                        'joline': np.random.randn(20).round(3) * 10})

        for idx in permutations(df['jim'].unique()):
            for i in range(3):
                verify_first_level(df, 'jim', idx[:i+1])

        i = [2,3,4,0,1,8,9,5,6,7,10,11,12,13,14,18,19,15,16,17]
        verify(df, 'joe', ['1st', '2nd', '3rd'], i)

        i = [0,1,2,3,4,10,11,12,5,6,7,8,9,15,16,17,18,19,13,14]
        verify(df, 'joe', ['3rd', '2nd', '1st'], i)

        i = [0,1,5,6,7,10,11,12,18,19,15,16,17]
        verify(df, 'joe', ['2nd', '3rd'], i)

        i = [0,1,2,3,4,10,11,12,8,9,15,16,17,13,14]
        verify(df, 'joe', ['3rd', '1st'], i)

    def test_getitem_ix_float_duplicates(self):
        df = pd.DataFrame(np.random.randn(3, 3),
                          index=[0.1, 0.2, 0.2], columns=list('abc'))
        expect = df.iloc[1:]
        tm.assert_frame_equal(df.loc[0.2], expect)
        tm.assert_frame_equal(df.ix[0.2], expect)

        expect = df.iloc[1:, 0]
        tm.assert_series_equal(df.loc[0.2, 'a'], expect)

        df.index = [1, 0.2, 0.2]
        expect = df.iloc[1:]
        tm.assert_frame_equal(df.loc[0.2], expect)
        tm.assert_frame_equal(df.ix[0.2], expect)

        expect = df.iloc[1:, 0]
        tm.assert_series_equal(df.loc[0.2, 'a'], expect)

        df = pd.DataFrame(np.random.randn(4, 3),
                          index=[1, 0.2, 0.2, 1], columns=list('abc'))
        expect = df.iloc[1:-1]
        tm.assert_frame_equal(df.loc[0.2], expect)
        tm.assert_frame_equal(df.ix[0.2], expect)

        expect = df.iloc[1:-1, 0]
        tm.assert_series_equal(df.loc[0.2, 'a'], expect)

        df.index = [0.1, 0.2, 2, 0.2]
        expect = df.iloc[[1, -1]]
        tm.assert_frame_equal(df.loc[0.2], expect)
        tm.assert_frame_equal(df.ix[0.2], expect)

        expect = df.iloc[[1, -1], 0]
        tm.assert_series_equal(df.loc[0.2, 'a'], expect)

    def test_setitem_with_sparse_value(self):
        # GH8131
        df = pd.DataFrame({'c_1':['a', 'b', 'c'], 'n_1': [1., 2., 3.]})
        sp_series = pd.Series([0, 0, 1]).to_sparse(fill_value=0)
        df['new_column'] = sp_series
        tm.assert_series_equal(df['new_column'], sp_series, check_names=False)

    def test_setitem_with_unaligned_sparse_value(self):
        df = pd.DataFrame({'c_1':['a', 'b', 'c'], 'n_1': [1., 2., 3.]})
        sp_series = (pd.Series([0, 0, 1], index=[2, 1, 0])
                     .to_sparse(fill_value=0))
        df['new_column'] = sp_series
        exp = pd.Series([1, 0, 0], name='new_column')
        tm.assert_series_equal(df['new_column'], exp)


_seriesd = tm.getSeriesData()
_tsd = tm.getTimeSeriesData()

_frame = DataFrame(_seriesd)
_frame2 = DataFrame(_seriesd, columns=['D', 'C', 'B', 'A'])
_intframe = DataFrame(dict((k, v.astype(int))
                           for k, v in compat.iteritems(_seriesd)))

_tsframe = DataFrame(_tsd)

_mixed_frame = _frame.copy()
_mixed_frame['foo'] = 'bar'


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


class TestDataFrame(tm.TestCase, CheckIndexing,
                    SafeForSparse):
    klass = DataFrame

    _multiprocess_can_split_ = True

    def setUp(self):
        import warnings

        self.frame = _frame.copy()
        self.frame2 = _frame2.copy()

        # force these all to int64 to avoid platform testing issues
        self.intframe = DataFrame(dict([ (c,s) for c,s in compat.iteritems(_intframe) ]), dtype = np.int64)
        self.tsframe = _tsframe.copy()
        self.mixed_frame = _mixed_frame.copy()
        self.mixed_float  = DataFrame({ 'A': _frame['A'].copy().astype('float32'),
                                        'B': _frame['B'].copy().astype('float32'),
                                        'C': _frame['C'].copy().astype('float16'),
                                        'D': _frame['D'].copy().astype('float64') })
        self.mixed_float2 = DataFrame({ 'A': _frame2['A'].copy().astype('float32'),
                                        'B': _frame2['B'].copy().astype('float32'),
                                        'C': _frame2['C'].copy().astype('float16'),
                                        'D': _frame2['D'].copy().astype('float64') })
        self.mixed_int    = DataFrame({ 'A': _intframe['A'].copy().astype('int32'),
                                        'B': np.ones(len(_intframe['B']),dtype='uint64'),
                                        'C': _intframe['C'].copy().astype('uint8'),
                                        'D': _intframe['D'].copy().astype('int64') })
        self.all_mixed    = DataFrame({'a': 1., 'b': 2, 'c': 'foo', 'float32' : np.array([1.]*10,dtype='float32'),
                                       'int32' : np.array([1]*10,dtype='int32'),
                                       }, index=np.arange(10))
        self.tzframe = DataFrame({'A' : date_range('20130101',periods=3),
                                  'B' : date_range('20130101',periods=3,tz='US/Eastern'),
                                  'C' : date_range('20130101',periods=3,tz='CET')})
        self.tzframe.iloc[1,1] = pd.NaT
        self.tzframe.iloc[1,2] = pd.NaT

        self.ts1 = tm.makeTimeSeries()
        self.ts2 = tm.makeTimeSeries()[5:]
        self.ts3 = tm.makeTimeSeries()[-5:]
        self.ts4 = tm.makeTimeSeries()[1:-1]

        self.ts_dict = {
            'col1': self.ts1,
            'col2': self.ts2,
            'col3': self.ts3,
            'col4': self.ts4,
        }
        self.empty = DataFrame({})

        arr = np.array([[1., 2., 3.],
                        [4., 5., 6.],
                        [7., 8., 9.]])

        self.simple = DataFrame(arr, columns=['one', 'two', 'three'],
                                index=['a', 'b', 'c'])

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
        assertRaisesRegexp(ValueError, 'No axis named', f._get_axis_number, None)

    def test_set_index(self):
        idx = Index(np.arange(len(self.mixed_frame)))

        # cache it
        _ = self.mixed_frame['foo']
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
        assert_frame_equal(result, expected, check_names=False) # TODO should set_index check_names ?

    def test_construction_with_categorical_index(self):

        ci = tm.makeCategoricalIndex(10)

        # with Categorical
        df = DataFrame({'A' : np.random.randn(10),
                        'B' : ci.values })
        idf = df.set_index('B')
        str(idf)
        tm.assert_index_equal(idf.index, ci, check_names=False)
        self.assertEqual(idf.index.name, 'B')

        # from a CategoricalIndex
        df = DataFrame({'A' : np.random.randn(10),
                        'B' : ci })
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

    def test_constructor(self):
        df = DataFrame()
        self.assertEqual(len(df.index), 0)

        df = DataFrame(data={})
        self.assertEqual(len(df.index), 0)

    def test_constructor_mixed(self):
        index, data = tm.getMixedTypeDict()

        indexed_frame = DataFrame(data, index=index)
        unindexed_frame = DataFrame(data)

        self.assertEqual(self.mixed_frame['foo'].dtype, np.object_)

    def test_constructor_cast_failure(self):
        foo = DataFrame({'a': ['a', 'b', 'c']}, dtype=np.float64)
        self.assertEqual(foo['a'].dtype, object)

        # GH 3010, constructing with odd arrays
        df = DataFrame(np.ones((4,2)))

        # this is ok
        df['foo'] = np.ones((4,2)).tolist()

        # this is not ok
        self.assertRaises(ValueError, df.__setitem__, tuple(['test']), np.ones((4,2)))

        # this is ok
        df['foo2'] = np.ones((4,2)).tolist()

    def test_constructor_dtype_copy(self):
        orig_df = DataFrame({
            'col1': [1.],
            'col2': [2.],
            'col3': [3.]})

        new_df = pd.DataFrame(orig_df, dtype=float, copy=True)

        new_df['col1'] = 200.
        self.assertEqual(orig_df['col1'][0], 1.)

    def test_constructor_dtype_nocast_view(self):
        df = DataFrame([[1, 2]])
        should_be_view = DataFrame(df, dtype=df[0].dtype)
        should_be_view[0][0] = 99
        self.assertEqual(df.values[0, 0], 99)

        should_be_view = DataFrame(df.values, dtype=df[0].dtype)
        should_be_view[0][0] = 97
        self.assertEqual(df.values[0, 0], 97)

    def test_constructor_dtype_list_data(self):
        df = DataFrame([[1, '2'],
                        [None, 'a']], dtype=object)
        self.assertIsNone(df.ix[1, 0])
        self.assertEqual(df.ix[0, 1], '2')

    def test_constructor_list_frames(self):

        # GH 3243
        result = DataFrame([DataFrame([])])
        self.assertEqual(result.shape, (1,0))

        result = DataFrame([DataFrame(dict(A = lrange(5)))])
        tm.assertIsInstance(result.iloc[0,0], DataFrame)

    def test_constructor_mixed_dtypes(self):

        def _make_mixed_dtypes_df(typ, ad = None):

            if typ == 'int':
                dtypes = MIXED_INT_DTYPES
                arrays = [ np.array(np.random.rand(10), dtype = d) for d in dtypes ]
            elif typ == 'float':
                dtypes = MIXED_FLOAT_DTYPES
                arrays = [ np.array(np.random.randint(10, size=10), dtype = d) for d in dtypes ]

            zipper = lzip(dtypes,arrays)
            for d,a in zipper:
                assert(a.dtype == d)
            if ad is None:
                ad = dict()
            ad.update(dict([ (d,a) for d,a in zipper ]))
            return DataFrame(ad)

        def _check_mixed_dtypes(df, dtypes = None):
            if dtypes is None:
                dtypes = MIXED_FLOAT_DTYPES + MIXED_INT_DTYPES
            for d in dtypes:
                if d in df:
                    assert(df.dtypes[d] == d)

        # mixed floating and integer coexinst in the same frame
        df     = _make_mixed_dtypes_df('float')
        _check_mixed_dtypes(df)

        # add lots of types
        df     = _make_mixed_dtypes_df('float', dict(A = 1, B = 'foo', C = 'bar'))
        _check_mixed_dtypes(df)

        # GH 622
        df     = _make_mixed_dtypes_df('int')
        _check_mixed_dtypes(df)

    def test_constructor_complex_dtypes(self):
        # GH10952
        a = np.random.rand(10).astype(np.complex64)
        b = np.random.rand(10).astype(np.complex128)

        df = DataFrame({'a': a, 'b': b})
        self.assertEqual(a.dtype, df.a.dtype)
        self.assertEqual(b.dtype, df.b.dtype)

    def test_constructor_rec(self):
        rec = self.frame.to_records(index=False)

        # Assigning causes segfault in NumPy < 1.5.1
        # rec.dtype.names = list(rec.dtype.names)[::-1]

        index = self.frame.index

        df = DataFrame(rec)
        self.assert_numpy_array_equal(df.columns, rec.dtype.names)

        df2 = DataFrame(rec, index=index)
        self.assert_numpy_array_equal(df2.columns, rec.dtype.names)
        self.assertTrue(df2.index.equals(index))

        rng = np.arange(len(rec))[::-1]
        df3 = DataFrame(rec, index=rng, columns=['C', 'B'])
        expected = DataFrame(rec, index=rng).reindex(columns=['C', 'B'])
        assert_frame_equal(df3, expected)

    def test_constructor_bool(self):
        df = DataFrame({0: np.ones(10, dtype=bool),
                        1: np.zeros(10, dtype=bool)})
        self.assertEqual(df.values.dtype, np.bool_)

    def test_constructor_overflow_int64(self):
        values = np.array([2 ** 64 - i for i in range(1, 10)],
                          dtype=np.uint64)

        result = DataFrame({'a': values})
        self.assertEqual(result['a'].dtype, object)

        # #2355
        data_scores = [(6311132704823138710, 273), (2685045978526272070, 23),
                       (8921811264899370420, 45), (long(17019687244989530680), 270),
                       (long(9930107427299601010), 273)]
        dtype = [('uid', 'u8'), ('score', 'u8')]
        data = np.zeros((len(data_scores),), dtype=dtype)
        data[:] = data_scores
        df_crawls = DataFrame(data)
        self.assertEqual(df_crawls['uid'].dtype, object)

    def test_constructor_ordereddict(self):
        import random
        nitems = 100
        nums = lrange(nitems)
        random.shuffle(nums)
        expected = ['A%d' % i for i in nums]
        df = DataFrame(OrderedDict(zip(expected, [[0]] * nitems)))
        self.assertEqual(expected, list(df.columns))

    def test_constructor_dict(self):
        frame = DataFrame({'col1': self.ts1,
                           'col2': self.ts2})

        tm.assert_dict_equal(self.ts1, frame['col1'], compare_keys=False)
        tm.assert_dict_equal(self.ts2, frame['col2'], compare_keys=False)

        frame = DataFrame({'col1': self.ts1,
                           'col2': self.ts2},
                          columns=['col2', 'col3', 'col4'])

        self.assertEqual(len(frame), len(self.ts2))
        self.assertNotIn('col1', frame)
        self.assertTrue(isnull(frame['col3']).all())

        # Corner cases
        self.assertEqual(len(DataFrame({})), 0)

        # mix dict and array, wrong size - no spec for which error should raise
        # first
        with tm.assertRaises(ValueError):
            DataFrame({'A': {'a': 'a', 'b': 'b'}, 'B': ['a', 'b', 'c']})

        # Length-one dict micro-optimization
        frame = DataFrame({'A': {'1': 1, '2': 2}})
        self.assert_numpy_array_equal(frame.index, ['1', '2'])

        # empty dict plus index
        idx = Index([0, 1, 2])
        frame = DataFrame({}, index=idx)
        self.assertIs(frame.index, idx)

        # empty with index and columns
        idx = Index([0, 1, 2])
        frame = DataFrame({}, index=idx, columns=idx)
        self.assertIs(frame.index, idx)
        self.assertIs(frame.columns, idx)
        self.assertEqual(len(frame._series), 3)

        # with dict of empty list and Series
        frame = DataFrame({'A': [], 'B': []}, columns=['A', 'B'])
        self.assertTrue(frame.index.equals(Index([])))

        # GH10856
        # dict with scalar values should raise error, even if columns passed
        with tm.assertRaises(ValueError):
            DataFrame({'a': 0.7})

        with tm.assertRaises(ValueError):
            DataFrame({'a': 0.7}, columns=['a'])

        with tm.assertRaises(ValueError):
            DataFrame({'a': 0.7}, columns=['b'])

    def test_constructor_multi_index(self):
        # GH 4078
        # construction error with mi and all-nan frame
        tuples = [(2, 3), (3, 3), (3, 3)]
        mi = MultiIndex.from_tuples(tuples)
        df = DataFrame(index=mi,columns=mi)
        self.assertTrue(pd.isnull(df).values.ravel().all())

        tuples = [(3, 3), (2, 3), (3, 3)]
        mi = MultiIndex.from_tuples(tuples)
        df = DataFrame(index=mi,columns=mi)
        self.assertTrue(pd.isnull(df).values.ravel().all())

    def test_constructor_error_msgs(self):
        msg = "Mixing dicts with non-Series may lead to ambiguous ordering."
        # mix dict and array, wrong size
        with assertRaisesRegexp(ValueError, msg):
            DataFrame({'A': {'a': 'a', 'b': 'b'},
                       'B': ['a', 'b', 'c']})

        # wrong size ndarray, GH 3105
        msg = "Shape of passed values is \(3, 4\), indices imply \(3, 3\)"
        with assertRaisesRegexp(ValueError, msg):
            DataFrame(np.arange(12).reshape((4, 3)),
                      columns=['foo', 'bar', 'baz'],
                      index=date_range('2000-01-01', periods=3))


        # higher dim raise exception
        with assertRaisesRegexp(ValueError, 'Must pass 2-d input'):
            DataFrame(np.zeros((3, 3, 3)), columns=['A', 'B', 'C'], index=[1])

        # wrong size axis labels
        with assertRaisesRegexp(ValueError, "Shape of passed values is \(3, 2\), indices imply \(3, 1\)"):
            DataFrame(np.random.rand(2,3), columns=['A', 'B', 'C'], index=[1])

        with assertRaisesRegexp(ValueError, "Shape of passed values is \(3, 2\), indices imply \(2, 2\)"):
            DataFrame(np.random.rand(2,3), columns=['A', 'B'], index=[1, 2])

        with assertRaisesRegexp(ValueError, 'If using all scalar values, you must pass an index'):
            DataFrame({'a': False, 'b': True})

    def test_constructor_with_embedded_frames(self):

        # embedded data frames
        df1 = DataFrame({'a':[1, 2, 3], 'b':[3, 4, 5]})
        df2 = DataFrame([df1, df1+10])

        df2.dtypes
        str(df2)

        result = df2.loc[0,0]
        assert_frame_equal(result,df1)

        result = df2.loc[1,0]
        assert_frame_equal(result,df1+10)

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

    def test_constructor_subclass_dict(self):
        # Test for passing dict subclass to constructor
        data = {'col1': tm.TestSubDict((x, 10.0 * x) for x in range(10)),
                'col2': tm.TestSubDict((x, 20.0 * x) for x in range(10))}
        df = DataFrame(data)
        refdf = DataFrame(dict((col, dict(compat.iteritems(val)))
                               for col, val in compat.iteritems(data)))
        assert_frame_equal(refdf, df)

        data = tm.TestSubDict(compat.iteritems(data))
        df = DataFrame(data)
        assert_frame_equal(refdf, df)

        # try with defaultdict
        from collections import defaultdict
        data = {}
        self.frame['B'][:10] = np.nan
        for k, v in compat.iteritems(self.frame):
            dct = defaultdict(dict)
            dct.update(v.to_dict())
            data[k] = dct
        frame = DataFrame(data)
        assert_frame_equal(self.frame.sort_index(), frame)

    def test_constructor_dict_block(self):
        expected = [[4., 3., 2., 1.]]
        df = DataFrame({'d': [4.], 'c': [3.], 'b': [2.], 'a': [1.]},
                       columns=['d', 'c', 'b', 'a'])
        assert_almost_equal(df.values, expected)

    def test_constructor_dict_cast(self):
        # cast float tests
        test_data = {
            'A': {'1': 1, '2': 2},
            'B': {'1': '1', '2': '2', '3': '3'},
        }
        frame = DataFrame(test_data, dtype=float)
        self.assertEqual(len(frame), 3)
        self.assertEqual(frame['B'].dtype, np.float64)
        self.assertEqual(frame['A'].dtype, np.float64)

        frame = DataFrame(test_data)
        self.assertEqual(len(frame), 3)
        self.assertEqual(frame['B'].dtype, np.object_)
        self.assertEqual(frame['A'].dtype, np.float64)

        # can't cast to float
        test_data = {
            'A': dict(zip(range(20), tm.makeStringIndex(20))),
            'B': dict(zip(range(15), randn(15)))
        }
        frame = DataFrame(test_data, dtype=float)
        self.assertEqual(len(frame), 20)
        self.assertEqual(frame['A'].dtype, np.object_)
        self.assertEqual(frame['B'].dtype, np.float64)

    def test_constructor_dict_dont_upcast(self):
        d = {'Col1': {'Row1': 'A String', 'Row2': np.nan}}
        df = DataFrame(d)
        tm.assertIsInstance(df['Col1']['Row2'], float)

        dm = DataFrame([[1, 2], ['a', 'b']], index=[1, 2], columns=[1, 2])
        tm.assertIsInstance(dm[1][1], int)

    def test_constructor_dict_of_tuples(self):
        # GH #1491
        data = {'a': (1, 2, 3), 'b': (4, 5, 6)}

        result = DataFrame(data)
        expected = DataFrame(dict((k, list(v)) for k, v in compat.iteritems(data)))
        assert_frame_equal(result, expected, check_dtype=False)

    def test_constructor_dict_multiindex(self):
        check = lambda result, expected: tm.assert_frame_equal(
            result, expected, check_dtype=True, check_index_type=True,
            check_column_type=True, check_names=True)
        d = {('a', 'a'): {('i', 'i'): 0, ('i', 'j'): 1, ('j', 'i'): 2},
             ('b', 'a'): {('i', 'i'): 6, ('i', 'j'): 5, ('j', 'i'): 4},
             ('b', 'c'): {('i', 'i'): 7, ('i', 'j'): 8, ('j', 'i'): 9}}
        _d = sorted(d.items())
        df = DataFrame(d)
        expected = DataFrame(
            [x[1] for x in _d],
            index=MultiIndex.from_tuples([x[0] for x in _d])).T
        expected.index = MultiIndex.from_tuples(expected.index)
        check(df, expected)

        d['z'] = {'y': 123., ('i', 'i'): 111, ('i', 'j'): 111, ('j', 'i'): 111}
        _d.insert(0, ('z', d['z']))
        expected = DataFrame(
            [x[1] for x in _d],
            index=Index([x[0] for x in _d], tupleize_cols=False)).T
        expected.index = Index(expected.index, tupleize_cols=False)
        df = DataFrame(d)
        df = df.reindex(columns=expected.columns, index=expected.index)
        check(df, expected)

    def test_constructor_dict_datetime64_index(self):
        # GH 10160
        dates_as_str = ['1984-02-19', '1988-11-06', '1989-12-03', '1990-03-15']

        def create_data(constructor):
            return dict((i, {constructor(s): 2*i}) for i, s in enumerate(dates_as_str))

        data_datetime64 = create_data(np.datetime64)
        data_datetime = create_data(lambda x: datetime.strptime(x, '%Y-%m-%d'))
        data_Timestamp = create_data(Timestamp)

        expected = DataFrame([{0: 0, 1: None, 2: None, 3: None},
                              {0: None, 1: 2, 2: None, 3: None},
                              {0: None, 1: None, 2: 4, 3: None},
                              {0: None, 1: None, 2: None, 3: 6}],
                             index=[Timestamp(dt) for dt in dates_as_str])

        result_datetime64 = DataFrame(data_datetime64)
        result_datetime = DataFrame(data_datetime)
        result_Timestamp = DataFrame(data_Timestamp)
        assert_frame_equal(result_datetime64, expected)
        assert_frame_equal(result_datetime, expected)
        assert_frame_equal(result_Timestamp, expected)

    def test_constructor_dict_timedelta64_index(self):
        # GH 10160
        td_as_int = [1, 2, 3, 4]

        def create_data(constructor):
            return dict((i, {constructor(s): 2*i}) for i, s in enumerate(td_as_int))

        data_timedelta64 = create_data(lambda x: np.timedelta64(x, 'D'))
        data_timedelta = create_data(lambda x: timedelta(days=x))
        data_Timedelta = create_data(lambda x: Timedelta(x, 'D'))

        expected = DataFrame([{0: 0, 1: None, 2: None, 3: None},
                              {0: None, 1: 2, 2: None, 3: None},
                              {0: None, 1: None, 2: 4, 3: None},
                              {0: None, 1: None, 2: None, 3: 6}],
                             index=[Timedelta(td, 'D') for td in td_as_int])

        result_timedelta64 = DataFrame(data_timedelta64)
        result_timedelta = DataFrame(data_timedelta)
        result_Timedelta = DataFrame(data_Timedelta)
        assert_frame_equal(result_timedelta64, expected)
        assert_frame_equal(result_timedelta, expected)
        assert_frame_equal(result_Timedelta, expected)

    def test_nested_dict_frame_constructor(self):
        rng = period_range('1/1/2000', periods=5)
        df = DataFrame(randn(10, 5), columns=rng)

        data = {}
        for col in df.columns:
            for row in df.index:
                data.setdefault(col, {})[row] = df.get_value(row, col)

        result = DataFrame(data, columns=rng)
        tm.assert_frame_equal(result, df)

        data = {}
        for col in df.columns:
            for row in df.index:
                data.setdefault(row, {})[col] = df.get_value(row, col)

        result = DataFrame(data, index=rng).T
        tm.assert_frame_equal(result, df)


    def _check_basic_constructor(self, empty):
        "mat: 2d matrix with shpae (3, 2) to input. empty - makes sized objects"
        mat = empty((2, 3), dtype=float)
        # 2-D input
        frame = DataFrame(mat, columns=['A', 'B', 'C'], index=[1, 2])

        self.assertEqual(len(frame.index), 2)
        self.assertEqual(len(frame.columns), 3)

        # 1-D input
        frame = DataFrame(empty((3,)), columns=['A'], index=[1, 2, 3])
        self.assertEqual(len(frame.index), 3)
        self.assertEqual(len(frame.columns), 1)


        # cast type
        frame = DataFrame(mat, columns=['A', 'B', 'C'],
                          index=[1, 2], dtype=np.int64)
        self.assertEqual(frame.values.dtype, np.int64)

        # wrong size axis labels
        msg = r'Shape of passed values is \(3, 2\), indices imply \(3, 1\)'
        with assertRaisesRegexp(ValueError, msg):
            DataFrame(mat, columns=['A', 'B', 'C'], index=[1])
        msg = r'Shape of passed values is \(3, 2\), indices imply \(2, 2\)'
        with assertRaisesRegexp(ValueError, msg):
            DataFrame(mat, columns=['A', 'B'], index=[1, 2])

        # higher dim raise exception
        with assertRaisesRegexp(ValueError, 'Must pass 2-d input'):
            DataFrame(empty((3, 3, 3)), columns=['A', 'B', 'C'],
                      index=[1])

        # automatic labeling
        frame = DataFrame(mat)
        self.assert_numpy_array_equal(frame.index, lrange(2))
        self.assert_numpy_array_equal(frame.columns, lrange(3))

        frame = DataFrame(mat, index=[1, 2])
        self.assert_numpy_array_equal(frame.columns, lrange(3))

        frame = DataFrame(mat, columns=['A', 'B', 'C'])
        self.assert_numpy_array_equal(frame.index, lrange(2))

        # 0-length axis
        frame = DataFrame(empty((0, 3)))
        self.assertEqual(len(frame.index), 0)

        frame = DataFrame(empty((3, 0)))
        self.assertEqual(len(frame.columns), 0)

    def test_constructor_ndarray(self):
        mat = np.zeros((2, 3), dtype=float)
        self._check_basic_constructor(np.ones)

        frame = DataFrame(['foo', 'bar'], index=[0, 1], columns=['A'])
        self.assertEqual(len(frame), 2)

    def test_constructor_maskedarray(self):
        self._check_basic_constructor(ma.masked_all)

        # Check non-masked values
        mat = ma.masked_all((2, 3), dtype=float)
        mat[0, 0] = 1.0
        mat[1, 2] = 2.0
        frame = DataFrame(mat, columns=['A', 'B', 'C'], index=[1, 2])
        self.assertEqual(1.0, frame['A'][1])
        self.assertEqual(2.0, frame['C'][2])

        # what is this even checking??
        mat = ma.masked_all((2, 3), dtype=float)
        frame = DataFrame(mat, columns=['A', 'B', 'C'], index=[1, 2])
        self.assertTrue(np.all(~np.asarray(frame == frame)))

    def test_constructor_maskedarray_nonfloat(self):
        # masked int promoted to float
        mat = ma.masked_all((2, 3), dtype=int)
        # 2-D input
        frame = DataFrame(mat, columns=['A', 'B', 'C'], index=[1, 2])

        self.assertEqual(len(frame.index), 2)
        self.assertEqual(len(frame.columns), 3)
        self.assertTrue(np.all(~np.asarray(frame == frame)))

        # cast type
        frame = DataFrame(mat, columns=['A', 'B', 'C'],
                          index=[1, 2], dtype=np.float64)
        self.assertEqual(frame.values.dtype, np.float64)

        # Check non-masked values
        mat2 = ma.copy(mat)
        mat2[0, 0] = 1
        mat2[1, 2] = 2
        frame = DataFrame(mat2, columns=['A', 'B', 'C'], index=[1, 2])
        self.assertEqual(1, frame['A'][1])
        self.assertEqual(2, frame['C'][2])

        # masked np.datetime64 stays (use lib.NaT as null)
        mat = ma.masked_all((2, 3), dtype='M8[ns]')
        # 2-D input
        frame = DataFrame(mat, columns=['A', 'B', 'C'], index=[1, 2])

        self.assertEqual(len(frame.index), 2)
        self.assertEqual(len(frame.columns), 3)
        self.assertTrue(isnull(frame).values.all())

        # cast type
        frame = DataFrame(mat, columns=['A', 'B', 'C'],
                          index=[1, 2], dtype=np.int64)
        self.assertEqual(frame.values.dtype, np.int64)

        # Check non-masked values
        mat2 = ma.copy(mat)
        mat2[0, 0] = 1
        mat2[1, 2] = 2
        frame = DataFrame(mat2, columns=['A', 'B', 'C'], index=[1, 2])
        self.assertEqual(1, frame['A'].view('i8')[1])
        self.assertEqual(2, frame['C'].view('i8')[2])

        # masked bool promoted to object
        mat = ma.masked_all((2, 3), dtype=bool)
        # 2-D input
        frame = DataFrame(mat, columns=['A', 'B', 'C'], index=[1, 2])

        self.assertEqual(len(frame.index), 2)
        self.assertEqual(len(frame.columns), 3)
        self.assertTrue(np.all(~np.asarray(frame == frame)))

        # cast type
        frame = DataFrame(mat, columns=['A', 'B', 'C'],
                          index=[1, 2], dtype=object)
        self.assertEqual(frame.values.dtype, object)

        # Check non-masked values
        mat2 = ma.copy(mat)
        mat2[0, 0] = True
        mat2[1, 2] = False
        frame = DataFrame(mat2, columns=['A', 'B', 'C'], index=[1, 2])
        self.assertEqual(True, frame['A'][1])
        self.assertEqual(False, frame['C'][2])

    def test_constructor_mrecarray(self):
        # Ensure mrecarray produces frame identical to dict of masked arrays
        # from GH3479

        assert_fr_equal = functools.partial(assert_frame_equal,
                                            check_index_type=True,
                                            check_column_type=True,
                                            check_frame_type=True)
        arrays = [
            ('float', np.array([1.5, 2.0])),
            ('int', np.array([1, 2])),
            ('str', np.array(['abc', 'def'])),
        ]
        for name, arr in arrays[:]:
            arrays.append(('masked1_' + name,
                           np.ma.masked_array(arr, mask=[False, True])))
        arrays.append(('masked_all', np.ma.masked_all((2,))))
        arrays.append(('masked_none',
                       np.ma.masked_array([1.0, 2.5], mask=False)))

        # call assert_frame_equal for all selections of 3 arrays
        for comb in itertools.combinations(arrays, 3):
            names, data = zip(*comb)
            mrecs = mrecords.fromarrays(data, names=names)

            # fill the comb
            comb = dict([ (k, v.filled()) if hasattr(v,'filled') else (k, v) for k, v in comb ])

            expected = DataFrame(comb,columns=names)
            result = DataFrame(mrecs)
            assert_fr_equal(result,expected)

            # specify columns
            expected = DataFrame(comb,columns=names[::-1])
            result = DataFrame(mrecs, columns=names[::-1])
            assert_fr_equal(result,expected)

            # specify index
            expected = DataFrame(comb,columns=names,index=[1,2])
            result = DataFrame(mrecs, index=[1,2])
            assert_fr_equal(result,expected)

    def test_constructor_corner(self):
        df = DataFrame(index=[])
        self.assertEqual(df.values.shape, (0, 0))

        # empty but with specified dtype
        df = DataFrame(index=lrange(10), columns=['a', 'b'], dtype=object)
        self.assertEqual(df.values.dtype, np.object_)

        # does not error but ends up float
        df = DataFrame(index=lrange(10), columns=['a', 'b'], dtype=int)
        self.assertEqual(df.values.dtype, np.object_)

        # #1783 empty dtype object
        df = DataFrame({}, columns=['foo', 'bar'])
        self.assertEqual(df.values.dtype, np.object_)

        df = DataFrame({'b': 1}, index=lrange(10), columns=list('abc'),
                       dtype=int)
        self.assertEqual(df.values.dtype, np.object_)


    def test_constructor_scalar_inference(self):
        data = {'int': 1, 'bool': True,
                'float': 3., 'complex': 4j, 'object': 'foo'}
        df = DataFrame(data, index=np.arange(10))

        self.assertEqual(df['int'].dtype, np.int64)
        self.assertEqual(df['bool'].dtype, np.bool_)
        self.assertEqual(df['float'].dtype, np.float64)
        self.assertEqual(df['complex'].dtype, np.complex128)
        self.assertEqual(df['object'].dtype, np.object_)

    def test_constructor_arrays_and_scalars(self):
        df = DataFrame({'a': randn(10), 'b': True})
        exp = DataFrame({'a': df['a'].values, 'b': [True] * 10})

        assert_frame_equal(df, exp)
        with tm.assertRaisesRegexp(ValueError, 'must pass an index'):
            DataFrame({'a': False, 'b': True})

    def test_constructor_DataFrame(self):
        df = DataFrame(self.frame)
        assert_frame_equal(df, self.frame)

        df_casted = DataFrame(self.frame, dtype=np.int64)
        self.assertEqual(df_casted.values.dtype, np.int64)

    def test_constructor_more(self):
        # used to be in test_matrix.py
        arr = randn(10)
        dm = DataFrame(arr, columns=['A'], index=np.arange(10))
        self.assertEqual(dm.values.ndim, 2)

        arr = randn(0)
        dm = DataFrame(arr)
        self.assertEqual(dm.values.ndim, 2)
        self.assertEqual(dm.values.ndim, 2)

        # no data specified
        dm = DataFrame(columns=['A', 'B'], index=np.arange(10))
        self.assertEqual(dm.values.shape, (10, 2))

        dm = DataFrame(columns=['A', 'B'])
        self.assertEqual(dm.values.shape, (0, 2))

        dm = DataFrame(index=np.arange(10))
        self.assertEqual(dm.values.shape, (10, 0))

        # corner, silly
        # TODO: Fix this Exception to be better...
        with assertRaisesRegexp(PandasError, 'constructor not properly called'):
            DataFrame((1, 2, 3))

        # can't cast
        mat = np.array(['foo', 'bar'], dtype=object).reshape(2, 1)
        with assertRaisesRegexp(ValueError, 'cast'):
            DataFrame(mat, index=[0, 1], columns=[0], dtype=float)

        dm = DataFrame(DataFrame(self.frame._series))
        tm.assert_frame_equal(dm, self.frame)

        # int cast
        dm = DataFrame({'A': np.ones(10, dtype=int),
                        'B': np.ones(10, dtype=np.float64)},
                       index=np.arange(10))

        self.assertEqual(len(dm.columns), 2)
        self.assertEqual(dm.values.dtype, np.float64)

    def test_constructor_empty_list(self):
        df = DataFrame([], index=[])
        expected = DataFrame(index=[])
        assert_frame_equal(df, expected)

        # GH 9939
        df = DataFrame([], columns=['A', 'B'])
        expected = DataFrame({}, columns=['A', 'B'])
        assert_frame_equal(df, expected)

        # Empty generator: list(empty_gen()) == []
        def empty_gen():
            return
            yield

        df = DataFrame(empty_gen(), columns=['A', 'B'])
        assert_frame_equal(df, expected)

    def test_constructor_list_of_lists(self):
        # GH #484
        l = [[1, 'a'], [2, 'b']]
        df = DataFrame(data=l, columns=["num", "str"])
        self.assertTrue(com.is_integer_dtype(df['num']))
        self.assertEqual(df['str'].dtype, np.object_)

        # GH 4851
        # list of 0-dim ndarrays
        expected = DataFrame({ 0: range(10) })
        data = [np.array(x) for x in range(10)]
        result = DataFrame(data)
        assert_frame_equal(result, expected)

    def test_constructor_sequence_like(self):
        # GH 3783
        # collections.Squence like
        import collections

        class DummyContainer(collections.Sequence):
            def __init__(self, lst):
                self._lst = lst
            def __getitem__(self, n):
                return self._lst.__getitem__(n)
            def __len__(self, n):
                return self._lst.__len__()

        l = [DummyContainer([1, 'a']), DummyContainer([2, 'b'])]
        columns = ["num", "str"]
        result = DataFrame(l, columns=columns)
        expected = DataFrame([[1,'a'],[2,'b']],columns=columns)
        assert_frame_equal(result, expected, check_dtype=False)

        # GH 4297
        # support Array
        import array
        result = DataFrame.from_items([('A', array.array('i', range(10)))])
        expected = DataFrame({ 'A' : list(range(10)) })
        assert_frame_equal(result, expected, check_dtype=False)

        expected = DataFrame([ list(range(10)), list(range(10)) ])
        result = DataFrame([ array.array('i', range(10)), array.array('i',range(10)) ])
        assert_frame_equal(result, expected, check_dtype=False)

    def test_constructor_iterator(self):

        expected = DataFrame([ list(range(10)), list(range(10)) ])
        result = DataFrame([ range(10), range(10) ])
        assert_frame_equal(result, expected)

    def test_constructor_generator(self):
        #related #2305

        gen1 = (i for i in range(10))
        gen2 = (i for i in range(10))

        expected = DataFrame([ list(range(10)), list(range(10)) ])
        result = DataFrame([ gen1, gen2 ])
        assert_frame_equal(result, expected)

        gen = ([ i, 'a'] for i in range(10))
        result = DataFrame(gen)
        expected = DataFrame({ 0 : range(10), 1 : 'a' })
        assert_frame_equal(result, expected, check_dtype=False)

    def test_constructor_list_of_dicts(self):
        data = [OrderedDict([['a', 1.5], ['b', 3], ['c', 4], ['d', 6]]),
                OrderedDict([['a', 1.5], ['b', 3], ['d', 6]]),
                OrderedDict([['a', 1.5], ['d', 6]]),
                OrderedDict(),
                OrderedDict([['a', 1.5], ['b', 3], ['c', 4]]),
                OrderedDict([['b', 3], ['c', 4], ['d', 6]])]

        result = DataFrame(data)
        expected = DataFrame.from_dict(dict(zip(range(len(data)), data)),
                                       orient='index')
        assert_frame_equal(result, expected.reindex(result.index))

        result = DataFrame([{}])
        expected = DataFrame(index=[0])
        assert_frame_equal(result, expected)

    def test_constructor_list_of_series(self):
        data = [OrderedDict([['a', 1.5], ['b', 3.0], ['c', 4.0]]),
                OrderedDict([['a', 1.5], ['b', 3.0], ['c', 6.0]])]
        sdict = OrderedDict(zip(['x', 'y'], data))
        idx = Index(['a', 'b', 'c'])

        # all named
        data2 = [Series([1.5, 3, 4], idx, dtype='O', name='x'),
                 Series([1.5, 3, 6], idx, name='y')]
        result = DataFrame(data2)
        expected = DataFrame.from_dict(sdict, orient='index')
        assert_frame_equal(result, expected)

        # some unnamed
        data2 = [Series([1.5, 3, 4], idx, dtype='O', name='x'),
                 Series([1.5, 3, 6], idx)]
        result = DataFrame(data2)

        sdict = OrderedDict(zip(['x', 'Unnamed 0'], data))
        expected = DataFrame.from_dict(sdict, orient='index')
        assert_frame_equal(result.sort_index(), expected)

        # none named
        data = [OrderedDict([['a', 1.5], ['b', 3], ['c', 4], ['d', 6]]),
                OrderedDict([['a', 1.5], ['b', 3], ['d', 6]]),
                OrderedDict([['a', 1.5], ['d', 6]]),
                OrderedDict(),
                OrderedDict([['a', 1.5], ['b', 3], ['c', 4]]),
                OrderedDict([['b', 3], ['c', 4], ['d', 6]])]
        data = [Series(d) for d in data]

        result = DataFrame(data)
        sdict = OrderedDict(zip(range(len(data)), data))
        expected = DataFrame.from_dict(sdict, orient='index')
        assert_frame_equal(result, expected.reindex(result.index))

        result2 = DataFrame(data, index=np.arange(6))
        assert_frame_equal(result, result2)

        result = DataFrame([Series({})])
        expected = DataFrame(index=[0])
        assert_frame_equal(result, expected)

        data = [OrderedDict([['a', 1.5], ['b', 3.0], ['c', 4.0]]),
                OrderedDict([['a', 1.5], ['b', 3.0], ['c', 6.0]])]
        sdict = OrderedDict(zip(range(len(data)), data))

        idx = Index(['a', 'b', 'c'])
        data2 = [Series([1.5, 3, 4], idx, dtype='O'),
                 Series([1.5, 3, 6], idx)]
        result = DataFrame(data2)
        expected = DataFrame.from_dict(sdict, orient='index')
        assert_frame_equal(result, expected)

    def test_constructor_list_of_derived_dicts(self):
        class CustomDict(dict):
            pass
        d = {'a': 1.5, 'b': 3}

        data_custom = [CustomDict(d)]
        data = [d]

        result_custom = DataFrame(data_custom)
        result = DataFrame(data)
        assert_frame_equal(result, result_custom)

    def test_constructor_ragged(self):
        data = {'A': randn(10),
                'B': randn(8)}
        with assertRaisesRegexp(ValueError, 'arrays must all be same length'):
            DataFrame(data)

    def test_constructor_scalar(self):
        idx = Index(lrange(3))
        df = DataFrame({"a": 0}, index=idx)
        expected = DataFrame({"a": [0, 0, 0]}, index=idx)
        assert_frame_equal(df, expected, check_dtype=False)

    def test_constructor_Series_copy_bug(self):
        df = DataFrame(self.frame['A'], index=self.frame.index, columns=['A'])
        df.copy()

    def test_constructor_mixed_dict_and_Series(self):
        data = {}
        data['A'] = {'foo': 1, 'bar': 2, 'baz': 3}
        data['B'] = Series([4, 3, 2, 1], index=['bar', 'qux', 'baz', 'foo'])

        result = DataFrame(data)
        self.assertTrue(result.index.is_monotonic)

        # ordering ambiguous, raise exception
        with assertRaisesRegexp(ValueError, 'ambiguous ordering'):
            DataFrame({'A': ['a', 'b'], 'B': {'a': 'a', 'b': 'b'}})

        # this is OK though
        result = DataFrame({'A': ['a', 'b'],
                            'B': Series(['a', 'b'], index=['a', 'b'])})
        expected = DataFrame({'A': ['a', 'b'], 'B': ['a', 'b']},
                             index=['a', 'b'])
        assert_frame_equal(result, expected)

    def test_constructor_tuples(self):
        result = DataFrame({'A': [(1, 2), (3, 4)]})
        expected = DataFrame({'A': Series([(1, 2), (3, 4)])})
        assert_frame_equal(result, expected)

    def test_constructor_orient(self):
        data_dict = self.mixed_frame.T._series
        recons = DataFrame.from_dict(data_dict, orient='index')
        expected = self.mixed_frame.sort_index()
        assert_frame_equal(recons, expected)

        # dict of sequence
        a = {'hi': [32, 3, 3],
             'there': [3, 5, 3]}
        rs = DataFrame.from_dict(a, orient='index')
        xp = DataFrame.from_dict(a).T.reindex(list(a.keys()))
        assert_frame_equal(rs, xp)

    def test_constructor_Series_named(self):
        a = Series([1, 2, 3], index=['a', 'b', 'c'], name='x')
        df = DataFrame(a)
        self.assertEqual(df.columns[0], 'x')
        self.assertTrue(df.index.equals(a.index))

        # ndarray like
        arr = np.random.randn(10)
        s = Series(arr,name='x')
        df = DataFrame(s)
        expected = DataFrame(dict(x = s))
        assert_frame_equal(df,expected)

        s = Series(arr,index=range(3,13))
        df = DataFrame(s)
        expected = DataFrame({ 0 : s })
        assert_frame_equal(df,expected)

        self.assertRaises(ValueError, DataFrame, s, columns=[1,2])

        # #2234
        a = Series([], name='x')
        df = DataFrame(a)
        self.assertEqual(df.columns[0], 'x')

        # series with name and w/o
        s1 = Series(arr,name='x')
        df = DataFrame([s1, arr]).T
        expected = DataFrame({ 'x' : s1, 'Unnamed 0' : arr },columns=['x','Unnamed 0'])
        assert_frame_equal(df,expected)

        # this is a bit non-intuitive here; the series collapse down to arrays
        df = DataFrame([arr, s1]).T
        expected = DataFrame({ 1 : s1, 0 : arr },columns=[0,1])
        assert_frame_equal(df,expected)

    def test_constructor_Series_differently_indexed(self):
        # name
        s1 = Series([1, 2, 3], index=['a', 'b', 'c'], name='x')

        # no name
        s2 = Series([1, 2, 3], index=['a', 'b', 'c'])

        other_index = Index(['a', 'b'])

        df1 = DataFrame(s1, index=other_index)
        exp1 = DataFrame(s1.reindex(other_index))
        self.assertEqual(df1.columns[0], 'x')
        assert_frame_equal(df1, exp1)

        df2 = DataFrame(s2, index=other_index)
        exp2 = DataFrame(s2.reindex(other_index))
        self.assertEqual(df2.columns[0], 0)
        self.assertTrue(df2.index.equals(other_index))
        assert_frame_equal(df2, exp2)

    def test_constructor_manager_resize(self):
        index = list(self.frame.index[:5])
        columns = list(self.frame.columns[:3])

        result = DataFrame(self.frame._data, index=index,
                           columns=columns)
        self.assert_numpy_array_equal(result.index, index)
        self.assert_numpy_array_equal(result.columns, columns)

    def test_constructor_from_items(self):
        items = [(c, self.frame[c]) for c in self.frame.columns]
        recons = DataFrame.from_items(items)
        assert_frame_equal(recons, self.frame)

        # pass some columns
        recons = DataFrame.from_items(items, columns=['C', 'B', 'A'])
        assert_frame_equal(recons, self.frame.ix[:, ['C', 'B', 'A']])

        # orient='index'

        row_items = [(idx, self.mixed_frame.xs(idx))
                     for idx in self.mixed_frame.index]

        recons = DataFrame.from_items(row_items,
                                      columns=self.mixed_frame.columns,
                                      orient='index')
        assert_frame_equal(recons, self.mixed_frame)
        self.assertEqual(recons['A'].dtype, np.float64)

        with tm.assertRaisesRegexp(TypeError,
                                   "Must pass columns with orient='index'"):
            DataFrame.from_items(row_items, orient='index')

        # orient='index', but thar be tuples
        arr = lib.list_to_object_array(
            [('bar', 'baz')] * len(self.mixed_frame))
        self.mixed_frame['foo'] = arr
        row_items = [(idx, list(self.mixed_frame.xs(idx)))
                     for idx in self.mixed_frame.index]
        recons = DataFrame.from_items(row_items,
                                      columns=self.mixed_frame.columns,
                                      orient='index')
        assert_frame_equal(recons, self.mixed_frame)
        tm.assertIsInstance(recons['foo'][0], tuple)

        rs = DataFrame.from_items([('A', [1, 2, 3]), ('B', [4, 5, 6])],
                                  orient='index', columns=['one', 'two', 'three'])
        xp = DataFrame([[1, 2, 3], [4, 5, 6]], index=['A', 'B'],
                       columns=['one', 'two', 'three'])
        assert_frame_equal(rs, xp)

    def test_constructor_mix_series_nonseries(self):
        df = DataFrame({'A': self.frame['A'],
                        'B': list(self.frame['B'])}, columns=['A', 'B'])
        assert_frame_equal(df, self.frame.ix[:, ['A', 'B']])

        with tm.assertRaisesRegexp(ValueError, 'does not match index length'):
            DataFrame({'A': self.frame['A'], 'B': list(self.frame['B'])[:-2]})

    def test_constructor_miscast_na_int_dtype(self):
        df = DataFrame([[np.nan, 1], [1, 0]], dtype=np.int64)
        expected = DataFrame([[np.nan, 1], [1, 0]])
        assert_frame_equal(df, expected)

    def test_constructor_iterator_failure(self):
        with assertRaisesRegexp(TypeError, 'iterator'):
            df = DataFrame(iter([1, 2, 3]))

    def test_constructor_column_duplicates(self):
        # it works! #2079
        df = DataFrame([[8, 5]], columns=['a', 'a'])
        edf = DataFrame([[8, 5]])
        edf.columns = ['a', 'a']

        assert_frame_equal(df, edf)

        idf = DataFrame.from_items(
            [('a', [8]), ('a', [5])], columns=['a', 'a'])
        assert_frame_equal(idf, edf)

        self.assertRaises(ValueError, DataFrame.from_items,
                          [('a', [8]), ('a', [5]), ('b', [6])],
                          columns=['b', 'a', 'a'])

    def test_constructor_empty_with_string_dtype(self):
        # GH 9428
        expected = DataFrame(index=[0, 1], columns=[0, 1], dtype=object)

        df = DataFrame(index=[0, 1], columns=[0, 1], dtype=str)
        assert_frame_equal(df, expected)
        df = DataFrame(index=[0, 1], columns=[0, 1], dtype=np.str_)
        assert_frame_equal(df, expected)
        df = DataFrame(index=[0, 1], columns=[0, 1], dtype=np.unicode_)
        assert_frame_equal(df, expected)
        df = DataFrame(index=[0, 1], columns=[0, 1], dtype='U5')
        assert_frame_equal(df, expected)


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

    def test_constructor_single_value(self):

        # expecting single value upcasting here
        df = DataFrame(0., index=[1, 2, 3], columns=['a', 'b', 'c'])
        assert_frame_equal(df, DataFrame(np.zeros(df.shape).astype('float64'), df.index,
                                         df.columns))

        df = DataFrame(0, index=[1, 2, 3], columns=['a', 'b', 'c'])
        assert_frame_equal(df, DataFrame(np.zeros(df.shape).astype('int64'), df.index,
                                         df.columns))


        df = DataFrame('a', index=[1, 2], columns=['a', 'c'])
        assert_frame_equal(df, DataFrame(np.array([['a', 'a'],
                                                   ['a', 'a']],
                                                  dtype=object),
                                         index=[1, 2],
                                         columns=['a', 'c']))

        self.assertRaises(com.PandasError, DataFrame, 'a', [1, 2])
        self.assertRaises(com.PandasError, DataFrame, 'a', columns=['a', 'c'])
        with tm.assertRaisesRegexp(TypeError, 'incompatible data and dtype'):
            DataFrame('a', [1, 2], ['a', 'c'], float)

    def test_constructor_with_datetimes(self):
        intname = np.dtype(np.int_).name
        floatname = np.dtype(np.float_).name
        datetime64name = np.dtype('M8[ns]').name
        objectname = np.dtype(np.object_).name

        # single item
        df = DataFrame({'A' : 1, 'B' : 'foo', 'C' : 'bar', 'D' : Timestamp("20010101"), 'E' : datetime(2001,1,2,0,0) },
                       index=np.arange(10))
        result = df.get_dtype_counts()
        expected = Series({'int64': 1, datetime64name: 2, objectname : 2})
        result.sort_index()
        expected.sort_index()
        assert_series_equal(result, expected)

        # check with ndarray construction ndim==0 (e.g. we are passing a ndim 0 ndarray with a dtype specified)
        df = DataFrame({'a': 1., 'b': 2, 'c': 'foo', floatname : np.array(1.,dtype=floatname),
                        intname : np.array(1,dtype=intname)}, index=np.arange(10))
        result = df.get_dtype_counts()
        expected = { objectname : 1 }
        if intname == 'int64':
            expected['int64'] = 2
        else:
            expected['int64'] = 1
            expected[intname] = 1
        if floatname == 'float64':
            expected['float64'] = 2
        else:
            expected['float64'] = 1
            expected[floatname] = 1

        result.sort_index()
        expected = Series(expected)
        expected.sort_index()
        assert_series_equal(result, expected)

        # check with ndarray construction ndim>0
        df = DataFrame({'a': 1., 'b': 2, 'c': 'foo', floatname : np.array([1.]*10,dtype=floatname),
                        intname : np.array([1]*10,dtype=intname)}, index=np.arange(10))
        result = df.get_dtype_counts()
        result.sort_index()
        assert_series_equal(result, expected)

        # GH 2809
        ind = date_range(start="2000-01-01", freq="D", periods=10)
        datetimes = [ts.to_pydatetime() for ts in ind]
        datetime_s = Series(datetimes)
        self.assertEqual(datetime_s.dtype, 'M8[ns]')
        df = DataFrame({'datetime_s':datetime_s})
        result = df.get_dtype_counts()
        expected = Series({ datetime64name : 1 })
        result.sort_index()
        expected.sort_index()
        assert_series_equal(result, expected)

        # GH 2810
        ind = date_range(start="2000-01-01", freq="D", periods=10)
        datetimes = [ts.to_pydatetime() for ts in ind]
        dates = [ts.date() for ts in ind]
        df = DataFrame({'datetimes': datetimes, 'dates':dates})
        result = df.get_dtype_counts()
        expected = Series({ datetime64name : 1, objectname : 1 })
        result.sort_index()
        expected.sort_index()
        assert_series_equal(result, expected)

        # GH 7594
        # don't coerce tz-aware
        import pytz
        tz = pytz.timezone('US/Eastern')
        dt = tz.localize(datetime(2012, 1, 1))

        df = DataFrame({'End Date': dt}, index=[0])
        self.assertEqual(df.iat[0,0],dt)
        assert_series_equal(df.dtypes,Series({'End Date' : 'datetime64[ns, US/Eastern]' }))

        df = DataFrame([{'End Date': dt}])
        self.assertEqual(df.iat[0,0],dt)
        assert_series_equal(df.dtypes,Series({'End Date' : 'datetime64[ns, US/Eastern]' }))

        # tz-aware (UTC and other tz's)
        # GH 8411
        dr = date_range('20130101',periods=3)
        df = DataFrame({ 'value' : dr})
        self.assertTrue(df.iat[0,0].tz is None)
        dr = date_range('20130101',periods=3,tz='UTC')
        df = DataFrame({ 'value' : dr})
        self.assertTrue(str(df.iat[0,0].tz) == 'UTC')
        dr = date_range('20130101',periods=3,tz='US/Eastern')
        df = DataFrame({ 'value' : dr})
        self.assertTrue(str(df.iat[0,0].tz) == 'US/Eastern')

        # GH 7822
        # preserver an index with a tz on dict construction
        i = date_range('1/1/2011', periods=5, freq='10s', tz = 'US/Eastern')

        expected = DataFrame( {'a' : i.to_series(keep_tz=True).reset_index(drop=True) })
        df = DataFrame()
        df['a'] = i
        assert_frame_equal(df, expected)

        df = DataFrame( {'a' : i } )
        assert_frame_equal(df, expected)

        # multiples
        i_no_tz = date_range('1/1/2011', periods=5, freq='10s')
        df = DataFrame( {'a' : i, 'b' :  i_no_tz } )
        expected = DataFrame( {'a' : i.to_series(keep_tz=True).reset_index(drop=True), 'b': i_no_tz })
        assert_frame_equal(df, expected)

    def test_constructor_with_datetime_tz(self):

        # 8260
        # support datetime64 with tz

        idx = Index(date_range('20130101',periods=3,tz='US/Eastern'),
                    name='foo')
        dr = date_range('20130110',periods=3)

        # construction
        df = DataFrame({'A' : idx, 'B' : dr})
        self.assertTrue(df['A'].dtype,'M8[ns, US/Eastern')
        self.assertTrue(df['A'].name == 'A')
        assert_series_equal(df['A'],Series(idx,name='A'))
        assert_series_equal(df['B'],Series(dr,name='B'))

        # construction from dict
        df2 = DataFrame(dict(A=Timestamp('20130102', tz='US/Eastern'), B=Timestamp('20130603', tz='CET')), index=range(5))
        assert_series_equal(df2.dtypes, Series(['datetime64[ns, US/Eastern]', 'datetime64[ns, CET]'], index=['A','B']))

        # dtypes
        tzframe = DataFrame({'A' : date_range('20130101',periods=3),
                             'B' : date_range('20130101',periods=3,tz='US/Eastern'),
                             'C' : date_range('20130101',periods=3,tz='CET')})
        tzframe.iloc[1,1] = pd.NaT
        tzframe.iloc[1,2] = pd.NaT
        result = tzframe.dtypes.sort_index()
        expected = Series([ np.dtype('datetime64[ns]'),
                            DatetimeTZDtype('datetime64[ns, US/Eastern]'),
                            DatetimeTZDtype('datetime64[ns, CET]') ],
                          ['A','B','C'])

        # concat
        df3 = pd.concat([df2.A.to_frame(),df2.B.to_frame()],axis=1)
        assert_frame_equal(df2, df3)

        # select_dtypes
        result = df3.select_dtypes(include=['datetime64[ns]'])
        expected = df3.reindex(columns=[])
        assert_frame_equal(result, expected)

        # this will select based on issubclass, and these are the same class
        result = df3.select_dtypes(include=['datetime64[ns, CET]'])
        expected = df3
        assert_frame_equal(result, expected)

        # from index
        idx2 = date_range('20130101',periods=3,tz='US/Eastern',name='foo')
        df2 = DataFrame(idx2)
        assert_series_equal(df2['foo'],Series(idx2,name='foo'))
        df2 = DataFrame(Series(idx2))
        assert_series_equal(df2['foo'],Series(idx2,name='foo'))

        idx2 = date_range('20130101',periods=3,tz='US/Eastern')
        df2 = DataFrame(idx2)
        assert_series_equal(df2[0],Series(idx2,name=0))
        df2 = DataFrame(Series(idx2))
        assert_series_equal(df2[0],Series(idx2,name=0))

        # interleave with object
        result = self.tzframe.assign(D = 'foo').values
        expected = np.array([[Timestamp('2013-01-01 00:00:00'),
                              Timestamp('2013-01-02 00:00:00'),
                              Timestamp('2013-01-03 00:00:00')],
                             [Timestamp('2013-01-01 00:00:00-0500', tz='US/Eastern'),
                              pd.NaT,
                              Timestamp('2013-01-03 00:00:00-0500', tz='US/Eastern')],
                             [Timestamp('2013-01-01 00:00:00+0100', tz='CET'),
                              pd.NaT,
                              Timestamp('2013-01-03 00:00:00+0100', tz='CET')],
                             ['foo','foo','foo']], dtype=object).T
        self.assert_numpy_array_equal(result, expected)

        # interleave with only datetime64[ns]
        result = self.tzframe.values
        expected = np.array([[Timestamp('2013-01-01 00:00:00'),
                              Timestamp('2013-01-02 00:00:00'),
                              Timestamp('2013-01-03 00:00:00')],
                             [Timestamp('2013-01-01 00:00:00-0500', tz='US/Eastern'),
                              pd.NaT,
                              Timestamp('2013-01-03 00:00:00-0500', tz='US/Eastern')],
                             [Timestamp('2013-01-01 00:00:00+0100', tz='CET'),
                              pd.NaT,
                              Timestamp('2013-01-03 00:00:00+0100', tz='CET')]], dtype=object).T
        self.assert_numpy_array_equal(result, expected)

        # astype
        expected = np.array([[Timestamp('2013-01-01 00:00:00'),
                              Timestamp('2013-01-02 00:00:00'),
                              Timestamp('2013-01-03 00:00:00')],
                             [Timestamp('2013-01-01 00:00:00-0500', tz='US/Eastern'),
                              pd.NaT,
                              Timestamp('2013-01-03 00:00:00-0500', tz='US/Eastern')],
                             [Timestamp('2013-01-01 00:00:00+0100', tz='CET'),
                              pd.NaT,
                              Timestamp('2013-01-03 00:00:00+0100', tz='CET')]], dtype=object).T
        result = self.tzframe.astype(object)
        assert_frame_equal(result, DataFrame(expected, index=self.tzframe.index, columns=self.tzframe.columns))

        result = self.tzframe.astype('datetime64[ns]')
        expected = DataFrame({'A' : date_range('20130101',periods=3),
                              'B' : date_range('20130101',periods=3,tz='US/Eastern').tz_convert('UTC').tz_localize(None),
                              'C' : date_range('20130101',periods=3,tz='CET').tz_convert('UTC').tz_localize(None)})
        expected.iloc[1,1] = pd.NaT
        expected.iloc[1,2] = pd.NaT
        assert_frame_equal(result, expected)

        # str formatting
        result = self.tzframe.astype(str)
        expected = np.array([['2013-01-01', '2013-01-01 00:00:00-05:00',
                              '2013-01-01 00:00:00+01:00'],
                             ['2013-01-02', 'NaT', 'NaT'],
                             ['2013-01-03', '2013-01-03 00:00:00-05:00',
                              '2013-01-03 00:00:00+01:00']], dtype=object)
        self.assert_numpy_array_equal(result, expected)

        result = str(self.tzframe)
        self.assertTrue('0 2013-01-01 2013-01-01 00:00:00-05:00 2013-01-01 00:00:00+01:00' in result)
        self.assertTrue('1 2013-01-02                       NaT                       NaT' in result)
        self.assertTrue('2 2013-01-03 2013-01-03 00:00:00-05:00 2013-01-03 00:00:00+01:00' in result)

        # setitem
        df['C'] = idx
        assert_series_equal(df['C'],Series(idx,name='C'))

        df['D'] = 'foo'
        df['D'] = idx
        assert_series_equal(df['D'],Series(idx,name='D'))
        del df['D']

        # assert that A & C are not sharing the same base (e.g. they
        # are copies)
        b1 = df._data.blocks[1]
        b2 = df._data.blocks[2]
        self.assertTrue(b1.values.equals(b2.values))
        self.assertFalse(id(b1.values.values.base) == id(b2.values.values.base))

        # with nan
        df2 = df.copy()
        df2.iloc[1,1] = pd.NaT
        df2.iloc[1,2] = pd.NaT
        result = df2['B']
        assert_series_equal(notnull(result), Series([True,False,True],name='B'))
        assert_series_equal(df2.dtypes, df.dtypes)

        # set/reset
        df = DataFrame({'A' : [0,1,2] }, index=idx)
        result = df.reset_index()
        self.assertTrue(result['foo'].dtype,'M8[ns, US/Eastern')

        result = result.set_index('foo')
        tm.assert_index_equal(df.index,idx)

    def test_constructor_for_list_with_dtypes(self):
        intname = np.dtype(np.int_).name
        floatname = np.dtype(np.float_).name
        datetime64name = np.dtype('M8[ns]').name
        objectname = np.dtype(np.object_).name

        # test list of lists/ndarrays
        df = DataFrame([np.arange(5) for x in range(5)])
        result = df.get_dtype_counts()
        expected = Series({'int64' : 5})

        df = DataFrame([np.array(np.arange(5),dtype='int32') for x in range(5)])
        result = df.get_dtype_counts()
        expected = Series({'int32' : 5})

        # overflow issue? (we always expecte int64 upcasting here)
        df = DataFrame({'a' : [2**31,2**31+1]})
        result = df.get_dtype_counts()
        expected = Series({'int64' : 1 })
        assert_series_equal(result, expected)

        # GH #2751 (construction with no index specified), make sure we cast to platform values
        df = DataFrame([1, 2])
        result = df.get_dtype_counts()
        expected = Series({'int64': 1 })
        assert_series_equal(result, expected)

        df = DataFrame([1.,2.])
        result = df.get_dtype_counts()
        expected = Series({'float64' : 1 })
        assert_series_equal(result, expected)

        df = DataFrame({'a' : [1, 2]})
        result = df.get_dtype_counts()
        expected = Series({'int64' : 1})
        assert_series_equal(result, expected)

        df = DataFrame({'a' : [1., 2.]})
        result = df.get_dtype_counts()
        expected = Series({'float64' : 1})
        assert_series_equal(result, expected)

        df = DataFrame({'a' : 1 }, index=lrange(3))
        result = df.get_dtype_counts()
        expected = Series({'int64': 1})
        assert_series_equal(result, expected)

        df = DataFrame({'a' : 1. }, index=lrange(3))
        result = df.get_dtype_counts()
        expected = Series({'float64': 1 })
        assert_series_equal(result, expected)

        # with object list
        df = DataFrame({'a':[1,2,4,7], 'b':[1.2, 2.3, 5.1, 6.3],
                        'c':list('abcd'), 'd':[datetime(2000,1,1) for i in range(4)],
                        'e' : [1.,2,4.,7]})
        result = df.get_dtype_counts()
        expected = Series({'int64': 1, 'float64' : 2, datetime64name: 1, objectname : 1})
        result.sort_index()
        expected.sort_index()
        assert_series_equal(result, expected)

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

        from datetime import datetime, timedelta
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

    def test_from_records_to_records(self):
        # from numpy documentation
        arr = np.zeros((2,), dtype=('i4,f4,a10'))
        arr[:] = [(1, 2., 'Hello'), (2, 3., "World")]

        frame = DataFrame.from_records(arr)

        index = np.arange(len(arr))[::-1]
        indexed_frame = DataFrame.from_records(arr, index=index)
        self.assert_numpy_array_equal(indexed_frame.index, index)

        # without names, it should go to last ditch
        arr2 = np.zeros((2,3))
        tm.assert_frame_equal(DataFrame.from_records(arr2), DataFrame(arr2))

        # wrong length
        msg = r'Shape of passed values is \(3, 2\), indices imply \(3, 1\)'
        with assertRaisesRegexp(ValueError, msg):
            DataFrame.from_records(arr, index=index[:-1])

        indexed_frame = DataFrame.from_records(arr, index='f1')

        # what to do?
        records = indexed_frame.to_records()
        self.assertEqual(len(records.dtype.names), 3)

        records = indexed_frame.to_records(index=False)
        self.assertEqual(len(records.dtype.names), 2)
        self.assertNotIn('index', records.dtype.names)

    def test_from_records_nones(self):
        tuples = [(1, 2, None, 3),
                  (1, 2, None, 3),
                  (None, 2, 5, 3)]

        df = DataFrame.from_records(tuples, columns=['a', 'b', 'c', 'd'])
        self.assertTrue(np.isnan(df['c'][0]))

    def test_from_records_iterator(self):
        arr = np.array([(1.0, 1.0, 2, 2), (3.0, 3.0, 4, 4), (5., 5., 6, 6), (7., 7., 8, 8)],
                       dtype=[('x', np.float64), ('u', np.float32), ('y', np.int64), ('z', np.int32) ])
        df = DataFrame.from_records(iter(arr), nrows=2)
        xp = DataFrame({'x': np.array([1.0, 3.0], dtype=np.float64),
                        'u': np.array([1.0, 3.0], dtype=np.float32),
                        'y': np.array([2, 4], dtype=np.int64),
                        'z': np.array([2, 4], dtype=np.int32)})
        assert_frame_equal(df.reindex_like(xp), xp)

        # no dtypes specified here, so just compare with the default
        arr = [(1.0, 2), (3.0, 4), (5., 6), (7., 8)]
        df = DataFrame.from_records(iter(arr), columns=['x', 'y'],
                                    nrows=2)
        assert_frame_equal(df, xp.reindex(columns=['x','y']), check_dtype=False)

    def test_from_records_tuples_generator(self):
        def tuple_generator(length):
            for i in range(length):
                letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
                yield (i, letters[i % len(letters)], i/length)

        columns_names = ['Integer', 'String', 'Float']
        columns = [[i[j] for i in tuple_generator(10)] for j in range(len(columns_names))]
        data = {'Integer': columns[0], 'String': columns[1], 'Float': columns[2]}
        expected = DataFrame(data, columns=columns_names)

        generator = tuple_generator(10)
        result = DataFrame.from_records(generator, columns=columns_names)
        assert_frame_equal(result, expected)

    def test_from_records_lists_generator(self):
        def list_generator(length):
            for i in range(length):
                letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
                yield [i, letters[i % len(letters)], i/length]

        columns_names = ['Integer', 'String', 'Float']
        columns = [[i[j] for i in list_generator(10)] for j in range(len(columns_names))]
        data = {'Integer': columns[0], 'String': columns[1], 'Float': columns[2]}
        expected = DataFrame(data, columns=columns_names)

        generator = list_generator(10)
        result = DataFrame.from_records(generator, columns=columns_names)
        assert_frame_equal(result, expected)

    def test_from_records_columns_not_modified(self):
        tuples = [(1, 2, 3),
                  (1, 2, 3),
                  (2, 5, 3)]

        columns = ['a', 'b', 'c']
        original_columns = list(columns)
        df = DataFrame.from_records(tuples, columns=columns, index='a')
        self.assertEqual(columns, original_columns)

    def test_from_records_decimal(self):
        from decimal import Decimal

        tuples = [(Decimal('1.5'),), (Decimal('2.5'),), (None,)]

        df = DataFrame.from_records(tuples, columns=['a'])
        self.assertEqual(df['a'].dtype, object)

        df = DataFrame.from_records(tuples, columns=['a'], coerce_float=True)
        self.assertEqual(df['a'].dtype, np.float64)
        self.assertTrue(np.isnan(df['a'].values[-1]))

    def test_from_records_duplicates(self):
        result = DataFrame.from_records([(1, 2, 3), (4, 5, 6)],
                                        columns=['a', 'b', 'a'])

        expected = DataFrame([(1, 2, 3), (4, 5, 6)],
                             columns=['a', 'b', 'a'])

        assert_frame_equal(result, expected)

    def test_from_records_set_index_name(self):
        def create_dict(order_id):
            return {'order_id': order_id, 'quantity': np.random.randint(1, 10),
                    'price': np.random.randint(1, 10)}
        documents = [create_dict(i) for i in range(10)]
        # demo missing data
        documents.append({'order_id': 10, 'quantity': 5})

        result = DataFrame.from_records(documents, index='order_id')
        self.assertEqual(result.index.name, 'order_id')

        # MultiIndex
        result = DataFrame.from_records(documents,
                                        index=['order_id', 'quantity'])
        self.assertEqual(result.index.names, ('order_id', 'quantity'))

    def test_from_records_misc_brokenness(self):
        # #2179

        data = {1: ['foo'], 2: ['bar']}

        result = DataFrame.from_records(data, columns=['a', 'b'])
        exp = DataFrame(data, columns=['a', 'b'])
        assert_frame_equal(result, exp)

        # overlap in index/index_names

        data = {'a': [1, 2, 3], 'b': [4, 5, 6]}

        result = DataFrame.from_records(data, index=['a', 'b', 'c'])
        exp = DataFrame(data, index=['a', 'b', 'c'])
        assert_frame_equal(result, exp)


        # GH 2623
        rows = []
        rows.append([datetime(2010, 1, 1), 1])
        rows.append([datetime(2010, 1, 2), 'hi']) # test col upconverts to obj
        df2_obj = DataFrame.from_records(rows, columns=['date', 'test'])
        results = df2_obj.get_dtype_counts()
        expected = Series({ 'datetime64[ns]' : 1, 'object' : 1 })

        rows = []
        rows.append([datetime(2010, 1, 1), 1])
        rows.append([datetime(2010, 1, 2), 1])
        df2_obj = DataFrame.from_records(rows, columns=['date', 'test'])
        results = df2_obj.get_dtype_counts()
        expected = Series({ 'datetime64[ns]' : 1, 'int64' : 1 })

    def test_from_records_empty(self):
        # 3562
        result = DataFrame.from_records([], columns=['a','b','c'])
        expected = DataFrame(columns=['a','b','c'])
        assert_frame_equal(result, expected)

        result = DataFrame.from_records([], columns=['a','b','b'])
        expected = DataFrame(columns=['a','b','b'])
        assert_frame_equal(result, expected)

    def test_from_records_empty_with_nonempty_fields_gh3682(self):
        a = np.array([(1, 2)], dtype=[('id', np.int64), ('value', np.int64)])
        df = DataFrame.from_records(a, index='id')
        assert_numpy_array_equal(df.index, Index([1], name='id'))
        self.assertEqual(df.index.name, 'id')
        assert_numpy_array_equal(df.columns, Index(['value']))

        b = np.array([], dtype=[('id', np.int64), ('value', np.int64)])
        df = DataFrame.from_records(b, index='id')
        assert_numpy_array_equal(df.index, Index([], name='id'))
        self.assertEqual(df.index.name, 'id')

    def test_from_records_with_datetimes(self):
        if sys.version < LooseVersion('2.7'):
            raise nose.SkipTest('rec arrays dont work properly with py2.6')

        # this may fail on certain platforms because of a numpy issue
        # related GH6140
        if not is_little_endian():
            raise nose.SkipTest("known failure of test on non-little endian")

        # construction with a null in a recarray
        # GH 6140
        expected = DataFrame({ 'EXPIRY'  : [datetime(2005, 3, 1, 0, 0), None ]})

        arrdata = [np.array([datetime(2005, 3, 1, 0, 0), None])]
        dtypes = [('EXPIRY', '<M8[ns]')]

        try:
            recarray = np.core.records.fromarrays(arrdata, dtype=dtypes)
        except (ValueError):
            raise nose.SkipTest("known failure of numpy rec array creation")

        result = DataFrame.from_records(recarray)
        assert_frame_equal(result,expected)

        # coercion should work too
        arrdata = [np.array([datetime(2005, 3, 1, 0, 0), None])]
        dtypes = [('EXPIRY', '<M8[m]')]
        recarray = np.core.records.fromarrays(arrdata, dtype=dtypes)
        result = DataFrame.from_records(recarray)
        assert_frame_equal(result,expected)

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
        tm.assert_frame_equal(df1.join(df2, how='left'), exp)
        tm.assert_frame_equal(df2.join(df1, how='right'), exp[['value2', 'value1']])

        exp_idx = pd.MultiIndex.from_product([['a', 'b'], ['x', 'y', 'z']],
                                             names=['first', 'second'])
        exp = pd.DataFrame([[0.471780, 10], [0.774908, 10], [0.563634, 10],
                            [-0.353756, 20], [0.368062, 20], [-1.721840, 20]],
                           index=exp_idx, columns=['value1', 'value2'])

        tm.assert_frame_equal(df1.join(df2, how='right'), exp)
        tm.assert_frame_equal(df2.join(df1, how='left'), exp[['value2', 'value1']])

    def test_from_records_sequencelike(self):
        df = DataFrame({'A' : np.array(np.random.randn(6), dtype = np.float64),
                        'A1': np.array(np.random.randn(6), dtype = np.float64),
                        'B' : np.array(np.arange(6), dtype = np.int64),
                        'C' : ['foo'] * 6,
                        'D' : np.array([True, False] * 3, dtype=bool),
                        'E' : np.array(np.random.randn(6), dtype = np.float32),
                        'E1': np.array(np.random.randn(6), dtype = np.float32),
                        'F' : np.array(np.arange(6), dtype = np.int32) })

        # this is actually tricky to create the recordlike arrays and have the dtypes be intact
        blocks = df.blocks
        tuples = []
        columns = []
        dtypes  = []
        for dtype, b in compat.iteritems(blocks):
            columns.extend(b.columns)
            dtypes.extend([ (c,np.dtype(dtype).descr[0][1]) for c in b.columns ])
        for i in range(len(df.index)):
            tup = []
            for _, b in compat.iteritems(blocks):
                tup.extend(b.iloc[i].values)
            tuples.append(tuple(tup))

        recarray  = np.array(tuples, dtype=dtypes).view(np.recarray)
        recarray2 = df.to_records()
        lists     = [list(x) for x in tuples]

        # tuples (lose the dtype info)
        result  = DataFrame.from_records(tuples,    columns=columns).reindex(columns=df.columns)

        # created recarray and with to_records recarray (have dtype info)
        result2 = DataFrame.from_records(recarray,  columns=columns).reindex(columns=df.columns)
        result3 = DataFrame.from_records(recarray2, columns=columns).reindex(columns=df.columns)

        # list of tupels (no dtype info)
        result4 = DataFrame.from_records(lists,     columns=columns).reindex(columns=df.columns)

        assert_frame_equal(result, df, check_dtype=False)
        assert_frame_equal(result2, df)
        assert_frame_equal(result3, df)
        assert_frame_equal(result4, df, check_dtype=False)

        # tuples is in the order of the columns
        result = DataFrame.from_records(tuples)
        self.assert_numpy_array_equal(result.columns, lrange(8))

        # test exclude parameter & we are casting the results here (as we don't have dtype info to recover)
        columns_to_test = [ columns.index('C'), columns.index('E1') ]

        exclude = list(set(range(8))-set(columns_to_test))
        result = DataFrame.from_records(tuples, exclude=exclude)
        result.columns = [ columns[i] for i in sorted(columns_to_test) ]
        assert_series_equal(result['C'], df['C'])
        assert_series_equal(result['E1'], df['E1'].astype('float64'))

        # empty case
        result = DataFrame.from_records([], columns=['foo', 'bar', 'baz'])
        self.assertEqual(len(result), 0)
        self.assert_numpy_array_equal(result.columns, ['foo', 'bar', 'baz'])

        result = DataFrame.from_records([])
        self.assertEqual(len(result), 0)
        self.assertEqual(len(result.columns), 0)

    def test_from_records_dictlike(self):

        # test the dict methods
        df = DataFrame({'A' : np.array(np.random.randn(6), dtype = np.float64),
                        'A1': np.array(np.random.randn(6), dtype = np.float64),
                        'B' : np.array(np.arange(6), dtype = np.int64),
                        'C' : ['foo'] * 6,
                        'D' : np.array([True, False] * 3, dtype=bool),
                        'E' : np.array(np.random.randn(6), dtype = np.float32),
                        'E1': np.array(np.random.randn(6), dtype = np.float32),
                        'F' : np.array(np.arange(6), dtype = np.int32) })

        # columns is in a different order here than the actual items iterated from the dict
        columns = []
        for dtype, b in compat.iteritems(df.blocks):
            columns.extend(b.columns)

        asdict    = dict((x, y) for x, y in compat.iteritems(df))
        asdict2   = dict((x, y.values) for x, y in compat.iteritems(df))

        # dict of series & dict of ndarrays (have dtype info)
        results = []
        results.append(DataFrame.from_records(asdict).reindex(columns=df.columns))
        results.append(DataFrame.from_records(asdict,    columns=columns).reindex(columns=df.columns))
        results.append(DataFrame.from_records(asdict2,   columns=columns).reindex(columns=df.columns))

        for r in results:
            assert_frame_equal(r, df)

    def test_from_records_with_index_data(self):
        df = DataFrame(np.random.randn(10, 3), columns=['A', 'B', 'C'])

        data = np.random.randn(10)
        df1 = DataFrame.from_records(df, index=data)
        assert(df1.index.equals(Index(data)))

    def test_from_records_bad_index_column(self):
        df = DataFrame(np.random.randn(10, 3), columns=['A', 'B', 'C'])

        # should pass
        df1 = DataFrame.from_records(df, index=['C'])
        assert(df1.index.equals(Index(df.C)))

        df1 = DataFrame.from_records(df, index='C')
        assert(df1.index.equals(Index(df.C)))

        # should fail
        self.assertRaises(ValueError, DataFrame.from_records, df, index=[2])
        self.assertRaises(KeyError, DataFrame.from_records, df, index=2)

    def test_from_records_non_tuple(self):
        class Record(object):

            def __init__(self, *args):
                self.args = args

            def __getitem__(self, i):
                return self.args[i]

            def __iter__(self):
                return iter(self.args)

        recs = [Record(1, 2, 3), Record(4, 5, 6), Record(7, 8, 9)]
        tups = lmap(tuple, recs)

        result = DataFrame.from_records(recs)
        expected = DataFrame.from_records(tups)
        assert_frame_equal(result, expected)

    def test_from_records_len0_with_columns(self):
        # #2633
        result = DataFrame.from_records([], index='foo',
                                        columns=['foo', 'bar'])

        self.assertTrue(np.array_equal(result.columns, ['bar']))
        self.assertEqual(len(result), 0)
        self.assertEqual(result.index.name, 'foo')

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
        assert_frame_equal(self.frame.head(0), self.frame)
        assert_frame_equal(self.frame.tail(0), self.frame)
        assert_frame_equal(self.frame.head(-1), self.frame[:-1])
        assert_frame_equal(self.frame.tail(-1), self.frame[1:])
        assert_frame_equal(self.frame.head(1), self.frame[:1])
        assert_frame_equal(self.frame.tail(1), self.frame[-1:])
        # with a float index
        df = self.frame.copy()
        df.index = np.arange(len(self.frame)) + 0.1
        assert_frame_equal(df.head(), df.iloc[:5])
        assert_frame_equal(df.tail(), df.iloc[-5:])
        assert_frame_equal(df.head(0), df)
        assert_frame_equal(df.tail(0), df)
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

    def test_len(self):
        self.assertEqual(len(self.frame), len(self.frame.index))

    def test_operators(self):
        garbage = random.random(4)
        colSeries = Series(garbage, index=np.array(self.frame.columns))

        idSum = self.frame + self.frame
        seriesSum = self.frame + colSeries

        for col, series in compat.iteritems(idSum):
            for idx, val in compat.iteritems(series):
                origVal = self.frame[col][idx] * 2
                if not np.isnan(val):
                    self.assertEqual(val, origVal)
                else:
                    self.assertTrue(np.isnan(origVal))

        for col, series in compat.iteritems(seriesSum):
            for idx, val in compat.iteritems(series):
                origVal = self.frame[col][idx] + colSeries[col]
                if not np.isnan(val):
                    self.assertEqual(val, origVal)
                else:
                    self.assertTrue(np.isnan(origVal))

        added = self.frame2 + self.frame2
        expected = self.frame2 * 2
        assert_frame_equal(added, expected)

        df = DataFrame({'a': ['a', None, 'b']})
        assert_frame_equal(df + df, DataFrame({'a': ['aa', np.nan, 'bb']}))

        # Test for issue #10181
        for dtype in ('float', 'int64'):
            frames = [
                DataFrame(dtype=dtype),
                DataFrame(columns=['A'], dtype=dtype),
                DataFrame(index=[0], dtype=dtype),
            ]
            for df in frames:
                self.assertTrue((df + df).equals(df))
                assert_frame_equal(df + df, df)

    def test_ops_np_scalar(self):
        vals, xs = np.random.rand(5, 3), [nan, 7, -23, 2.718, -3.14, np.inf]
        f = lambda x: DataFrame(x, index=list('ABCDE'),
                columns=['jim', 'joe', 'jolie'])

        df = f(vals)

        for x in xs:
            assert_frame_equal(df / np.array(x), f(vals / x))
            assert_frame_equal(np.array(x) * df, f(vals * x))
            assert_frame_equal(df + np.array(x), f(vals + x))
            assert_frame_equal(np.array(x) - df, f(x - vals))

    def test_operators_boolean(self):

        # GH 5808
        # empty frames, non-mixed dtype

        result = DataFrame(index=[1]) & DataFrame(index=[1])
        assert_frame_equal(result,DataFrame(index=[1]))

        result = DataFrame(index=[1]) | DataFrame(index=[1])
        assert_frame_equal(result,DataFrame(index=[1]))

        result = DataFrame(index=[1]) & DataFrame(index=[1,2])
        assert_frame_equal(result,DataFrame(index=[1,2]))

        result = DataFrame(index=[1],columns=['A']) & DataFrame(index=[1],columns=['A'])
        assert_frame_equal(result,DataFrame(index=[1],columns=['A']))

        result = DataFrame(True,index=[1],columns=['A']) & DataFrame(True,index=[1],columns=['A'])
        assert_frame_equal(result,DataFrame(True,index=[1],columns=['A']))

        result = DataFrame(True,index=[1],columns=['A']) | DataFrame(True,index=[1],columns=['A'])
        assert_frame_equal(result,DataFrame(True,index=[1],columns=['A']))

        # boolean ops
        result = DataFrame(1,index=[1],columns=['A']) | DataFrame(True,index=[1],columns=['A'])
        assert_frame_equal(result,DataFrame(1,index=[1],columns=['A']))

        def f():
            DataFrame(1.0,index=[1],columns=['A']) | DataFrame(True,index=[1],columns=['A'])
        self.assertRaises(TypeError, f)

        def f():
            DataFrame('foo',index=[1],columns=['A']) | DataFrame(True,index=[1],columns=['A'])
        self.assertRaises(TypeError, f)

    def test_operators_none_as_na(self):
        df = DataFrame({"col1": [2, 5.0, 123, None],
                        "col2": [1, 2, 3, 4]}, dtype=object)

        ops = [operator.add, operator.sub, operator.mul, operator.truediv]

        # since filling converts dtypes from object, changed expected to be object
        for op in ops:
            filled = df.fillna(np.nan)
            result = op(df, 3)
            expected = op(filled, 3).astype(object)
            expected[com.isnull(expected)] = None
            assert_frame_equal(result, expected)

            result = op(df, df)
            expected = op(filled, filled).astype(object)
            expected[com.isnull(expected)] = None
            assert_frame_equal(result, expected)

            result = op(df, df.fillna(7))
            assert_frame_equal(result, expected)

            result = op(df.fillna(7), df)
            assert_frame_equal(result, expected, check_dtype=False)

    def test_comparison_invalid(self):

        def check(df,df2):

            for (x, y) in [(df,df2),(df2,df)]:
                self.assertRaises(TypeError, lambda : x == y)
                self.assertRaises(TypeError, lambda : x != y)
                self.assertRaises(TypeError, lambda : x >= y)
                self.assertRaises(TypeError, lambda : x > y)
                self.assertRaises(TypeError, lambda : x < y)
                self.assertRaises(TypeError, lambda : x <= y)

        # GH4968
        # invalid date/int comparisons
        df = DataFrame(np.random.randint(10, size=(10, 1)), columns=['a'])
        df['dates'] = date_range('20010101', periods=len(df))

        df2 = df.copy()
        df2['dates'] = df['a']
        check(df,df2)

        df = DataFrame(np.random.randint(10, size=(10, 2)), columns=['a', 'b'])
        df2 = DataFrame({'a': date_range('20010101', periods=len(df)), 'b': date_range('20100101', periods=len(df))})
        check(df,df2)

    def test_timestamp_compare(self):
        # make sure we can compare Timestamps on the right AND left hand side
        # GH4982
        df = DataFrame({'dates1': date_range('20010101', periods=10),
                        'dates2': date_range('20010102', periods=10),
                        'intcol': np.random.randint(1000000000, size=10),
                        'floatcol': np.random.randn(10),
                        'stringcol': list(tm.rands(10))})
        df.loc[np.random.rand(len(df)) > 0.5, 'dates2'] = pd.NaT
        ops = {'gt': 'lt', 'lt': 'gt', 'ge': 'le', 'le': 'ge', 'eq': 'eq',
               'ne': 'ne'}
        for left, right in ops.items():
            left_f = getattr(operator, left)
            right_f = getattr(operator, right)

            # no nats
            expected = left_f(df, Timestamp('20010109'))
            result = right_f(Timestamp('20010109'), df)
            tm.assert_frame_equal(result, expected)

            # nats
            expected = left_f(df, Timestamp('nat'))
            result = right_f(Timestamp('nat'), df)
            tm.assert_frame_equal(result, expected)

    def test_modulo(self):

        # GH3590, modulo as ints
        p = DataFrame({ 'first' : [3,4,5,8], 'second' : [0,0,0,3] })

        ### this is technically wrong as the integer portion is coerced to float ###
        expected = DataFrame({ 'first' : Series([0,0,0,0],dtype='float64'), 'second' : Series([np.nan,np.nan,np.nan,0]) })
        result = p % p
        assert_frame_equal(result,expected)

        # numpy has a slightly different (wrong) treatement
        result2 = DataFrame(p.values % p.values,index=p.index,columns=p.columns,dtype='float64')
        result2.iloc[0:3,1] = np.nan
        assert_frame_equal(result2,expected)

        result = p % 0
        expected = DataFrame(np.nan,index=p.index,columns=p.columns)
        assert_frame_equal(result,expected)

        # numpy has a slightly different (wrong) treatement
        result2 = DataFrame(p.values.astype('float64') % 0,index=p.index,columns=p.columns)
        assert_frame_equal(result2,expected)

        # not commutative with series
        p = DataFrame(np.random.randn(10, 5))
        s = p[0]
        res = s % p
        res2 = p % s
        self.assertFalse(np.array_equal(res.fillna(0), res2.fillna(0)))

    def test_div(self):

        # integer div, but deal with the 0's (GH 9144)
        p = DataFrame({ 'first' : [3,4,5,8], 'second' : [0,0,0,3] })
        result = p / p

        expected = DataFrame({'first': Series([1.0, 1.0, 1.0, 1.0]),
                              'second': Series([nan, nan, nan, 1])})
        assert_frame_equal(result,expected)

        result2 = DataFrame(p.values.astype('float') / p.values, index=p.index,
                            columns=p.columns)
        assert_frame_equal(result2,expected)

        result = p / 0
        expected = DataFrame(inf, index=p.index, columns=p.columns)
        expected.iloc[0:3, 1] = nan
        assert_frame_equal(result,expected)

        # numpy has a slightly different (wrong) treatement
        result2 = DataFrame(p.values.astype('float64') / 0, index=p.index,
                            columns=p.columns)
        assert_frame_equal(result2,expected)

        p = DataFrame(np.random.randn(10, 5))
        s = p[0]
        res = s / p
        res2 = p / s
        self.assertFalse(np.array_equal(res.fillna(0), res2.fillna(0)))

    def test_logical_operators(self):

        def _check_bin_op(op):
            result = op(df1, df2)
            expected = DataFrame(op(df1.values, df2.values), index=df1.index,
                                 columns=df1.columns)
            self.assertEqual(result.values.dtype, np.bool_)
            assert_frame_equal(result, expected)

        def _check_unary_op(op):
            result = op(df1)
            expected = DataFrame(op(df1.values), index=df1.index,
                                 columns=df1.columns)
            self.assertEqual(result.values.dtype, np.bool_)
            assert_frame_equal(result, expected)

        df1 = {'a': {'a': True, 'b': False, 'c': False, 'd': True, 'e': True},
               'b': {'a': False, 'b': True, 'c': False,
                     'd': False, 'e': False},
               'c': {'a': False, 'b': False, 'c': True,
                     'd': False, 'e': False},
               'd': {'a': True, 'b': False, 'c': False, 'd': True, 'e': True},
               'e': {'a': True, 'b': False, 'c': False, 'd': True, 'e': True}}

        df2 = {'a': {'a': True, 'b': False, 'c': True, 'd': False, 'e': False},
               'b': {'a': False, 'b': True, 'c': False,
                     'd': False, 'e': False},
               'c': {'a': True, 'b': False, 'c': True, 'd': False, 'e': False},
               'd': {'a': False, 'b': False, 'c': False,
                     'd': True, 'e': False},
               'e': {'a': False, 'b': False, 'c': False,
                     'd': False, 'e': True}}

        df1 = DataFrame(df1)
        df2 = DataFrame(df2)

        _check_bin_op(operator.and_)
        _check_bin_op(operator.or_)
        _check_bin_op(operator.xor)

        # operator.neg is deprecated in numpy >= 1.9
        _check_unary_op(operator.inv)

    def test_logical_typeerror(self):
        if not compat.PY3:
            self.assertRaises(TypeError, self.frame.__eq__, 'foo')
            self.assertRaises(TypeError, self.frame.__lt__, 'foo')
            self.assertRaises(TypeError, self.frame.__gt__, 'foo')
            self.assertRaises(TypeError, self.frame.__ne__, 'foo')
        else:
            raise nose.SkipTest('test_logical_typeerror not tested on PY3')

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
                self.assertRaises(ValueError, f)


        df = DataFrame([[1,2,3],[4,5,6]], index=[1,np.nan])
        check(df)

        df = DataFrame([[1,2,3],[4,5,6]], columns=[1.1,2.2,np.nan])
        check(df)

        df = DataFrame([[0,1,2,3],[4,5,6,7]], columns=[np.nan,1.1,2.2,np.nan])
        check(df)

        df = DataFrame([[0.0,1,2,3.0],[4,5,6,7]], columns=[np.nan,1.1,2.2,np.nan])
        check(df)

    def test_logical_with_nas(self):
        d = DataFrame({'a': [np.nan, False], 'b': [True, True]})

        # GH4947
        # bool comparisons should return bool
        result = d['a'] | d['b']
        expected = Series([False, True])
        assert_series_equal(result, expected)

        # GH4604, automatic casting here
        result = d['a'].fillna(False) | d['b']
        expected = Series([True, True])
        assert_series_equal(result, expected)

        result = d['a'].fillna(False,downcast=False) | d['b']
        expected = Series([True, True])
        assert_series_equal(result, expected)

    def test_neg(self):
        # what to do?
        assert_frame_equal(-self.frame, -1 * self.frame)

    def test_invert(self):
        assert_frame_equal(-(self.frame < 0), ~(self.frame < 0))

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

    def test_arith_flex_frame(self):
        ops = ['add', 'sub', 'mul', 'div', 'truediv', 'pow', 'floordiv', 'mod']
        if not compat.PY3:
            aliases = {}
        else:
            aliases = {'div': 'truediv'}

        for op in ops:
            try:
                alias = aliases.get(op, op)
                f = getattr(operator, alias)
                result = getattr(self.frame, op)(2 * self.frame)
                exp = f(self.frame, 2 * self.frame)
                assert_frame_equal(result, exp)

                # vs mix float
                result = getattr(self.mixed_float, op)(2 * self.mixed_float)
                exp = f(self.mixed_float, 2 * self.mixed_float)
                assert_frame_equal(result, exp)
                _check_mixed_float(result, dtype = dict(C = None))

                # vs mix int
                if op in ['add','sub','mul']:
                    result = getattr(self.mixed_int, op)(2 + self.mixed_int)
                    exp = f(self.mixed_int, 2 + self.mixed_int)

                    # overflow in the uint
                    dtype = None
                    if op in ['sub']:
                        dtype = dict(B = 'object', C = None)
                    elif op in ['add','mul']:
                        dtype = dict(C = None)
                    assert_frame_equal(result, exp)
                    _check_mixed_int(result, dtype = dtype)

                    # rops
                    r_f = lambda x, y: f(y, x)
                    result = getattr(self.frame, 'r' + op)(2 * self.frame)
                    exp = r_f(self.frame, 2 * self.frame)
                    assert_frame_equal(result, exp)

                    # vs mix float
                    result = getattr(self.mixed_float, op)(2 * self.mixed_float)
                    exp = f(self.mixed_float, 2 * self.mixed_float)
                    assert_frame_equal(result, exp)
                    _check_mixed_float(result, dtype = dict(C = None))

                    result = getattr(self.intframe, op)(2 * self.intframe)
                    exp = f(self.intframe, 2 * self.intframe)
                    assert_frame_equal(result, exp)

                    # vs mix int
                    if op in ['add','sub','mul']:
                        result = getattr(self.mixed_int, op)(2 + self.mixed_int)
                        exp = f(self.mixed_int, 2 + self.mixed_int)

                        # overflow in the uint
                        dtype = None
                        if op in ['sub']:
                            dtype = dict(B = 'object', C = None)
                        elif op in ['add','mul']:
                            dtype = dict(C = None)
                        assert_frame_equal(result, exp)
                        _check_mixed_int(result, dtype = dtype)
            except:
                com.pprint_thing("Failing operation %r" % op)
                raise

            # ndim >= 3
            ndim_5 = np.ones(self.frame.shape + (3, 4, 5))
            with assertRaisesRegexp(ValueError, 'shape'):
                f(self.frame, ndim_5)

            with assertRaisesRegexp(ValueError, 'shape'):
                getattr(self.frame, op)(ndim_5)


        # res_add = self.frame.add(self.frame)
        # res_sub = self.frame.sub(self.frame)
        # res_mul = self.frame.mul(self.frame)
        # res_div = self.frame.div(2 * self.frame)

        # assert_frame_equal(res_add, self.frame + self.frame)
        # assert_frame_equal(res_sub, self.frame - self.frame)
        # assert_frame_equal(res_mul, self.frame * self.frame)
        # assert_frame_equal(res_div, self.frame / (2 * self.frame))

        const_add = self.frame.add(1)
        assert_frame_equal(const_add, self.frame + 1)

        # corner cases
        result = self.frame.add(self.frame[:0])
        assert_frame_equal(result, self.frame * np.nan)

        result = self.frame[:0].add(self.frame)
        assert_frame_equal(result, self.frame * np.nan)
        with assertRaisesRegexp(NotImplementedError, 'fill_value'):
            self.frame.add(self.frame.iloc[0], fill_value=3)
        with assertRaisesRegexp(NotImplementedError, 'fill_value'):
            self.frame.add(self.frame.iloc[0], axis='index', fill_value=3)

    def test_binary_ops_align(self):

        # test aligning binary ops

        # GH 6681
        index=MultiIndex.from_product([list('abc'),
                                       ['one','two','three'],
                                       [1,2,3]],
                                      names=['first','second','third'])

        df = DataFrame(np.arange(27*3).reshape(27,3),
                       index=index,
                       columns=['value1','value2','value3']).sortlevel()

        idx = pd.IndexSlice
        for op in ['add','sub','mul','div','truediv']:
            opa = getattr(operator,op,None)
            if opa is None:
                continue

            x = Series([ 1.0, 10.0, 100.0], [1,2,3])
            result = getattr(df,op)(x,level='third',axis=0)

            expected = pd.concat([ opa(df.loc[idx[:,:,i],:],v) for i, v in x.iteritems() ]).sortlevel()
            assert_frame_equal(result, expected)

            x = Series([ 1.0, 10.0], ['two','three'])
            result = getattr(df,op)(x,level='second',axis=0)

            expected = pd.concat([ opa(df.loc[idx[:,i],:],v) for i, v in x.iteritems() ]).reindex_like(df).sortlevel()
            assert_frame_equal(result, expected)

        ## GH9463 (alignment level of dataframe with series)

        midx = MultiIndex.from_product([['A', 'B'],['a', 'b']])
        df = DataFrame(np.ones((2,4), dtype='int64'), columns=midx)
        s = pd.Series({'a':1, 'b':2})

        df2 = df.copy()
        df2.columns.names = ['lvl0', 'lvl1']
        s2 = s.copy()
        s2.index.name = 'lvl1'

        # different cases of integer/string level names:
        res1 = df.mul(s, axis=1, level=1)
        res2 = df.mul(s2, axis=1, level=1)
        res3 = df2.mul(s, axis=1, level=1)
        res4 = df2.mul(s2, axis=1, level=1)
        res5 = df2.mul(s, axis=1, level='lvl1')
        res6 = df2.mul(s2, axis=1, level='lvl1')

        exp = DataFrame(np.array([[1, 2, 1, 2], [1, 2, 1, 2]], dtype='int64'),
                        columns=midx)

        for res in [res1, res2]:
            assert_frame_equal(res, exp)

        exp.columns.names = ['lvl0', 'lvl1']
        for res in [res3, res4, res5, res6]:
            assert_frame_equal(res, exp)

    def test_arith_mixed(self):

        left = DataFrame({'A': ['a', 'b', 'c'],
                          'B': [1, 2, 3]})

        result = left + left
        expected = DataFrame({'A': ['aa', 'bb', 'cc'],
                              'B': [2, 4, 6]})
        assert_frame_equal(result, expected)

    def test_arith_getitem_commute(self):
        df = DataFrame({'A': [1.1, 3.3], 'B': [2.5, -3.9]})

        self._test_op(df, operator.add)
        self._test_op(df, operator.sub)
        self._test_op(df, operator.mul)
        self._test_op(df, operator.truediv)
        self._test_op(df, operator.floordiv)
        self._test_op(df, operator.pow)

        self._test_op(df, lambda x, y: y + x)
        self._test_op(df, lambda x, y: y - x)
        self._test_op(df, lambda x, y: y * x)
        self._test_op(df, lambda x, y: y / x)
        self._test_op(df, lambda x, y: y ** x)

        self._test_op(df, lambda x, y: x + y)
        self._test_op(df, lambda x, y: x - y)
        self._test_op(df, lambda x, y: x * y)
        self._test_op(df, lambda x, y: x / y)
        self._test_op(df, lambda x, y: x ** y)

    @staticmethod
    def _test_op(df, op):
        result = op(df, 1)

        if not df.columns.is_unique:
            raise ValueError("Only unique columns supported by this test")

        for col in result.columns:
            assert_series_equal(result[col], op(df[col], 1))

    def test_bool_flex_frame(self):
        data = np.random.randn(5, 3)
        other_data = np.random.randn(5, 3)
        df = DataFrame(data)
        other = DataFrame(other_data)
        ndim_5 = np.ones(df.shape + (1, 3))

        # Unaligned
        def _check_unaligned_frame(meth, op, df, other):
            part_o = other.ix[3:, 1:].copy()
            rs = meth(part_o)
            xp = op(df, part_o.reindex(index=df.index, columns=df.columns))
            assert_frame_equal(rs, xp)

        # DataFrame
        self.assertTrue(df.eq(df).values.all())
        self.assertFalse(df.ne(df).values.any())
        for op in ['eq', 'ne', 'gt', 'lt', 'ge', 'le']:
            f = getattr(df, op)
            o = getattr(operator, op)
            # No NAs
            assert_frame_equal(f(other), o(df, other))
            _check_unaligned_frame(f, o, df, other)
            # ndarray
            assert_frame_equal(f(other.values), o(df, other.values))
            # scalar
            assert_frame_equal(f(0), o(df, 0))
            # NAs
            assert_frame_equal(f(np.nan), o(df, np.nan))
            with assertRaisesRegexp(ValueError, 'shape'):
                f(ndim_5)

        # Series
        def _test_seq(df, idx_ser, col_ser):
            idx_eq = df.eq(idx_ser, axis=0)
            col_eq = df.eq(col_ser)
            idx_ne = df.ne(idx_ser, axis=0)
            col_ne = df.ne(col_ser)
            assert_frame_equal(col_eq, df == Series(col_ser))
            assert_frame_equal(col_eq, -col_ne)
            assert_frame_equal(idx_eq, -idx_ne)
            assert_frame_equal(idx_eq, df.T.eq(idx_ser).T)
            assert_frame_equal(col_eq, df.eq(list(col_ser)))
            assert_frame_equal(idx_eq, df.eq(Series(idx_ser), axis=0))
            assert_frame_equal(idx_eq, df.eq(list(idx_ser), axis=0))

            idx_gt = df.gt(idx_ser, axis=0)
            col_gt = df.gt(col_ser)
            idx_le = df.le(idx_ser, axis=0)
            col_le = df.le(col_ser)

            assert_frame_equal(col_gt, df > Series(col_ser))
            assert_frame_equal(col_gt, -col_le)
            assert_frame_equal(idx_gt, -idx_le)
            assert_frame_equal(idx_gt, df.T.gt(idx_ser).T)

            idx_ge = df.ge(idx_ser, axis=0)
            col_ge = df.ge(col_ser)
            idx_lt = df.lt(idx_ser, axis=0)
            col_lt = df.lt(col_ser)
            assert_frame_equal(col_ge, df >= Series(col_ser))
            assert_frame_equal(col_ge, -col_lt)
            assert_frame_equal(idx_ge, -idx_lt)
            assert_frame_equal(idx_ge, df.T.ge(idx_ser).T)

        idx_ser = Series(np.random.randn(5))
        col_ser = Series(np.random.randn(3))
        _test_seq(df, idx_ser, col_ser)


        # list/tuple
        _test_seq(df, idx_ser.values, col_ser.values)

        # NA
        df.ix[0, 0] = np.nan
        rs = df.eq(df)
        self.assertFalse(rs.ix[0, 0])
        rs = df.ne(df)
        self.assertTrue(rs.ix[0, 0])
        rs = df.gt(df)
        self.assertFalse(rs.ix[0, 0])
        rs = df.lt(df)
        self.assertFalse(rs.ix[0, 0])
        rs = df.ge(df)
        self.assertFalse(rs.ix[0, 0])
        rs = df.le(df)
        self.assertFalse(rs.ix[0, 0])



        # complex
        arr = np.array([np.nan, 1, 6, np.nan])
        arr2 = np.array([2j, np.nan, 7, None])
        df = DataFrame({'a': arr})
        df2 = DataFrame({'a': arr2})
        rs = df.gt(df2)
        self.assertFalse(rs.values.any())
        rs = df.ne(df2)
        self.assertTrue(rs.values.all())

        arr3 = np.array([2j, np.nan, None])
        df3 = DataFrame({'a': arr3})
        rs = df3.gt(2j)
        self.assertFalse(rs.values.any())

        # corner, dtype=object
        df1 = DataFrame({'col': ['foo', np.nan, 'bar']})
        df2 = DataFrame({'col': ['foo', datetime.now(), 'bar']})
        result = df1.ne(df2)
        exp = DataFrame({'col': [False, True, False]})
        assert_frame_equal(result, exp)

    def test_arith_flex_series(self):
        df = self.simple

        row = df.xs('a')
        col = df['two']
        # after arithmetic refactor, add truediv here
        ops = ['add', 'sub', 'mul', 'mod']
        for op in ops:
            f = getattr(df, op)
            op = getattr(operator, op)
            assert_frame_equal(f(row), op(df, row))
            assert_frame_equal(f(col, axis=0), op(df.T, col).T)

        # special case for some reason
        assert_frame_equal(df.add(row, axis=None), df + row)

        # cases which will be refactored after big arithmetic refactor
        assert_frame_equal(df.div(row), df / row)
        assert_frame_equal(df.div(col, axis=0), (df.T / col).T)

        # broadcasting issue in GH7325
        df = DataFrame(np.arange(3*2).reshape((3,2)),dtype='int64')
        expected = DataFrame([[nan, inf], [1.0, 1.5], [1.0, 1.25]])
        result = df.div(df[0],axis='index')
        assert_frame_equal(result,expected)

        df = DataFrame(np.arange(3*2).reshape((3,2)),dtype='float64')
        expected = DataFrame([[np.nan,np.inf],[1.0,1.5],[1.0,1.25]])
        result = df.div(df[0],axis='index')
        assert_frame_equal(result,expected)

    def test_arith_non_pandas_object(self):
        df = self.simple

        val1 = df.xs('a').values
        added = DataFrame(df.values + val1, index=df.index, columns=df.columns)
        assert_frame_equal(df + val1, added)

        added = DataFrame((df.values.T + val1).T,
                          index=df.index, columns=df.columns)
        assert_frame_equal(df.add(val1, axis=0), added)

        val2 = list(df['two'])

        added = DataFrame(df.values + val2, index=df.index, columns=df.columns)
        assert_frame_equal(df + val2, added)

        added = DataFrame((df.values.T + val2).T, index=df.index,
                          columns=df.columns)
        assert_frame_equal(df.add(val2, axis='index'), added)

        val3 = np.random.rand(*df.shape)
        added = DataFrame(df.values + val3, index=df.index, columns=df.columns)
        assert_frame_equal(df.add(val3), added)

    def test_combineFrame(self):
        frame_copy = self.frame.reindex(self.frame.index[::2])

        del frame_copy['D']
        frame_copy['C'][:5] = nan

        added = self.frame + frame_copy
        tm.assert_dict_equal(added['A'].valid(),
                             self.frame['A'] * 2,
                             compare_keys=False)

        self.assertTrue(np.isnan(added['C'].reindex(frame_copy.index)[:5]).all())

        # assert(False)

        self.assertTrue(np.isnan(added['D']).all())

        self_added = self.frame + self.frame
        self.assertTrue(self_added.index.equals(self.frame.index))

        added_rev = frame_copy + self.frame
        self.assertTrue(np.isnan(added['D']).all())

        # corner cases

        # empty
        plus_empty = self.frame + self.empty
        self.assertTrue(np.isnan(plus_empty.values).all())

        empty_plus = self.empty + self.frame
        self.assertTrue(np.isnan(empty_plus.values).all())

        empty_empty = self.empty + self.empty
        self.assertTrue(empty_empty.empty)

        # out of order
        reverse = self.frame.reindex(columns=self.frame.columns[::-1])

        assert_frame_equal(reverse + self.frame, self.frame * 2)

        # mix vs float64, upcast
        added = self.frame + self.mixed_float
        _check_mixed_float(added, dtype = 'float64')
        added = self.mixed_float + self.frame
        _check_mixed_float(added, dtype = 'float64')

        # mix vs mix
        added = self.mixed_float + self.mixed_float2
        _check_mixed_float(added, dtype = dict(C = None))
        added = self.mixed_float2 + self.mixed_float
        _check_mixed_float(added, dtype = dict(C = None))

        # with int
        added = self.frame + self.mixed_int
        _check_mixed_float(added, dtype = 'float64')

    def test_combineSeries(self):

        # Series
        series = self.frame.xs(self.frame.index[0])

        added = self.frame + series

        for key, s in compat.iteritems(added):
            assert_series_equal(s, self.frame[key] + series[key])

        larger_series = series.to_dict()
        larger_series['E'] = 1
        larger_series = Series(larger_series)
        larger_added = self.frame + larger_series

        for key, s in compat.iteritems(self.frame):
            assert_series_equal(larger_added[key], s + series[key])
        self.assertIn('E', larger_added)
        self.assertTrue(np.isnan(larger_added['E']).all())

        # vs mix (upcast) as needed
        added = self.mixed_float + series
        _check_mixed_float(added, dtype = 'float64')
        added = self.mixed_float + series.astype('float32')
        _check_mixed_float(added, dtype = dict(C = None))
        added = self.mixed_float + series.astype('float16')
        _check_mixed_float(added, dtype = dict(C = None))

        #### these raise with numexpr.....as we are adding an int64 to an uint64....weird
        # vs int
        #added = self.mixed_int + (100*series).astype('int64')
        #_check_mixed_int(added, dtype = dict(A = 'int64', B = 'float64', C = 'int64', D = 'int64'))
        #added = self.mixed_int + (100*series).astype('int32')
        #_check_mixed_int(added, dtype = dict(A = 'int32', B = 'float64', C = 'int32', D = 'int64'))


        # TimeSeries
        ts = self.tsframe['A']

        # 10890
        # we no longer allow auto timeseries broadcasting
        # and require explict broadcasting
        added = self.tsframe.add(ts, axis='index')

        for key, col in compat.iteritems(self.tsframe):
            result = col + ts
            assert_series_equal(added[key], result, check_names=False)
            self.assertEqual(added[key].name, key)
            if col.name == ts.name:
                self.assertEqual(result.name, 'A')
            else:
                self.assertTrue(result.name is None)

        smaller_frame = self.tsframe[:-5]
        smaller_added = smaller_frame.add(ts, axis='index')

        self.assertTrue(smaller_added.index.equals(self.tsframe.index))

        smaller_ts = ts[:-5]
        smaller_added2 = self.tsframe.add(smaller_ts, axis='index')
        assert_frame_equal(smaller_added, smaller_added2)

        # length 0, result is all-nan
        result = self.tsframe.add(ts[:0], axis='index')
        expected = DataFrame(np.nan,index=self.tsframe.index,columns=self.tsframe.columns)
        assert_frame_equal(result, expected)

        # Frame is all-nan
        result = self.tsframe[:0].add(ts, axis='index')
        expected = DataFrame(np.nan,index=self.tsframe.index,columns=self.tsframe.columns)
        assert_frame_equal(result, expected)

        # empty but with non-empty index
        frame = self.tsframe[:1].reindex(columns=[])
        result = frame.mul(ts,axis='index')
        self.assertEqual(len(result), len(ts))

    def test_combineFunc(self):
        result = self.frame * 2
        self.assert_numpy_array_equal(result.values, self.frame.values * 2)

        # vs mix
        result = self.mixed_float * 2
        for c, s in compat.iteritems(result):
            self.assert_numpy_array_equal(s.values, self.mixed_float[c].values * 2)
        _check_mixed_float(result, dtype = dict(C = None))

        result = self.empty * 2
        self.assertIs(result.index, self.empty.index)
        self.assertEqual(len(result.columns), 0)

    def test_comparisons(self):
        df1 = tm.makeTimeDataFrame()
        df2 = tm.makeTimeDataFrame()

        row = self.simple.xs('a')
        ndim_5 = np.ones(df1.shape + (1, 1, 1))

        def test_comp(func):
            result = func(df1, df2)
            self.assert_numpy_array_equal(result.values,
                                          func(df1.values, df2.values))
            with assertRaisesRegexp(ValueError, 'Wrong number of dimensions'):
                func(df1, ndim_5)

            result2 = func(self.simple, row)
            self.assert_numpy_array_equal(result2.values,
                                          func(self.simple.values, row.values))

            result3 = func(self.frame, 0)
            self.assert_numpy_array_equal(result3.values,
                                          func(self.frame.values, 0))


            with assertRaisesRegexp(ValueError, 'Can only compare '
                                    'identically-labeled DataFrame'):
                func(self.simple, self.simple[:2])

        test_comp(operator.eq)
        test_comp(operator.ne)
        test_comp(operator.lt)
        test_comp(operator.gt)
        test_comp(operator.ge)
        test_comp(operator.le)

    def test_string_comparison(self):
        df = DataFrame([{"a": 1, "b": "foo"}, {"a": 2, "b": "bar"}])
        mask_a = df.a > 1
        assert_frame_equal(df[mask_a], df.ix[1:1, :])
        assert_frame_equal(df[-mask_a], df.ix[0:0, :])

        mask_b = df.b == "foo"
        assert_frame_equal(df[mask_b], df.ix[0:0, :])
        assert_frame_equal(df[-mask_b], df.ix[1:1, :])

    def test_float_none_comparison(self):
        df = DataFrame(np.random.randn(8, 3), index=lrange(8),
                       columns=['A', 'B', 'C'])

        self.assertRaises(TypeError, df.__eq__, None)

    def test_boolean_comparison(self):

        # GH 4576
        # boolean comparisons with a tuple/list give unexpected results
        df = DataFrame(np.arange(6).reshape((3,2)))
        b = np.array([2, 2])
        b_r = np.atleast_2d([2,2])
        b_c = b_r.T
        l = (2,2,2)
        tup = tuple(l)

        # gt
        expected = DataFrame([[False,False],[False,True],[True,True]])
        result = df>b
        assert_frame_equal(result,expected)

        result = df.values>b
        assert_numpy_array_equal(result,expected.values)

        result = df>l
        assert_frame_equal(result,expected)

        result = df>tup
        assert_frame_equal(result,expected)

        result = df>b_r
        assert_frame_equal(result,expected)

        result = df.values>b_r
        assert_numpy_array_equal(result,expected.values)

        self.assertRaises(ValueError, df.__gt__, b_c)
        self.assertRaises(ValueError, df.values.__gt__, b_c)

        # ==
        expected = DataFrame([[False,False],[True,False],[False,False]])
        result = df == b
        assert_frame_equal(result,expected)

        result = df==l
        assert_frame_equal(result,expected)

        result = df==tup
        assert_frame_equal(result,expected)

        result = df == b_r
        assert_frame_equal(result,expected)

        result = df.values == b_r
        assert_numpy_array_equal(result,expected.values)

        self.assertRaises(ValueError, lambda : df == b_c)
        self.assertFalse((df.values == b_c))

        # with alignment
        df = DataFrame(np.arange(6).reshape((3,2)),columns=list('AB'),index=list('abc'))
        expected.index=df.index
        expected.columns=df.columns

        result = df==l
        assert_frame_equal(result,expected)

        result = df==tup
        assert_frame_equal(result,expected)

        # not shape compatible
        self.assertRaises(ValueError, lambda : df == (2,2))
        self.assertRaises(ValueError, lambda : df == [2,2])

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

    def test_to_csv_from_csv(self):

        pname = '__tmp_to_csv_from_csv__'
        with ensure_clean(pname) as path:

            self.frame['A'][:5] = nan

            self.frame.to_csv(path)
            self.frame.to_csv(path, columns=['A', 'B'])
            self.frame.to_csv(path, header=False)
            self.frame.to_csv(path, index=False)

            # test roundtrip
            self.tsframe.to_csv(path)
            recons = DataFrame.from_csv(path)

            assert_frame_equal(self.tsframe, recons)

            self.tsframe.to_csv(path, index_label='index')
            recons = DataFrame.from_csv(path, index_col=None)
            assert(len(recons.columns) == len(self.tsframe.columns) + 1)

            # no index
            self.tsframe.to_csv(path, index=False)
            recons = DataFrame.from_csv(path, index_col=None)
            assert_almost_equal(self.tsframe.values, recons.values)

            # corner case
            dm = DataFrame({'s1': Series(lrange(3), lrange(3)),
                            's2': Series(lrange(2), lrange(2))})
            dm.to_csv(path)
            recons = DataFrame.from_csv(path)
            assert_frame_equal(dm, recons)

        with ensure_clean(pname) as path:

            # duplicate index
            df = DataFrame(np.random.randn(3, 3), index=['a', 'a', 'b'],
                           columns=['x', 'y', 'z'])
            df.to_csv(path)
            result = DataFrame.from_csv(path)
            assert_frame_equal(result, df)

            midx = MultiIndex.from_tuples([('A', 1, 2), ('A', 1, 2), ('B', 1, 2)])
            df = DataFrame(np.random.randn(3, 3), index=midx,
                           columns=['x', 'y', 'z'])
            df.to_csv(path)
            result = DataFrame.from_csv(path, index_col=[0, 1, 2],
                                        parse_dates=False)
            assert_frame_equal(result, df, check_names=False)  # TODO from_csv names index ['Unnamed: 1', 'Unnamed: 2'] should it ?

            # column aliases
            col_aliases = Index(['AA', 'X', 'Y', 'Z'])
            self.frame2.to_csv(path, header=col_aliases)
            rs = DataFrame.from_csv(path)
            xp = self.frame2.copy()
            xp.columns = col_aliases

            assert_frame_equal(xp, rs)

            self.assertRaises(ValueError, self.frame2.to_csv, path,
                              header=['AA', 'X'])

        with ensure_clean(pname) as path:
            df1 = DataFrame(np.random.randn(3, 1))
            df2 = DataFrame(np.random.randn(3, 1))

            df1.to_csv(path)
            df2.to_csv(path,mode='a',header=False)
            xp = pd.concat([df1,df2])
            rs = pd.read_csv(path,index_col=0)
            rs.columns = lmap(int,rs.columns)
            xp.columns = lmap(int,xp.columns)
            assert_frame_equal(xp,rs)

        with ensure_clean() as path:
            # GH 10833 (TimedeltaIndex formatting)
            dt = pd.Timedelta(seconds=1)
            df = pd.DataFrame({'dt_data': [i*dt for i in range(3)]},
                              index=pd.Index([i*dt for i in range(3)],
                                             name='dt_index'))
            df.to_csv(path)

            result = pd.read_csv(path, index_col='dt_index')
            result.index = pd.to_timedelta(result.index)
            # TODO: remove renaming when GH 10875 is solved
            result.index = result.index.rename('dt_index')
            result['dt_data'] = pd.to_timedelta(result['dt_data'])

            assert_frame_equal(df, result, check_index_type=True)

        # tz, 8260
        with ensure_clean(pname) as path:

            self.tzframe.to_csv(path)
            result = pd.read_csv(path, index_col=0, parse_dates=['A'])

            converter = lambda c: pd.to_datetime(result[c]).dt.tz_localize('UTC').dt.tz_convert(self.tzframe[c].dt.tz)
            result['B'] = converter('B')
            result['C'] = converter('C')
            assert_frame_equal(result, self.tzframe)

    def test_to_csv_cols_reordering(self):
        # GH3454
        import pandas as pd

        def _check_df(df,cols=None):
            with ensure_clean() as path:
                df.to_csv(path,columns = cols,engine='python')
                rs_p = pd.read_csv(path,index_col=0)
                df.to_csv(path,columns = cols,chunksize=chunksize)
                rs_c = pd.read_csv(path,index_col=0)

            if cols:
                df = df[cols]
            assert (rs_c.columns==rs_p.columns).all()
            assert_frame_equal(df,rs_c,check_names=False)

        chunksize=5
        N = int(chunksize*2.5)

        df= mkdf(N, 3)
        cs = df.columns
        cols = [cs[2],cs[0]]
        _check_df(df,cols)

    def test_to_csv_legacy_raises_on_dupe_cols(self):
        df= mkdf(10, 3)
        df.columns = ['a','a','b']
        with ensure_clean() as path:
            self.assertRaises(NotImplementedError,df.to_csv,path,engine='python')

    def test_to_csv_new_dupe_cols(self):
        import pandas as pd
        def _check_df(df,cols=None):
            with ensure_clean() as path:
                df.to_csv(path,columns = cols,chunksize=chunksize)
                rs_c = pd.read_csv(path,index_col=0)

                # we wrote them in a different order
                # so compare them in that order
                if cols is not None:

                    if df.columns.is_unique:
                        rs_c.columns = cols
                    else:
                        indexer, missing = df.columns.get_indexer_non_unique(cols)
                        rs_c.columns = df.columns.take(indexer)

                    for c in cols:
                        obj_df = df[c]
                        obj_rs = rs_c[c]
                        if isinstance(obj_df,Series):
                            assert_series_equal(obj_df,obj_rs)
                        else:
                            assert_frame_equal(obj_df,obj_rs,check_names=False)

                # wrote in the same order
                else:
                    rs_c.columns = df.columns
                    assert_frame_equal(df,rs_c,check_names=False)

        chunksize=5
        N = int(chunksize*2.5)

        # dupe cols
        df= mkdf(N, 3)
        df.columns = ['a','a','b']
        _check_df(df,None)

        # dupe cols with selection
        cols = ['b','a']
        _check_df(df,cols)

    @slow
    def test_to_csv_moar(self):
        path = '__tmp_to_csv_moar__'

        def _do_test(df,path,r_dtype=None,c_dtype=None,rnlvl=None,cnlvl=None,
                     dupe_col=False):

            kwargs = dict(parse_dates=False)
            if cnlvl:
                if rnlvl is not None:
                    kwargs['index_col'] = lrange(rnlvl)
                kwargs['header'] = lrange(cnlvl)
                with ensure_clean(path) as path:
                    df.to_csv(path,encoding='utf8',chunksize=chunksize,tupleize_cols=False)
                    recons = DataFrame.from_csv(path,tupleize_cols=False,**kwargs)
            else:
                kwargs['header'] = 0
                with ensure_clean(path) as path:
                    df.to_csv(path,encoding='utf8',chunksize=chunksize)
                    recons = DataFrame.from_csv(path,**kwargs)

            def _to_uni(x):
                if not isinstance(x, compat.text_type):
                    return x.decode('utf8')
                return x
            if dupe_col:
                # read_Csv disambiguates the columns by
                # labeling them dupe.1,dupe.2, etc'. monkey patch columns
                recons.columns = df.columns
            if rnlvl and not cnlvl:
                delta_lvl = [recons.iloc[:, i].values for i in range(rnlvl-1)]
                ix=MultiIndex.from_arrays([list(recons.index)]+delta_lvl)
                recons.index = ix
                recons = recons.iloc[:,rnlvl-1:]

            type_map = dict(i='i',f='f',s='O',u='O',dt='O',p='O')
            if r_dtype:
                if r_dtype == 'u': # unicode
                    r_dtype='O'
                    recons.index = np.array(lmap(_to_uni,recons.index),
                                            dtype=r_dtype)
                    df.index = np.array(lmap(_to_uni,df.index),dtype=r_dtype)
                elif r_dtype == 'dt': # unicode
                    r_dtype='O'
                    recons.index = np.array(lmap(Timestamp,recons.index),
                                            dtype=r_dtype)
                    df.index = np.array(lmap(Timestamp,df.index),dtype=r_dtype)
                elif r_dtype == 'p':
                    r_dtype='O'
                    recons.index = np.array(list(map(Timestamp,
                                                     recons.index.to_datetime())),
                                            dtype=r_dtype)
                    df.index = np.array(list(map(Timestamp,
                                                 df.index.to_datetime())),
                                        dtype=r_dtype)
                else:
                    r_dtype= type_map.get(r_dtype)
                    recons.index = np.array(recons.index,dtype=r_dtype )
                    df.index = np.array(df.index,dtype=r_dtype )
            if c_dtype:
                if c_dtype == 'u':
                    c_dtype='O'
                    recons.columns = np.array(lmap(_to_uni,recons.columns),
                                              dtype=c_dtype)
                    df.columns = np.array(lmap(_to_uni,df.columns),dtype=c_dtype )
                elif c_dtype == 'dt':
                    c_dtype='O'
                    recons.columns = np.array(lmap(Timestamp,recons.columns),
                                                dtype=c_dtype )
                    df.columns = np.array(lmap(Timestamp,df.columns),dtype=c_dtype)
                elif c_dtype == 'p':
                    c_dtype='O'
                    recons.columns = np.array(lmap(Timestamp,recons.columns.to_datetime()),
                                              dtype=c_dtype)
                    df.columns = np.array(lmap(Timestamp,df.columns.to_datetime()),dtype=c_dtype )
                else:
                    c_dtype= type_map.get(c_dtype)
                    recons.columns = np.array(recons.columns,dtype=c_dtype )
                    df.columns = np.array(df.columns,dtype=c_dtype )

            assert_frame_equal(df,recons,check_names=False,check_less_precise=True)

        N = 100
        chunksize=1000

        # GH3437
        from pandas import NaT
        def make_dtnat_arr(n,nnat=None):
            if nnat is None:
                nnat= int(n*0.1) # 10%
            s=list(date_range('2000',freq='5min',periods=n))
            if nnat:
                for i in np.random.randint(0,len(s),nnat):
                    s[i] = NaT
                i = np.random.randint(100)
                s[-i] = NaT
                s[i] = NaT
            return s

        # N=35000
        s1=make_dtnat_arr(chunksize+5)
        s2=make_dtnat_arr(chunksize+5,0)
        path = '1.csv'

        # s3=make_dtnjat_arr(chunksize+5,0)
        with ensure_clean('.csv') as pth:
            df=DataFrame(dict(a=s1,b=s2))
            df.to_csv(pth,chunksize=chunksize)
            recons = DataFrame.from_csv(pth)._convert(datetime=True,
                                                             coerce=True)
            assert_frame_equal(df, recons,check_names=False,check_less_precise=True)

        for ncols in [4]:
            base = int((chunksize// ncols or 1) or 1)
            for nrows in [2,10,N-1,N,N+1,N+2,2*N-2,2*N-1,2*N,2*N+1,2*N+2,
                  base-1,base,base+1]:
                _do_test(mkdf(nrows, ncols,r_idx_type='dt',
                              c_idx_type='s'),path, 'dt','s')


        for ncols in [4]:
            base = int((chunksize// ncols or 1) or 1)
            for nrows in [2,10,N-1,N,N+1,N+2,2*N-2,2*N-1,2*N,2*N+1,2*N+2,
                  base-1,base,base+1]:
                _do_test(mkdf(nrows, ncols,r_idx_type='dt',
                              c_idx_type='s'),path, 'dt','s')
                pass

        for r_idx_type,c_idx_type  in [('i','i'),('s','s'),('u','dt'),('p','p')]:
            for ncols in [1,2,3,4]:
                base = int((chunksize// ncols or 1) or 1)
                for nrows in [2,10,N-1,N,N+1,N+2,2*N-2,2*N-1,2*N,2*N+1,2*N+2,
                      base-1,base,base+1]:
                    _do_test(mkdf(nrows, ncols,r_idx_type=r_idx_type,
                                  c_idx_type=c_idx_type),path,r_idx_type,c_idx_type)

        for ncols in [1,2,3,4]:
            base = int((chunksize// ncols or 1) or 1)
            for nrows in [10,N-2,N-1,N,N+1,N+2,2*N-2,2*N-1,2*N,2*N+1,2*N+2,
                      base-1,base,base+1]:
                _do_test(mkdf(nrows, ncols),path)

        for nrows in [10,N-2,N-1,N,N+1,N+2]:
            df = mkdf(nrows, 3)
            cols = list(df.columns)
            cols[:2] = ["dupe","dupe"]
            cols[-2:] = ["dupe","dupe"]
            ix = list(df.index)
            ix[:2] = ["rdupe","rdupe"]
            ix[-2:] = ["rdupe","rdupe"]
            df.index=ix
            df.columns=cols
            _do_test(df,path,dupe_col=True)


        _do_test(DataFrame(index=lrange(10)),path)
        _do_test(mkdf(chunksize//2+1, 2,r_idx_nlevels=2),path,rnlvl=2)
        for ncols in [2,3,4]:
            base = int(chunksize//ncols)
            for nrows in [10,N-2,N-1,N,N+1,N+2,2*N-2,2*N-1,2*N,2*N+1,2*N+2,
                      base-1,base,base+1]:
                _do_test(mkdf(nrows, ncols,r_idx_nlevels=2),path,rnlvl=2)
                _do_test(mkdf(nrows, ncols,c_idx_nlevels=2),path,cnlvl=2)
                _do_test(mkdf(nrows, ncols,r_idx_nlevels=2,c_idx_nlevels=2),
                         path,rnlvl=2,cnlvl=2)

    def test_to_csv_from_csv_w_some_infs(self):

        # test roundtrip with inf, -inf, nan, as full columns and mix
        self.frame['G'] = np.nan
        f = lambda x: [np.inf, np.nan][np.random.rand() < .5]
        self.frame['H'] = self.frame.index.map(f)

        with ensure_clean() as path:
            self.frame.to_csv(path)
            recons = DataFrame.from_csv(path)

            assert_frame_equal(self.frame, recons, check_names=False)  # TODO to_csv drops column name
            assert_frame_equal(np.isinf(self.frame), np.isinf(recons), check_names=False)

    def test_to_csv_from_csv_w_all_infs(self):

        # test roundtrip with inf, -inf, nan, as full columns and mix
        self.frame['E'] = np.inf
        self.frame['F'] = -np.inf

        with ensure_clean() as path:
            self.frame.to_csv(path)
            recons = DataFrame.from_csv(path)

            assert_frame_equal(self.frame, recons, check_names=False)  # TODO to_csv drops column name
            assert_frame_equal(np.isinf(self.frame), np.isinf(recons), check_names=False)

    def test_to_csv_no_index(self):
        # GH 3624, after appending columns, to_csv fails
        pname = '__tmp_to_csv_no_index__'
        with ensure_clean(pname) as path:
            df = DataFrame({'c1':[1,2,3], 'c2':[4,5,6]})
            df.to_csv(path, index=False)
            result = read_csv(path)
            assert_frame_equal(df,result)
            df['c3'] = Series([7,8,9],dtype='int64')
            df.to_csv(path, index=False)
            result = read_csv(path)
            assert_frame_equal(df,result)

    def test_to_csv_headers(self):
        # GH6186, the presence or absence of `index` incorrectly
        # causes to_csv to have different header semantics.
        pname = '__tmp_to_csv_headers__'
        from_df = DataFrame([[1, 2], [3, 4]], columns=['A', 'B'])
        to_df  = DataFrame([[1, 2], [3, 4]], columns=['X', 'Y'])
        with ensure_clean(pname) as path:
            from_df.to_csv(path, header=['X', 'Y'])
            recons = DataFrame.from_csv(path)
            assert_frame_equal(to_df, recons)

            from_df.to_csv(path, index=False, header=['X', 'Y'])
            recons = DataFrame.from_csv(path)
            recons.reset_index(inplace=True)
            assert_frame_equal(to_df, recons)

    def test_to_csv_multiindex(self):

        pname = '__tmp_to_csv_multiindex__'
        frame = self.frame
        old_index = frame.index
        arrays = np.arange(len(old_index) * 2).reshape(2, -1)
        new_index = MultiIndex.from_arrays(arrays, names=['first', 'second'])
        frame.index = new_index

        with ensure_clean(pname) as path:

            frame.to_csv(path, header=False)
            frame.to_csv(path, columns=['A', 'B'])

            # round trip
            frame.to_csv(path)
            df = DataFrame.from_csv(path, index_col=[0, 1], parse_dates=False)

            assert_frame_equal(frame, df, check_names=False)  # TODO to_csv drops column name
            self.assertEqual(frame.index.names, df.index.names)
            self.frame.index = old_index  # needed if setUP becomes a classmethod

            # try multiindex with dates
            tsframe = self.tsframe
            old_index = tsframe.index
            new_index = [old_index, np.arange(len(old_index))]
            tsframe.index = MultiIndex.from_arrays(new_index)

            tsframe.to_csv(path, index_label=['time', 'foo'])
            recons = DataFrame.from_csv(path, index_col=[0, 1])
            assert_frame_equal(tsframe, recons, check_names=False)  # TODO to_csv drops column name

            # do not load index
            tsframe.to_csv(path)
            recons = DataFrame.from_csv(path, index_col=None)
            np.testing.assert_equal(len(recons.columns), len(tsframe.columns) + 2)

            # no index
            tsframe.to_csv(path, index=False)
            recons = DataFrame.from_csv(path, index_col=None)
            assert_almost_equal(recons.values, self.tsframe.values)
            self.tsframe.index = old_index  # needed if setUP becomes classmethod

        with ensure_clean(pname) as path:
            # GH3571, GH1651, GH3141

            def _make_frame(names=None):
                if names is True:
                    names = ['first','second']
                return DataFrame(np.random.randint(0,10,size=(3,3)),
                                 columns=MultiIndex.from_tuples([('bah', 'foo'),
                                                                 ('bah', 'bar'),
                                                                 ('ban', 'baz')],
                                                                names=names),
                                 dtype='int64')

            # column & index are multi-index
            df = mkdf(5,3,r_idx_nlevels=2,c_idx_nlevels=4)
            df.to_csv(path,tupleize_cols=False)
            result = read_csv(path,header=[0,1,2,3],index_col=[0,1],tupleize_cols=False)
            assert_frame_equal(df,result)

            # column is mi
            df = mkdf(5,3,r_idx_nlevels=1,c_idx_nlevels=4)
            df.to_csv(path,tupleize_cols=False)
            result = read_csv(path,header=[0,1,2,3],index_col=0,tupleize_cols=False)
            assert_frame_equal(df,result)

            # dup column names?
            df = mkdf(5,3,r_idx_nlevels=3,c_idx_nlevels=4)
            df.to_csv(path,tupleize_cols=False)
            result = read_csv(path,header=[0,1,2,3],index_col=[0,1,2],tupleize_cols=False)
            assert_frame_equal(df,result)

            # writing with no index
            df = _make_frame()
            df.to_csv(path,tupleize_cols=False,index=False)
            result = read_csv(path,header=[0,1],tupleize_cols=False)
            assert_frame_equal(df,result)

            # we lose the names here
            df = _make_frame(True)
            df.to_csv(path,tupleize_cols=False,index=False)
            result = read_csv(path,header=[0,1],tupleize_cols=False)
            self.assertTrue(all([ x is None for x in result.columns.names ]))
            result.columns.names = df.columns.names
            assert_frame_equal(df,result)

            # tupleize_cols=True and index=False
            df = _make_frame(True)
            df.to_csv(path,tupleize_cols=True,index=False)
            result = read_csv(path,header=0,tupleize_cols=True,index_col=None)
            result.columns = df.columns
            assert_frame_equal(df,result)

            # whatsnew example
            df = _make_frame()
            df.to_csv(path,tupleize_cols=False)
            result = read_csv(path,header=[0,1],index_col=[0],tupleize_cols=False)
            assert_frame_equal(df,result)

            df = _make_frame(True)
            df.to_csv(path,tupleize_cols=False)
            result = read_csv(path,header=[0,1],index_col=[0],tupleize_cols=False)
            assert_frame_equal(df,result)

            # column & index are multi-index (compatibility)
            df = mkdf(5,3,r_idx_nlevels=2,c_idx_nlevels=4)
            df.to_csv(path,tupleize_cols=True)
            result = read_csv(path,header=0,index_col=[0,1],tupleize_cols=True)
            result.columns = df.columns
            assert_frame_equal(df,result)

            # invalid options
            df = _make_frame(True)
            df.to_csv(path,tupleize_cols=False)

            # catch invalid headers
            with assertRaisesRegexp(CParserError, 'Passed header=\[0,1,2\] are too many rows for this multi_index of columns'):
                read_csv(path,tupleize_cols=False,header=lrange(3),index_col=0)

            with assertRaisesRegexp(CParserError, 'Passed header=\[0,1,2,3,4,5,6\], len of 7, but only 6 lines in file'):
                read_csv(path,tupleize_cols=False,header=lrange(7),index_col=0)

            for i in [4,5,6]:
                with tm.assertRaises(CParserError):
                    read_csv(path, tupleize_cols=False, header=lrange(i), index_col=0)

            # write with cols
            with assertRaisesRegexp(TypeError, 'cannot specify cols with a MultiIndex'):
                df.to_csv(path, tupleize_cols=False, columns=['foo', 'bar'])

        with ensure_clean(pname) as path:
            # empty
            tsframe[:0].to_csv(path)
            recons = DataFrame.from_csv(path)
            exp = tsframe[:0]
            exp.index = []

            self.assertTrue(recons.columns.equals(exp.columns))
            self.assertEqual(len(recons), 0)

    def test_to_csv_float32_nanrep(self):
        df = DataFrame(np.random.randn(1, 4).astype(np.float32))
        df[1] = np.nan

        with ensure_clean('__tmp_to_csv_float32_nanrep__.csv') as path:
            df.to_csv(path, na_rep=999)

            with open(path) as f:
                lines = f.readlines()
                self.assertEqual(lines[1].split(',')[2], '999')

    def test_to_csv_withcommas(self):

        # Commas inside fields should be correctly escaped when saving as CSV.
        df = DataFrame({'A': [1, 2, 3], 'B': ['5,6', '7,8', '9,0']})

        with ensure_clean('__tmp_to_csv_withcommas__.csv') as path:
            df.to_csv(path)
            df2 = DataFrame.from_csv(path)
            assert_frame_equal(df2, df)

    def test_to_csv_mixed(self):

        def create_cols(name):
            return [ "%s%03d" % (name,i) for i in range(5) ]

        df_float  = DataFrame(np.random.randn(100, 5),dtype='float64',columns=create_cols('float'))
        df_int    = DataFrame(np.random.randn(100, 5),dtype='int64',columns=create_cols('int'))
        df_bool   = DataFrame(True,index=df_float.index,columns=create_cols('bool'))
        df_object = DataFrame('foo',index=df_float.index,columns=create_cols('object'))
        df_dt     = DataFrame(Timestamp('20010101'),index=df_float.index,columns=create_cols('date'))

        # add in some nans
        df_float.ix[30:50,1:3] = np.nan

        #### this is a bug in read_csv right now ####
        #df_dt.ix[30:50,1:3] = np.nan

        df        = pd.concat([ df_float, df_int, df_bool, df_object, df_dt ], axis=1)

        # dtype
        dtypes = dict()
        for n,dtype in [('float',np.float64),('int',np.int64),('bool',np.bool),('object',np.object)]:
            for c in create_cols(n):
                dtypes[c] = dtype

        with ensure_clean() as filename:
            df.to_csv(filename)
            rs = read_csv(filename, index_col=0, dtype=dtypes, parse_dates=create_cols('date'))
            assert_frame_equal(rs, df)

    def test_to_csv_dups_cols(self):

        df        = DataFrame(np.random.randn(1000, 30),columns=lrange(15)+lrange(15),dtype='float64')

        with ensure_clean() as filename:
            df.to_csv(filename) # single dtype, fine
            result = read_csv(filename,index_col=0)
            result.columns = df.columns
            assert_frame_equal(result,df)

        df_float  = DataFrame(np.random.randn(1000, 3),dtype='float64')
        df_int    = DataFrame(np.random.randn(1000, 3),dtype='int64')
        df_bool   = DataFrame(True,index=df_float.index,columns=lrange(3))
        df_object = DataFrame('foo',index=df_float.index,columns=lrange(3))
        df_dt     = DataFrame(Timestamp('20010101'),index=df_float.index,columns=lrange(3))
        df        = pd.concat([ df_float, df_int, df_bool, df_object, df_dt ], axis=1, ignore_index=True)

        cols = []
        for i in range(5):
            cols.extend([0,1,2])
        df.columns = cols

        from pandas import to_datetime
        with ensure_clean() as filename:
            df.to_csv(filename)
            result = read_csv(filename,index_col=0)

            # date cols
            for i in ['0.4','1.4','2.4']:
                result[i] = to_datetime(result[i])

            result.columns = df.columns
            assert_frame_equal(result,df)

        # GH3457
        from pandas.util.testing import makeCustomDataframe as mkdf

        N=10
        df= mkdf(N, 3)
        df.columns = ['a','a','b']

        with ensure_clean() as filename:
            df.to_csv(filename)

            # read_csv will rename the dups columns
            result = read_csv(filename,index_col=0)
            result = result.rename(columns={ 'a.1' : 'a' })
            assert_frame_equal(result,df)

    def test_to_csv_chunking(self):

        aa=DataFrame({'A':lrange(100000)})
        aa['B'] = aa.A + 1.0
        aa['C'] = aa.A + 2.0
        aa['D'] = aa.A + 3.0

        for chunksize in [10000,50000,100000]:
            with ensure_clean() as filename:
                aa.to_csv(filename,chunksize=chunksize)
                rs = read_csv(filename,index_col=0)
                assert_frame_equal(rs, aa)

    def test_to_csv_wide_frame_formatting(self):
        # Issue #8621
        df = DataFrame(np.random.randn(1, 100010), columns=None, index=None)
        with ensure_clean() as filename:
            df.to_csv(filename, header=False, index=False)
            rs = read_csv(filename, header=None)
            assert_frame_equal(rs, df)

    def test_to_csv_bug(self):
        f1 = StringIO('a,1.0\nb,2.0')
        df = DataFrame.from_csv(f1, header=None)
        newdf = DataFrame({'t': df[df.columns[0]]})

        with ensure_clean() as path:
            newdf.to_csv(path)

            recons = read_csv(path, index_col=0)
            assert_frame_equal(recons, newdf, check_names=False)  # don't check_names as t != 1

    def test_to_csv_unicode(self):

        df = DataFrame({u('c/\u03c3'): [1, 2, 3]})
        with ensure_clean() as path:

            df.to_csv(path, encoding='UTF-8')
            df2 = read_csv(path, index_col=0, encoding='UTF-8')
            assert_frame_equal(df, df2)

            df.to_csv(path, encoding='UTF-8', index=False)
            df2 = read_csv(path, index_col=None, encoding='UTF-8')
            assert_frame_equal(df, df2)

    def test_to_csv_unicode_index_col(self):
        buf = StringIO('')
        df = DataFrame(
            [[u("\u05d0"), "d2", "d3", "d4"], ["a1", "a2", "a3", "a4"]],
            columns=[u("\u05d0"),
                     u("\u05d1"), u("\u05d2"), u("\u05d3")],
            index=[u("\u05d0"), u("\u05d1")])

        df.to_csv(buf, encoding='UTF-8')
        buf.seek(0)

        df2 = read_csv(buf, index_col=0, encoding='UTF-8')
        assert_frame_equal(df, df2)

    def test_to_csv_stringio(self):
        buf = StringIO()
        self.frame.to_csv(buf)
        buf.seek(0)
        recons = read_csv(buf, index_col=0)
        assert_frame_equal(recons, self.frame, check_names=False)  # TODO to_csv drops column name

    def test_to_csv_float_format(self):

        df = DataFrame([[0.123456, 0.234567, 0.567567],
                        [12.32112, 123123.2, 321321.2]],
                       index=['A', 'B'], columns=['X', 'Y', 'Z'])

        with ensure_clean() as filename:

            df.to_csv(filename, float_format='%.2f')

            rs = read_csv(filename, index_col=0)
            xp = DataFrame([[0.12, 0.23, 0.57],
                            [12.32, 123123.20, 321321.20]],
                           index=['A', 'B'], columns=['X', 'Y', 'Z'])
            assert_frame_equal(rs, xp)

    def test_to_csv_quoting(self):
        df = DataFrame({'A': [1, 2, 3], 'B': ['foo', 'bar', 'baz']})

        buf = StringIO()
        df.to_csv(buf, index=False, quoting=csv.QUOTE_NONNUMERIC)

        result = buf.getvalue()
        expected = ('"A","B"\n'
                    '1,"foo"\n'
                    '2,"bar"\n'
                    '3,"baz"\n')

        self.assertEqual(result, expected)

        # quoting windows line terminators, presents with encoding?
        # #3503
        text = 'a,b,c\n1,"test \r\n",3\n'
        df = pd.read_csv(StringIO(text))
        buf = StringIO()
        df.to_csv(buf, encoding='utf-8', index=False)
        self.assertEqual(buf.getvalue(), text)

    def test_to_csv_unicodewriter_quoting(self):
        df = DataFrame({'A': [1, 2, 3], 'B': ['foo', 'bar', 'baz']})

        buf = StringIO()
        df.to_csv(buf, index=False, quoting=csv.QUOTE_NONNUMERIC,
                  encoding='utf-8')

        result = buf.getvalue()
        expected = ('"A","B"\n'
                    '1,"foo"\n'
                    '2,"bar"\n'
                    '3,"baz"\n')

        self.assertEqual(result, expected)

    def test_to_csv_quote_none(self):
        # GH4328
        df = DataFrame({'A': ['hello', '{"hello"}']})
        for encoding in (None, 'utf-8'):
            buf = StringIO()
            df.to_csv(buf, quoting=csv.QUOTE_NONE,
                      encoding=encoding, index=False)
            result = buf.getvalue()
            expected = 'A\nhello\n{"hello"}\n'
            self.assertEqual(result, expected)

    def test_to_csv_index_no_leading_comma(self):
        df = DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]},
                       index=['one', 'two', 'three'])

        buf = StringIO()
        df.to_csv(buf, index_label=False)
        expected = ('A,B\n'
                    'one,1,4\n'
                    'two,2,5\n'
                    'three,3,6\n')
        self.assertEqual(buf.getvalue(), expected)

    def test_to_csv_line_terminators(self):
        df = DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]},
                       index=['one', 'two', 'three'])

        buf = StringIO()
        df.to_csv(buf, line_terminator='\r\n')
        expected = (',A,B\r\n'
                    'one,1,4\r\n'
                    'two,2,5\r\n'
                    'three,3,6\r\n')
        self.assertEqual(buf.getvalue(), expected)

        buf = StringIO()
        df.to_csv(buf)  # The default line terminator remains \n
        expected = (',A,B\n'
                    'one,1,4\n'
                    'two,2,5\n'
                    'three,3,6\n')
        self.assertEqual(buf.getvalue(), expected)

    def test_to_csv_from_csv_categorical(self):

        # CSV with categoricals should result in the same output as when one would add a "normal"
        # Series/DataFrame.
        s = Series(pd.Categorical(['a', 'b', 'b', 'a', 'a', 'c', 'c', 'c']))
        s2 = Series(['a', 'b', 'b', 'a', 'a', 'c', 'c', 'c'])
        res = StringIO()
        s.to_csv(res)
        exp = StringIO()
        s2.to_csv(exp)
        self.assertEqual(res.getvalue(), exp.getvalue())

        df = DataFrame({"s":s})
        df2 = DataFrame({"s":s2})
        res = StringIO()
        df.to_csv(res)
        exp = StringIO()
        df2.to_csv(exp)
        self.assertEqual(res.getvalue(), exp.getvalue())

    def test_to_csv_path_is_none(self):
        # GH 8215
        # Make sure we return string for consistency with
        # Series.to_csv()
        csv_str = self.frame.to_csv(path=None)
        self.assertIsInstance(csv_str, str)
        recons = pd.read_csv(StringIO(csv_str), index_col=0)
        assert_frame_equal(self.frame, recons)

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
        exp_size = len(dtypes) * n * 8  # cols * rows * bytes
        self.assertEqual(df_size, exp_size)
        # Ensure number of cols in memory_usage is the same as df
        size_df = np.size(df.columns.values)  # index=False; default
        self.assertEqual(size_df, np.size(df.memory_usage()))

        # test for validity
        DataFrame(1,index=['a'],columns=['A']).memory_usage(index=True)
        DataFrame(1,index=['a'],columns=['A']).index.nbytes
        DataFrame(1,index=pd.MultiIndex.from_product([['a'],range(1000)]),columns=['A']).index.nbytes
        DataFrame(1,index=pd.MultiIndex.from_product([['a'],range(1000)]),columns=['A']).index.values.nbytes
        DataFrame(1,index=pd.MultiIndex.from_product([['a'],range(1000)]),columns=['A']).memory_usage(index=True)
        DataFrame(1,index=pd.MultiIndex.from_product([['a'],range(1000)]),columns=['A']).index.nbytes
        DataFrame(1,index=pd.MultiIndex.from_product([['a'],range(1000)]),columns=['A']).index.values.nbytes

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

    def _check_method(self, method='pearson', check_minp=False):
        if not check_minp:
            correls = self.frame.corr(method=method)
            exp = self.frame['A'].corr(self.frame['C'], method=method)
            assert_almost_equal(correls['A']['C'], exp)
        else:
            result = self.frame.corr(min_periods=len(self.frame) - 8)
            expected = self.frame.corr()
            expected.ix['A', 'B'] = expected.ix['B', 'A'] = nan

    def test_corr_pearson(self):
        tm._skip_if_no_scipy()
        self.frame['A'][:5] = nan
        self.frame['B'][5:10] = nan

        self._check_method('pearson')

    def test_corr_kendall(self):
        tm._skip_if_no_scipy()
        self.frame['A'][:5] = nan
        self.frame['B'][5:10] = nan

        self._check_method('kendall')

    def test_corr_spearman(self):
        tm._skip_if_no_scipy()
        self.frame['A'][:5] = nan
        self.frame['B'][5:10] = nan

        self._check_method('spearman')

    def test_corr_non_numeric(self):
        tm._skip_if_no_scipy()
        self.frame['A'][:5] = nan
        self.frame['B'][5:10] = nan

        # exclude non-numeric types
        result = self.mixed_frame.corr()
        expected = self.mixed_frame.ix[:, ['A', 'B', 'C', 'D']].corr()
        assert_frame_equal(result, expected)

    def test_corr_nooverlap(self):
        tm._skip_if_no_scipy()

        # nothing in common
        for meth in ['pearson', 'kendall', 'spearman']:
            df = DataFrame({'A': [1, 1.5, 1, np.nan, np.nan, np.nan],
                            'B': [np.nan, np.nan, np.nan, 1, 1.5, 1]})
            rs = df.corr(meth)
            self.assertTrue(isnull(rs.ix['A', 'B']))
            self.assertTrue(isnull(rs.ix['B', 'A']))
            self.assertEqual(rs.ix['A', 'A'], 1)
            self.assertEqual(rs.ix['B', 'B'], 1)

    def test_corr_constant(self):
        tm._skip_if_no_scipy()

        # constant --> all NA

        for meth in ['pearson', 'spearman']:
            df = DataFrame({'A': [1, 1, 1, np.nan, np.nan, np.nan],
                            'B': [np.nan, np.nan, np.nan, 1, 1, 1]})
            rs = df.corr(meth)
            self.assertTrue(isnull(rs.values).all())

    def test_corr_int(self):
        # dtypes other than float64 #1761
        df3 = DataFrame({"a": [1, 2, 3, 4], "b": [1, 2, 3, 4]})

        # it works!
        df3.cov()
        df3.corr()

    def test_cov(self):
        # min_periods no NAs (corner case)
        expected = self.frame.cov()
        result = self.frame.cov(min_periods=len(self.frame))

        assert_frame_equal(expected, result)

        result = self.frame.cov(min_periods=len(self.frame) + 1)
        self.assertTrue(isnull(result.values).all())

        # with NAs
        frame = self.frame.copy()
        frame['A'][:5] = nan
        frame['B'][5:10] = nan
        result = self.frame.cov(min_periods=len(self.frame) - 8)
        expected = self.frame.cov()
        expected.ix['A', 'B'] = np.nan
        expected.ix['B', 'A'] = np.nan

        # regular
        self.frame['A'][:5] = nan
        self.frame['B'][:10] = nan
        cov = self.frame.cov()

        assert_almost_equal(cov['A']['C'],
                            self.frame['A'].cov(self.frame['C']))

        # exclude non-numeric types
        result = self.mixed_frame.cov()
        expected = self.mixed_frame.ix[:, ['A', 'B', 'C', 'D']].cov()
        assert_frame_equal(result, expected)

        # Single column frame
        df = DataFrame(np.linspace(0.0,1.0,10))
        result = df.cov()
        expected = DataFrame(np.cov(df.values.T).reshape((1,1)),
                             index=df.columns,columns=df.columns)
        assert_frame_equal(result, expected)
        df.ix[0] = np.nan
        result = df.cov()
        expected = DataFrame(np.cov(df.values[1:].T).reshape((1,1)),
                             index=df.columns,columns=df.columns)
        assert_frame_equal(result, expected)

    def test_corrwith(self):
        a = self.tsframe
        noise = Series(randn(len(a)), index=a.index)

        b = self.tsframe + noise

        # make sure order does not matter
        b = b.reindex(columns=b.columns[::-1], index=b.index[::-1][10:])
        del b['B']

        colcorr = a.corrwith(b, axis=0)
        assert_almost_equal(colcorr['A'], a['A'].corr(b['A']))

        rowcorr = a.corrwith(b, axis=1)
        assert_series_equal(rowcorr, a.T.corrwith(b.T, axis=0))

        dropped = a.corrwith(b, axis=0, drop=True)
        assert_almost_equal(dropped['A'], a['A'].corr(b['A']))
        self.assertNotIn('B', dropped)

        dropped = a.corrwith(b, axis=1, drop=True)
        self.assertNotIn(a.index[-1], dropped.index)

        # non time-series data
        index = ['a', 'b', 'c', 'd', 'e']
        columns = ['one', 'two', 'three', 'four']
        df1 = DataFrame(randn(5, 4), index=index, columns=columns)
        df2 = DataFrame(randn(4, 4), index=index[:4], columns=columns)
        correls = df1.corrwith(df2, axis=1)
        for row in index[:4]:
            assert_almost_equal(correls[row], df1.ix[row].corr(df2.ix[row]))

    def test_corrwith_with_objects(self):
        df1 = tm.makeTimeDataFrame()
        df2 = tm.makeTimeDataFrame()
        cols = ['A', 'B', 'C', 'D']

        df1['obj'] = 'foo'
        df2['obj'] = 'bar'

        result = df1.corrwith(df2)
        expected = df1.ix[:, cols].corrwith(df2.ix[:, cols])
        assert_series_equal(result, expected)

        result = df1.corrwith(df2, axis=1)
        expected = df1.ix[:, cols].corrwith(df2.ix[:, cols], axis=1)
        assert_series_equal(result, expected)

    def test_corrwith_series(self):
        result = self.tsframe.corrwith(self.tsframe['A'])
        expected = self.tsframe.apply(self.tsframe['A'].corr)

        assert_series_equal(result, expected)

    def test_corrwith_matches_corrcoef(self):
        df1 = DataFrame(np.arange(10000), columns=['a'])
        df2 = DataFrame(np.arange(10000)**2, columns=['a'])
        c1 = df1.corrwith(df2)['a']
        c2 = np.corrcoef(df1['a'],df2['a'])[0][1]

        assert_almost_equal(c1, c2)
        self.assertTrue(c1 < 1)

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

    def test_replace_inplace(self):
        self.tsframe['A'][:5] = nan
        self.tsframe['A'][-5:] = nan

        tsframe = self.tsframe.copy()
        tsframe.replace(nan, 0, inplace=True)
        assert_frame_equal(tsframe, self.tsframe.fillna(0))

        self.assertRaises(TypeError, self.tsframe.replace, nan, inplace=True)
        self.assertRaises(TypeError, self.tsframe.replace, nan)

        # mixed type
        self.mixed_frame.ix[5:20,'foo'] = nan
        self.mixed_frame.ix[-10:,'A'] = nan

        result = self.mixed_frame.replace(np.nan, 0)
        expected = self.mixed_frame.fillna(value=0)
        assert_frame_equal(result, expected)

        tsframe = self.tsframe.copy()
        tsframe.replace([nan], [0], inplace=True)
        assert_frame_equal(tsframe, self.tsframe.fillna(0))

    def test_regex_replace_scalar(self):
        obj = {'a': list('ab..'), 'b': list('efgh')}
        dfobj = DataFrame(obj)
        mix = {'a': lrange(4), 'b': list('ab..')}
        dfmix = DataFrame(mix)

        ### simplest cases
        ## regex -> value
        # obj frame
        res = dfobj.replace(r'\s*\.\s*', nan, regex=True)
        assert_frame_equal(dfobj, res.fillna('.'))

        # mixed
        res = dfmix.replace(r'\s*\.\s*', nan, regex=True)
        assert_frame_equal(dfmix, res.fillna('.'))

        ## regex -> regex
        # obj frame
        res = dfobj.replace(r'\s*(\.)\s*', r'\1\1\1', regex=True)
        objc = obj.copy()
        objc['a'] = ['a', 'b', '...', '...']
        expec = DataFrame(objc)
        assert_frame_equal(res, expec)

        # with mixed
        res = dfmix.replace(r'\s*(\.)\s*', r'\1\1\1', regex=True)
        mixc = mix.copy()
        mixc['b'] = ['a', 'b', '...', '...']
        expec = DataFrame(mixc)
        assert_frame_equal(res, expec)

        # everything with compiled regexs as well
        res = dfobj.replace(re.compile(r'\s*\.\s*'), nan, regex=True)
        assert_frame_equal(dfobj, res.fillna('.'))

        # mixed
        res = dfmix.replace(re.compile(r'\s*\.\s*'), nan, regex=True)
        assert_frame_equal(dfmix, res.fillna('.'))

        ## regex -> regex
        # obj frame
        res = dfobj.replace(re.compile(r'\s*(\.)\s*'), r'\1\1\1')
        objc = obj.copy()
        objc['a'] = ['a', 'b', '...', '...']
        expec = DataFrame(objc)
        assert_frame_equal(res, expec)

        # with mixed
        res = dfmix.replace(re.compile(r'\s*(\.)\s*'), r'\1\1\1')
        mixc = mix.copy()
        mixc['b'] = ['a', 'b', '...', '...']
        expec = DataFrame(mixc)
        assert_frame_equal(res, expec)

        res = dfmix.replace(regex=re.compile(r'\s*(\.)\s*'), value=r'\1\1\1')
        mixc = mix.copy()
        mixc['b'] = ['a', 'b', '...', '...']
        expec = DataFrame(mixc)
        assert_frame_equal(res, expec)

        res = dfmix.replace(regex=r'\s*(\.)\s*', value=r'\1\1\1')
        mixc = mix.copy()
        mixc['b'] = ['a', 'b', '...', '...']
        expec = DataFrame(mixc)
        assert_frame_equal(res, expec)

    def test_regex_replace_scalar_inplace(self):
        obj = {'a': list('ab..'), 'b': list('efgh')}
        dfobj = DataFrame(obj)
        mix = {'a': lrange(4), 'b': list('ab..')}
        dfmix = DataFrame(mix)

        ### simplest cases
        ## regex -> value
        # obj frame
        res = dfobj.copy()
        res.replace(r'\s*\.\s*', nan, regex=True, inplace=True)
        assert_frame_equal(dfobj, res.fillna('.'))

        # mixed
        res = dfmix.copy()
        res.replace(r'\s*\.\s*', nan, regex=True, inplace=True)
        assert_frame_equal(dfmix, res.fillna('.'))

        ## regex -> regex
        # obj frame
        res = dfobj.copy()
        res.replace(r'\s*(\.)\s*', r'\1\1\1', regex=True, inplace=True)
        objc = obj.copy()
        objc['a'] = ['a', 'b', '...', '...']
        expec = DataFrame(objc)
        assert_frame_equal(res, expec)

        # with mixed
        res = dfmix.copy()
        res.replace(r'\s*(\.)\s*', r'\1\1\1', regex=True, inplace=True)
        mixc = mix.copy()
        mixc['b'] = ['a', 'b', '...', '...']
        expec = DataFrame(mixc)
        assert_frame_equal(res, expec)

        # everything with compiled regexs as well
        res = dfobj.copy()
        res.replace(re.compile(r'\s*\.\s*'), nan, regex=True, inplace=True)
        assert_frame_equal(dfobj, res.fillna('.'))

        # mixed
        res = dfmix.copy()
        res.replace(re.compile(r'\s*\.\s*'), nan, regex=True, inplace=True)
        assert_frame_equal(dfmix, res.fillna('.'))

        ## regex -> regex
        # obj frame
        res = dfobj.copy()
        res.replace(re.compile(r'\s*(\.)\s*'), r'\1\1\1', regex=True,
                    inplace=True)
        objc = obj.copy()
        objc['a'] = ['a', 'b', '...', '...']
        expec = DataFrame(objc)
        assert_frame_equal(res, expec)

        # with mixed
        res = dfmix.copy()
        res.replace(re.compile(r'\s*(\.)\s*'), r'\1\1\1', regex=True,
                    inplace=True)
        mixc = mix.copy()
        mixc['b'] = ['a', 'b', '...', '...']
        expec = DataFrame(mixc)
        assert_frame_equal(res, expec)

        res = dfobj.copy()
        res.replace(regex=r'\s*\.\s*', value=nan, inplace=True)
        assert_frame_equal(dfobj, res.fillna('.'))

        # mixed
        res = dfmix.copy()
        res.replace(regex=r'\s*\.\s*', value=nan, inplace=True)
        assert_frame_equal(dfmix, res.fillna('.'))

        ## regex -> regex
        # obj frame
        res = dfobj.copy()
        res.replace(regex=r'\s*(\.)\s*', value=r'\1\1\1', inplace=True)
        objc = obj.copy()
        objc['a'] = ['a', 'b', '...', '...']
        expec = DataFrame(objc)
        assert_frame_equal(res, expec)

        # with mixed
        res = dfmix.copy()
        res.replace(regex=r'\s*(\.)\s*', value=r'\1\1\1', inplace=True)
        mixc = mix.copy()
        mixc['b'] = ['a', 'b', '...', '...']
        expec = DataFrame(mixc)
        assert_frame_equal(res, expec)

        # everything with compiled regexs as well
        res = dfobj.copy()
        res.replace(regex=re.compile(r'\s*\.\s*'), value=nan, inplace=True)
        assert_frame_equal(dfobj, res.fillna('.'))

        # mixed
        res = dfmix.copy()
        res.replace(regex=re.compile(r'\s*\.\s*'), value=nan, inplace=True)
        assert_frame_equal(dfmix, res.fillna('.'))

        ## regex -> regex
        # obj frame
        res = dfobj.copy()
        res.replace(regex=re.compile(r'\s*(\.)\s*'), value=r'\1\1\1',
                    inplace=True)
        objc = obj.copy()
        objc['a'] = ['a', 'b', '...', '...']
        expec = DataFrame(objc)
        assert_frame_equal(res, expec)

        # with mixed
        res = dfmix.copy()
        res.replace(regex=re.compile(r'\s*(\.)\s*'), value=r'\1\1\1',
                    inplace=True)
        mixc = mix.copy()
        mixc['b'] = ['a', 'b', '...', '...']
        expec = DataFrame(mixc)
        assert_frame_equal(res, expec)

    def test_regex_replace_list_obj(self):
        obj = {'a': list('ab..'), 'b': list('efgh'), 'c': list('helo')}
        dfobj = DataFrame(obj)

        ## lists of regexes and values
        # list of [re1, re2, ..., reN] -> [v1, v2, ..., vN]
        to_replace_res = [r'\s*\.\s*', r'e|f|g']
        values = [nan, 'crap']
        res = dfobj.replace(to_replace_res, values, regex=True)
        expec = DataFrame({'a': ['a', 'b', nan, nan], 'b': ['crap'] * 3 +
                           ['h'], 'c': ['h', 'crap', 'l', 'o']})
        assert_frame_equal(res, expec)

        # list of [re1, re2, ..., reN] -> [re1, re2, .., reN]
        to_replace_res = [r'\s*(\.)\s*', r'(e|f|g)']
        values = [r'\1\1', r'\1_crap']
        res = dfobj.replace(to_replace_res, values, regex=True)
        expec = DataFrame({'a': ['a', 'b', '..', '..'], 'b': ['e_crap',
                                                              'f_crap',
                                                              'g_crap', 'h'],
                           'c': ['h', 'e_crap', 'l', 'o']})

        assert_frame_equal(res, expec)

        # list of [re1, re2, ..., reN] -> [(re1 or v1), (re2 or v2), ..., (reN
        # or vN)]
        to_replace_res = [r'\s*(\.)\s*', r'e']
        values = [r'\1\1', r'crap']
        res = dfobj.replace(to_replace_res, values, regex=True)
        expec = DataFrame({'a': ['a', 'b', '..', '..'], 'b': ['crap', 'f', 'g',
                                                              'h'],
                           'c': ['h', 'crap', 'l', 'o']})
        assert_frame_equal(res, expec)

        to_replace_res = [r'\s*(\.)\s*', r'e']
        values = [r'\1\1', r'crap']
        res = dfobj.replace(value=values, regex=to_replace_res)
        expec = DataFrame({'a': ['a', 'b', '..', '..'], 'b': ['crap', 'f', 'g',
                                                              'h'],
                           'c': ['h', 'crap', 'l', 'o']})
        assert_frame_equal(res, expec)

    def test_regex_replace_list_obj_inplace(self):
        ### same as above with inplace=True
        ## lists of regexes and values
        obj = {'a': list('ab..'), 'b': list('efgh'), 'c': list('helo')}
        dfobj = DataFrame(obj)

        ## lists of regexes and values
        # list of [re1, re2, ..., reN] -> [v1, v2, ..., vN]
        to_replace_res = [r'\s*\.\s*', r'e|f|g']
        values = [nan, 'crap']
        res = dfobj.copy()
        res.replace(to_replace_res, values, inplace=True, regex=True)
        expec = DataFrame({'a': ['a', 'b', nan, nan], 'b': ['crap'] * 3 +
                           ['h'], 'c': ['h', 'crap', 'l', 'o']})
        assert_frame_equal(res, expec)

        # list of [re1, re2, ..., reN] -> [re1, re2, .., reN]
        to_replace_res = [r'\s*(\.)\s*', r'(e|f|g)']
        values = [r'\1\1', r'\1_crap']
        res = dfobj.copy()
        res.replace(to_replace_res, values, inplace=True, regex=True)
        expec = DataFrame({'a': ['a', 'b', '..', '..'], 'b': ['e_crap',
                                                              'f_crap',
                                                              'g_crap', 'h'],
                           'c': ['h', 'e_crap', 'l', 'o']})

        assert_frame_equal(res, expec)

        # list of [re1, re2, ..., reN] -> [(re1 or v1), (re2 or v2), ..., (reN
        # or vN)]
        to_replace_res = [r'\s*(\.)\s*', r'e']
        values = [r'\1\1', r'crap']
        res = dfobj.copy()
        res.replace(to_replace_res, values, inplace=True, regex=True)
        expec = DataFrame({'a': ['a', 'b', '..', '..'], 'b': ['crap', 'f', 'g',
                                                              'h'],
                           'c': ['h', 'crap', 'l', 'o']})
        assert_frame_equal(res, expec)

        to_replace_res = [r'\s*(\.)\s*', r'e']
        values = [r'\1\1', r'crap']
        res = dfobj.copy()
        res.replace(value=values, regex=to_replace_res, inplace=True)
        expec = DataFrame({'a': ['a', 'b', '..', '..'], 'b': ['crap', 'f', 'g',
                                                              'h'],
                           'c': ['h', 'crap', 'l', 'o']})
        assert_frame_equal(res, expec)

    def test_regex_replace_list_mixed(self):
        ## mixed frame to make sure this doesn't break things
        mix = {'a': lrange(4), 'b': list('ab..')}
        dfmix = DataFrame(mix)

        ## lists of regexes and values
        # list of [re1, re2, ..., reN] -> [v1, v2, ..., vN]
        to_replace_res = [r'\s*\.\s*', r'a']
        values = [nan, 'crap']
        mix2 = {'a': lrange(4), 'b': list('ab..'), 'c': list('halo')}
        dfmix2 = DataFrame(mix2)
        res = dfmix2.replace(to_replace_res, values, regex=True)
        expec = DataFrame({'a': mix2['a'], 'b': ['crap', 'b', nan, nan],
                           'c': ['h', 'crap', 'l', 'o']})
        assert_frame_equal(res, expec)

        # list of [re1, re2, ..., reN] -> [re1, re2, .., reN]
        to_replace_res = [r'\s*(\.)\s*', r'(a|b)']
        values = [r'\1\1', r'\1_crap']
        res = dfmix.replace(to_replace_res, values, regex=True)
        expec = DataFrame({'a': mix['a'], 'b': ['a_crap', 'b_crap', '..',
                                                '..']})

        assert_frame_equal(res, expec)

        # list of [re1, re2, ..., reN] -> [(re1 or v1), (re2 or v2), ..., (reN
        # or vN)]
        to_replace_res = [r'\s*(\.)\s*', r'a', r'(b)']
        values = [r'\1\1', r'crap', r'\1_crap']
        res = dfmix.replace(to_replace_res, values, regex=True)
        expec = DataFrame({'a': mix['a'], 'b': ['crap', 'b_crap', '..', '..']})
        assert_frame_equal(res, expec)

        to_replace_res = [r'\s*(\.)\s*', r'a', r'(b)']
        values = [r'\1\1', r'crap', r'\1_crap']
        res = dfmix.replace(regex=to_replace_res, value=values)
        expec = DataFrame({'a': mix['a'], 'b': ['crap', 'b_crap', '..', '..']})
        assert_frame_equal(res, expec)

    def test_regex_replace_list_mixed_inplace(self):
        mix = {'a': lrange(4), 'b': list('ab..')}
        dfmix = DataFrame(mix)
        # the same inplace
        ## lists of regexes and values
        # list of [re1, re2, ..., reN] -> [v1, v2, ..., vN]
        to_replace_res = [r'\s*\.\s*', r'a']
        values = [nan, 'crap']
        res = dfmix.copy()
        res.replace(to_replace_res, values, inplace=True, regex=True)
        expec = DataFrame({'a': mix['a'], 'b': ['crap', 'b', nan, nan]})
        assert_frame_equal(res, expec)

        # list of [re1, re2, ..., reN] -> [re1, re2, .., reN]
        to_replace_res = [r'\s*(\.)\s*', r'(a|b)']
        values = [r'\1\1', r'\1_crap']
        res = dfmix.copy()
        res.replace(to_replace_res, values, inplace=True, regex=True)
        expec = DataFrame({'a': mix['a'], 'b': ['a_crap', 'b_crap', '..',
                                                '..']})

        assert_frame_equal(res, expec)

        # list of [re1, re2, ..., reN] -> [(re1 or v1), (re2 or v2), ..., (reN
        # or vN)]
        to_replace_res = [r'\s*(\.)\s*', r'a', r'(b)']
        values = [r'\1\1', r'crap', r'\1_crap']
        res = dfmix.copy()
        res.replace(to_replace_res, values, inplace=True, regex=True)
        expec = DataFrame({'a': mix['a'], 'b': ['crap', 'b_crap', '..', '..']})
        assert_frame_equal(res, expec)

        to_replace_res = [r'\s*(\.)\s*', r'a', r'(b)']
        values = [r'\1\1', r'crap', r'\1_crap']
        res = dfmix.copy()
        res.replace(regex=to_replace_res, value=values, inplace=True)
        expec = DataFrame({'a': mix['a'], 'b': ['crap', 'b_crap', '..', '..']})
        assert_frame_equal(res, expec)

    def test_regex_replace_dict_mixed(self):
        mix = {'a': lrange(4), 'b': list('ab..'), 'c': ['a', 'b', nan, 'd']}
        dfmix = DataFrame(mix)

        ## dicts
        # single dict {re1: v1}, search the whole frame
        # need test for this...

        # list of dicts {re1: v1, re2: v2, ..., re3: v3}, search the whole
        # frame
        res = dfmix.replace({'b': r'\s*\.\s*'}, {'b': nan}, regex=True)
        res2 = dfmix.copy()
        res2.replace({'b': r'\s*\.\s*'}, {'b': nan}, inplace=True, regex=True)
        expec = DataFrame({'a': mix['a'], 'b': ['a', 'b', nan, nan], 'c':
                           mix['c']})
        assert_frame_equal(res, expec)
        assert_frame_equal(res2, expec)

        # list of dicts {re1: re11, re2: re12, ..., reN: re1N}, search the
        # whole frame
        res = dfmix.replace({'b': r'\s*(\.)\s*'}, {'b': r'\1ty'}, regex=True)
        res2 = dfmix.copy()
        res2.replace({'b': r'\s*(\.)\s*'}, {'b': r'\1ty'}, inplace=True,
                     regex=True)
        expec = DataFrame({'a': mix['a'], 'b': ['a', 'b', '.ty', '.ty'], 'c':
                           mix['c']})
        assert_frame_equal(res, expec)
        assert_frame_equal(res2, expec)

        res = dfmix.replace(regex={'b': r'\s*(\.)\s*'}, value={'b': r'\1ty'})
        res2 = dfmix.copy()
        res2.replace(regex={'b': r'\s*(\.)\s*'}, value={'b': r'\1ty'},
                     inplace=True)
        expec = DataFrame({'a': mix['a'], 'b': ['a', 'b', '.ty', '.ty'], 'c':
                           mix['c']})
        assert_frame_equal(res, expec)
        assert_frame_equal(res2, expec)

        # scalar -> dict
        # to_replace regex, {value: value}
        expec = DataFrame({'a': mix['a'], 'b': [nan, 'b', '.', '.'], 'c':
                           mix['c']})
        res = dfmix.replace('a', {'b': nan}, regex=True)
        res2 = dfmix.copy()
        res2.replace('a', {'b': nan}, regex=True, inplace=True)
        assert_frame_equal(res, expec)
        assert_frame_equal(res2, expec)

        res = dfmix.replace('a', {'b': nan}, regex=True)
        res2 = dfmix.copy()
        res2.replace(regex='a', value={'b': nan}, inplace=True)
        expec = DataFrame({'a': mix['a'], 'b': [nan, 'b', '.', '.'], 'c':
                           mix['c']})
        assert_frame_equal(res, expec)
        assert_frame_equal(res2, expec)

    def test_regex_replace_dict_nested(self):
        # nested dicts will not work until this is implemented for Series
        mix = {'a': lrange(4), 'b': list('ab..'), 'c': ['a', 'b', nan, 'd']}
        dfmix = DataFrame(mix)
        res = dfmix.replace({'b': {r'\s*\.\s*': nan}}, regex=True)
        res2 = dfmix.copy()
        res4 = dfmix.copy()
        res2.replace({'b': {r'\s*\.\s*': nan}}, inplace=True, regex=True)
        res3 = dfmix.replace(regex={'b': {r'\s*\.\s*': nan}})
        res4.replace(regex={'b': {r'\s*\.\s*': nan}}, inplace=True)
        expec = DataFrame({'a': mix['a'], 'b': ['a', 'b', nan, nan], 'c':
                           mix['c']})
        assert_frame_equal(res, expec)
        assert_frame_equal(res2, expec)
        assert_frame_equal(res3, expec)
        assert_frame_equal(res4, expec)

    def test_regex_replace_dict_nested_gh4115(self):
        df = pd.DataFrame({'Type':['Q','T','Q','Q','T'], 'tmp':2})
        expected = DataFrame({'Type': [0,1,0,0,1], 'tmp': 2})
        assert_frame_equal(df.replace({'Type': {'Q':0,'T':1}}), expected)

    def test_regex_replace_list_to_scalar(self):
        mix = {'a': lrange(4), 'b': list('ab..'), 'c': ['a', 'b', nan, 'd']}
        df = DataFrame(mix)
        res = df.replace([r'\s*\.\s*', 'a|b'], nan, regex=True)
        res2 = df.copy()
        res3 = df.copy()
        res2.replace([r'\s*\.\s*', 'a|b'], nan, regex=True, inplace=True)
        res3.replace(regex=[r'\s*\.\s*', 'a|b'], value=nan, inplace=True)
        expec = DataFrame({'a': mix['a'], 'b': np.array([nan] * 4),
                           'c': [nan, nan, nan, 'd']})
        assert_frame_equal(res, expec)
        assert_frame_equal(res2, expec)
        assert_frame_equal(res3, expec)

    def test_regex_replace_str_to_numeric(self):
        # what happens when you try to replace a numeric value with a regex?
        mix = {'a': lrange(4), 'b': list('ab..'), 'c': ['a', 'b', nan, 'd']}
        df = DataFrame(mix)
        res = df.replace(r'\s*\.\s*', 0, regex=True)
        res2 = df.copy()
        res2.replace(r'\s*\.\s*', 0, inplace=True, regex=True)
        res3 = df.copy()
        res3.replace(regex=r'\s*\.\s*', value=0, inplace=True)
        expec = DataFrame({'a': mix['a'], 'b': ['a', 'b', 0, 0], 'c':
                           mix['c']})
        assert_frame_equal(res, expec)
        assert_frame_equal(res2, expec)
        assert_frame_equal(res3, expec)

    def test_regex_replace_regex_list_to_numeric(self):
        mix = {'a': lrange(4), 'b': list('ab..'), 'c': ['a', 'b', nan, 'd']}
        df = DataFrame(mix)
        res = df.replace([r'\s*\.\s*', 'b'], 0, regex=True)
        res2 = df.copy()
        res2.replace([r'\s*\.\s*', 'b'], 0, regex=True, inplace=True)
        res3 = df.copy()
        res3.replace(regex=[r'\s*\.\s*', 'b'], value=0, inplace=True)
        expec = DataFrame({'a': mix['a'], 'b': ['a', 0, 0, 0], 'c': ['a', 0,
                                                                     nan,
                                                                     'd']})
        assert_frame_equal(res, expec)
        assert_frame_equal(res2, expec)
        assert_frame_equal(res3, expec)

    def test_regex_replace_series_of_regexes(self):
        mix = {'a': lrange(4), 'b': list('ab..'), 'c': ['a', 'b', nan, 'd']}
        df = DataFrame(mix)
        s1 = Series({'b': r'\s*\.\s*'})
        s2 = Series({'b': nan})
        res = df.replace(s1, s2, regex=True)
        res2 = df.copy()
        res2.replace(s1, s2, inplace=True, regex=True)
        res3 = df.copy()
        res3.replace(regex=s1, value=s2, inplace=True)
        expec = DataFrame({'a': mix['a'], 'b': ['a', 'b', nan, nan], 'c':
                           mix['c']})
        assert_frame_equal(res, expec)
        assert_frame_equal(res2, expec)
        assert_frame_equal(res3, expec)

    def test_regex_replace_numeric_to_object_conversion(self):
        mix = {'a': lrange(4), 'b': list('ab..'), 'c': ['a', 'b', nan, 'd']}
        df = DataFrame(mix)
        res = df.replace(0, 'a')
        expec = DataFrame({'a': ['a', 1, 2, 3], 'b': mix['b'], 'c': mix['c']})
        assert_frame_equal(res, expec)
        self.assertEqual(res.a.dtype, np.object_)

    def test_replace_regex_metachar(self):
        metachars = '[]', '()', '\d', '\w', '\s'

        for metachar in metachars:
            df = DataFrame({'a': [metachar, 'else']})
            result = df.replace({'a': {metachar: 'paren'}})
            expected = DataFrame({'a': ['paren', 'else']})
            tm.assert_frame_equal(result, expected)

    def test_replace(self):
        self.tsframe['A'][:5] = nan
        self.tsframe['A'][-5:] = nan

        zero_filled = self.tsframe.replace(nan, -1e8)
        assert_frame_equal(zero_filled, self.tsframe.fillna(-1e8))
        assert_frame_equal(zero_filled.replace(-1e8, nan), self.tsframe)

        self.tsframe['A'][:5] = nan
        self.tsframe['A'][-5:] = nan
        self.tsframe['B'][:5] = -1e8

        # empty
        df = DataFrame(index=['a', 'b'])
        assert_frame_equal(df, df.replace(5, 7))

    def test_replace_list(self):
        obj = {'a': list('ab..'), 'b': list('efgh'), 'c': list('helo')}
        dfobj = DataFrame(obj)

        ## lists of regexes and values
        # list of [v1, v2, ..., vN] -> [v1, v2, ..., vN]
        to_replace_res = [r'.', r'e']
        values = [nan, 'crap']
        res = dfobj.replace(to_replace_res, values)
        expec = DataFrame({'a': ['a', 'b', nan, nan],
                           'b': ['crap', 'f', 'g', 'h'], 'c': ['h', 'crap',
                                                               'l', 'o']})
        assert_frame_equal(res, expec)

        # list of [v1, v2, ..., vN] -> [v1, v2, .., vN]
        to_replace_res = [r'.', r'f']
        values = [r'..', r'crap']
        res = dfobj.replace(to_replace_res, values)
        expec = DataFrame({'a': ['a', 'b', '..', '..'], 'b': ['e', 'crap', 'g',
                                                              'h'],
                           'c': ['h', 'e', 'l', 'o']})

        assert_frame_equal(res, expec)

    def test_replace_series_dict(self):
        # from GH 3064
        df = DataFrame({'zero': {'a': 0.0, 'b': 1}, 'one': {'a': 2.0, 'b': 0}})
        result = df.replace(0, {'zero': 0.5, 'one': 1.0})
        expected = DataFrame({'zero': {'a': 0.5, 'b': 1}, 'one': {'a': 2.0, 'b': 1.0}})
        assert_frame_equal(result, expected)

        result = df.replace(0, df.mean())
        assert_frame_equal(result, expected)

        # series to series/dict
        df = DataFrame({'zero': {'a': 0.0, 'b': 1}, 'one': {'a': 2.0, 'b': 0}})
        s = Series({'zero': 0.0, 'one': 2.0})
        result = df.replace(s, {'zero': 0.5, 'one': 1.0})
        expected = DataFrame({'zero': {'a': 0.5, 'b': 1}, 'one': {'a': 1.0, 'b': 0.0}})
        assert_frame_equal(result, expected)

        result = df.replace(s, df.mean())
        assert_frame_equal(result, expected)

    def test_replace_convert(self):
        # gh 3907
        df = DataFrame([['foo', 'bar', 'bah'], ['bar', 'foo', 'bah']])
        m = {'foo': 1, 'bar': 2, 'bah': 3}
        rep = df.replace(m)
        expec = Series([ np.int64] * 3)
        res = rep.dtypes
        assert_series_equal(expec, res)

    def test_replace_mixed(self):
        self.mixed_frame.ix[5:20,'foo'] = nan
        self.mixed_frame.ix[-10:,'A'] = nan

        result = self.mixed_frame.replace(np.nan, -18)
        expected = self.mixed_frame.fillna(value=-18)
        assert_frame_equal(result, expected)
        assert_frame_equal(result.replace(-18, nan), self.mixed_frame)

        result = self.mixed_frame.replace(np.nan, -1e8)
        expected = self.mixed_frame.fillna(value=-1e8)
        assert_frame_equal(result, expected)
        assert_frame_equal(result.replace(-1e8, nan), self.mixed_frame)

        # int block upcasting
        df = DataFrame({ 'A' : Series([1.0,2.0],dtype='float64'), 'B' : Series([0,1],dtype='int64') })
        expected = DataFrame({ 'A' : Series([1.0,2.0],dtype='float64'), 'B' : Series([0.5,1],dtype='float64') })
        result = df.replace(0, 0.5)
        assert_frame_equal(result,expected)

        df.replace(0, 0.5, inplace=True)
        assert_frame_equal(df,expected)

        # int block splitting
        df = DataFrame({ 'A' : Series([1.0,2.0],dtype='float64'), 'B' : Series([0,1],dtype='int64'), 'C' : Series([1,2],dtype='int64') })
        expected = DataFrame({ 'A' : Series([1.0,2.0],dtype='float64'), 'B' : Series([0.5,1],dtype='float64'), 'C' : Series([1,2],dtype='int64') })
        result = df.replace(0, 0.5)
        assert_frame_equal(result,expected)

        # to object block upcasting
        df = DataFrame({ 'A' : Series([1.0,2.0],dtype='float64'), 'B' : Series([0,1],dtype='int64') })
        expected = DataFrame({ 'A' : Series([1,'foo'],dtype='object'), 'B' : Series([0,1],dtype='int64') })
        result = df.replace(2, 'foo')
        assert_frame_equal(result,expected)

        expected = DataFrame({ 'A' : Series(['foo','bar'],dtype='object'), 'B' : Series([0,'foo'],dtype='object') })
        result = df.replace([1,2], ['foo','bar'])
        assert_frame_equal(result,expected)

        # test case from
        from pandas.util.testing import makeCustomDataframe as mkdf
        df = DataFrame({'A' : Series([3,0],dtype='int64'), 'B' : Series([0,3],dtype='int64') })
        result = df.replace(3, df.mean().to_dict())
        expected = df.copy().astype('float64')
        m = df.mean()
        expected.iloc[0,0] = m[0]
        expected.iloc[1,1] = m[1]
        assert_frame_equal(result,expected)

    def test_replace_simple_nested_dict(self):
        df = DataFrame({'col': range(1, 5)})
        expected = DataFrame({'col': ['a', 2, 3, 'b']})

        result = df.replace({'col': {1: 'a', 4: 'b'}})
        tm.assert_frame_equal(expected, result)

        # in this case, should be the same as the not nested version
        result = df.replace({1: 'a', 4: 'b'})
        tm.assert_frame_equal(expected, result)

    def test_replace_simple_nested_dict_with_nonexistent_value(self):
        df = DataFrame({'col': range(1, 5)})
        expected = DataFrame({'col': ['a', 2, 3, 'b']})

        result = df.replace({-1: '-', 1: 'a', 4: 'b'})
        tm.assert_frame_equal(expected, result)

        result = df.replace({'col': {-1: '-', 1: 'a', 4: 'b'}})
        tm.assert_frame_equal(expected, result)

    def test_interpolate(self):
        pass

    def test_replace_value_is_none(self):
        self.assertRaises(TypeError, self.tsframe.replace, nan)
        orig_value = self.tsframe.iloc[0, 0]
        orig2 = self.tsframe.iloc[1, 0]

        self.tsframe.iloc[0, 0] = nan
        self.tsframe.iloc[1, 0] = 1

        result = self.tsframe.replace(to_replace={nan: 0})
        expected = self.tsframe.T.replace(to_replace={nan: 0}).T
        assert_frame_equal(result, expected)

        result = self.tsframe.replace(to_replace={nan: 0, 1: -1e8})
        tsframe = self.tsframe.copy()
        tsframe.iloc[0, 0] = 0
        tsframe.iloc[1, 0] = -1e8
        expected = tsframe
        assert_frame_equal(expected, result)
        self.tsframe.iloc[0, 0] = orig_value
        self.tsframe.iloc[1, 0] = orig2

    def test_replace_for_new_dtypes(self):

        # dtypes
        tsframe = self.tsframe.copy().astype(np.float32)
        tsframe['A'][:5] = nan
        tsframe['A'][-5:] = nan

        zero_filled = tsframe.replace(nan, -1e8)
        assert_frame_equal(zero_filled, tsframe.fillna(-1e8))
        assert_frame_equal(zero_filled.replace(-1e8, nan), tsframe)

        tsframe['A'][:5] = nan
        tsframe['A'][-5:] = nan
        tsframe['B'][:5] = -1e8

        b = tsframe['B']
        b[b == -1e8] = nan
        tsframe['B'] = b
        result = tsframe.fillna(method='bfill')
        assert_frame_equal(result, tsframe.fillna(method='bfill'))

    def test_replace_dtypes(self):
        # int
        df = DataFrame({'ints': [1, 2, 3]})
        result = df.replace(1, 0)
        expected = DataFrame({'ints': [0, 2, 3]})
        assert_frame_equal(result, expected)

        df = DataFrame({'ints': [1, 2, 3]}, dtype=np.int32)
        result = df.replace(1, 0)
        expected = DataFrame({'ints': [0, 2, 3]}, dtype=np.int32)
        assert_frame_equal(result, expected)

        df = DataFrame({'ints': [1, 2, 3]}, dtype=np.int16)
        result = df.replace(1, 0)
        expected = DataFrame({'ints': [0, 2, 3]}, dtype=np.int16)
        assert_frame_equal(result, expected)

        # bools
        df = DataFrame({'bools': [True, False, True]})
        result = df.replace(False, True)
        self.assertTrue(result.values.all())

        # complex blocks
        df = DataFrame({'complex': [1j, 2j, 3j]})
        result = df.replace(1j, 0j)
        expected = DataFrame({'complex': [0j, 2j, 3j]})
        assert_frame_equal(result, expected)

        # datetime blocks
        prev = datetime.today()
        now = datetime.today()
        df = DataFrame({'datetime64': Index([prev, now, prev])})
        result = df.replace(prev, now)
        expected = DataFrame({'datetime64': Index([now] * 3)})
        assert_frame_equal(result, expected)

    def test_replace_input_formats(self):
        # both dicts
        to_rep = {'A': np.nan, 'B': 0, 'C': ''}
        values = {'A': 0, 'B': -1, 'C': 'missing'}
        df = DataFrame({'A': [np.nan, 0, np.inf], 'B': [0, 2, 5],
                        'C': ['', 'asdf', 'fd']})
        filled = df.replace(to_rep, values)
        expected = {}
        for k, v in compat.iteritems(df):
            expected[k] = v.replace(to_rep[k], values[k])
        assert_frame_equal(filled, DataFrame(expected))

        result = df.replace([0, 2, 5], [5, 2, 0])
        expected = DataFrame({'A': [np.nan, 5, np.inf], 'B': [5, 2, 0],
                              'C': ['', 'asdf', 'fd']})
        assert_frame_equal(result, expected)

        # dict to scalar
        filled = df.replace(to_rep, 0)
        expected = {}
        for k, v in compat.iteritems(df):
            expected[k] = v.replace(to_rep[k], 0)
        assert_frame_equal(filled, DataFrame(expected))

        self.assertRaises(TypeError, df.replace, to_rep, [np.nan, 0, ''])

        # scalar to dict
        values = {'A': 0, 'B': -1, 'C': 'missing'}
        df = DataFrame({'A': [np.nan, 0, np.nan], 'B': [0, 2, 5],
                        'C': ['', 'asdf', 'fd']})
        filled = df.replace(np.nan, values)
        expected = {}
        for k, v in compat.iteritems(df):
            expected[k] = v.replace(np.nan, values[k])
        assert_frame_equal(filled, DataFrame(expected))

        # list to list
        to_rep = [np.nan, 0, '']
        values = [-2, -1, 'missing']
        result = df.replace(to_rep, values)
        expected = df.copy()
        for i in range(len(to_rep)):
            expected.replace(to_rep[i], values[i], inplace=True)
        assert_frame_equal(result, expected)

        self.assertRaises(ValueError, df.replace, to_rep, values[1:])

        # list to scalar
        to_rep = [np.nan, 0, '']
        result = df.replace(to_rep, -1)
        expected = df.copy()
        for i in range(len(to_rep)):
            expected.replace(to_rep[i], -1, inplace=True)
        assert_frame_equal(result, expected)

    def test_replace_limit(self):
        pass

    def test_replace_dict_no_regex(self):
        answer = Series({0: 'Strongly Agree', 1: 'Agree', 2: 'Neutral', 3:
                         'Disagree', 4: 'Strongly Disagree'})
        weights = {'Agree': 4, 'Disagree': 2, 'Neutral': 3, 'Strongly Agree':
                   5, 'Strongly Disagree': 1}
        expected = Series({0: 5, 1: 4, 2: 3, 3: 2, 4: 1})
        result = answer.replace(weights)
        tm.assert_series_equal(result, expected)

    def test_replace_series_no_regex(self):
        answer = Series({0: 'Strongly Agree', 1: 'Agree', 2: 'Neutral', 3:
                         'Disagree', 4: 'Strongly Disagree'})
        weights = Series({'Agree': 4, 'Disagree': 2, 'Neutral': 3,
                          'Strongly Agree': 5, 'Strongly Disagree': 1})
        expected = Series({0: 5, 1: 4, 2: 3, 3: 2, 4: 1})
        result = answer.replace(weights)
        tm.assert_series_equal(result, expected)

    def test_replace_dict_tuple_list_ordering_remains_the_same(self):
        df = DataFrame(dict(A=[nan, 1]))
        res1 = df.replace(to_replace={nan: 0, 1: -1e8})
        res2 = df.replace(to_replace=(1, nan), value=[-1e8, 0])
        res3 = df.replace(to_replace=[1, nan], value=[-1e8, 0])

        expected = DataFrame({'A': [0, -1e8]})
        tm.assert_frame_equal(res1, res2)
        tm.assert_frame_equal(res2, res3)
        tm.assert_frame_equal(res3, expected)

    def test_replace_doesnt_replace_without_regex(self):
        from pandas.compat import StringIO
        raw = """fol T_opp T_Dir T_Enh
        0    1     0     0    vo
        1    2    vr     0     0
        2    2     0     0     0
        3    3     0    bt     0"""
        df = read_csv(StringIO(raw), sep=r'\s+')
        res = df.replace({'\D': 1})
        tm.assert_frame_equal(df, res)

    def test_replace_bool_with_string(self):
        df = DataFrame({'a': [True, False], 'b': list('ab')})
        result = df.replace(True, 'a')
        expected = DataFrame({'a': ['a', False], 'b': df.b})
        tm.assert_frame_equal(result, expected)

    def test_replace_pure_bool_with_string_no_op(self):
        df = DataFrame(np.random.rand(2, 2) > 0.5)
        result = df.replace('asdf', 'fdsa')
        tm.assert_frame_equal(df, result)

    def test_replace_bool_with_bool(self):
        df = DataFrame(np.random.rand(2, 2) > 0.5)
        result = df.replace(False, True)
        expected = DataFrame(np.ones((2, 2), dtype=bool))
        tm.assert_frame_equal(result, expected)

    def test_replace_with_dict_with_bool_keys(self):
        df = DataFrame({0: [True, False], 1: [False, True]})
        with tm.assertRaisesRegexp(TypeError, 'Cannot compare types .+'):
            df.replace({'asdf': 'asdb', True: 'yes'})

    def test_replace_truthy(self):
        df = DataFrame({'a': [True, True]})
        r = df.replace([np.inf, -np.inf], np.nan)
        e = df
        tm.assert_frame_equal(r, e)

    def test_replace_int_to_int_chain(self):
        df = DataFrame({'a': lrange(1, 5)})
        with tm.assertRaisesRegexp(ValueError, "Replacement not allowed .+"):
            df.replace({'a': dict(zip(range(1, 5), range(2, 6)))})

    def test_replace_str_to_str_chain(self):
        a = np.arange(1, 5)
        astr = a.astype(str)
        bstr = np.arange(2, 6).astype(str)
        df = DataFrame({'a': astr})
        with tm.assertRaisesRegexp(ValueError, "Replacement not allowed .+"):
            df.replace({'a': dict(zip(astr, bstr))})

    def test_replace_swapping_bug(self):
        df = pd.DataFrame({'a': [True, False, True]})
        res = df.replace({'a': {True: 'Y', False: 'N'}})
        expect = pd.DataFrame({'a': ['Y', 'N', 'Y']})
        tm.assert_frame_equal(res, expect)

        df = pd.DataFrame({'a': [0, 1, 0]})
        res = df.replace({'a': {0: 'Y', 1: 'N'}})
        expect = pd.DataFrame({'a': ['Y', 'N', 'Y']})
        tm.assert_frame_equal(res, expect)

    def test_replace_period(self):
        d = {'fname':
             {'out_augmented_AUG_2011.json': pd.Period(year=2011, month=8, freq='M'),
              'out_augmented_JAN_2011.json': pd.Period(year=2011, month=1, freq='M'),
              'out_augmented_MAY_2012.json': pd.Period(year=2012, month=5, freq='M'),
              'out_augmented_SUBSIDY_WEEK.json': pd.Period(year=2011, month=4, freq='M'),
              'out_augmented_AUG_2012.json': pd.Period(year=2012, month=8, freq='M'),
              'out_augmented_MAY_2011.json': pd.Period(year=2011, month=5, freq='M'),
              'out_augmented_SEP_2013.json': pd.Period(year=2013, month=9, freq='M')}}

        df = pd.DataFrame(['out_augmented_AUG_2012.json',
                           'out_augmented_SEP_2013.json',
                           'out_augmented_SUBSIDY_WEEK.json',
                           'out_augmented_MAY_2012.json',
                           'out_augmented_MAY_2011.json',
                           'out_augmented_AUG_2011.json',
                           'out_augmented_JAN_2011.json'], columns=['fname'])
        tm.assert_equal(set(df.fname.values), set(d['fname'].keys()))
        expected = DataFrame({'fname': [d['fname'][k]
                                        for k in df.fname.values]})
        result = df.replace(d)
        tm.assert_frame_equal(result, expected)

    def test_replace_datetime(self):
        d = {'fname':
             {'out_augmented_AUG_2011.json': pd.Timestamp('2011-08'),
              'out_augmented_JAN_2011.json': pd.Timestamp('2011-01'),
              'out_augmented_MAY_2012.json': pd.Timestamp('2012-05'),
              'out_augmented_SUBSIDY_WEEK.json': pd.Timestamp('2011-04'),
              'out_augmented_AUG_2012.json': pd.Timestamp('2012-08'),
              'out_augmented_MAY_2011.json': pd.Timestamp('2011-05'),
              'out_augmented_SEP_2013.json': pd.Timestamp('2013-09')}}

        df = pd.DataFrame(['out_augmented_AUG_2012.json',
                           'out_augmented_SEP_2013.json',
                           'out_augmented_SUBSIDY_WEEK.json',
                           'out_augmented_MAY_2012.json',
                           'out_augmented_MAY_2011.json',
                           'out_augmented_AUG_2011.json',
                           'out_augmented_JAN_2011.json'], columns=['fname'])
        tm.assert_equal(set(df.fname.values), set(d['fname'].keys()))
        expected = DataFrame({'fname': [d['fname'][k]
                                        for k in df.fname.values]})
        result = df.replace(d)
        tm.assert_frame_equal(result, expected)

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

        # no columns but index
        df = DataFrame(index=['a', 'b', 'c'])
        result = df.xs('a')
        expected = Series([], name='a')
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

    def test_pivot(self):
        data = {
            'index': ['A', 'B', 'C', 'C', 'B', 'A'],
            'columns': ['One', 'One', 'One', 'Two', 'Two', 'Two'],
            'values': [1., 2., 3., 3., 2., 1.]
        }

        frame = DataFrame(data)
        pivoted = frame.pivot(
            index='index', columns='columns', values='values')

        expected = DataFrame({
            'One': {'A': 1., 'B': 2., 'C': 3.},
            'Two': {'A': 1., 'B': 2., 'C': 3.}
        })
        expected.index.name, expected.columns.name = 'index', 'columns'

        assert_frame_equal(pivoted, expected)

        # name tracking
        self.assertEqual(pivoted.index.name, 'index')
        self.assertEqual(pivoted.columns.name, 'columns')

        # don't specify values
        pivoted = frame.pivot(index='index', columns='columns')
        self.assertEqual(pivoted.index.name, 'index')
        self.assertEqual(pivoted.columns.names, (None, 'columns'))

        # pivot multiple columns
        wp = tm.makePanel()
        lp = wp.to_frame()
        df = lp.reset_index()
        assert_frame_equal(df.pivot('major', 'minor'), lp.unstack())

    def test_pivot_duplicates(self):
        data = DataFrame({'a': ['bar', 'bar', 'foo', 'foo', 'foo'],
                          'b': ['one', 'two', 'one', 'one', 'two'],
                          'c': [1., 2., 3., 3., 4.]})
        with assertRaisesRegexp(ValueError, 'duplicate entries'):
            data.pivot('a', 'b', 'c')

    def test_pivot_empty(self):
        df = DataFrame({}, columns=['a', 'b', 'c'])
        result = df.pivot('a', 'b', 'c')
        expected = DataFrame({})
        assert_frame_equal(result, expected, check_names=False)

    def test_pivot_integer_bug(self):
        df = DataFrame(data=[("A", "1", "A1"), ("B", "2", "B2")])

        result = df.pivot(index=1, columns=0, values=2)
        repr(result)
        self.assert_numpy_array_equal(result.columns, ['A', 'B'])

    def test_pivot_index_none(self):
        # gh-3962
        data = {
            'index': ['A', 'B', 'C', 'C', 'B', 'A'],
            'columns': ['One', 'One', 'One', 'Two', 'Two', 'Two'],
            'values': [1., 2., 3., 3., 2., 1.]
        }

        frame = DataFrame(data).set_index('index')
        result = frame.pivot(columns='columns', values='values')
        expected = DataFrame({
            'One': {'A': 1., 'B': 2., 'C': 3.},
            'Two': {'A': 1., 'B': 2., 'C': 3.}
        })

        expected.index.name, expected.columns.name = 'index', 'columns'
        assert_frame_equal(result, expected)

        # omit values
        result = frame.pivot(columns='columns')

        expected.columns = pd.MultiIndex.from_tuples([('values', 'One'),
                                                     ('values', 'Two')],
                                                     names=[None, 'columns'])
        expected.index.name = 'index'
        assert_frame_equal(result, expected, check_names=False)
        self.assertEqual(result.index.name, 'index',)
        self.assertEqual(result.columns.names, (None, 'columns'))
        expected.columns = expected.columns.droplevel(0)

        data = {
            'index': range(7),
            'columns': ['One', 'One', 'One', 'Two', 'Two', 'Two'],
            'values': [1., 2., 3., 3., 2., 1.]
        }

        result = frame.pivot(columns='columns', values='values')

        expected.columns.name = 'columns'
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
                index=[2, np.nan, 1, 5], columns=['joe', 'jim'])

        i, j = [np.nan, 5, 5, np.nan, 1, 2, np.nan], [1, 3, 3, 1, 2, 0, 1]
        tm.assert_frame_equal(df.reindex(i), df.iloc[j])

        df.index = df.index.astype('object')
        tm.assert_frame_equal(df.reindex(i), df.iloc[j])

        # GH10388
        df = pd.DataFrame({'other':['a', 'b', np.nan, 'c'],
                           'date':['2015-03-22', np.nan, '2012-01-08', np.nan],
                           'amount':[2, 3, 4, 5]})

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
        tm.assert_frame_equal(expl, res1l)
        tm.assert_frame_equal(expl, res2r)
        expr = pd.DataFrame([0, 0, 1, 1, np.nan, np.nan] * 2, index=midx)
        tm.assert_frame_equal(expr, res1r)
        tm.assert_frame_equal(expr, res2l)

        res1l, res1r = df1.align(df2, join='right')
        res2l, res2r = df2.align(df1, join='left')

        exp_idx = pd.MultiIndex.from_product([range(2), range(2), range(2)],
                                             names=('a', 'b', 'c'))
        expl = pd.DataFrame([0, 1, 2, 3, 6, 7, 8, 9], index=exp_idx)
        tm.assert_frame_equal(expl, res1l)
        tm.assert_frame_equal(expl, res2r)
        expr = pd.DataFrame([0, 0, 1, 1] * 2, index=exp_idx)
        tm.assert_frame_equal(expr, res1r)
        tm.assert_frame_equal(expr, res2l)

    def test_where(self):
        default_frame = DataFrame(np.random.randn(5, 3),columns=['A','B','C'])

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
        assert_series_equal(result, Series([]))

        empty_with_cols = DataFrame(columns=['a', 'b', 'c'])
        result = empty_with_cols.apply(x.append, axis=1, reduce=False)
        assert_frame_equal(result, empty_with_cols)
        result = empty_with_cols.apply(x.append, axis=1, reduce=True)
        assert_series_equal(result, Series([]))

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
        expected = Series(np.nan, index=[])
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
        expected = DataFrame(0., index=[0, 1, 2], columns=[1, 2])
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
        import random
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

    def test_combineAdd(self):

        with tm.assert_produces_warning(FutureWarning):
            # trivial
            comb = self.frame.combineAdd(self.frame)
            assert_frame_equal(comb, self.frame * 2)

            # more rigorous
            a = DataFrame([[1., nan, nan, 2., nan]],
                          columns=np.arange(5))
            b = DataFrame([[2., 3., nan, 2., 6., nan]],
                          columns=np.arange(6))
            expected = DataFrame([[3., 3., nan, 4., 6., nan]],
                                 columns=np.arange(6))

            result = a.combineAdd(b)
            assert_frame_equal(result, expected)
            result2 = a.T.combineAdd(b.T)
            assert_frame_equal(result2, expected.T)

            expected2 = a.combine(b, operator.add, fill_value=0.)
            assert_frame_equal(expected, expected2)

            # corner cases
            comb = self.frame.combineAdd(self.empty)
            assert_frame_equal(comb, self.frame)

            comb = self.empty.combineAdd(self.frame)
            assert_frame_equal(comb, self.frame)

            # integer corner case
            df1 = DataFrame({'x': [5]})
            df2 = DataFrame({'x': [1]})
            df3 = DataFrame({'x': [6]})
            comb = df1.combineAdd(df2)
            assert_frame_equal(comb, df3)

            # mixed type GH2191
            df1 = DataFrame({'A': [1, 2], 'B': [3, 4]})
            df2 = DataFrame({'A': [1, 2], 'C': [5, 6]})
            rs = df1.combineAdd(df2)
            xp = DataFrame({'A': [2, 4], 'B': [3, 4.], 'C': [5, 6.]})
            assert_frame_equal(xp, rs)

        # TODO: test integer fill corner?

    def test_combineMult(self):

        with tm.assert_produces_warning(FutureWarning):
            # trivial
            comb = self.frame.combineMult(self.frame)

            assert_frame_equal(comb, self.frame ** 2)

            # corner cases
            comb = self.frame.combineMult(self.empty)
            assert_frame_equal(comb, self.frame)

            comb = self.empty.combineMult(self.frame)
            assert_frame_equal(comb, self.frame)

    def test_combine_generic(self):
        df1 = self.frame
        df2 = self.frame.ix[:-5, ['A', 'B', 'C']]

        combined = df1.combine(df2, np.add)
        combined2 = df2.combine(df1, np.add)
        self.assertTrue(combined['D'].isnull().all())
        self.assertTrue(combined2['D'].isnull().all())

        chunk = combined.ix[:-5, ['A', 'B', 'C']]
        chunk2 = combined2.ix[:-5, ['A', 'B', 'C']]

        exp = self.frame.ix[:-5, ['A', 'B', 'C']].reindex_like(chunk) * 2
        assert_frame_equal(chunk, exp)
        assert_frame_equal(chunk2, exp)

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
            self.assertTrue((clipped_df.values[mask] == df.values[mask]).all() == True)

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

    def test_bool_describe_in_mixed_frame(self):
        df = DataFrame({
            'string_data': ['a', 'b', 'c', 'd', 'e'],
            'bool_data': [True, True, False, False, False],
            'int_data': [10, 20, 30, 40, 50],
        })

        # Boolean data and integer data is included in .describe() output, string data isn't
        self.assert_numpy_array_equal(df.describe().columns, ['bool_data', 'int_data'])

        bool_describe = df.describe()['bool_data']

        # Both the min and the max values should stay booleans
        self.assertEqual(bool_describe['min'].dtype, np.bool_)
        self.assertEqual(bool_describe['max'].dtype, np.bool_)

        self.assertFalse(bool_describe['min'])
        self.assertTrue(bool_describe['max'])

        # For numeric operations, like mean or median, the values True/False are cast to
        # the integer values 1 and 0
        assert_almost_equal(bool_describe['mean'], 0.4)
        assert_almost_equal(bool_describe['50%'], 0)

    def test_reduce_mixed_frame(self):
        # GH 6806
        df = DataFrame({
            'bool_data': [True, True, False, False, False],
            'int_data': [10, 20, 30, 40, 50],
            'string_data': ['a', 'b', 'c', 'd', 'e'],
        })
        df.reindex(columns=['bool_data', 'int_data', 'string_data'])
        test = df.sum(axis=0)
        assert_almost_equal(test.values, [2, 150, 'abcde'])
        assert_series_equal(test, df.T.sum(axis=1))

    def test_count(self):
        f = lambda s: notnull(s).sum()
        self._check_stat_op('count', f,
                            has_skipna=False,
                            has_numeric_only=True,
                            check_dtype=False,
                            check_dates=True)

        # corner case
        frame = DataFrame()
        ct1 = frame.count(1)
        tm.assertIsInstance(ct1, Series)

        ct2 = frame.count(0)
        tm.assertIsInstance(ct2, Series)

        # GH #423
        df = DataFrame(index=lrange(10))
        result = df.count(1)
        expected = Series(0, index=df.index)
        assert_series_equal(result, expected)

        df = DataFrame(columns=lrange(10))
        result = df.count(0)
        expected = Series(0, index=df.columns)
        assert_series_equal(result, expected)

        df = DataFrame()
        result = df.count()
        expected = Series(0, index=[])
        assert_series_equal(result, expected)

    def test_sum(self):
        self._check_stat_op('sum', np.sum, has_numeric_only=True)

        # mixed types (with upcasting happening)
        self._check_stat_op('sum', np.sum, frame=self.mixed_float.astype('float32'),
                            has_numeric_only=True, check_dtype=False, check_less_precise=True)

    def test_stat_operators_attempt_obj_array(self):
        data = {
            'a': [-0.00049987540199591344, -0.0016467257772919831,
                  0.00067695870775883013],
            'b': [-0, -0, 0.0],
            'c': [0.00031111847529610595, 0.0014902627951905339,
                  -0.00094099200035979691]
        }
        df1 = DataFrame(data, index=['foo', 'bar', 'baz'],
                        dtype='O')
        methods = ['sum', 'mean', 'prod', 'var', 'std', 'skew', 'min', 'max']

        # GH #676
        df2 = DataFrame({0: [np.nan, 2], 1: [np.nan, 3],
                        2: [np.nan, 4]}, dtype=object)

        for df in [df1, df2]:
            for meth in methods:
                self.assertEqual(df.values.dtype, np.object_)
                result = getattr(df, meth)(1)
                expected = getattr(df.astype('f8'), meth)(1)

                if not tm._incompat_bottleneck_version(meth):
                    assert_series_equal(result, expected)

    def test_mean(self):
        self._check_stat_op('mean', np.mean, check_dates=True)

    def test_product(self):
        self._check_stat_op('product', np.prod)

    def test_median(self):
        def wrapper(x):
            if isnull(x).any():
                return np.nan
            return np.median(x)

        self._check_stat_op('median', wrapper, check_dates=True)

    def test_min(self):
        self._check_stat_op('min', np.min, check_dates=True)
        self._check_stat_op('min', np.min, frame=self.intframe)

    def test_cummin(self):
        self.tsframe.ix[5:10, 0] = nan
        self.tsframe.ix[10:15, 1] = nan
        self.tsframe.ix[15:, 2] = nan

        # axis = 0
        cummin = self.tsframe.cummin()
        expected = self.tsframe.apply(Series.cummin)
        assert_frame_equal(cummin, expected)

        # axis = 1
        cummin = self.tsframe.cummin(axis=1)
        expected = self.tsframe.apply(Series.cummin, axis=1)
        assert_frame_equal(cummin, expected)

        # works
        df = DataFrame({'A': np.arange(20)}, index=np.arange(20))
        result = df.cummin()

        # fix issue
        cummin_xs = self.tsframe.cummin(axis=1)
        self.assertEqual(np.shape(cummin_xs), np.shape(self.tsframe))

    def test_cummax(self):
        self.tsframe.ix[5:10, 0] = nan
        self.tsframe.ix[10:15, 1] = nan
        self.tsframe.ix[15:, 2] = nan

        # axis = 0
        cummax = self.tsframe.cummax()
        expected = self.tsframe.apply(Series.cummax)
        assert_frame_equal(cummax, expected)

        # axis = 1
        cummax = self.tsframe.cummax(axis=1)
        expected = self.tsframe.apply(Series.cummax, axis=1)
        assert_frame_equal(cummax, expected)

        # works
        df = DataFrame({'A': np.arange(20)}, index=np.arange(20))
        result = df.cummax()

        # fix issue
        cummax_xs = self.tsframe.cummax(axis=1)
        self.assertEqual(np.shape(cummax_xs), np.shape(self.tsframe))

    def test_max(self):
        self._check_stat_op('max', np.max, check_dates=True)
        self._check_stat_op('max', np.max, frame=self.intframe)

    def test_mad(self):
        f = lambda x: np.abs(x - x.mean()).mean()
        self._check_stat_op('mad', f)

    def test_var_std(self):
        alt = lambda x: np.var(x, ddof=1)
        self._check_stat_op('var', alt)

        alt = lambda x: np.std(x, ddof=1)
        self._check_stat_op('std', alt)

        result = self.tsframe.std(ddof=4)
        expected = self.tsframe.apply(lambda x: x.std(ddof=4))
        assert_almost_equal(result, expected)

        result = self.tsframe.var(ddof=4)
        expected = self.tsframe.apply(lambda x: x.var(ddof=4))
        assert_almost_equal(result, expected)

        arr = np.repeat(np.random.random((1, 1000)), 1000, 0)
        result = nanops.nanvar(arr, axis=0)
        self.assertFalse((result < 0).any())
        if nanops._USE_BOTTLENECK:
            nanops._USE_BOTTLENECK = False
            result = nanops.nanvar(arr, axis=0)
            self.assertFalse((result < 0).any())
            nanops._USE_BOTTLENECK = True

    def test_numeric_only_flag(self):
        # GH #9201
        methods = ['sem', 'var', 'std']
        df1 = DataFrame(np.random.randn(5, 3), columns=['foo', 'bar', 'baz'])
        # set one entry to a number in str format
        df1.ix[0, 'foo'] = '100'

        df2 = DataFrame(np.random.randn(5, 3), columns=['foo', 'bar', 'baz'])
        # set one entry to a non-number str
        df2.ix[0, 'foo'] = 'a'

        for meth in methods:
            result = getattr(df1, meth)(axis=1, numeric_only=True)
            expected = getattr(df1[['bar', 'baz']], meth)(axis=1)
            assert_series_equal(expected, result)

            result = getattr(df2, meth)(axis=1, numeric_only=True)
            expected = getattr(df2[['bar', 'baz']], meth)(axis=1)
            assert_series_equal(expected, result)

            # df1 has all numbers, df2 has a letter inside
            self.assertRaises(TypeError, lambda : getattr(df1, meth)(axis=1, numeric_only=False))
            self.assertRaises(TypeError, lambda : getattr(df2, meth)(axis=1, numeric_only=False))

    def test_sem(self):
        alt = lambda x: np.std(x, ddof=1)/np.sqrt(len(x))
        self._check_stat_op('sem', alt)

        result = self.tsframe.sem(ddof=4)
        expected = self.tsframe.apply(lambda x: x.std(ddof=4)/np.sqrt(len(x)))
        assert_almost_equal(result, expected)

        arr = np.repeat(np.random.random((1, 1000)), 1000, 0)
        result = nanops.nansem(arr, axis=0)
        self.assertFalse((result < 0).any())
        if nanops._USE_BOTTLENECK:
            nanops._USE_BOTTLENECK = False
            result = nanops.nansem(arr, axis=0)
            self.assertFalse((result < 0).any())
            nanops._USE_BOTTLENECK = True

    def test_skew(self):
        tm._skip_if_no_scipy()
        from scipy.stats import skew

        def alt(x):
            if len(x) < 3:
                return np.nan
            return skew(x, bias=False)

        self._check_stat_op('skew', alt)

    def test_kurt(self):
        tm._skip_if_no_scipy()

        from scipy.stats import kurtosis

        def alt(x):
            if len(x) < 4:
                return np.nan
            return kurtosis(x, bias=False)

        self._check_stat_op('kurt', alt)

        index = MultiIndex(levels=[['bar'], ['one', 'two', 'three'], [0, 1]],
                           labels=[[0, 0, 0, 0, 0, 0],
                                   [0, 1, 2, 0, 1, 2],
                                   [0, 1, 0, 1, 0, 1]])
        df = DataFrame(np.random.randn(6, 3), index=index)

        kurt = df.kurt()
        kurt2 = df.kurt(level=0).xs('bar')
        assert_series_equal(kurt, kurt2, check_names=False)
        self.assertTrue(kurt.name is None)
        self.assertEqual(kurt2.name, 'bar')

    def _check_stat_op(self, name, alternative, frame=None, has_skipna=True,
                       has_numeric_only=False, check_dtype=True, check_dates=False,
                       check_less_precise=False):
        if frame is None:
            frame = self.frame
            # set some NAs
            frame.ix[5:10] = np.nan
            frame.ix[15:20, -2:] = np.nan

        f = getattr(frame, name)

        if check_dates:
            df = DataFrame({'b': date_range('1/1/2001', periods=2)})
            _f = getattr(df, name)
            result = _f()
            self.assertIsInstance(result, Series)

            df['a'] = lrange(len(df))
            result = getattr(df, name)()
            self.assertIsInstance(result, Series)
            self.assertTrue(len(result))

        if has_skipna:
            def skipna_wrapper(x):
                nona = x.dropna()
                if len(nona) == 0:
                    return np.nan
                return alternative(nona)

            def wrapper(x):
                return alternative(x.values)

            result0 = f(axis=0, skipna=False)
            result1 = f(axis=1, skipna=False)
            assert_series_equal(result0, frame.apply(wrapper),
                                check_dtype=check_dtype,
                                check_less_precise=check_less_precise)
            assert_series_equal(result1, frame.apply(wrapper, axis=1),
                                check_dtype=False,
                                check_less_precise=check_less_precise)  # HACK: win32
        else:
            skipna_wrapper = alternative
            wrapper = alternative

        result0 = f(axis=0)
        result1 = f(axis=1)
        assert_series_equal(result0, frame.apply(skipna_wrapper),
                            check_dtype=check_dtype,
                            check_less_precise=check_less_precise)
        if not tm._incompat_bottleneck_version(name):
            assert_series_equal(result1, frame.apply(skipna_wrapper, axis=1),
                                check_dtype=False,
                                check_less_precise=check_less_precise)

        # check dtypes
        if check_dtype:
            lcd_dtype = frame.values.dtype
            self.assertEqual(lcd_dtype, result0.dtype)
            self.assertEqual(lcd_dtype, result1.dtype)

        # result = f(axis=1)
        # comp = frame.apply(alternative, axis=1).reindex(result.index)
        # assert_series_equal(result, comp)

        # bad axis
        assertRaisesRegexp(ValueError, 'No axis named 2', f, axis=2)
        # make sure works on mixed-type frame
        getattr(self.mixed_frame, name)(axis=0)
        getattr(self.mixed_frame, name)(axis=1)

        if has_numeric_only:
            getattr(self.mixed_frame, name)(axis=0, numeric_only=True)
            getattr(self.mixed_frame, name)(axis=1, numeric_only=True)
            getattr(self.frame, name)(axis=0, numeric_only=False)
            getattr(self.frame, name)(axis=1, numeric_only=False)

        # all NA case
        if has_skipna:
            all_na = self.frame * np.NaN
            r0 = getattr(all_na, name)(axis=0)
            r1 = getattr(all_na, name)(axis=1)
            if not tm._incompat_bottleneck_version(name):
                self.assertTrue(np.isnan(r0).all())
                self.assertTrue(np.isnan(r1).all())

    def test_mode(self):
        df = pd.DataFrame({"A": [12, 12, 11, 12, 19, 11],
                           "B": [10, 10, 10, np.nan, 3, 4],
                           "C": [8, 8, 8, 9, 9, 9],
                           "D": np.arange(6,dtype='int64'),
                           "E": [8, 8, 1, 1, 3, 3]})
        assert_frame_equal(df[["A"]].mode(),
                           pd.DataFrame({"A": [12]}))
        expected = pd.Series([], dtype='int64', name='D').to_frame()
        assert_frame_equal(df[["D"]].mode(), expected)
        expected = pd.Series([1, 3, 8], dtype='int64', name='E').to_frame()
        assert_frame_equal(df[["E"]].mode(), expected)
        assert_frame_equal(df[["A", "B"]].mode(),
                           pd.DataFrame({"A": [12], "B": [10.]}))
        assert_frame_equal(df.mode(),
                           pd.DataFrame({"A": [12, np.nan, np.nan],
                                         "B": [10, np.nan, np.nan],
                                         "C": [8, 9, np.nan],
                                         "D": [np.nan, np.nan, np.nan],
                                         "E": [1, 3, 8]}))

        # outputs in sorted order
        df["C"] = list(reversed(df["C"]))
        com.pprint_thing(df["C"])
        com.pprint_thing(df["C"].mode())
        a, b = (df[["A", "B", "C"]].mode(),
                           pd.DataFrame({"A": [12, np.nan],
                                         "B": [10, np.nan],
                                         "C": [8, 9]}))
        com.pprint_thing(a)
        com.pprint_thing(b)
        assert_frame_equal(a, b)
        # should work with heterogeneous types
        df = pd.DataFrame({"A": np.arange(6,dtype='int64'),
                           "B": pd.date_range('2011', periods=6),
                           "C": list('abcdef')})
        exp = pd.DataFrame({"A": pd.Series([], dtype=df["A"].dtype),
                            "B": pd.Series([], dtype=df["B"].dtype),
                            "C": pd.Series([], dtype=df["C"].dtype)})
        assert_frame_equal(df.mode(), exp)

        # and also when not empty
        df.loc[1, "A"] = 0
        df.loc[4, "B"] = df.loc[3, "B"]
        df.loc[5, "C"] = 'e'
        exp = pd.DataFrame({"A": pd.Series([0], dtype=df["A"].dtype),
                            "B": pd.Series([df.loc[3, "B"]], dtype=df["B"].dtype),
                            "C": pd.Series(['e'], dtype=df["C"].dtype)})

        assert_frame_equal(df.mode(), exp)

    def test_sum_corner(self):
        axis0 = self.empty.sum(0)
        axis1 = self.empty.sum(1)
        tm.assertIsInstance(axis0, Series)
        tm.assertIsInstance(axis1, Series)
        self.assertEqual(len(axis0), 0)
        self.assertEqual(len(axis1), 0)

    def test_sum_object(self):
        values = self.frame.values.astype(int)
        frame = DataFrame(values, index=self.frame.index,
                          columns=self.frame.columns)
        deltas = frame * timedelta(1)
        deltas.sum()

    def test_sum_bool(self):
        # ensure this works, bug report
        bools = np.isnan(self.frame)
        bools.sum(1)
        bools.sum(0)

    def test_mean_corner(self):
        # unit test when have object data
        the_mean = self.mixed_frame.mean(axis=0)
        the_sum = self.mixed_frame.sum(axis=0, numeric_only=True)
        self.assertTrue(the_sum.index.equals(the_mean.index))
        self.assertTrue(len(the_mean.index) < len(self.mixed_frame.columns))

        # xs sum mixed type, just want to know it works...
        the_mean = self.mixed_frame.mean(axis=1)
        the_sum = self.mixed_frame.sum(axis=1, numeric_only=True)
        self.assertTrue(the_sum.index.equals(the_mean.index))

        # take mean of boolean column
        self.frame['bool'] = self.frame['A'] > 0
        means = self.frame.mean(0)
        self.assertEqual(means['bool'], self.frame['bool'].values.mean())

    def test_stats_mixed_type(self):
        # don't blow up
        self.mixed_frame.std(1)
        self.mixed_frame.var(1)
        self.mixed_frame.mean(1)
        self.mixed_frame.skew(1)

    def test_median_corner(self):
        def wrapper(x):
            if isnull(x).any():
                return np.nan
            return np.median(x)

        self._check_stat_op('median', wrapper, frame=self.intframe,
                            check_dtype=False, check_dates=True)

    def test_quantile(self):
        from numpy import percentile

        q = self.tsframe.quantile(0.1, axis=0)
        self.assertEqual(q['A'], percentile(self.tsframe['A'], 10))
        q = self.tsframe.quantile(0.9, axis=1)
        q = self.intframe.quantile(0.1)
        self.assertEqual(q['A'], percentile(self.intframe['A'], 10))

        # test degenerate case
        q = DataFrame({'x': [], 'y': []}).quantile(0.1, axis=0)
        assert(np.isnan(q['x']) and np.isnan(q['y']))

        # non-numeric exclusion
        df = DataFrame({'col1':['A','A','B','B'], 'col2':[1,2,3,4]})
        rs = df.quantile(0.5)
        xp = df.median()
        assert_series_equal(rs, xp)

        # axis
        df = DataFrame({"A": [1, 2, 3], "B": [2, 3, 4]}, index=[1, 2, 3])
        result = df.quantile(.5, axis=1)
        expected = Series([1.5, 2.5, 3.5], index=[1, 2, 3])
        assert_series_equal(result, expected)

        result = df.quantile([.5, .75], axis=1)
        expected = DataFrame({1: [1.5, 1.75], 2: [2.5, 2.75],
                              3: [3.5, 3.75]}, index=[0.5, 0.75])
        assert_frame_equal(result, expected, check_index_type=True)

        # We may want to break API in the future to change this
        # so that we exclude non-numeric along the same axis
        # See GH #7312
        df = DataFrame([[1, 2, 3],
                        ['a', 'b', 4]])
        result = df.quantile(.5, axis=1)
        expected = Series([3., 4.], index=[0, 1])
        assert_series_equal(result, expected)

    def test_quantile_axis_parameter(self):
        # GH 9543/9544
        from numpy import percentile

        df = DataFrame({"A": [1, 2, 3], "B": [2, 3, 4]}, index=[1, 2, 3])

        result = df.quantile(.5, axis=0)

        expected = Series([2., 3.], index=["A", "B"])
        assert_series_equal(result, expected)

        expected = df.quantile(.5, axis="index")
        assert_series_equal(result, expected)

        result = df.quantile(.5, axis=1)

        expected = Series([1.5, 2.5, 3.5], index=[1, 2, 3])
        assert_series_equal(result, expected)

        result = df.quantile(.5, axis="columns")
        assert_series_equal(result, expected)

        self.assertRaises(ValueError, df.quantile, 0.1, axis=-1)
        self.assertRaises(ValueError, df.quantile, 0.1, axis="column")

    def test_quantile_multi(self):
        df = DataFrame([[1, 1, 1], [2, 2, 2], [3, 3, 3]],
                       columns=['a', 'b', 'c'])
        result = df.quantile([.25, .5])
        expected = DataFrame([[1.5, 1.5, 1.5], [2., 2., 2.]],
                             index=[.25, .5], columns=['a', 'b', 'c'])
        assert_frame_equal(result, expected)

        # axis = 1
        result = df.quantile([.25, .5], axis=1)
        expected = DataFrame([[1.5, 1.5, 1.5], [2., 2., 2.]],
                             index=[.25, .5], columns=[0, 1, 2])

        # empty
        result = DataFrame({'x': [], 'y': []}).quantile([0.1, .9], axis=0)
        expected = DataFrame({'x': [np.nan, np.nan], 'y': [np.nan, np.nan]},
                             index=[.1, .9])
        assert_frame_equal(result, expected)

    def test_quantile_datetime(self):
        df = DataFrame({'a': pd.to_datetime(['2010', '2011']), 'b': [0, 5]})

        # exclude datetime
        result = df.quantile(.5)
        expected = Series([2.5], index=['b'])

        # datetime
        result = df.quantile(.5, numeric_only=False)
        expected = Series([Timestamp('2010-07-02 12:00:00'), 2.5],
                          index=['a', 'b'])
        assert_series_equal(result, expected)

        # datetime w/ multi
        result = df.quantile([.5], numeric_only=False)
        expected = DataFrame([[Timestamp('2010-07-02 12:00:00'), 2.5]],
                             index=[.5], columns=['a', 'b'])
        assert_frame_equal(result, expected)

        # axis = 1
        df['c'] = pd.to_datetime(['2011', '2012'])
        result = df[['a', 'c']].quantile(.5, axis=1, numeric_only=False)
        expected = Series([Timestamp('2010-07-02 12:00:00'),
                           Timestamp('2011-07-02 12:00:00')],
                          index=[0, 1])
        assert_series_equal(result, expected)

        result = df[['a', 'c']].quantile([.5], axis=1, numeric_only=False)
        expected = DataFrame([[Timestamp('2010-07-02 12:00:00'),
                               Timestamp('2011-07-02 12:00:00')]],
                             index=[0.5], columns=[0, 1])
        assert_frame_equal(result, expected)

    def test_quantile_invalid(self):
        msg = 'percentiles should all be in the interval \\[0, 1\\]'
        for invalid in [-1, 2, [0.5, -1], [0.5, 2]]:
            with tm.assertRaisesRegexp(ValueError, msg):
                self.tsframe.quantile(invalid)

    def test_cumsum(self):
        self.tsframe.ix[5:10, 0] = nan
        self.tsframe.ix[10:15, 1] = nan
        self.tsframe.ix[15:, 2] = nan

        # axis = 0
        cumsum = self.tsframe.cumsum()
        expected = self.tsframe.apply(Series.cumsum)
        assert_frame_equal(cumsum, expected)

        # axis = 1
        cumsum = self.tsframe.cumsum(axis=1)
        expected = self.tsframe.apply(Series.cumsum, axis=1)
        assert_frame_equal(cumsum, expected)

        # works
        df = DataFrame({'A': np.arange(20)}, index=np.arange(20))
        result = df.cumsum()

        # fix issue
        cumsum_xs = self.tsframe.cumsum(axis=1)
        self.assertEqual(np.shape(cumsum_xs), np.shape(self.tsframe))

    def test_cumprod(self):
        self.tsframe.ix[5:10, 0] = nan
        self.tsframe.ix[10:15, 1] = nan
        self.tsframe.ix[15:, 2] = nan

        # axis = 0
        cumprod = self.tsframe.cumprod()
        expected = self.tsframe.apply(Series.cumprod)
        assert_frame_equal(cumprod, expected)

        # axis = 1
        cumprod = self.tsframe.cumprod(axis=1)
        expected = self.tsframe.apply(Series.cumprod, axis=1)
        assert_frame_equal(cumprod, expected)

        # fix issue
        cumprod_xs = self.tsframe.cumprod(axis=1)
        self.assertEqual(np.shape(cumprod_xs), np.shape(self.tsframe))

        # ints
        df = self.tsframe.fillna(0).astype(int)
        df.cumprod(0)
        df.cumprod(1)

        # ints32
        df = self.tsframe.fillna(0).astype(np.int32)
        df.cumprod(0)
        df.cumprod(1)

    def test_rank(self):
        tm._skip_if_no_scipy()
        from scipy.stats import rankdata

        self.frame['A'][::2] = np.nan
        self.frame['B'][::3] = np.nan
        self.frame['C'][::4] = np.nan
        self.frame['D'][::5] = np.nan

        ranks0 = self.frame.rank()
        ranks1 = self.frame.rank(1)
        mask = np.isnan(self.frame.values)

        fvals = self.frame.fillna(np.inf).values

        exp0 = np.apply_along_axis(rankdata, 0, fvals)
        exp0[mask] = np.nan

        exp1 = np.apply_along_axis(rankdata, 1, fvals)
        exp1[mask] = np.nan

        assert_almost_equal(ranks0.values, exp0)
        assert_almost_equal(ranks1.values, exp1)

        # integers
        df = DataFrame(np.random.randint(0, 5, size=40).reshape((10, 4)))

        result = df.rank()
        exp = df.astype(float).rank()
        assert_frame_equal(result, exp)

        result = df.rank(1)
        exp = df.astype(float).rank(1)
        assert_frame_equal(result, exp)

    def test_rank2(self):
        from datetime import datetime
        df = DataFrame([[1, 3, 2], [1, 2, 3]])
        expected = DataFrame([[1.0, 3.0, 2.0], [1, 2, 3]]) / 3.0
        result = df.rank(1, pct=True)
        assert_frame_equal(result, expected)

        df = DataFrame([[1, 3, 2], [1, 2, 3]])
        expected = df.rank(0) / 2.0
        result = df.rank(0, pct=True)
        assert_frame_equal(result, expected)



        df = DataFrame([['b', 'c', 'a'], ['a', 'c', 'b']])
        expected = DataFrame([[2.0, 3.0, 1.0], [1, 3, 2]])
        result = df.rank(1, numeric_only=False)
        assert_frame_equal(result, expected)


        expected = DataFrame([[2.0, 1.5, 1.0], [1, 1.5, 2]])
        result = df.rank(0, numeric_only=False)
        assert_frame_equal(result, expected)

        df = DataFrame([['b', np.nan, 'a'], ['a', 'c', 'b']])
        expected = DataFrame([[2.0, nan, 1.0], [1.0, 3.0, 2.0]])
        result = df.rank(1, numeric_only=False)
        assert_frame_equal(result, expected)

        expected = DataFrame([[2.0, nan, 1.0], [1.0, 1.0, 2.0]])
        result = df.rank(0, numeric_only=False)
        assert_frame_equal(result, expected)

        # f7u12, this does not work without extensive workaround
        data = [[datetime(2001, 1, 5), nan, datetime(2001, 1, 2)],
                [datetime(2000, 1, 2), datetime(2000, 1, 3),
                 datetime(2000, 1, 1)]]
        df = DataFrame(data)

        # check the rank
        expected = DataFrame([[2., nan, 1.],
                              [2., 3., 1.]])
        result = df.rank(1, numeric_only=False)
        assert_frame_equal(result, expected)

        # mixed-type frames
        self.mixed_frame['datetime'] = datetime.now()
        self.mixed_frame['timedelta'] = timedelta(days=1,seconds=1)

        result = self.mixed_frame.rank(1)
        expected = self.mixed_frame.rank(1, numeric_only=True)
        assert_frame_equal(result, expected)

        df = DataFrame({"a":[1e-20, -5, 1e-20+1e-40, 10, 1e60, 1e80, 1e-30]})
        exp = DataFrame({"a":[ 3.5,  1. ,  3.5,  5. ,  6. ,  7. ,  2. ]})
        assert_frame_equal(df.rank(), exp)

    def test_rank_na_option(self):
        tm._skip_if_no_scipy()
        from scipy.stats import rankdata

        self.frame['A'][::2] = np.nan
        self.frame['B'][::3] = np.nan
        self.frame['C'][::4] = np.nan
        self.frame['D'][::5] = np.nan

        # bottom
        ranks0 = self.frame.rank(na_option='bottom')
        ranks1 = self.frame.rank(1, na_option='bottom')

        fvals = self.frame.fillna(np.inf).values

        exp0 = np.apply_along_axis(rankdata, 0, fvals)
        exp1 = np.apply_along_axis(rankdata, 1, fvals)

        assert_almost_equal(ranks0.values, exp0)
        assert_almost_equal(ranks1.values, exp1)

        # top
        ranks0 = self.frame.rank(na_option='top')
        ranks1 = self.frame.rank(1, na_option='top')

        fval0 = self.frame.fillna((self.frame.min() - 1).to_dict()).values
        fval1 = self.frame.T
        fval1 = fval1.fillna((fval1.min() - 1).to_dict()).T
        fval1 = fval1.fillna(np.inf).values

        exp0 = np.apply_along_axis(rankdata, 0, fval0)
        exp1 = np.apply_along_axis(rankdata, 1, fval1)

        assert_almost_equal(ranks0.values, exp0)
        assert_almost_equal(ranks1.values, exp1)

        # descending

        # bottom
        ranks0 = self.frame.rank(na_option='top', ascending=False)
        ranks1 = self.frame.rank(1, na_option='top', ascending=False)

        fvals = self.frame.fillna(np.inf).values

        exp0 = np.apply_along_axis(rankdata, 0, -fvals)
        exp1 = np.apply_along_axis(rankdata, 1, -fvals)

        assert_almost_equal(ranks0.values, exp0)
        assert_almost_equal(ranks1.values, exp1)

        # descending

        # top
        ranks0 = self.frame.rank(na_option='bottom', ascending=False)
        ranks1 = self.frame.rank(1, na_option='bottom', ascending=False)

        fval0 = self.frame.fillna((self.frame.min() - 1).to_dict()).values
        fval1 = self.frame.T
        fval1 = fval1.fillna((fval1.min() - 1).to_dict()).T
        fval1 = fval1.fillna(np.inf).values

        exp0 = np.apply_along_axis(rankdata, 0, -fval0)
        exp1 = np.apply_along_axis(rankdata, 1, -fval1)

        assert_almost_equal(ranks0.values, exp0)
        assert_almost_equal(ranks1.values, exp1)

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

    def test_count_objects(self):
        dm = DataFrame(self.mixed_frame._series)
        df = DataFrame(self.mixed_frame._series)

        tm.assert_series_equal(dm.count(), df.count())
        tm.assert_series_equal(dm.count(1), df.count(1))

    def test_cumsum_corner(self):
        dm = DataFrame(np.arange(20).reshape(4, 5),
                       index=lrange(4), columns=lrange(5))
        result = dm.cumsum()

    #----------------------------------------------------------------------
    # Stacking / unstacking

    def test_stack_unstack(self):
        stacked = self.frame.stack()
        stacked_df = DataFrame({'foo': stacked, 'bar': stacked})

        unstacked = stacked.unstack()
        unstacked_df = stacked_df.unstack()

        assert_frame_equal(unstacked, self.frame)
        assert_frame_equal(unstacked_df['bar'], self.frame)

        unstacked_cols = stacked.unstack(0)
        unstacked_cols_df = stacked_df.unstack(0)
        assert_frame_equal(unstacked_cols.T, self.frame)
        assert_frame_equal(unstacked_cols_df['bar'].T, self.frame)

    def test_stack_ints(self):
        df = DataFrame(
             np.random.randn(30, 27),
             columns=MultiIndex.from_tuples(
                list(itertools.product(range(3), repeat=3))
            )
        )
        assert_frame_equal(
            df.stack(level=[1, 2]),
            df.stack(level=1).stack(level=1)
        )
        assert_frame_equal(
            df.stack(level=[-2, -1]),
            df.stack(level=1).stack(level=1)
        )

        df_named = df.copy()
        df_named.columns.set_names(range(3), inplace=True)
        assert_frame_equal(
            df_named.stack(level=[1, 2]),
            df_named.stack(level=1).stack(level=1)
        )

    def test_stack_mixed_levels(self):
        columns = MultiIndex.from_tuples(
            [('A', 'cat', 'long'), ('B', 'cat', 'long'),
             ('A', 'dog', 'short'), ('B', 'dog', 'short')],
            names=['exp', 'animal', 'hair_length']
        )
        df = DataFrame(randn(4, 4), columns=columns)

        animal_hair_stacked = df.stack(level=['animal', 'hair_length'])
        exp_hair_stacked = df.stack(level=['exp', 'hair_length'])

        # GH #8584: Need to check that stacking works when a number
        # is passed that is both a level name and in the range of
        # the level numbers
        df2 = df.copy()
        df2.columns.names = ['exp', 'animal', 1]
        assert_frame_equal(df2.stack(level=['animal', 1]),
                           animal_hair_stacked,  check_names=False)
        assert_frame_equal(df2.stack(level=['exp', 1]),
                           exp_hair_stacked,  check_names=False)

        # When mixed types are passed and the ints are not level
        # names, raise
        self.assertRaises(ValueError, df2.stack, level=['animal', 0])

        # GH #8584: Having 0 in the level names could raise a
        # strange error about lexsort depth
        df3 = df.copy()
        df3.columns.names = ['exp', 'animal', 0]
        assert_frame_equal(df3.stack(level=['animal', 0]),
                           animal_hair_stacked, check_names=False)

    def test_stack_int_level_names(self):
        columns = MultiIndex.from_tuples(
            [('A', 'cat', 'long'), ('B', 'cat', 'long'),
             ('A', 'dog', 'short'), ('B', 'dog', 'short')],
            names=['exp', 'animal', 'hair_length']
        )
        df = DataFrame(randn(4, 4), columns=columns)

        exp_animal_stacked = df.stack(level=['exp', 'animal'])
        animal_hair_stacked = df.stack(level=['animal', 'hair_length'])
        exp_hair_stacked = df.stack(level=['exp', 'hair_length'])

        df2 = df.copy()
        df2.columns.names = [0, 1, 2]
        assert_frame_equal(df2.stack(level=[1, 2]), animal_hair_stacked,
                           check_names=False )
        assert_frame_equal(df2.stack(level=[0, 1]), exp_animal_stacked,
                           check_names=False)
        assert_frame_equal(df2.stack(level=[0, 2]), exp_hair_stacked,
                           check_names=False)

        # Out-of-order int column names
        df3 = df.copy()
        df3.columns.names = [2, 0, 1]
        assert_frame_equal(df3.stack(level=[0, 1]), animal_hair_stacked,
                           check_names=False)
        assert_frame_equal(df3.stack(level=[2, 0]), exp_animal_stacked,
                           check_names=False)
        assert_frame_equal(df3.stack(level=[2, 1]), exp_hair_stacked,
                           check_names=False)


    def test_unstack_bool(self):
        df = DataFrame([False, False],
                       index=MultiIndex.from_arrays([['a', 'b'], ['c', 'l']]),
                       columns=['col'])
        rs = df.unstack()
        xp = DataFrame(np.array([[False, np.nan], [np.nan, False]],
                                dtype=object),
                       index=['a', 'b'],
                       columns=MultiIndex.from_arrays([['col', 'col'],
                                                       ['c', 'l']]))
        assert_frame_equal(rs, xp)

    def test_unstack_level_binding(self):
        # GH9856
        mi = pd.MultiIndex(
                levels=[[u('foo'), u('bar')], [u('one'), u('two')],
                        [u('a'), u('b')]],
                labels=[[0, 0, 1, 1], [0, 1, 0, 1], [1, 0, 1, 0]],
                names=[u('first'), u('second'), u('third')])
        s = pd.Series(0, index=mi)
        result = s.unstack([1, 2]).stack(0)

        expected_mi = pd.MultiIndex(
                        levels=[['foo', 'bar'], ['one', 'two']],
                        labels=[[0, 0, 1, 1], [0, 1, 0, 1]],
                        names=['first', 'second'])

        expected = pd.DataFrame(np.array([[np.nan, 0],
                                          [0, np.nan],
                                          [np.nan, 0],
                                          [0, np.nan]],
                                         dtype=np.float64),
                                index=expected_mi,
                                columns=pd.Index(['a', 'b'], name='third'))

        assert_frame_equal(result, expected)

    def test_unstack_to_series(self):
        # check reversibility
        data = self.frame.unstack()

        self.assertTrue(isinstance(data, Series))
        undo = data.unstack().T
        assert_frame_equal(undo, self.frame)

        # check NA handling
        data = DataFrame({'x': [1, 2, np.NaN], 'y': [3.0, 4, np.NaN]})
        data.index = Index(['a', 'b', 'c'])
        result = data.unstack()

        midx = MultiIndex(levels=[['x', 'y'], ['a', 'b', 'c']],
                          labels=[[0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 1, 2]])
        expected = Series([1, 2, np.NaN, 3, 4, np.NaN], index=midx)

        assert_series_equal(result, expected)

        # check composability of unstack
        old_data = data.copy()
        for _ in range(4):
            data = data.unstack()
        assert_frame_equal(old_data, data)

    def test_unstack_dtypes(self):

        # GH 2929
        rows = [[1, 1, 3, 4],
                [1, 2, 3, 4],
                [2, 1, 3, 4],
                [2, 2, 3, 4]]

        df = DataFrame(rows, columns=list('ABCD'))
        result = df.get_dtype_counts()
        expected = Series({'int64' : 4})
        assert_series_equal(result, expected)

        # single dtype
        df2 = df.set_index(['A','B'])
        df3 = df2.unstack('B')
        result = df3.get_dtype_counts()
        expected = Series({'int64' : 4})
        assert_series_equal(result, expected)

        # mixed
        df2 = df.set_index(['A','B'])
        df2['C'] = 3.
        df3 = df2.unstack('B')
        result = df3.get_dtype_counts()
        expected = Series({'int64' : 2, 'float64' : 2})
        assert_series_equal(result, expected)

        df2['D'] = 'foo'
        df3 = df2.unstack('B')
        result = df3.get_dtype_counts()
        expected = Series({'float64' : 2, 'object' : 2})
        assert_series_equal(result, expected)

        # GH7405
        for c, d in (np.zeros(5), np.zeros(5)), \
                    (np.arange(5, dtype='f8'), np.arange(5, 10, dtype='f8')):

            df = DataFrame({'A': ['a']*5, 'C':c, 'D':d,
                            'B':pd.date_range('2012-01-01', periods=5)})

            right = df.iloc[:3].copy(deep=True)

            df = df.set_index(['A', 'B'])
            df['D'] = df['D'].astype('int64')

            left = df.iloc[:3].unstack(0)
            right = right.set_index(['A', 'B']).unstack(0)
            right[('D', 'a')] = right[('D', 'a')].astype('int64')

            self.assertEqual(left.shape, (3, 2))
            tm.assert_frame_equal(left, right)

    def test_unstack_non_unique_index_names(self):
        idx = MultiIndex.from_tuples([('a', 'b'), ('c', 'd')],
                                     names=['c1', 'c1'])
        df = DataFrame([1, 2], index=idx)
        with tm.assertRaises(ValueError):
            df.unstack('c1')

        with tm.assertRaises(ValueError):
            df.T.stack('c1')

    def test_unstack_nan_index(self):  # GH7466
        cast = lambda val: '{0:1}'.format('' if val != val else val)
        nan = np.nan

        def verify(df):
            mk_list = lambda a: list(a) if isinstance(a, tuple) else [a]
            rows, cols = df.notnull().values.nonzero()
            for i, j in zip(rows, cols):
                left = sorted(df.iloc[i, j].split('.'))
                right = mk_list(df.index[i]) + mk_list(df.columns[j])
                right = sorted(list(map(cast, right)))
                self.assertEqual(left, right)

        df = DataFrame({'jim':['a', 'b', nan, 'd'],
                        'joe':['w', 'x', 'y', 'z'],
                        'jolie':['a.w', 'b.x', ' .y', 'd.z']})

        left  = df.set_index(['jim', 'joe']).unstack()['jolie']
        right = df.set_index(['joe', 'jim']).unstack()['jolie'].T
        assert_frame_equal(left, right)

        for idx in permutations(df.columns[:2]):
            mi = df.set_index(list(idx))
            for lev in range(2):
                udf = mi.unstack(level=lev)
                self.assertEqual(udf.notnull().values.sum(), len(df))
                verify(udf['jolie'])

        df = DataFrame({'1st':['d'] * 3 + [nan] * 5 + ['a'] * 2 +
                              ['c'] * 3 + ['e'] * 2 + ['b'] * 5,
                        '2nd':['y'] * 2 + ['w'] * 3 + [nan] * 3 +
                              ['z'] * 4 + [nan] * 3 + ['x'] * 3 + [nan] * 2,
                        '3rd':[67,39,53,72,57,80,31,18,11,30,59,
                               50,62,59,76,52,14,53,60,51]})

        df['4th'], df['5th'] = \
                df.apply(lambda r: '.'.join(map(cast, r)), axis=1), \
                df.apply(lambda r: '.'.join(map(cast, r.iloc[::-1])), axis=1)

        for idx in permutations(['1st', '2nd', '3rd']):
            mi = df.set_index(list(idx))
            for lev in range(3):
                udf = mi.unstack(level=lev)
                self.assertEqual(udf.notnull().values.sum(), 2 * len(df))
                for col in ['4th', '5th']:
                    verify(udf[col])

        # GH7403
        df = pd.DataFrame({'A': list('aaaabbbb'),'B':range(8), 'C':range(8)})
        df.iloc[3, 1] = np.NaN
        left = df.set_index(['A', 'B']).unstack(0)

        vals = [[3, 0, 1, 2, nan, nan, nan, nan],
                [nan, nan, nan, nan, 4, 5, 6, 7]]
        vals = list(map(list, zip(*vals)))
        idx = Index([nan, 0, 1, 2, 4, 5, 6, 7], name='B')
        cols = MultiIndex(levels=[['C'], ['a', 'b']],
                          labels=[[0, 0], [0, 1]],
                          names=[None, 'A'])

        right = DataFrame(vals, columns=cols, index=idx)
        assert_frame_equal(left, right)

        df = DataFrame({'A': list('aaaabbbb'), 'B':list(range(4))*2,
                        'C':range(8)})
        df.iloc[2,1] = np.NaN
        left = df.set_index(['A', 'B']).unstack(0)

        vals = [[2, nan], [0, 4], [1, 5], [nan, 6], [3, 7]]
        cols = MultiIndex(levels=[['C'], ['a', 'b']],
                          labels=[[0, 0], [0, 1]],
                          names=[None, 'A'])
        idx = Index([nan, 0, 1, 2, 3], name='B')
        right = DataFrame(vals, columns=cols, index=idx)
        assert_frame_equal(left, right)

        df = pd.DataFrame({'A': list('aaaabbbb'),'B':list(range(4))*2,
                           'C':range(8)})
        df.iloc[3,1] = np.NaN
        left = df.set_index(['A', 'B']).unstack(0)

        vals = [[3, nan], [0, 4], [1, 5], [2, 6], [nan, 7]]
        cols = MultiIndex(levels=[['C'], ['a', 'b']],
                          labels=[[0, 0], [0, 1]],
                          names=[None, 'A'])
        idx = Index([nan, 0, 1, 2, 3], name='B')
        right = DataFrame(vals, columns=cols, index=idx)
        assert_frame_equal(left, right)

        # GH7401
        df = pd.DataFrame({'A': list('aaaaabbbbb'), 'C':np.arange(10),
            'B':date_range('2012-01-01', periods=5).tolist()*2 })

        df.iloc[3,1] = np.NaN
        left = df.set_index(['A', 'B']).unstack()

        vals = np.array([[3, 0, 1, 2, nan, 4], [nan, 5, 6, 7, 8, 9]])
        idx = Index(['a', 'b'], name='A')
        cols = MultiIndex(levels=[['C'], date_range('2012-01-01', periods=5)],
                          labels=[[0, 0, 0, 0, 0, 0], [-1, 0, 1, 2, 3, 4]],
                          names=[None, 'B'])

        right = DataFrame(vals, columns=cols, index=idx)
        assert_frame_equal(left, right)

        # GH4862
        vals = [['Hg', nan, nan, 680585148],
                ['U', 0.0, nan, 680585148],
                ['Pb', 7.07e-06, nan, 680585148],
                ['Sn', 2.3614e-05, 0.0133, 680607017],
                ['Ag', 0.0, 0.0133, 680607017],
                ['Hg', -0.00015, 0.0133, 680607017]]
        df = DataFrame(vals, columns=['agent', 'change', 'dosage', 's_id'],
                index=[17263, 17264, 17265, 17266, 17267, 17268])

        left = df.copy().set_index(['s_id','dosage','agent']).unstack()

        vals = [[nan, nan, 7.07e-06, nan, 0.0],
                [0.0, -0.00015, nan, 2.3614e-05, nan]]

        idx = MultiIndex(levels=[[680585148, 680607017], [0.0133]],
                         labels=[[0, 1], [-1, 0]],
                         names=['s_id', 'dosage'])

        cols = MultiIndex(levels=[['change'], ['Ag', 'Hg', 'Pb', 'Sn', 'U']],
                          labels=[[0, 0, 0, 0, 0], [0, 1, 2, 3, 4]],
                          names=[None, 'agent'])

        right = DataFrame(vals, columns=cols, index=idx)
        assert_frame_equal(left, right)

        left = df.ix[17264:].copy().set_index(['s_id','dosage','agent'])
        assert_frame_equal(left.unstack(), right)

        # GH9497 - multiple unstack with nulls
        df = DataFrame({'1st':[1, 2, 1, 2, 1, 2],
                        '2nd':pd.date_range('2014-02-01', periods=6, freq='D'),
                        'jim':100 + np.arange(6),
                        'joe':(np.random.randn(6) * 10).round(2)})

        df['3rd'] = df['2nd'] - pd.Timestamp('2014-02-02')
        df.loc[1, '2nd'] = df.loc[3, '2nd'] = nan
        df.loc[1, '3rd'] = df.loc[4, '3rd'] = nan

        left = df.set_index(['1st', '2nd', '3rd']).unstack(['2nd', '3rd'])
        self.assertEqual(left.notnull().values.sum(), 2 * len(df))

        for col in ['jim', 'joe']:
           for _, r in df.iterrows():
               key = r['1st'], (col, r['2nd'], r['3rd'])
               self.assertEqual(r[col], left.loc[key])

    def test_stack_datetime_column_multiIndex(self):
        # GH 8039
        t = datetime(2014, 1, 1)
        df = DataFrame([1, 2, 3, 4], columns=MultiIndex.from_tuples([(t, 'A', 'B')]))
        result = df.stack()

        eidx = MultiIndex.from_product([(0, 1, 2, 3), ('B',)])
        ecols = MultiIndex.from_tuples([(t, 'A')])
        expected = DataFrame([1, 2, 3, 4], index=eidx, columns=ecols)
        assert_frame_equal(result, expected)

    def test_stack_partial_multiIndex(self):
        # GH 8844
        def _test_stack_with_multiindex(multiindex):
            df = DataFrame(np.arange(3 * len(multiindex)).reshape(3, len(multiindex)),
                           columns=multiindex)
            for level in (-1, 0, 1, [0, 1], [1, 0]):
                result = df.stack(level=level, dropna=False)

                if isinstance(level, int):
                    # Stacking a single level should not make any all-NaN rows,
                    # so df.stack(level=level, dropna=False) should be the same
                    # as df.stack(level=level, dropna=True).
                    expected = df.stack(level=level, dropna=True)
                    if isinstance(expected, Series):
                        assert_series_equal(result, expected)
                    else:
                        assert_frame_equal(result, expected)

                df.columns = MultiIndex.from_tuples(df.columns.get_values(),
                                                    names=df.columns.names)
                expected = df.stack(level=level, dropna=False)
                if isinstance(expected, Series):
                    assert_series_equal(result, expected)
                else:
                    assert_frame_equal(result, expected)

        full_multiindex = MultiIndex.from_tuples([('B', 'x'), ('B', 'z'),
                                                  ('A', 'y'),
                                                  ('C', 'x'), ('C', 'u')],
                                                 names=['Upper', 'Lower'])
        for multiindex_columns in ([0, 1, 2, 3, 4],
                                   [0, 1, 2, 3], [0, 1, 2, 4],
                                   [0, 1, 2], [1, 2, 3], [2, 3, 4],
                                   [0, 1], [0, 2], [0, 3],
                                   [0], [2], [4]):
            _test_stack_with_multiindex(full_multiindex[multiindex_columns])
            if len(multiindex_columns) > 1:
                multiindex_columns.reverse()
                _test_stack_with_multiindex(full_multiindex[multiindex_columns])

        df = DataFrame(np.arange(6).reshape(2, 3), columns=full_multiindex[[0, 1, 3]])
        result = df.stack(dropna=False)
        expected = DataFrame([[0, 2], [1, nan], [3, 5], [4, nan]],
                             index=MultiIndex(levels=[[0, 1], ['u', 'x', 'y', 'z']],
                                              labels=[[0, 0, 1, 1], [1, 3, 1, 3]],
                                              names=[None, 'Lower']),
                             columns=Index(['B', 'C'], name='Upper'),
                             dtype=df.dtypes[0])
        assert_frame_equal(result, expected)

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
        """
        in 0.14 this will return a view if possible
        a copy otherwise, but this is numpy dependent
        """

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

    def test_sum_bools(self):
        df = DataFrame(index=lrange(1), columns=lrange(10))
        bools = isnull(df)
        self.assertEqual(bools.sum(axis=1)[0], 10)

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

    def test_idxmin(self):
        frame = self.frame
        frame.ix[5:10] = np.nan
        frame.ix[15:20, -2:] = np.nan
        for skipna in [True, False]:
            for axis in [0, 1]:
                for df in [frame, self.intframe]:
                    result = df.idxmin(axis=axis, skipna=skipna)
                    expected = df.apply(
                        Series.idxmin, axis=axis, skipna=skipna)
                    assert_series_equal(result, expected)

        self.assertRaises(ValueError, frame.idxmin, axis=2)

    def test_idxmax(self):
        frame = self.frame
        frame.ix[5:10] = np.nan
        frame.ix[15:20, -2:] = np.nan
        for skipna in [True, False]:
            for axis in [0, 1]:
                for df in [frame, self.intframe]:
                    result = df.idxmax(axis=axis, skipna=skipna)
                    expected = df.apply(
                        Series.idxmax, axis=axis, skipna=skipna)
                    assert_series_equal(result, expected)

        self.assertRaises(ValueError, frame.idxmax, axis=2)

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

    def test_any_all(self):

        self._check_bool_op('any', np.any, has_skipna=True, has_bool_only=True)
        self._check_bool_op('all', np.all, has_skipna=True, has_bool_only=True)

        df = DataFrame(randn(10, 4)) > 0
        df.any(1)
        df.all(1)
        df.any(1, bool_only=True)
        df.all(1, bool_only=True)

        # skip pathological failure cases
        # class CantNonzero(object):

        #     def __nonzero__(self):
        #         raise ValueError

        # df[4] = CantNonzero()

        # it works!
        # df.any(1)
        # df.all(1)
        # df.any(1, bool_only=True)
        # df.all(1, bool_only=True)

        # df[4][4] = np.nan
        # df.any(1)
        # df.all(1)
        # df.any(1, bool_only=True)
        # df.all(1, bool_only=True)

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

    def _check_bool_op(self, name, alternative, frame=None, has_skipna=True,
                       has_bool_only=False):
        if frame is None:
            frame = self.frame > 0
            # set some NAs
            frame = DataFrame(frame.values.astype(object), frame.index,
                              frame.columns)
            frame.ix[5:10] = np.nan
            frame.ix[15:20, -2:] = np.nan

        f = getattr(frame, name)

        if has_skipna:
            def skipna_wrapper(x):
                nona = x.dropna().values
                return alternative(nona)

            def wrapper(x):
                return alternative(x.values)

            result0 = f(axis=0, skipna=False)
            result1 = f(axis=1, skipna=False)
            assert_series_equal(result0, frame.apply(wrapper))
            assert_series_equal(result1, frame.apply(wrapper, axis=1),
                                check_dtype=False)  # HACK: win32
        else:
            skipna_wrapper = alternative
            wrapper = alternative

        result0 = f(axis=0)
        result1 = f(axis=1)
        assert_series_equal(result0, frame.apply(skipna_wrapper))
        assert_series_equal(result1, frame.apply(skipna_wrapper, axis=1),
                            check_dtype=False)

        # result = f(axis=1)
        # comp = frame.apply(alternative, axis=1).reindex(result.index)
        # assert_series_equal(result, comp)

        # bad axis
        self.assertRaises(ValueError, f, axis=2)

        # make sure works on mixed-type frame
        mixed = self.mixed_frame
        mixed['_bool_'] = np.random.randn(len(mixed)) > 0
        getattr(mixed, name)(axis=0)
        getattr(mixed, name)(axis=1)

        class NonzeroFail:

            def __nonzero__(self):
                raise ValueError

        mixed['_nonzero_fail_'] = NonzeroFail()

        if has_bool_only:
            getattr(mixed, name)(axis=0, bool_only=True)
            getattr(mixed, name)(axis=1, bool_only=True)
            getattr(frame, name)(axis=0, bool_only=False)
            getattr(frame, name)(axis=1, bool_only=False)

        # all NA case
        if has_skipna:
            all_na = frame * np.NaN
            r0 = getattr(all_na, name)(axis=0)
            r1 = getattr(all_na, name)(axis=1)
            if name == 'any':
                self.assertFalse(r0.any())
                self.assertFalse(r1.any())
            else:
                self.assertTrue(r0.all())
                self.assertTrue(r1.all())

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

    def test_to_csv_date_format(self):
        from pandas import to_datetime
        pname = '__tmp_to_csv_date_format__'
        with ensure_clean(pname) as path:
            for engine in [None, 'python']:
                dt_index = self.tsframe.index
                datetime_frame = DataFrame({'A': dt_index, 'B': dt_index.shift(1)}, index=dt_index)

                datetime_frame.to_csv(path, date_format='%Y%m%d', engine=engine)
                # Check that the data was put in the specified format
                test = read_csv(path, index_col=0)

                datetime_frame_int = datetime_frame.applymap(lambda x: int(x.strftime('%Y%m%d')))
                datetime_frame_int.index = datetime_frame_int.index.map(lambda x: int(x.strftime('%Y%m%d')))

                assert_frame_equal(test, datetime_frame_int)

                datetime_frame.to_csv(path, date_format='%Y-%m-%d', engine=engine)
                # Check that the data was put in the specified format
                test = read_csv(path, index_col=0)
                datetime_frame_str = datetime_frame.applymap(lambda x: x.strftime('%Y-%m-%d'))
                datetime_frame_str.index = datetime_frame_str.index.map(lambda x: x.strftime('%Y-%m-%d'))

                assert_frame_equal(test, datetime_frame_str)

                # Check that columns get converted
                datetime_frame_columns = datetime_frame.T

                datetime_frame_columns.to_csv(path, date_format='%Y%m%d', engine=engine)

                test = read_csv(path, index_col=0)

                datetime_frame_columns = datetime_frame_columns.applymap(lambda x: int(x.strftime('%Y%m%d')))
                # Columns don't get converted to ints by read_csv
                datetime_frame_columns.columns = datetime_frame_columns.columns.map(lambda x: x.strftime('%Y%m%d'))

                assert_frame_equal(test, datetime_frame_columns)

                # test NaTs
                nat_index = to_datetime(['NaT'] * 10 + ['2000-01-01', '1/1/2000', '1-1-2000'])
                nat_frame = DataFrame({'A': nat_index}, index=nat_index)

                nat_frame.to_csv(path, date_format='%Y-%m-%d', engine=engine)

                test = read_csv(path, parse_dates=[0, 1], index_col=0)

                assert_frame_equal(test, nat_frame)

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
        tm.assert_frame_equal(ri, ei)

        ri = df.select_dtypes(include=[np.number,'category'])
        ei = df[['b', 'c', 'd', 'f']]
        tm.assert_frame_equal(ri, ei)

    def test_select_dtypes_exclude(self):
        df = DataFrame({'a': list('abc'),
                        'b': list(range(1, 4)),
                        'c': np.arange(3, 6).astype('u1'),
                        'd': np.arange(4.0, 7.0, dtype='float64'),
                        'e': [True, False, True]})
        re = df.select_dtypes(exclude=[np.number])
        ee = df[['a', 'e']]
        tm.assert_frame_equal(re, ee)

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
        tm.assert_frame_equal(r, e)

        exclude = 'datetime',
        include = 'bool', 'int64', 'int32'
        r = df.select_dtypes(include=include, exclude=exclude)
        e = df[['b', 'e']]
        tm.assert_frame_equal(r, e)

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
        tm.assert_frame_equal(r, e)

        r = df.select_dtypes(include=['i8', 'O', 'timedelta64[ns]'])
        e = df[['a', 'b', 'g']]
        tm.assert_frame_equal(r, e)

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

    def test_nlargest(self):
        # GH10393
        from string import ascii_lowercase
        df = pd.DataFrame({'a': np.random.permutation(10),
                           'b': list(ascii_lowercase[:10])})
        result = df.nlargest(5, 'a')
        expected = df.sort_values('a', ascending=False).head(5)
        tm.assert_frame_equal(result, expected)

    def test_nlargest_multiple_columns(self):
        from string import ascii_lowercase
        df = pd.DataFrame({'a': np.random.permutation(10),
                           'b': list(ascii_lowercase[:10]),
                           'c': np.random.permutation(10).astype('float64')})
        result = df.nlargest(5, ['a', 'b'])
        expected = df.sort_values(['a', 'b'], ascending=False).head(5)
        tm.assert_frame_equal(result, expected)

    def test_nsmallest(self):
        from string import ascii_lowercase
        df = pd.DataFrame({'a': np.random.permutation(10),
                           'b': list(ascii_lowercase[:10])})
        result = df.nsmallest(5, 'a')
        expected = df.sort_values('a').head(5)
        tm.assert_frame_equal(result, expected)

    def test_nsmallest_multiple_columns(self):
        from string import ascii_lowercase
        df = pd.DataFrame({'a': np.random.permutation(10),
                           'b': list(ascii_lowercase[:10]),
                           'c': np.random.permutation(10).astype('float64')})
        result = df.nsmallest(5, ['a', 'c'])
        expected = df.sort_values(['a', 'c']).head(5)
        tm.assert_frame_equal(result, expected)

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


def skip_if_no_ne(engine='numexpr'):
    if engine == 'numexpr':
        try:
            import numexpr as ne
        except ImportError:
            raise nose.SkipTest("cannot query engine numexpr when numexpr not "
                                "installed")


def skip_if_no_pandas_parser(parser):
    if parser != 'pandas':
        raise nose.SkipTest("cannot evaluate with parser {0!r}".format(parser))


class TestDataFrameQueryWithMultiIndex(object):
    def check_query_with_named_multiindex(self, parser, engine):
        tm.skip_if_no_ne(engine)
        a = tm.choice(['red', 'green'], size=10)
        b = tm.choice(['eggs', 'ham'], size=10)
        index = MultiIndex.from_arrays([a, b], names=['color', 'food'])
        df = DataFrame(randn(10, 2), index=index)
        ind = Series(df.index.get_level_values('color').values, index=index,
                     name='color')

        # equality
        res1 = df.query('color == "red"', parser=parser, engine=engine)
        res2 = df.query('"red" == color', parser=parser, engine=engine)
        exp = df[ind == 'red']
        assert_frame_equal(res1, exp)
        assert_frame_equal(res2, exp)

        # inequality
        res1 = df.query('color != "red"', parser=parser, engine=engine)
        res2 = df.query('"red" != color', parser=parser, engine=engine)
        exp = df[ind != 'red']
        assert_frame_equal(res1, exp)
        assert_frame_equal(res2, exp)

        # list equality (really just set membership)
        res1 = df.query('color == ["red"]', parser=parser, engine=engine)
        res2 = df.query('["red"] == color', parser=parser, engine=engine)
        exp = df[ind.isin(['red'])]
        assert_frame_equal(res1, exp)
        assert_frame_equal(res2, exp)

        res1 = df.query('color != ["red"]', parser=parser, engine=engine)
        res2 = df.query('["red"] != color', parser=parser, engine=engine)
        exp = df[~ind.isin(['red'])]
        assert_frame_equal(res1, exp)
        assert_frame_equal(res2, exp)

        # in/not in ops
        res1 = df.query('["red"] in color', parser=parser, engine=engine)
        res2 = df.query('"red" in color', parser=parser, engine=engine)
        exp = df[ind.isin(['red'])]
        assert_frame_equal(res1, exp)
        assert_frame_equal(res2, exp)

        res1 = df.query('["red"] not in color', parser=parser, engine=engine)
        res2 = df.query('"red" not in color', parser=parser, engine=engine)
        exp = df[~ind.isin(['red'])]
        assert_frame_equal(res1, exp)
        assert_frame_equal(res2, exp)

    def test_query_with_named_multiindex(self):
        for parser, engine in product(['pandas'], ENGINES):
            yield self.check_query_with_named_multiindex, parser, engine

    def check_query_with_unnamed_multiindex(self, parser, engine):
        tm.skip_if_no_ne(engine)
        a = tm.choice(['red', 'green'], size=10)
        b = tm.choice(['eggs', 'ham'], size=10)
        index = MultiIndex.from_arrays([a, b])
        df = DataFrame(randn(10, 2), index=index)
        ind = Series(df.index.get_level_values(0).values, index=index)

        res1 = df.query('ilevel_0 == "red"', parser=parser, engine=engine)
        res2 = df.query('"red" == ilevel_0', parser=parser, engine=engine)
        exp = df[ind == 'red']
        assert_frame_equal(res1, exp)
        assert_frame_equal(res2, exp)

        # inequality
        res1 = df.query('ilevel_0 != "red"', parser=parser, engine=engine)
        res2 = df.query('"red" != ilevel_0', parser=parser, engine=engine)
        exp = df[ind != 'red']
        assert_frame_equal(res1, exp)
        assert_frame_equal(res2, exp)

        # list equality (really just set membership)
        res1 = df.query('ilevel_0 == ["red"]', parser=parser, engine=engine)
        res2 = df.query('["red"] == ilevel_0', parser=parser, engine=engine)
        exp = df[ind.isin(['red'])]
        assert_frame_equal(res1, exp)
        assert_frame_equal(res2, exp)

        res1 = df.query('ilevel_0 != ["red"]', parser=parser, engine=engine)
        res2 = df.query('["red"] != ilevel_0', parser=parser, engine=engine)
        exp = df[~ind.isin(['red'])]
        assert_frame_equal(res1, exp)
        assert_frame_equal(res2, exp)

        # in/not in ops
        res1 = df.query('["red"] in ilevel_0', parser=parser, engine=engine)
        res2 = df.query('"red" in ilevel_0', parser=parser, engine=engine)
        exp = df[ind.isin(['red'])]
        assert_frame_equal(res1, exp)
        assert_frame_equal(res2, exp)

        res1 = df.query('["red"] not in ilevel_0', parser=parser, engine=engine)
        res2 = df.query('"red" not in ilevel_0', parser=parser, engine=engine)
        exp = df[~ind.isin(['red'])]
        assert_frame_equal(res1, exp)
        assert_frame_equal(res2, exp)

        #### LEVEL 1 ####
        ind = Series(df.index.get_level_values(1).values, index=index)
        res1 = df.query('ilevel_1 == "eggs"', parser=parser, engine=engine)
        res2 = df.query('"eggs" == ilevel_1', parser=parser, engine=engine)
        exp = df[ind == 'eggs']
        assert_frame_equal(res1, exp)
        assert_frame_equal(res2, exp)

        # inequality
        res1 = df.query('ilevel_1 != "eggs"', parser=parser, engine=engine)
        res2 = df.query('"eggs" != ilevel_1', parser=parser, engine=engine)
        exp = df[ind != 'eggs']
        assert_frame_equal(res1, exp)
        assert_frame_equal(res2, exp)

        # list equality (really just set membership)
        res1 = df.query('ilevel_1 == ["eggs"]', parser=parser, engine=engine)
        res2 = df.query('["eggs"] == ilevel_1', parser=parser, engine=engine)
        exp = df[ind.isin(['eggs'])]
        assert_frame_equal(res1, exp)
        assert_frame_equal(res2, exp)

        res1 = df.query('ilevel_1 != ["eggs"]', parser=parser, engine=engine)
        res2 = df.query('["eggs"] != ilevel_1', parser=parser, engine=engine)
        exp = df[~ind.isin(['eggs'])]
        assert_frame_equal(res1, exp)
        assert_frame_equal(res2, exp)

        # in/not in ops
        res1 = df.query('["eggs"] in ilevel_1', parser=parser, engine=engine)
        res2 = df.query('"eggs" in ilevel_1', parser=parser, engine=engine)
        exp = df[ind.isin(['eggs'])]
        assert_frame_equal(res1, exp)
        assert_frame_equal(res2, exp)

        res1 = df.query('["eggs"] not in ilevel_1', parser=parser, engine=engine)
        res2 = df.query('"eggs" not in ilevel_1', parser=parser, engine=engine)
        exp = df[~ind.isin(['eggs'])]
        assert_frame_equal(res1, exp)
        assert_frame_equal(res2, exp)

    def test_query_with_unnamed_multiindex(self):
        for parser, engine in product(['pandas'], ENGINES):
            yield self.check_query_with_unnamed_multiindex, parser, engine

    def check_query_with_partially_named_multiindex(self, parser, engine):
        tm.skip_if_no_ne(engine)
        a = tm.choice(['red', 'green'], size=10)
        b = np.arange(10)
        index = MultiIndex.from_arrays([a, b])
        index.names = [None, 'rating']
        df = DataFrame(randn(10, 2), index=index)
        res = df.query('rating == 1', parser=parser, engine=engine)
        ind = Series(df.index.get_level_values('rating').values, index=index,
                     name='rating')
        exp = df[ind == 1]
        assert_frame_equal(res, exp)

        res = df.query('rating != 1', parser=parser, engine=engine)
        ind = Series(df.index.get_level_values('rating').values, index=index,
                     name='rating')
        exp = df[ind != 1]
        assert_frame_equal(res, exp)

        res = df.query('ilevel_0 == "red"', parser=parser, engine=engine)
        ind = Series(df.index.get_level_values(0).values, index=index)
        exp = df[ind == "red"]
        assert_frame_equal(res, exp)

        res = df.query('ilevel_0 != "red"', parser=parser, engine=engine)
        ind = Series(df.index.get_level_values(0).values, index=index)
        exp = df[ind != "red"]
        assert_frame_equal(res, exp)

    def test_query_with_partially_named_multiindex(self):
        for parser, engine in product(['pandas'], ENGINES):
            yield self.check_query_with_partially_named_multiindex, parser, engine

    def test_query_multiindex_get_index_resolvers(self):
        for parser, engine in product(['pandas'], ENGINES):
            yield self.check_query_multiindex_get_index_resolvers, parser, engine

    def check_query_multiindex_get_index_resolvers(self, parser, engine):
        df = mkdf(10, 3, r_idx_nlevels=2, r_idx_names=['spam', 'eggs'])
        resolvers = df._get_index_resolvers()

        def to_series(mi, level):
            level_values = mi.get_level_values(level)
            s = level_values.to_series()
            s.index = mi
            return s

        col_series = df.columns.to_series()
        expected = {'index': df.index,
                    'columns': col_series,
                    'spam': to_series(df.index, 'spam'),
                    'eggs': to_series(df.index, 'eggs'),
                    'C0': col_series}
        for k, v in resolvers.items():
            if isinstance(v, Index):
                assert v.is_(expected[k])
            elif isinstance(v, Series):
                tm.assert_series_equal(v, expected[k])
            else:
                raise AssertionError("object must be a Series or Index")

    def test_raise_on_panel_with_multiindex(self):
        for parser, engine in product(PARSERS, ENGINES):
            yield self.check_raise_on_panel_with_multiindex, parser, engine

    def check_raise_on_panel_with_multiindex(self, parser, engine):
        tm.skip_if_no_ne()
        p = tm.makePanel(7)
        p.items = tm.makeCustomIndex(len(p.items), nlevels=2)
        with tm.assertRaises(NotImplementedError):
            pd.eval('p + 1', parser=parser, engine=engine)

    def test_raise_on_panel4d_with_multiindex(self):
        for parser, engine in product(PARSERS, ENGINES):
            yield self.check_raise_on_panel4d_with_multiindex, parser, engine

    def check_raise_on_panel4d_with_multiindex(self, parser, engine):
        tm.skip_if_no_ne()
        p4d = tm.makePanel4D(7)
        p4d.items = tm.makeCustomIndex(len(p4d.items), nlevels=2)
        with tm.assertRaises(NotImplementedError):
            pd.eval('p4d + 1', parser=parser, engine=engine)


class TestDataFrameQueryNumExprPandas(tm.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestDataFrameQueryNumExprPandas, cls).setUpClass()
        cls.engine = 'numexpr'
        cls.parser = 'pandas'
        tm.skip_if_no_ne(cls.engine)

    @classmethod
    def tearDownClass(cls):
        super(TestDataFrameQueryNumExprPandas, cls).tearDownClass()
        del cls.engine, cls.parser

    def test_date_query_with_attribute_access(self):
        engine, parser = self.engine, self.parser
        skip_if_no_pandas_parser(parser)
        df = DataFrame(randn(5, 3))
        df['dates1'] = date_range('1/1/2012', periods=5)
        df['dates2'] = date_range('1/1/2013', periods=5)
        df['dates3'] = date_range('1/1/2014', periods=5)
        res = df.query('@df.dates1 < 20130101 < @df.dates3', engine=engine,
                       parser=parser)
        expec = df[(df.dates1 < '20130101') & ('20130101' < df.dates3)]
        assert_frame_equal(res, expec)

    def test_date_query_no_attribute_access(self):
        engine, parser = self.engine, self.parser
        df = DataFrame(randn(5, 3))
        df['dates1'] = date_range('1/1/2012', periods=5)
        df['dates2'] = date_range('1/1/2013', periods=5)
        df['dates3'] = date_range('1/1/2014', periods=5)
        res = df.query('dates1 < 20130101 < dates3', engine=engine,
                       parser=parser)
        expec = df[(df.dates1 < '20130101') & ('20130101' < df.dates3)]
        tm.assert_frame_equal(res, expec)

    def test_date_query_with_NaT(self):
        engine, parser = self.engine, self.parser
        n = 10
        df = DataFrame(randn(n, 3))
        df['dates1'] = date_range('1/1/2012', periods=n)
        df['dates2'] = date_range('1/1/2013', periods=n)
        df['dates3'] = date_range('1/1/2014', periods=n)
        df.loc[np.random.rand(n) > 0.5, 'dates1'] = pd.NaT
        df.loc[np.random.rand(n) > 0.5, 'dates3'] = pd.NaT
        res = df.query('dates1 < 20130101 < dates3', engine=engine,
                       parser=parser)
        expec = df[(df.dates1 < '20130101') & ('20130101' < df.dates3)]
        assert_frame_equal(res, expec)

    def test_date_index_query(self):
        engine, parser = self.engine, self.parser
        n = 10
        df = DataFrame(randn(n, 3))
        df['dates1'] = date_range('1/1/2012', periods=n)
        df['dates3'] = date_range('1/1/2014', periods=n)
        df.set_index('dates1', inplace=True, drop=True)
        res = df.query('index < 20130101 < dates3', engine=engine,
                       parser=parser)
        expec = df[(df.index < '20130101') & ('20130101' < df.dates3)]
        assert_frame_equal(res, expec)

    def test_date_index_query_with_NaT(self):
        engine, parser = self.engine, self.parser
        n = 10
        df = DataFrame(randn(n, 3))
        df['dates1'] = date_range('1/1/2012', periods=n)
        df['dates3'] = date_range('1/1/2014', periods=n)
        df.iloc[0, 0] = pd.NaT
        df.set_index('dates1', inplace=True, drop=True)
        res = df.query('index < 20130101 < dates3', engine=engine,
                       parser=parser)
        expec = df[(df.index < '20130101') & ('20130101' < df.dates3)]
        assert_frame_equal(res, expec)

    def test_date_index_query_with_NaT_duplicates(self):
        engine, parser = self.engine, self.parser
        n = 10
        d = {}
        d['dates1'] = date_range('1/1/2012', periods=n)
        d['dates3'] = date_range('1/1/2014', periods=n)
        df = DataFrame(d)
        df.loc[np.random.rand(n) > 0.5, 'dates1'] = pd.NaT
        df.set_index('dates1', inplace=True, drop=True)
        res = df.query('index < 20130101 < dates3', engine=engine, parser=parser)
        expec = df[(df.index.to_series() < '20130101') & ('20130101' < df.dates3)]
        assert_frame_equal(res, expec)

    def test_date_query_with_non_date(self):
        engine, parser = self.engine, self.parser

        n = 10
        df = DataFrame({'dates': date_range('1/1/2012', periods=n),
                        'nondate': np.arange(n)})

        ops = '==', '!=', '<', '>', '<=', '>='

        for op in ops:
            with tm.assertRaises(TypeError):
                df.query('dates %s nondate' % op, parser=parser, engine=engine)

    def test_query_syntax_error(self):
        engine, parser = self.engine, self.parser
        df = DataFrame({"i": lrange(10), "+": lrange(3, 13),
                        "r": lrange(4, 14)})
        with tm.assertRaises(SyntaxError):
            df.query('i - +', engine=engine, parser=parser)

    def test_query_scope(self):
        from pandas.computation.ops import UndefinedVariableError
        engine, parser = self.engine, self.parser
        skip_if_no_pandas_parser(parser)

        df = DataFrame(np.random.randn(20, 2), columns=list('ab'))

        a, b = 1, 2
        res = df.query('a > b', engine=engine, parser=parser)
        expected = df[df.a > df.b]
        tm.assert_frame_equal(res, expected)

        res = df.query('@a > b', engine=engine, parser=parser)
        expected = df[a > df.b]
        tm.assert_frame_equal(res, expected)

        # no local variable c
        with tm.assertRaises(UndefinedVariableError):
            df.query('@a > b > @c', engine=engine, parser=parser)

        # no column named 'c'
        with tm.assertRaises(UndefinedVariableError):
            df.query('@a > b > c', engine=engine, parser=parser)

    def test_query_doesnt_pickup_local(self):
        from pandas.computation.ops import UndefinedVariableError

        engine, parser = self.engine, self.parser
        n = m = 10
        df = DataFrame(np.random.randint(m, size=(n, 3)), columns=list('abc'))

        from numpy import sin

        # we don't pick up the local 'sin'
        with tm.assertRaises(UndefinedVariableError):
            df.query('sin > 5', engine=engine, parser=parser)

    def test_query_builtin(self):
        from pandas.computation.engines import NumExprClobberingError
        engine, parser = self.engine, self.parser

        n = m = 10
        df = DataFrame(np.random.randint(m, size=(n, 3)), columns=list('abc'))

        df.index.name = 'sin'
        with tm.assertRaisesRegexp(NumExprClobberingError,
                                  'Variables in expression.+'):
            df.query('sin > 5', engine=engine, parser=parser)

    def test_query(self):
        engine, parser = self.engine, self.parser
        df = DataFrame(np.random.randn(10, 3), columns=['a', 'b', 'c'])

        assert_frame_equal(df.query('a < b', engine=engine, parser=parser),
                           df[df.a < df.b])
        assert_frame_equal(df.query('a + b > b * c', engine=engine,
                                    parser=parser),
                           df[df.a + df.b > df.b * df.c])

    def test_query_index_with_name(self):
        engine, parser = self.engine, self.parser
        df = DataFrame(np.random.randint(10, size=(10, 3)),
                       index=Index(range(10), name='blob'),
                       columns=['a', 'b', 'c'])
        res = df.query('(blob < 5) & (a < b)', engine=engine, parser=parser)
        expec = df[(df.index < 5) & (df.a < df.b)]
        assert_frame_equal(res, expec)

        res = df.query('blob < b', engine=engine, parser=parser)
        expec = df[df.index < df.b]

        assert_frame_equal(res, expec)

    def test_query_index_without_name(self):
        engine, parser = self.engine, self.parser
        df = DataFrame(np.random.randint(10, size=(10, 3)),
                       index=range(10), columns=['a', 'b', 'c'])

        # "index" should refer to the index
        res = df.query('index < b', engine=engine, parser=parser)
        expec = df[df.index < df.b]
        assert_frame_equal(res, expec)

        # test against a scalar
        res = df.query('index < 5', engine=engine, parser=parser)
        expec = df[df.index < 5]
        assert_frame_equal(res, expec)

    def test_nested_scope(self):
        engine = self.engine
        parser = self.parser

        skip_if_no_pandas_parser(parser)

        df = DataFrame(np.random.randn(5, 3))
        df2 = DataFrame(np.random.randn(5, 3))
        expected = df[(df > 0) & (df2 > 0)]

        result = df.query('(@df > 0) & (@df2 > 0)', engine=engine, parser=parser)
        assert_frame_equal(result, expected)

        result = pd.eval('df[df > 0 and df2 > 0]', engine=engine,
                         parser=parser)
        assert_frame_equal(result, expected)

        result = pd.eval('df[df > 0 and df2 > 0 and df[df > 0] > 0]',
                         engine=engine, parser=parser)
        expected = df[(df > 0) & (df2 > 0) & (df[df > 0] > 0)]
        assert_frame_equal(result, expected)

        result = pd.eval('df[(df>0) & (df2>0)]', engine=engine, parser=parser)
        expected = df.query('(@df>0) & (@df2>0)', engine=engine, parser=parser)
        assert_frame_equal(result, expected)

    def test_nested_raises_on_local_self_reference(self):
        from pandas.computation.ops import UndefinedVariableError

        df = DataFrame(np.random.randn(5, 3))

        # can't reference ourself b/c we're a local so @ is necessary
        with tm.assertRaises(UndefinedVariableError):
            df.query('df > 0', engine=self.engine, parser=self.parser)

    def test_local_syntax(self):
        skip_if_no_pandas_parser(self.parser)

        engine, parser = self.engine, self.parser
        df = DataFrame(randn(100, 10), columns=list('abcdefghij'))
        b = 1
        expect = df[df.a < b]
        result = df.query('a < @b', engine=engine, parser=parser)
        assert_frame_equal(result, expect)

        expect = df[df.a < df.b]
        result = df.query('a < b', engine=engine, parser=parser)
        assert_frame_equal(result, expect)

    def test_chained_cmp_and_in(self):
        skip_if_no_pandas_parser(self.parser)
        engine, parser = self.engine, self.parser
        cols = list('abc')
        df = DataFrame(randn(100, len(cols)), columns=cols)
        res = df.query('a < b < c and a not in b not in c', engine=engine,
                       parser=parser)
        ind = (df.a < df.b) & (df.b < df.c) & ~df.b.isin(df.a) & ~df.c.isin(df.b)
        expec = df[ind]
        assert_frame_equal(res, expec)

    def test_local_variable_with_in(self):
        engine, parser = self.engine, self.parser
        skip_if_no_pandas_parser(parser)
        a = Series(np.random.randint(3, size=15), name='a')
        b = Series(np.random.randint(10, size=15), name='b')
        df = DataFrame({'a': a, 'b': b})

        expected = df.loc[(df.b - 1).isin(a)]
        result = df.query('b - 1 in a', engine=engine, parser=parser)
        tm.assert_frame_equal(expected, result)

        b = Series(np.random.randint(10, size=15), name='b')
        expected = df.loc[(b - 1).isin(a)]
        result = df.query('@b - 1 in a', engine=engine, parser=parser)
        tm.assert_frame_equal(expected, result)

    def test_at_inside_string(self):
        engine, parser = self.engine, self.parser
        skip_if_no_pandas_parser(parser)
        c = 1
        df = DataFrame({'a': ['a', 'a', 'b', 'b', '@c', '@c']})
        result = df.query('a == "@c"', engine=engine, parser=parser)
        expected = df[df.a == "@c"]
        tm.assert_frame_equal(result, expected)

    def test_query_undefined_local(self):
        from pandas.computation.ops import UndefinedVariableError
        engine, parser = self.engine, self.parser
        skip_if_no_pandas_parser(parser)
        df = DataFrame(np.random.rand(10, 2), columns=list('ab'))
        with tm.assertRaisesRegexp(UndefinedVariableError,
                                   "local variable 'c' is not defined"):
            df.query('a == @c', engine=engine, parser=parser)

    def test_index_resolvers_come_after_columns_with_the_same_name(self):
        n = 1
        a = np.r_[20:101:20]

        df = DataFrame({'index': a, 'b': np.random.randn(a.size)})
        df.index.name = 'index'
        result = df.query('index > 5', engine=self.engine, parser=self.parser)
        expected = df[df['index'] > 5]
        tm.assert_frame_equal(result, expected)

        df = DataFrame({'index': a, 'b': np.random.randn(a.size)})
        result = df.query('ilevel_0 > 5', engine=self.engine, parser=self.parser)
        expected = df.loc[df.index[df.index > 5]]
        tm.assert_frame_equal(result, expected)

        df = DataFrame({'a': a, 'b': np.random.randn(a.size)})
        df.index.name = 'a'
        result = df.query('a > 5', engine=self.engine, parser=self.parser)
        expected = df[df.a > 5]
        tm.assert_frame_equal(result, expected)

        result = df.query('index > 5', engine=self.engine, parser=self.parser)
        expected = df.loc[df.index[df.index > 5]]
        tm.assert_frame_equal(result, expected)

    def test_inf(self):
        n = 10
        df = DataFrame({'a': np.random.rand(n), 'b': np.random.rand(n)})
        df.loc[::2, 0] = np.inf
        ops = '==', '!='
        d = dict(zip(ops, (operator.eq, operator.ne)))
        for op, f in d.items():
            q = 'a %s inf' % op
            expected = df[f(df.a, np.inf)]
            result = df.query(q, engine=self.engine, parser=self.parser)
            tm.assert_frame_equal(result, expected)


class TestDataFrameQueryNumExprPython(TestDataFrameQueryNumExprPandas):

    @classmethod
    def setUpClass(cls):
        super(TestDataFrameQueryNumExprPython, cls).setUpClass()
        cls.engine = 'numexpr'
        cls.parser = 'python'
        tm.skip_if_no_ne(cls.engine)
        cls.frame = _frame.copy()

    def test_date_query_no_attribute_access(self):
        engine, parser = self.engine, self.parser
        df = DataFrame(randn(5, 3))
        df['dates1'] = date_range('1/1/2012', periods=5)
        df['dates2'] = date_range('1/1/2013', periods=5)
        df['dates3'] = date_range('1/1/2014', periods=5)
        res = df.query('(dates1 < 20130101) & (20130101 < dates3)',
                       engine=engine, parser=parser)
        expec = df[(df.dates1 < '20130101') & ('20130101' < df.dates3)]
        tm.assert_frame_equal(res, expec)
    def test_date_query_with_NaT(self):
        engine, parser = self.engine, self.parser
        n = 10
        df = DataFrame(randn(n, 3))
        df['dates1'] = date_range('1/1/2012', periods=n)
        df['dates2'] = date_range('1/1/2013', periods=n)
        df['dates3'] = date_range('1/1/2014', periods=n)
        df.loc[np.random.rand(n) > 0.5, 'dates1'] = pd.NaT
        df.loc[np.random.rand(n) > 0.5, 'dates3'] = pd.NaT
        res = df.query('(dates1 < 20130101) & (20130101 < dates3)',
                       engine=engine, parser=parser)
        expec = df[(df.dates1 < '20130101') & ('20130101' < df.dates3)]
        assert_frame_equal(res, expec)

    def test_date_index_query(self):
        engine, parser = self.engine, self.parser
        n = 10
        df = DataFrame(randn(n, 3))
        df['dates1'] = date_range('1/1/2012', periods=n)
        df['dates3'] = date_range('1/1/2014', periods=n)
        df.set_index('dates1', inplace=True, drop=True)
        res = df.query('(index < 20130101) & (20130101 < dates3)',
                       engine=engine, parser=parser)
        expec = df[(df.index < '20130101') & ('20130101' < df.dates3)]
        assert_frame_equal(res, expec)

    def test_date_index_query_with_NaT(self):
        engine, parser = self.engine, self.parser
        n = 10
        df = DataFrame(randn(n, 3))
        df['dates1'] = date_range('1/1/2012', periods=n)
        df['dates3'] = date_range('1/1/2014', periods=n)
        df.iloc[0, 0] = pd.NaT
        df.set_index('dates1', inplace=True, drop=True)
        res = df.query('(index < 20130101) & (20130101 < dates3)',
                       engine=engine, parser=parser)
        expec = df[(df.index < '20130101') & ('20130101' < df.dates3)]
        assert_frame_equal(res, expec)

    def test_date_index_query_with_NaT_duplicates(self):
        engine, parser = self.engine, self.parser
        n = 10
        df = DataFrame(randn(n, 3))
        df['dates1'] = date_range('1/1/2012', periods=n)
        df['dates3'] = date_range('1/1/2014', periods=n)
        df.loc[np.random.rand(n) > 0.5, 'dates1'] = pd.NaT
        df.set_index('dates1', inplace=True, drop=True)
        with tm.assertRaises(NotImplementedError):
            df.query('index < 20130101 < dates3', engine=engine, parser=parser)

    def test_nested_scope(self):
        from pandas.computation.ops import UndefinedVariableError
        engine = self.engine
        parser = self.parser
        # smoke test
        x = 1
        result = pd.eval('x + 1', engine=engine, parser=parser)
        self.assertEqual(result, 2)

        df  = DataFrame(np.random.randn(5, 3))
        df2 = DataFrame(np.random.randn(5, 3))

        # don't have the pandas parser
        with tm.assertRaises(SyntaxError):
            df.query('(@df>0) & (@df2>0)', engine=engine, parser=parser)

        with tm.assertRaises(UndefinedVariableError):
            df.query('(df>0) & (df2>0)', engine=engine, parser=parser)

        expected = df[(df > 0) & (df2 > 0)]
        result = pd.eval('df[(df > 0) & (df2 > 0)]', engine=engine,
                         parser=parser)
        tm.assert_frame_equal(expected, result)

        expected = df[(df > 0) & (df2 > 0) & (df[df > 0] > 0)]
        result = pd.eval('df[(df > 0) & (df2 > 0) & (df[df > 0] > 0)]',
                         engine=engine, parser=parser)
        tm.assert_frame_equal(expected, result)


class TestDataFrameQueryPythonPandas(TestDataFrameQueryNumExprPandas):

    @classmethod
    def setUpClass(cls):
        super(TestDataFrameQueryPythonPandas, cls).setUpClass()
        cls.engine = 'python'
        cls.parser = 'pandas'
        cls.frame = _frame.copy()

    def test_query_builtin(self):
        from pandas.computation.engines import NumExprClobberingError
        engine, parser = self.engine, self.parser

        n = m = 10
        df = DataFrame(np.random.randint(m, size=(n, 3)), columns=list('abc'))

        df.index.name = 'sin'
        expected = df[df.index > 5]
        result = df.query('sin > 5', engine=engine, parser=parser)
        tm.assert_frame_equal(expected, result)


class TestDataFrameQueryPythonPython(TestDataFrameQueryNumExprPython):

    @classmethod
    def setUpClass(cls):
        super(TestDataFrameQueryPythonPython, cls).setUpClass()
        cls.engine = cls.parser = 'python'
        cls.frame = _frame.copy()

    def test_query_builtin(self):
        from pandas.computation.engines import NumExprClobberingError
        engine, parser = self.engine, self.parser

        n = m = 10
        df = DataFrame(np.random.randint(m, size=(n, 3)), columns=list('abc'))

        df.index.name = 'sin'
        expected = df[df.index > 5]
        result = df.query('sin > 5', engine=engine, parser=parser)
        tm.assert_frame_equal(expected, result)


PARSERS = 'python', 'pandas'
ENGINES = 'python', 'numexpr'


class TestDataFrameQueryStrings(object):
    def check_str_query_method(self, parser, engine):
        tm.skip_if_no_ne(engine)
        df = DataFrame(randn(10, 1), columns=['b'])
        df['strings'] = Series(list('aabbccddee'))
        expect = df[df.strings == 'a']

        if parser != 'pandas':
            col = 'strings'
            lst = '"a"'

            lhs = [col] * 2 + [lst] * 2
            rhs = lhs[::-1]

            eq, ne = '==', '!='
            ops = 2 * ([eq] + [ne])

            for lhs, op, rhs in zip(lhs, ops, rhs):
                ex = '{lhs} {op} {rhs}'.format(lhs=lhs, op=op, rhs=rhs)
                assertRaises(NotImplementedError, df.query, ex, engine=engine,
                             parser=parser, local_dict={'strings': df.strings})
        else:
            res = df.query('"a" == strings', engine=engine, parser=parser)
            assert_frame_equal(res, expect)

            res = df.query('strings == "a"', engine=engine, parser=parser)
            assert_frame_equal(res, expect)
            assert_frame_equal(res, df[df.strings.isin(['a'])])

            expect = df[df.strings != 'a']
            res = df.query('strings != "a"', engine=engine, parser=parser)
            assert_frame_equal(res, expect)

            res = df.query('"a" != strings', engine=engine, parser=parser)
            assert_frame_equal(res, expect)
            assert_frame_equal(res, df[~df.strings.isin(['a'])])

    def test_str_query_method(self):
        for parser, engine in product(PARSERS, ENGINES):
            yield self.check_str_query_method, parser, engine

    def test_str_list_query_method(self):
        for parser, engine in product(PARSERS, ENGINES):
            yield self.check_str_list_query_method, parser, engine

    def check_str_list_query_method(self, parser, engine):
        tm.skip_if_no_ne(engine)
        df = DataFrame(randn(10, 1), columns=['b'])
        df['strings'] = Series(list('aabbccddee'))
        expect = df[df.strings.isin(['a', 'b'])]

        if parser != 'pandas':
            col = 'strings'
            lst = '["a", "b"]'

            lhs = [col] * 2 + [lst] * 2
            rhs = lhs[::-1]

            eq, ne = '==', '!='
            ops = 2 * ([eq] + [ne])

            for lhs, op, rhs in zip(lhs, ops, rhs):
                ex = '{lhs} {op} {rhs}'.format(lhs=lhs, op=op, rhs=rhs)
                with tm.assertRaises(NotImplementedError):
                    df.query(ex, engine=engine, parser=parser)
        else:
            res = df.query('strings == ["a", "b"]', engine=engine,
                           parser=parser)
            assert_frame_equal(res, expect)

            res = df.query('["a", "b"] == strings', engine=engine,
                           parser=parser)
            assert_frame_equal(res, expect)

            expect = df[~df.strings.isin(['a', 'b'])]

            res = df.query('strings != ["a", "b"]', engine=engine,
                           parser=parser)
            assert_frame_equal(res, expect)

            res = df.query('["a", "b"] != strings', engine=engine,
                           parser=parser)
            assert_frame_equal(res, expect)

    def check_query_with_string_columns(self, parser, engine):
        tm.skip_if_no_ne(engine)
        df = DataFrame({'a': list('aaaabbbbcccc'),
                        'b': list('aabbccddeeff'),
                        'c': np.random.randint(5, size=12),
                        'd': np.random.randint(9, size=12)})
        if parser == 'pandas':
            res = df.query('a in b', parser=parser, engine=engine)
            expec = df[df.a.isin(df.b)]
            assert_frame_equal(res, expec)

            res = df.query('a in b and c < d', parser=parser, engine=engine)
            expec = df[df.a.isin(df.b) & (df.c < df.d)]
            assert_frame_equal(res, expec)
        else:
            with assertRaises(NotImplementedError):
                df.query('a in b', parser=parser, engine=engine)

            with assertRaises(NotImplementedError):
                df.query('a in b and c < d', parser=parser, engine=engine)

    def test_query_with_string_columns(self):
        for parser, engine in product(PARSERS, ENGINES):
            yield self.check_query_with_string_columns, parser, engine

    def check_object_array_eq_ne(self, parser, engine):
        tm.skip_if_no_ne(engine)
        df = DataFrame({'a': list('aaaabbbbcccc'),
                        'b': list('aabbccddeeff'),
                        'c': np.random.randint(5, size=12),
                        'd': np.random.randint(9, size=12)})
        res = df.query('a == b', parser=parser, engine=engine)
        exp = df[df.a == df.b]
        assert_frame_equal(res, exp)

        res = df.query('a != b', parser=parser, engine=engine)
        exp = df[df.a != df.b]
        assert_frame_equal(res, exp)

    def test_object_array_eq_ne(self):
        for parser, engine in product(PARSERS, ENGINES):
            yield self.check_object_array_eq_ne, parser, engine

    def check_query_with_nested_strings(self, parser, engine):
        tm.skip_if_no_ne(engine)
        skip_if_no_pandas_parser(parser)
        from pandas.compat import StringIO
        raw = """id          event          timestamp
        1   "page 1 load"   1/1/2014 0:00:01
        1   "page 1 exit"   1/1/2014 0:00:31
        2   "page 2 load"   1/1/2014 0:01:01
        2   "page 2 exit"   1/1/2014 0:01:31
        3   "page 3 load"   1/1/2014 0:02:01
        3   "page 3 exit"   1/1/2014 0:02:31
        4   "page 1 load"   2/1/2014 1:00:01
        4   "page 1 exit"   2/1/2014 1:00:31
        5   "page 2 load"   2/1/2014 1:01:01
        5   "page 2 exit"   2/1/2014 1:01:31
        6   "page 3 load"   2/1/2014 1:02:01
        6   "page 3 exit"   2/1/2014 1:02:31
        """
        df = pd.read_csv(StringIO(raw), sep=r'\s{2,}', engine='python',
                         parse_dates=['timestamp'])
        expected = df[df.event == '"page 1 load"']
        res = df.query("""'"page 1 load"' in event""", parser=parser,
                       engine=engine)
        tm.assert_frame_equal(expected, res)

    def test_query_with_nested_string(self):
        for parser, engine in product(PARSERS, ENGINES):
            yield self.check_query_with_nested_strings, parser, engine

    def check_query_with_nested_special_character(self, parser, engine):
        skip_if_no_pandas_parser(parser)
        tm.skip_if_no_ne(engine)
        df = DataFrame({'a': ['a', 'b', 'test & test'],
                        'b': [1, 2, 3]})
        res = df.query('a == "test & test"', parser=parser, engine=engine)
        expec = df[df.a == 'test & test']
        tm.assert_frame_equal(res, expec)

    def test_query_with_nested_special_character(self):
        for parser, engine in product(PARSERS, ENGINES):
            yield self.check_query_with_nested_special_character, parser, engine

    def check_query_lex_compare_strings(self, parser, engine):
        tm.skip_if_no_ne(engine=engine)
        import operator as opr

        a = Series(tm.choice(list('abcde'), 20))
        b = Series(np.arange(a.size))
        df = DataFrame({'X': a, 'Y': b})

        ops = {'<': opr.lt, '>': opr.gt, '<=': opr.le, '>=': opr.ge}

        for op, func in ops.items():
            res = df.query('X %s "d"' % op, engine=engine, parser=parser)
            expected = df[func(df.X, 'd')]
            assert_frame_equal(res, expected)

    def test_query_lex_compare_strings(self):
        for parser, engine in product(PARSERS, ENGINES):
            yield self.check_query_lex_compare_strings, parser, engine

    def check_query_single_element_booleans(self, parser, engine):
        tm.skip_if_no_ne(engine)
        columns = 'bid', 'bidsize', 'ask', 'asksize'
        data = np.random.randint(2, size=(1, len(columns))).astype(bool)
        df = DataFrame(data, columns=columns)
        res = df.query('bid & ask', engine=engine, parser=parser)
        expected = df[df.bid & df.ask]
        assert_frame_equal(res, expected)

    def test_query_single_element_booleans(self):
        for parser, engine in product(PARSERS, ENGINES):
            yield self.check_query_single_element_booleans, parser, engine

    def check_query_string_scalar_variable(self, parser, engine):
        tm.skip_if_no_ne(engine)
        df = pd.DataFrame({'Symbol': ['BUD US', 'BUD US', 'IBM US', 'IBM US'],
                           'Price': [109.70, 109.72, 183.30, 183.35]})
        e = df[df.Symbol == 'BUD US']
        symb = 'BUD US'
        r = df.query('Symbol == @symb', parser=parser, engine=engine)
        tm.assert_frame_equal(e, r)

    def test_query_string_scalar_variable(self):
        for parser, engine in product(['pandas'], ENGINES):
            yield self.check_query_string_scalar_variable, parser, engine


class TestDataFrameEvalNumExprPandas(tm.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestDataFrameEvalNumExprPandas, cls).setUpClass()
        cls.engine = 'numexpr'
        cls.parser = 'pandas'
        tm.skip_if_no_ne()

    def setUp(self):
        self.frame = DataFrame(randn(10, 3), columns=list('abc'))

    def tearDown(self):
        del self.frame

    def test_simple_expr(self):
        res = self.frame.eval('a + b', engine=self.engine, parser=self.parser)
        expect = self.frame.a + self.frame.b
        assert_series_equal(res, expect)

    def test_bool_arith_expr(self):
        res = self.frame.eval('a[a < 1] + b', engine=self.engine,
                              parser=self.parser)
        expect = self.frame.a[self.frame.a < 1] + self.frame.b
        assert_series_equal(res, expect)

    def test_invalid_type_for_operator_raises(self):
        df = DataFrame({'a': [1, 2], 'b': ['c', 'd']})
        ops = '+', '-', '*', '/'
        for op in ops:
            with tm.assertRaisesRegexp(TypeError,
                                       "unsupported operand type\(s\) for "
                                       ".+: '.+' and '.+'"):
                df.eval('a {0} b'.format(op), engine=self.engine,
                        parser=self.parser)


class TestDataFrameEvalNumExprPython(TestDataFrameEvalNumExprPandas):

    @classmethod
    def setUpClass(cls):
        super(TestDataFrameEvalNumExprPython, cls).setUpClass()
        cls.engine = 'numexpr'
        cls.parser = 'python'
        tm.skip_if_no_ne(cls.engine)


class TestDataFrameEvalPythonPandas(TestDataFrameEvalNumExprPandas):

    @classmethod
    def setUpClass(cls):
        super(TestDataFrameEvalPythonPandas, cls).setUpClass()
        cls.engine = 'python'
        cls.parser = 'pandas'


class TestDataFrameEvalPythonPython(TestDataFrameEvalNumExprPython):

    @classmethod
    def setUpClass(cls):
        super(TestDataFrameEvalPythonPython, cls).tearDownClass()
        cls.engine = cls.parser = 'python'


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
