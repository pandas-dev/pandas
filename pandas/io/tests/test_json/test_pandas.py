
# pylint: disable-msg=W0612,E1101
from copy import deepcopy
from datetime import datetime, timedelta
from StringIO import StringIO
import cPickle as pickle
import operator
import os
import unittest

import nose
import numpy as np

from pandas import Series, DataFrame, DatetimeIndex, Timestamp
import pandas as pd
read_json = pd.read_json

from pandas.util.testing import (assert_almost_equal, assert_frame_equal,
                                 assert_series_equal, network,
                                 ensure_clean)
import pandas.util.testing as tm
from numpy.testing.decorators import slow

_seriesd = tm.getSeriesData()
_tsd = tm.getTimeSeriesData()

_frame = DataFrame(_seriesd)
_frame2 = DataFrame(_seriesd, columns=['D', 'C', 'B', 'A'])
_intframe = DataFrame(dict((k, v.astype(int))
                           for k, v in _seriesd.iteritems()))

_tsframe = DataFrame(_tsd)

_mixed_frame = _frame.copy()

class TestPandasObjects(unittest.TestCase):

    def setUp(self):
        self.ts = tm.makeTimeSeries()
        self.ts.name = 'ts'

        self.series = tm.makeStringSeries()
        self.series.name = 'series'

        self.objSeries = tm.makeObjectSeries()
        self.objSeries.name = 'objects'

        self.empty_series = Series([], index=[])
        self.empty_frame = DataFrame({})

        self.frame = _frame.copy()
        self.frame2 = _frame2.copy()
        self.intframe = _intframe.copy()
        self.tsframe = _tsframe.copy()
        self.mixed_frame = _mixed_frame.copy()

    def test_frame_from_json_to_json(self):

        def _check_orient(df, orient, dtype=None, numpy=True):
            df = df.sort()
            dfjson = df.to_json(orient=orient)
            unser = read_json(dfjson, orient=orient, dtype=dtype,
                              numpy=numpy)
            unser = unser.sort()
            if df.index.dtype.type == np.datetime64:
                unser.index = DatetimeIndex(unser.index.values.astype('i8'))
            if orient == "records":
                # index is not captured in this orientation
                assert_almost_equal(df.values, unser.values)
                self.assert_(df.columns.equals(unser.columns))
            elif orient == "values":
                # index and cols are not captured in this orientation
                assert_almost_equal(df.values, unser.values)
            elif orient == "split":
                # index and col labels might not be strings
                unser.index = [str(i) for i in unser.index]
                unser.columns = [str(i) for i in unser.columns]
                unser = unser.sort()
                assert_almost_equal(df.values, unser.values)
            else:
                assert_frame_equal(df, unser)

        def _check_all_orients(df, dtype=None):
            _check_orient(df, "columns", dtype=dtype)
            _check_orient(df, "records", dtype=dtype)
            _check_orient(df, "split", dtype=dtype)
            _check_orient(df, "index", dtype=dtype)
            _check_orient(df, "values", dtype=dtype)

            _check_orient(df, "columns", dtype=dtype, numpy=False)
            _check_orient(df, "records", dtype=dtype, numpy=False)
            _check_orient(df, "split", dtype=dtype, numpy=False)
            _check_orient(df, "index", dtype=dtype, numpy=False)
            _check_orient(df, "values", dtype=dtype, numpy=False)

        # basic
        _check_all_orients(self.frame)
        self.assertEqual(self.frame.to_json(),
                         self.frame.to_json(orient="columns"))

        _check_all_orients(self.intframe, dtype=self.intframe.values.dtype)

        # big one
        # index and columns are strings as all unserialised JSON object keys
        # are assumed to be strings
        biggie = DataFrame(np.zeros((200, 4)),
                           columns=[str(i) for i in range(4)],
                           index=[str(i) for i in range(200)])
        _check_all_orients(biggie)

        # dtypes
        _check_all_orients(DataFrame(biggie, dtype=np.float64),
                           dtype=np.float64)
        _check_all_orients(DataFrame(biggie, dtype=np.int), dtype=np.int)
        _check_all_orients(DataFrame(biggie, dtype='<U3'), dtype='<U3')

        # empty
        _check_all_orients(self.empty_frame)

        # time series data
        _check_all_orients(self.tsframe)

        # mixed data
        index = pd.Index(['a', 'b', 'c', 'd', 'e'])
        data = {
            'A': [0., 1., 2., 3., 4.],
            'B': [0., 1., 0., 1., 0.],
            'C': ['foo1', 'foo2', 'foo3', 'foo4', 'foo5'],
            'D': [True, False, True, False, True]
        }
        df = DataFrame(data=data, index=index)
        _check_orient(df, "split")
        _check_orient(df, "records")
        _check_orient(df, "values")
        _check_orient(df, "columns")
        # index oriented is problematic as it is read back in in a transposed
        # state, so the columns are interpreted as having mixed data and
        # given object dtypes.
        # force everything to have object dtype beforehand
        _check_orient(df.transpose().transpose(), "index")

    def test_frame_from_json_bad_data(self):
        self.assertRaises(ValueError, read_json, StringIO('{"key":b:a:d}'))

        # too few indices
        json = StringIO('{"columns":["A","B"],'
                        '"index":["2","3"],'
                        '"data":[[1.0,"1"],[2.0,"2"],[null,"3"]]}"')
        self.assertRaises(ValueError, read_json, json,
                          orient="split")

        # too many columns
        json = StringIO('{"columns":["A","B","C"],'
                        '"index":["1","2","3"],'
                        '"data":[[1.0,"1"],[2.0,"2"],[null,"3"]]}"')
        self.assertRaises(AssertionError, read_json, json,
                          orient="split")

        # bad key
        json = StringIO('{"badkey":["A","B"],'
                        '"index":["2","3"],'
                        '"data":[[1.0,"1"],[2.0,"2"],[null,"3"]]}"')
        self.assertRaises(TypeError, read_json, json,
                          orient="split")

    def test_frame_from_json_nones(self):
        df = DataFrame([[1, 2], [4, 5, 6]])
        unser = read_json(df.to_json())
        self.assert_(np.isnan(unser['2'][0]))

        df = DataFrame([['1', '2'], ['4', '5', '6']])
        unser = read_json(df.to_json())
        self.assert_(unser['2'][0] is None)

        unser = read_json(df.to_json(), numpy=False)
        self.assert_(unser['2'][0] is None)

        # infinities get mapped to nulls which get mapped to NaNs during
        # deserialisation
        df = DataFrame([[1, 2], [4, 5, 6]])
        df[2][0] = np.inf
        unser = read_json(df.to_json())
        self.assert_(np.isnan(unser['2'][0]))

        df[2][0] = np.NINF
        unser = read_json(df.to_json())
        self.assert_(np.isnan(unser['2'][0]))

    def test_frame_to_json_except(self):
        df = DataFrame([1, 2, 3])
        self.assertRaises(ValueError, df.to_json, orient="garbage")

    def test_series_from_json_to_json(self):

        def _check_orient(series, orient, dtype=None, numpy=True):
            series = series.sort_index()
            unser = read_json(series.to_json(orient=orient), typ='series',
                              orient=orient, numpy=numpy, dtype=dtype)
            unser = unser.sort_index()
            if series.index.dtype.type == np.datetime64:
                unser.index = DatetimeIndex(unser.index.values.astype('i8'))
            if orient == "records" or orient == "values":
                assert_almost_equal(series.values, unser.values)
            else:
                try:
                    assert_series_equal(series, unser)
                except:
                    raise
                if orient == "split":
                    self.assert_(series.name == unser.name)

        def _check_all_orients(series, dtype=None):
            _check_orient(series, "columns", dtype=dtype)
            _check_orient(series, "records", dtype=dtype)
            _check_orient(series, "split", dtype=dtype)
            _check_orient(series, "index", dtype=dtype)
            _check_orient(series, "values", dtype=dtype)

            _check_orient(series, "columns", dtype=dtype, numpy=False)
            _check_orient(series, "records", dtype=dtype, numpy=False)
            _check_orient(series, "split", dtype=dtype, numpy=False)
            _check_orient(series, "index", dtype=dtype, numpy=False)
            _check_orient(series, "values", dtype=dtype, numpy=False)

        # basic
        _check_all_orients(self.series)
        self.assertEqual(self.series.to_json(),
                         self.series.to_json(orient="index"))

        objSeries = Series([str(d) for d in self.objSeries],
                           index=self.objSeries.index,
                           name=self.objSeries.name)
        _check_all_orients(objSeries)
        _check_all_orients(self.empty_series)
        _check_all_orients(self.ts)

        # dtype
        s = Series(range(6), index=['a','b','c','d','e','f'])
        _check_all_orients(Series(s, dtype=np.float64), dtype=np.float64)
        _check_all_orients(Series(s, dtype=np.int), dtype=np.int)

    def test_series_to_json_except(self):
        s = Series([1, 2, 3])
        self.assertRaises(ValueError, s.to_json, orient="garbage")

    def test_typ(self):

        s = Series(range(6), index=['a','b','c','d','e','f'])
        result = read_json(s.to_json(),typ=None)
        assert_series_equal(result,s)

    def test_reconstruction_index(self):

        df = DataFrame([[1, 2, 3], [4, 5, 6]])
        result = read_json(df.to_json())

        # the index is serialized as strings....correct?
        #assert_frame_equal(result,df)

    def test_path(self):
        with ensure_clean('test.json') as path:

            for df in [ self.frame, self.frame2, self.intframe, self.tsframe, self.mixed_frame ]:
                df.to_json(path)
                read_json(path)

    def test_axis_dates(self):

        # frame
        json = self.tsframe.to_json()
        result = read_json(json)
        assert_frame_equal(result,self.tsframe)

        # series
        json = self.ts.to_json()
        result = read_json(json,typ='series')
        assert_series_equal(result,self.ts)

    def test_parse_dates(self):

        # frame
        df = self.tsframe.copy()
        df['date'] = Timestamp('20130101')

        json = df.to_json()
        result = read_json(json,parse_dates=True)
        assert_frame_equal(result,df)

        df['foo'] = 1.
        json = df.to_json()
        result = read_json(json,parse_dates=True)
        assert_frame_equal(result,df)

        # series
        ts = Series(Timestamp('20130101'),index=self.ts.index)
        json = ts.to_json()
        result = read_json(json,typ='series',parse_dates=True)
        assert_series_equal(result,ts)

    def test_date_format(self):
        
        df = self.tsframe.copy()
        df['date'] = Timestamp('20130101')
        df_orig = df.copy()

        json = df.to_json(date_format='iso')
        result = read_json(json,parse_dates=True)
        assert_frame_equal(result,df_orig)

        # make sure that we did in fact copy
        assert_frame_equal(df,df_orig)

        ts = Series(Timestamp('20130101'),index=self.ts.index)
        json = ts.to_json(date_format='iso')
        result = read_json(json,typ='series',parse_dates=True)
        assert_series_equal(result,ts)

    def test_weird_nested_json(self):

        # this used to core dump the parser
        s = r'''{
        "status": "success",
        "data": {
        "posts": [
            {
            "id": 1,
            "title": "A blog post",
            "body": "Some useful content"
            },
            {
            "id": 2,
            "title": "Another blog post",
            "body": "More content"
            }
        ]
    }
}'''

        read_json(s)

    @network
    @slow
    def test_round_trip_exception_(self):
		# GH 3867

        df = pd.read_csv('https://raw.github.com/hayd/lahman2012/master/csvs/Teams.csv')
        s = df.to_json()
        result = pd.read_json(s)
        result.index = result.index.astype(int)
        result = result.reindex(columns=df.columns,index=df.index)
        assert_frame_equal(result,df)

    @network
    @slow
    def test_url(self):
        import urllib2
        try:

            url = 'https://api.github.com/repos/pydata/pandas/issues?per_page=5'
            result = read_json(url,parse_dates=True)
            for c in ['created_at','closed_at','updated_at']:
                self.assert_(result[c].dtype == 'datetime64[ns]')
            
            url = 'http://search.twitter.com/search.json?q=pandas%20python'
            result = read_json(url)
            
        except urllib2.URLError:
            raise nose.SkipTest
