# pylint: disable-msg=E1101,W0612

from datetime import datetime, timedelta
import operator
import nose

import numpy as np
from numpy import nan
import pandas as pd

from pandas import (Index, Series, DataFrame, Panel,
                    isnull, notnull,date_range, _np_version_under1p7)
from pandas.core.index import Index, MultiIndex
from pandas.tseries.index import Timestamp, DatetimeIndex

import pandas.core.common as com

from pandas.compat import StringIO, lrange, range, zip, u, OrderedDict, long
from pandas import compat
from pandas.util.testing import (assert_series_equal,
                                 assert_frame_equal,
                                 assert_panel_equal,
                                 assert_almost_equal,
                                 ensure_clean)
import pandas.util.testing as tm


def _skip_if_no_scipy():
    try:
        import scipy.interpolate
    except ImportError:
        raise nose.SkipTest('scipy.interpolate missing')


def _skip_if_no_pchip():
    try:
        from scipy.interpolate import pchip_interpolate
    except ImportError:
        raise nose.SkipTest('scipy.interpolate.pchip missing')

#------------------------------------------------------------------------------
# Generic types test cases


class Generic(object):

    _multiprocess_can_split_ = True

    def setUp(self):
        import warnings
        warnings.filterwarnings(action='ignore', category=FutureWarning)

    @property
    def _ndim(self):
        return self._typ._AXIS_LEN

    def _axes(self):
        """ return the axes for my object typ """
        return self._typ._AXIS_ORDERS

    def _construct(self, shape, value=None, dtype=None, **kwargs):
        """ construct an object for the given shape
            if value is specified use that if its a scalar
            if value is an array, repeat it as needed """

        if isinstance(shape,int):
            shape = tuple([shape] * self._ndim)
        if value is not None:
            if np.isscalar(value):
                if value == 'empty':
                    arr = None

                    # remove the info axis
                    kwargs.pop(self._typ._info_axis_name,None)
                else:
                    arr = np.empty(shape,dtype=dtype)
                    arr.fill(value)
            else:
                fshape = np.prod(shape)
                arr = value.ravel()
                new_shape = fshape/arr.shape[0]
                if fshape % arr.shape[0] != 0:
                    raise Exception("invalid value passed in _construct")

                arr = np.repeat(arr,new_shape).reshape(shape)
        else:
            arr = np.random.randn(*shape)
        return self._typ(arr,dtype=dtype,**kwargs)

    def _compare(self, result, expected):
        self._comparator(result,expected)

    def test_rename(self):

        # single axis
        for axis in self._axes():
            kwargs = { axis : list('ABCD') }
            obj = self._construct(4,**kwargs)

            # no values passed
            #self.assertRaises(Exception, o.rename(str.lower))

            # rename a single axis
            result = obj.rename(**{ axis : str.lower })
            expected = obj.copy()
            setattr(expected,axis,list('abcd'))
            self._compare(result, expected)

        # multiple axes at once

    def test_get_numeric_data(self):

        n = 4
        kwargs = { }
        for i in range(self._ndim):
            kwargs[self._typ._AXIS_NAMES[i]] = list(range(n))

        # get the numeric data
        o = self._construct(n,**kwargs)
        result = o._get_numeric_data()
        self._compare(result, o)

        # non-inclusion
        result = o._get_bool_data()
        expected = self._construct(n,value='empty',**kwargs)
        self._compare(result,expected)

        # get the bool data
        arr = np.array([True,True,False,True])
        o = self._construct(n,value=arr,**kwargs)
        result = o._get_numeric_data()
        self._compare(result, o)

        # _get_numeric_data is includes _get_bool_data, so can't test for non-inclusion

    def test_nonzero(self):

        # GH 4633
        # look at the boolean/nonzero behavior for objects
        obj = self._construct(shape=4)
        self.assertRaises(ValueError, lambda : bool(obj == 0))
        self.assertRaises(ValueError, lambda : bool(obj == 1))
        self.assertRaises(ValueError, lambda : bool(obj))

        obj = self._construct(shape=4,value=1)
        self.assertRaises(ValueError, lambda : bool(obj == 0))
        self.assertRaises(ValueError, lambda : bool(obj == 1))
        self.assertRaises(ValueError, lambda : bool(obj))

        obj = self._construct(shape=4,value=np.nan)
        self.assertRaises(ValueError, lambda : bool(obj == 0))
        self.assertRaises(ValueError, lambda : bool(obj == 1))
        self.assertRaises(ValueError, lambda : bool(obj))

        # empty
        obj = self._construct(shape=0)
        self.assertRaises(ValueError, lambda : bool(obj))

        # invalid behaviors

        obj1 = self._construct(shape=4,value=1)
        obj2 = self._construct(shape=4,value=1)

        def f():
            if obj1:
                print("this works and shouldn't")
        self.assertRaises(ValueError, f)
        self.assertRaises(ValueError, lambda : obj1 and obj2)
        self.assertRaises(ValueError, lambda : obj1 or obj2)
        self.assertRaises(ValueError, lambda : not obj1)

    def test_numpy_1_7_compat_numeric_methods(self):
        if _np_version_under1p7:
            raise nose.SkipTest("numpy < 1.7")

        # GH 4435
        # numpy in 1.7 tries to pass addtional arguments to pandas functions

        o = self._construct(shape=4)
        for op in ['min','max','max','var','std','prod','sum','cumsum','cumprod',
                   'median','skew','kurt','compound','cummax','cummin','all','any']:
            f = getattr(np,op,None)
            if f is not None:
                f(o)

    def test_downcast(self):
        # test close downcasting

        o = self._construct(shape=4, value=9, dtype=np.int64)
        result = o.copy()
        result._data = o._data.downcast(dtypes='infer')
        self._compare(result, o)

        o = self._construct(shape=4, value=9.)
        expected = o.astype(np.int64)
        result = o.copy()
        result._data = o._data.downcast(dtypes='infer')
        self._compare(result, expected)

        o = self._construct(shape=4, value=9.5)
        result = o.copy()
        result._data = o._data.downcast(dtypes='infer')
        self._compare(result, o)

        # are close
        o = self._construct(shape=4, value=9.000000000005)
        result = o.copy()
        result._data = o._data.downcast(dtypes='infer')
        expected = o.astype(np.int64)
        self._compare(result, expected)

    def test_constructor_compound_dtypes(self):
        # GH 5191
        # compound dtypes should raise not-implementederror

        def f(dtype):
            return self._construct(shape=3, dtype=dtype)

        self.assertRaises(NotImplementedError, f, [("A","datetime64[h]"), ("B","str"), ("C","int32")])

        # these work (though results may be unexpected)
        f('int64')
        f('float64')
        f('M8[ns]')

    def check_metadata(self, x, y=None):
        for m in x._metadata:
            v = getattr(x,m,None)
            if y is None:
                self.assert_(v is None)
            else:
                self.assert_(v == getattr(y,m,None))

    def test_metadata_propagation(self):
        # check that the metadata matches up on the resulting ops

        o = self._construct(shape=3)
        o.name = 'foo'
        o2 = self._construct(shape=3)
        o2.name = 'bar'

        # TODO
        # Once panel can do non-trivial combine operations
        # (currently there is an a raise in the Panel arith_ops to prevent
        # this, though it actually does work)
        # can remove all of these try: except: blocks on the actual operations


        # ----------
        # preserving
        # ----------

        # simple ops with scalars
        for op in [ '__add__','__sub__','__truediv__','__mul__' ]:
            result = getattr(o,op)(1)
            self.check_metadata(o,result)

        # ops with like
        for op in [ '__add__','__sub__','__truediv__','__mul__' ]:
            try:
                result = getattr(o,op)(o)
                self.check_metadata(o,result)
            except (ValueError, AttributeError):
                pass

        # simple boolean
        for op in [ '__eq__','__le__', '__ge__' ]:
            v1 = getattr(o,op)(o)
            self.check_metadata(o,v1)

            try:
                self.check_metadata(o, v1 & v1)
            except (ValueError):
                pass

            try:
                self.check_metadata(o, v1 | v1)
            except (ValueError):
                pass

        # combine_first
        try:
            result = o.combine_first(o2)
            self.check_metadata(o,result)
        except (AttributeError):
            pass

        # ---------------------------
        # non-preserving (by default)
        # ---------------------------

        # add non-like
        try:
            result = o + o2
            self.check_metadata(result)
        except (ValueError, AttributeError):
            pass

        # simple boolean
        for op in [ '__eq__','__le__', '__ge__' ]:

            # this is a name matching op
            v1 = getattr(o,op)(o)

            v2 = getattr(o,op)(o2)
            self.check_metadata(v2)

            try:
                self.check_metadata(v1 & v2)
            except (ValueError):
                pass

            try:
                self.check_metadata(v1 | v2)
            except (ValueError):
                pass

    def test_head_tail(self):
        # GH5370

        o = self._construct(shape=10)

        # check all index types
        for index in [ tm.makeFloatIndex, tm.makeIntIndex,
                       tm.makeStringIndex, tm.makeUnicodeIndex,
                       tm.makeDateIndex, tm.makePeriodIndex ]:
            axis = o._get_axis_name(0)
            setattr(o,axis,index(len(getattr(o,axis))))

            # Panel + dims
            try:
                o.head()
            except (NotImplementedError):
                raise nose.SkipTest('not implemented on {0}'.format(o.__class__.__name__))

            self._compare(o.head(), o.iloc[:5])
            self._compare(o.tail(), o.iloc[-5:])

            # 0-len
            self._compare(o.head(0), o.iloc[:])
            self._compare(o.tail(0), o.iloc[0:])

            # bounded
            self._compare(o.head(len(o)+1), o)
            self._compare(o.tail(len(o)+1), o)

            # neg index
            self._compare(o.head(-3), o.head(7))
            self._compare(o.tail(-3), o.tail(7))

class TestSeries(tm.TestCase, Generic):
    _typ = Series
    _comparator = lambda self, x, y: assert_series_equal(x,y)

    def setUp(self):
        self.ts = tm.makeTimeSeries()  # Was at top level in test_series
        self.ts.name = 'ts'

        self.series = tm.makeStringSeries()
        self.series.name = 'series'

    def test_rename_mi(self):
        s = Series([11,21,31],
                   index=MultiIndex.from_tuples([("A",x) for x in ["a","B","c"]]))
        result = s.rename(str.lower)

    def test_get_numeric_data_preserve_dtype(self):

        # get the numeric data
        o = Series([1,2,3])
        result = o._get_numeric_data()
        self._compare(result, o)

        o = Series([1,'2',3.])
        result = o._get_numeric_data()
        expected = Series([],dtype=object)
        self._compare(result, expected)

        o = Series([True,False,True])
        result = o._get_numeric_data()
        self._compare(result, o)

        o = Series([True,False,True])
        result = o._get_bool_data()
        self._compare(result, o)

        o = Series(date_range('20130101',periods=3))
        result = o._get_numeric_data()
        expected = Series([],dtype='M8[ns]')
        self._compare(result, expected)

    def test_nonzero_single_element(self):

        # allow single item via bool method
        s = Series([True])
        self.assert_(s.bool() is True)

        s = Series([False])
        self.assert_(s.bool() is False)

        # single item nan to raise
        for s in [ Series([np.nan]), Series([pd.NaT]), Series([True]), Series([False]) ]:
            self.assertRaises(ValueError, lambda : bool(s))

        for s in [ Series([np.nan]), Series([pd.NaT])]:
            self.assertRaises(ValueError, lambda : s.bool())

        # multiple bool are still an error
        for s in [Series([True,True]), Series([False, False])]:
            self.assertRaises(ValueError, lambda : bool(s))
            self.assertRaises(ValueError, lambda : s.bool())

        # single non-bool are an error
        for s in [Series([1]), Series([0]),
                  Series(['a']), Series([0.0])]:
                self.assertRaises(ValueError, lambda : bool(s))
                self.assertRaises(ValueError, lambda : s.bool())

    def test_metadata_propagation_indiv(self):
        # check that the metadata matches up on the resulting ops

        o = Series(range(3),range(3))
        o.name = 'foo'
        o2 = Series(range(3),range(3))
        o2.name = 'bar'

        result = o.T
        self.check_metadata(o,result)

        # resample
        ts = Series(np.random.rand(1000),
                    index=date_range('20130101',periods=1000,freq='s'),
                    name='foo')
        result = ts.resample('1T')
        self.check_metadata(ts,result)

        result = ts.resample('1T',how='min')
        self.check_metadata(ts,result)

        result = ts.resample('1T',how=lambda x: x.sum())
        self.check_metadata(ts,result)

    def test_interpolate(self):
        ts = Series(np.arange(len(self.ts), dtype=float), self.ts.index)

        ts_copy = ts.copy()
        ts_copy[5:10] = np.NaN

        linear_interp = ts_copy.interpolate(method='linear')
        self.assert_(np.array_equal(linear_interp, ts))

        ord_ts = Series([d.toordinal() for d in self.ts.index],
                        index=self.ts.index).astype(float)

        ord_ts_copy = ord_ts.copy()
        ord_ts_copy[5:10] = np.NaN

        time_interp = ord_ts_copy.interpolate(method='time')
        self.assert_(np.array_equal(time_interp, ord_ts))

        # try time interpolation on a non-TimeSeries
        self.assertRaises(ValueError, self.series.interpolate, method='time')

    def test_interp_regression(self):
        _skip_if_no_scipy()
        _skip_if_no_pchip()

        ser = Series(np.sort(np.random.uniform(size=100)))

        # interpolate at new_index
        new_index = ser.index + Index([49.25, 49.5, 49.75, 50.25, 50.5, 50.75])
        interp_s = ser.reindex(new_index).interpolate(method='pchip')
        # does not blow up, GH5977
        interp_s[49:51]

    def test_interpolate_corners(self):
        s = Series([np.nan, np.nan])
        assert_series_equal(s.interpolate(), s)

        s = Series([]).interpolate()
        assert_series_equal(s.interpolate(), s)

        _skip_if_no_scipy()
        s = Series([np.nan, np.nan])
        assert_series_equal(s.interpolate(method='polynomial', order=1), s)

        s = Series([]).interpolate()
        assert_series_equal(s.interpolate(method='polynomial', order=1), s)

    def test_interpolate_index_values(self):
        s = Series(np.nan, index=np.sort(np.random.rand(30)))
        s[::3] = np.random.randn(10)

        vals = s.index.values.astype(float)

        result = s.interpolate(method='values')

        expected = s.copy()
        bad = isnull(expected.values)
        good = -bad
        expected = Series(
            np.interp(vals[bad], vals[good], s.values[good]), index=s.index[bad])

        assert_series_equal(result[bad], expected)

    def test_interpolate_non_ts(self):
        s = Series([1, 3, np.nan, np.nan, np.nan, 11])
        with tm.assertRaises(ValueError):
            s.interpolate(method='time')

    # New interpolation tests
    def test_nan_interpolate(self):
        s = Series([0, 1, np.nan, 3])
        result = s.interpolate()
        expected = Series([0, 1, 2, 3])
        assert_series_equal(result, expected)

        _skip_if_no_scipy()
        result = s.interpolate(method='polynomial', order=1)
        assert_series_equal(result, expected)

    def test_nan_irregular_index(self):
        s = Series([1, 2, np.nan, 4], index=[1, 3, 5, 9])
        result = s.interpolate()
        expected = Series([1, 2, 3, 4], index=[1, 3, 5, 9])
        assert_series_equal(result, expected)

    def test_nan_str_index(self):
        s = Series([0, 1, 2, np.nan], index=list('abcd'))
        result = s.interpolate()
        expected = Series([0, 1, 2, 2], index=list('abcd'))
        assert_series_equal(result, expected)

    def test_interp_quad(self):
        _skip_if_no_scipy()
        sq = Series([1, 4, np.nan, 16], index=[1, 2, 3, 4])
        result = sq.interpolate(method='quadratic')
        expected = Series([1, 4, 9, 16], index=[1, 2, 3, 4])
        assert_series_equal(result, expected)

    def test_interp_scipy_basic(self):
        _skip_if_no_scipy()
        s = Series([1, 3, np.nan, 12, np.nan, 25])
        # slinear
        expected = Series([1., 3., 7.5, 12., 18.5, 25.])
        result = s.interpolate(method='slinear')
        assert_series_equal(result, expected)
        # nearest
        expected = Series([1, 3, 3, 12, 12, 25])
        result = s.interpolate(method='nearest')
        assert_series_equal(result, expected)
        # zero
        expected = Series([1, 3, 3, 12, 12, 25])
        result = s.interpolate(method='zero')
        assert_series_equal(result, expected)
        # quadratic
        expected = Series([1, 3., 6.769231, 12., 18.230769, 25.])
        result = s.interpolate(method='quadratic')
        assert_series_equal(result, expected)
        # cubic
        expected = Series([1., 3., 6.8, 12., 18.2, 25.])
        result = s.interpolate(method='cubic')
        assert_series_equal(result, expected)

    def test_interp_limit(self):
        s = Series([1, 3, np.nan, np.nan, np.nan, 11])
        expected = Series([1., 3., 5., 7., np.nan, 11.])
        result = s.interpolate(method='linear', limit=2)
        assert_series_equal(result, expected)

    def test_interp_all_good(self):
        # scipy
        _skip_if_no_scipy()
        s = Series([1, 2, 3])
        result = s.interpolate(method='polynomial', order=1)
        assert_series_equal(result, s)

        # non-scipy
        result = s.interpolate()
        assert_series_equal(result, s)

    def test_interp_multiIndex(self):
        idx = MultiIndex.from_tuples([(0, 'a'), (1, 'b'), (2, 'c')])
        s = Series([1, 2, np.nan], index=idx)

        expected = s.copy()
        expected.loc[2] = 2
        expected = expected.astype(np.int64)
        result = s.interpolate()
        assert_series_equal(result, expected)

        _skip_if_no_scipy()
        with tm.assertRaises(ValueError):
            s.interpolate(method='polynomial', order=1)

    def test_interp_nonmono_raise(self):
        _skip_if_no_scipy()
        s = pd.Series([1, 2, 3], index=[0, 2, 1])
        with tm.assertRaises(ValueError):
            s.interpolate(method='krogh')

    def test_interp_datetime64(self):
        _skip_if_no_scipy()
        df = Series([1, np.nan, 3], index=date_range('1/1/2000', periods=3))
        result = df.interpolate(method='nearest')
        expected = Series([1, 1, 3], index=date_range('1/1/2000', periods=3))
        assert_series_equal(result, expected)

class TestDataFrame(tm.TestCase, Generic):
    _typ = DataFrame
    _comparator = lambda self, x, y: assert_frame_equal(x,y)

    def test_rename_mi(self):
        df = DataFrame([11,21,31],
                       index=MultiIndex.from_tuples([("A",x) for x in ["a","B","c"]]))
        result = df.rename(str.lower)

    def test_nonzero_single_element(self):

        # allow single item via bool method
        df = DataFrame([[True]])
        self.assert_(df.bool() is True)

        df = DataFrame([[False]])
        self.assert_(df.bool() is False)

        df = DataFrame([[False, False]])
        self.assertRaises(ValueError, lambda : df.bool())
        self.assertRaises(ValueError, lambda : bool(df))

    def test_get_numeric_data_preserve_dtype(self):

        # get the numeric data
        o = DataFrame({'A': [1, '2', 3.]})
        result = o._get_numeric_data()
        expected = DataFrame(index=[0, 1, 2], dtype=object)
        self._compare(result, expected)

    def test_interp_basic(self):
        df = DataFrame({'A': [1, 2, np.nan, 4], 'B': [1, 4, 9, np.nan],
                        'C': [1, 2, 3, 5], 'D': list('abcd')})
        expected = DataFrame({'A': [1, 2, 3, 4], 'B': [1, 4, 9, 9],
                              'C': [1, 2, 3, 5], 'D': list('abcd')})
        result = df.interpolate()
        assert_frame_equal(result, expected)

        result = df.set_index('C').interpolate()
        expected = df.set_index('C')
        expected.A.loc[3] = 3
        expected.B.loc[5] = 9
        expected[['A', 'B']] = expected[['A', 'B']].astype(np.int64)

        assert_frame_equal(result, expected)

    def test_interp_bad_method(self):
        df = DataFrame({'A': [1, 2, np.nan, 4], 'B': [1, 4, 9, np.nan],
                        'C': [1, 2, 3, 5], 'D': list('abcd')})
        with tm.assertRaises(ValueError):
            df.interpolate(method='not_a_method')

    def test_interp_combo(self):
        df = DataFrame({'A': [1., 2., np.nan, 4.], 'B': [1, 4, 9, np.nan],
                        'C': [1, 2, 3, 5], 'D': list('abcd')})

        result = df['A'].interpolate()
        expected = Series([1, 2, 3, 4])
        assert_series_equal(result, expected)

    def test_interp_nan_idx(self):
        df = DataFrame({'A': [1, 2, np.nan, 4], 'B': [np.nan, 2, 3, 4]})
        df = df.set_index('A')
        with tm.assertRaises(NotImplementedError):
            df.interpolate(method='values')

    def test_interp_various(self):
        _skip_if_no_scipy()
        df = DataFrame({'A': [1, 2, np.nan, 4, 5, np.nan, 7],
                        'C': [1, 2, 3, 5, 8, 13, 21]})
        df = df.set_index('C')
        expected = df.copy()
        result = df.interpolate(method='polynomial', order=1)

        expected.A.loc[3] = 2.66666667
        expected.A.loc[13] = 5.76923076
        assert_frame_equal(result, expected)

        result = df.interpolate(method='cubic')
        expected.A.loc[3] = 2.81621174
        expected.A.loc[13] = 5.64146581
        assert_frame_equal(result, expected)

        result = df.interpolate(method='nearest')
        expected.A.loc[3] = 2
        expected.A.loc[13] = 5
        assert_frame_equal(result, expected, check_dtype=False)

        result = df.interpolate(method='quadratic')
        expected.A.loc[3] = 2.82533638
        expected.A.loc[13] = 6.02817974
        assert_frame_equal(result, expected)

        result = df.interpolate(method='slinear')
        expected.A.loc[3] = 2.66666667
        expected.A.loc[13] = 5.76923077
        assert_frame_equal(result, expected)

        result = df.interpolate(method='zero')
        expected.A.loc[3] = 2.
        expected.A.loc[13] = 5
        assert_frame_equal(result, expected, check_dtype=False)

        result = df.interpolate(method='quadratic')
        expected.A.loc[3] = 2.82533638
        expected.A.loc[13] = 6.02817974
        assert_frame_equal(result, expected)

    def test_interp_alt_scipy(self):
        _skip_if_no_scipy()
        df = DataFrame({'A': [1, 2, np.nan, 4, 5, np.nan, 7],
                        'C': [1, 2, 3, 5, 8, 13, 21]})
        result = df.interpolate(method='barycentric')
        expected = df.copy()
        expected['A'].iloc[2] = 3
        expected['A'].iloc[5] = 6
        assert_frame_equal(result, expected.astype(np.int64))

        result = df.interpolate(method='krogh')
        expectedk = df.copy()
        expectedk['A'].iloc[2] = 3
        expectedk['A'].iloc[5] = 6
        expectedk['A'] = expected['A'].astype(np.int64)
        assert_frame_equal(result, expectedk)

        _skip_if_no_pchip()
        result = df.interpolate(method='pchip')
        expected['A'].iloc[2] = 3
        expected['A'].iloc[5] = 6.125
        assert_frame_equal(result, expected)

    def test_interp_rowwise(self):
        df = DataFrame({0: [1, 2, np.nan, 4],
                        1: [2, 3, 4, np.nan],
                        2: [np.nan, 4, 5, 6],
                        3: [4, np.nan, 6, 7],
                        4: [1, 2, 3, 4]})
        result = df.interpolate(axis=1)
        expected = df.copy()
        expected[1].loc[3] = 5
        expected[2].loc[0] = 3
        expected[3].loc[1] = 3
        expected[4] = expected[4].astype(np.float64)
        assert_frame_equal(result, expected)

        # scipy route
        _skip_if_no_scipy()
        result = df.interpolate(axis=1, method='values')
        assert_frame_equal(result, expected)

        result = df.interpolate(axis=0)
        expected = df.interpolate()
        assert_frame_equal(result, expected)

    def test_rowwise_alt(self):
        df = DataFrame({0: [0, .5, 1., np.nan, 4, 8, np.nan, np.nan, 64],
                        1: [1, 2, 3, 4, 3, 2, 1, 0, -1]})
        df.interpolate(axis=0)

    def test_interp_leading_nans(self):
        df = DataFrame({"A": [np.nan, np.nan, .5, .25, 0],
                        "B": [np.nan, -3, -3.5, np.nan, -4]})
        result = df.interpolate()
        expected = df.copy()
        expected['B'].loc[3] = -3.75
        assert_frame_equal(result, expected)

        _skip_if_no_scipy()
        result = df.interpolate(method='polynomial', order=1)
        assert_frame_equal(result, expected)

    def test_interp_raise_on_only_mixed(self):
        df = DataFrame({'A': [1, 2, np.nan, 4], 'B': ['a', 'b', 'c', 'd'],
                        'C': [np.nan, 2, 5, 7], 'D': [np.nan, np.nan, 9, 9],
                        'E': [1, 2, 3, 4]})
        with tm.assertRaises(TypeError):
            df.interpolate(axis=1)

    def test_no_order(self):
        _skip_if_no_scipy()
        s = Series([0, 1, np.nan, 3])
        with tm.assertRaises(ValueError):
            s.interpolate(method='polynomial')
        with tm.assertRaises(ValueError):
            s.interpolate(method='spline')

    def test_spline(self):
        _skip_if_no_scipy()
        s = Series([1, 2, np.nan, 4, 5, np.nan, 7])
        result = s.interpolate(method='spline', order=1)
        expected = Series([1, 2, 3, 4, 5, 6, 7])
        assert_series_equal(result, expected)

    def test_metadata_propagation_indiv(self):

        # groupby
        df = DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                              'foo', 'bar', 'foo', 'foo'],
                        'B': ['one', 'one', 'two', 'three',
                              'two', 'two', 'one', 'three'],
                        'C': np.random.randn(8),
                        'D': np.random.randn(8)})
        result = df.groupby('A').sum()
        self.check_metadata(df,result)

        # resample
        df = DataFrame(np.random.randn(1000,2),
                       index=date_range('20130101',periods=1000,freq='s'))
        result = df.resample('1T')
        self.check_metadata(df,result)

class TestPanel(tm.TestCase, Generic):
    _typ = Panel
    _comparator = lambda self, x, y: assert_panel_equal(x, y)


class TestNDFrame(tm.TestCase):
    # tests that don't fit elsewhere

    def test_squeeze(self):
        # noop
        for s in [ tm.makeFloatSeries(), tm.makeStringSeries(), tm.makeObjectSeries() ]:
            tm.assert_series_equal(s.squeeze(),s)
        for df in [ tm.makeTimeDataFrame() ]:
            tm.assert_frame_equal(df.squeeze(),df)
        for p in [ tm.makePanel() ]:
            tm.assert_panel_equal(p.squeeze(),p)
        for p4d in [ tm.makePanel4D() ]:
            tm.assert_panel4d_equal(p4d.squeeze(),p4d)

        # squeezing
        df = tm.makeTimeDataFrame().reindex(columns=['A'])
        tm.assert_series_equal(df.squeeze(),df['A'])

        p = tm.makePanel().reindex(items=['ItemA'])
        tm.assert_frame_equal(p.squeeze(),p['ItemA'])

        p = tm.makePanel().reindex(items=['ItemA'],minor_axis=['A'])
        tm.assert_series_equal(p.squeeze(),p.ix['ItemA',:,'A'])

        p4d = tm.makePanel4D().reindex(labels=['label1'])
        tm.assert_panel_equal(p4d.squeeze(),p4d['label1'])

        p4d = tm.makePanel4D().reindex(labels=['label1'],items=['ItemA'])
        tm.assert_frame_equal(p4d.squeeze(),p4d.ix['label1','ItemA'])

    def test_equals(self):
        s1 = pd.Series([1, 2, 3], index=[0, 2, 1])
        s2 = s1.copy()
        self.assert_(s1.equals(s2))

        s1[1] = 99
        self.assert_(not s1.equals(s2))

        # NaNs compare as equal
        s1 = pd.Series([1, np.nan, 3, np.nan], index=[0, 2, 1, 3])
        s2 = s1.copy()
        self.assert_(s1.equals(s2))

        s2[0] = 9.9
        self.assert_(not s1.equals(s2))
        
        idx = MultiIndex.from_tuples([(0, 'a'), (1, 'b'), (2, 'c')])
        s1 = Series([1, 2, np.nan], index=idx)
        s2 = s1.copy()
        self.assert_(s1.equals(s2))

        # Add object dtype column with nans
        index = np.random.random(10)
        df1 = DataFrame(np.random.random(10,), index=index, columns=['floats'])
        df1['text'] = 'the sky is so blue. we could use more chocolate.'.split()
        df1['start'] = date_range('2000-1-1', periods=10, freq='T')
        df1['end'] = date_range('2000-1-1', periods=10, freq='D')
        df1['diff'] = df1['end'] - df1['start']
        df1['bool'] = (np.arange(10) % 3 == 0)
        df1.ix[::2] = nan
        df2 = df1.copy()
        self.assert_(df1['text'].equals(df2['text']))
        self.assert_(df1['start'].equals(df2['start']))
        self.assert_(df1['end'].equals(df2['end']))
        self.assert_(df1['diff'].equals(df2['diff']))
        self.assert_(df1['bool'].equals(df2['bool']))
        self.assert_(df1.equals(df2))
        self.assert_(not df1.equals(object))

        # different dtype
        different = df1.copy()
        different['floats'] = different['floats'].astype('float32')
        self.assert_(not df1.equals(different)) 

        # different index
        different_index = -index
        different = df2.set_index(different_index)
        self.assert_(not df1.equals(different))        

        # different columns
        different = df2.copy()
        different.columns = df2.columns[::-1]
        self.assert_(not df1.equals(different))        

        # DatetimeIndex
        index = pd.date_range('2000-1-1', periods=10, freq='T')
        df1 = df1.set_index(index)
        df2 = df1.copy()
        self.assert_(df1.equals(df2))

        # MultiIndex
        df3 = df1.set_index(['text'], append=True)
        df2 = df1.set_index(['text'], append=True)
        self.assert_(df3.equals(df2))

        df2 = df1.set_index(['floats'], append=True)
        self.assert_(not df3.equals(df2))

        # NaN in index
        df3 = df1.set_index(['floats'], append=True)
        df2 = df1.set_index(['floats'], append=True)
        self.assert_(df3.equals(df2))


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
