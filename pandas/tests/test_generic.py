# -*- coding: utf-8 -*-
# pylint: disable-msg=E1101,W0612

from datetime import datetime, timedelta
import nose
import numpy as np
from numpy import nan
import pandas as pd

from pandas import (Index, Series, DataFrame, Panel,
                    isnull, notnull, date_range, period_range)
from pandas.core.index import Index, MultiIndex

import pandas.core.common as com

from pandas.compat import StringIO, lrange, range, zip, u, OrderedDict, long
from pandas import compat
from pandas.util.testing import (assert_series_equal,
                                 assert_frame_equal,
                                 assert_panel_equal,
                                 assert_almost_equal,
                                 assert_equal,
                                 ensure_clean)
import pandas.util.testing as tm


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
        pass

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

    def test_get_default(self):

        # GH 7725
        d0 = "a", "b", "c", "d"
        d1 = np.arange(4, dtype='int64')
        others = "e", 10

        for data, index in ((d0, d1), (d1, d0)):
            s = Series(data, index=index)
            for i,d in zip(index, data):
                self.assertEqual(s.get(i), d)
                self.assertEqual(s.get(i, d), d)
                self.assertEqual(s.get(i, "z"), d)
                for other in others:
                    self.assertEqual(s.get(other, "z"), "z")
                    self.assertEqual(s.get(other, other), other)

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
                com.pprint_thing("this works and shouldn't")
        self.assertRaises(ValueError, f)
        self.assertRaises(ValueError, lambda : obj1 and obj2)
        self.assertRaises(ValueError, lambda : obj1 or obj2)
        self.assertRaises(ValueError, lambda : not obj1)

    def test_numpy_1_7_compat_numeric_methods(self):
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
                self.assertIsNone(v)
            else:
                self.assertEqual(v, getattr(y,m,None))

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

    def test_sample(self):
        # Fixes issue: 2419

        o = self._construct(shape=10)

        ###
        # Check behavior of random_state argument
        ###

        # Check for stability when receives seed or random state -- run 10 times.
        for test in range(10):
            seed = np.random.randint(0,100)
            self._compare(o.sample(n=4, random_state=seed), o.sample(n=4, random_state=seed))
            self._compare(o.sample(frac=0.7,random_state=seed), o.sample(frac=0.7, random_state=seed))

            self._compare(o.sample(n=4, random_state=np.random.RandomState(test)),
                          o.sample(n=4, random_state=np.random.RandomState(test)))

            self._compare(o.sample(frac=0.7,random_state=np.random.RandomState(test)),
                          o.sample(frac=0.7, random_state=np.random.RandomState(test)))


        # Check for error when random_state argument invalid.
        with tm.assertRaises(ValueError):
            o.sample(random_state='astring!')

        ###
        # Check behavior of `frac` and `N`
        ###

        # Giving both frac and N throws error
        with tm.assertRaises(ValueError):
            o.sample(n=3, frac=0.3)

        # Check that raises right error for negative lengths
        with tm.assertRaises(ValueError):
            o.sample(n=-3)
        with tm.assertRaises(ValueError):
            o.sample(frac=-0.3)

        # Make sure float values of `n` give error
        with tm.assertRaises(ValueError):
            o.sample(n= 3.2)

        # Check lengths are right
        self.assertTrue(len(o.sample(n=4) == 4))
        self.assertTrue(len(o.sample(frac=0.34) == 3))
        self.assertTrue(len(o.sample(frac=0.36) == 4))

        ###
        # Check weights
        ###

        # Weight length must be right
        with tm.assertRaises(ValueError):
            o.sample(n=3, weights=[0,1])

        with tm.assertRaises(ValueError):
            bad_weights = [0.5]*11
            o.sample(n=3, weights=bad_weights)

        with tm.assertRaises(ValueError):
            bad_weight_series = Series([0,0,0.2])
            o.sample(n=4, weights=bad_weight_series)

        # Check won't accept negative weights
        with tm.assertRaises(ValueError):
            bad_weights = [-0.1]*10
            o.sample(n=3, weights=bad_weights)

        # Check inf and -inf throw errors:
        with tm.assertRaises(ValueError):
            weights_with_inf = [0.1]*10
            weights_with_inf[0] = np.inf
            o.sample(n=3, weights=weights_with_inf)

        with tm.assertRaises(ValueError):
            weights_with_ninf = [0.1]*10
            weights_with_ninf[0] =  -np.inf
            o.sample(n=3, weights=weights_with_ninf)

        # All zeros raises errors
        zero_weights = [0]*10
        with tm.assertRaises(ValueError):
            o.sample(n=3, weights=zero_weights)

        # All missing weights
        nan_weights = [np.nan]*10
        with tm.assertRaises(ValueError):
            o.sample(n=3, weights=nan_weights)


        # A few dataframe test with degenerate weights.
        easy_weight_list = [0]*10
        easy_weight_list[5] = 1

        df = pd.DataFrame({'col1':range(10,20),
                           'col2':range(20,30),
                           'colString': ['a']*10,
                           'easyweights':easy_weight_list})
        sample1 = df.sample(n=1, weights='easyweights')
        assert_frame_equal(sample1, df.iloc[5:6])

        # Ensure proper error if string given as weight for Series, panel, or
        # DataFrame with axis = 1.
        s = Series(range(10))
        with tm.assertRaises(ValueError):
            s.sample(n=3, weights='weight_column')

        panel = pd.Panel(items = [0,1,2], major_axis = [2,3,4], minor_axis = [3,4,5])
        with tm.assertRaises(ValueError):
            panel.sample(n=1, weights='weight_column')

        with tm.assertRaises(ValueError):
            df.sample(n=1, weights='weight_column', axis = 1)

        # Check weighting key error
        with tm.assertRaises(KeyError):
            df.sample(n=3, weights='not_a_real_column_name')

         # Check np.nan are replaced by zeros.
        weights_with_nan = [np.nan]*10
        weights_with_nan[5] = 0.5
        self._compare(o.sample(n=1, axis=0, weights=weights_with_nan), o.iloc[5:6])

        # Check None are also replaced by zeros.
        weights_with_None = [None]*10
        weights_with_None[5] = 0.5
        self._compare(o.sample(n=1, axis=0, weights=weights_with_None), o.iloc[5:6])

        # Check that re-normalizes weights that don't sum to one.
        weights_less_than_1 = [0]*10
        weights_less_than_1[0] = 0.5
        tm.assert_frame_equal(df.sample(n=1, weights=weights_less_than_1), df.iloc[:1])


        ###
        # Test axis argument
        ###

        # Test axis argument
        df = pd.DataFrame({'col1':range(10), 'col2':['a']*10})
        second_column_weight = [0,1]
        assert_frame_equal(df.sample(n=1, axis=1, weights=second_column_weight), df[['col2']])

        # Different axis arg types
        assert_frame_equal(df.sample(n=1, axis='columns', weights=second_column_weight),
                           df[['col2']])

        weight = [0]*10
        weight[5] = 0.5
        assert_frame_equal(df.sample(n=1, axis='rows', weights=weight),
                           df.iloc[5:6])
        assert_frame_equal(df.sample(n=1, axis='index', weights=weight),
                           df.iloc[5:6])

        # Check out of range axis values
        with tm.assertRaises(ValueError):
            df.sample(n=1, axis=2)

        with tm.assertRaises(ValueError):
            df.sample(n=1, axis='not_a_name')

        with tm.assertRaises(ValueError):
            s = pd.Series(range(10))
            s.sample(n=1, axis=1)

        # Test weight length compared to correct axis
        with tm.assertRaises(ValueError):
            df.sample(n=1, axis=1, weights=[0.5]*10)

        # Check weights with axis = 1
        easy_weight_list = [0]*3
        easy_weight_list[2] = 1

        df = pd.DataFrame({'col1':range(10,20),
                           'col2':range(20,30),
                           'colString': ['a']*10})
        sample1 = df.sample(n=1, axis=1, weights=easy_weight_list)
        assert_frame_equal(sample1, df[['colString']])

        # Test default axes
        p = pd.Panel(items = ['a','b','c'], major_axis=[2,4,6], minor_axis=[1,3,5])
        assert_panel_equal(p.sample(n=3, random_state=42), p.sample(n=3, axis=1, random_state=42))
        assert_frame_equal(df.sample(n=3, random_state=42), df.sample(n=3, axis=0, random_state=42))

        # Test that function aligns weights with frame
        df = DataFrame({'col1':[5,6,7], 'col2':['a','b','c'], }, index = [9,5,3])
        s = Series([1,0,0], index=[3,5,9])
        assert_frame_equal(df.loc[[3]], df.sample(1, weights=s))

        # Weights have index values to be dropped because not in
        # sampled DataFrame
        s2 = Series([0.001,0,10000], index=[3,5,10])
        assert_frame_equal(df.loc[[3]], df.sample(1, weights=s2))

        # Weights have empty values to be filed with zeros
        s3 = Series([0.01,0], index=[3,5])
        assert_frame_equal(df.loc[[3]], df.sample(1, weights=s3))

        # No overlap in weight and sampled DataFrame indices
        s4 = Series([1,0], index=[1,2])
        with tm.assertRaises(ValueError):
            df.sample(1, weights=s4)


    def test_size_compat(self):
        # GH8846
        # size property should be defined

        o = self._construct(shape=10)
        self.assertTrue(o.size == np.prod(o.shape))
        self.assertTrue(o.size == 10**len(o.axes))

    def test_split_compat(self):
        # xref GH8846
        o = self._construct(shape=10)
        self.assertTrue(len(np.array_split(o,5)) == 5)
        self.assertTrue(len(np.array_split(o,2)) == 2)

    def test_unexpected_keyword(self):  # GH8597
        from pandas.util.testing import assertRaisesRegexp

        df = DataFrame(np.random.randn(5, 2), columns=['jim', 'joe'])
        ca = pd.Categorical([0, 0, 2, 2, 3, np.nan])
        ts = df['joe'].copy()
        ts[2] = np.nan

        with assertRaisesRegexp(TypeError, 'unexpected keyword'):
            df.drop('joe', axis=1, in_place=True)

        with assertRaisesRegexp(TypeError, 'unexpected keyword'):
            df.reindex([1, 0], inplace=True)

        with assertRaisesRegexp(TypeError, 'unexpected keyword'):
            ca.fillna(0, inplace=True)

        with assertRaisesRegexp(TypeError, 'unexpected keyword'):
            ts.fillna(0, in_place=True)

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
        o = Series([1, 2, 3])
        result = o._get_numeric_data()
        self._compare(result, o)

        o = Series([1, '2', 3.])
        result = o._get_numeric_data()
        expected = Series([], dtype=object, index=pd.Index([], dtype=object))
        self._compare(result, expected)

        o = Series([True, False, True])
        result = o._get_numeric_data()
        self._compare(result, o)

        o = Series([True, False, True])
        result = o._get_bool_data()
        self._compare(result, o)

        o = Series(date_range('20130101',periods=3))
        result = o._get_numeric_data()
        expected = Series([],dtype='M8[ns]', index=pd.Index([], dtype=object))
        self._compare(result, expected)

    def test_nonzero_single_element(self):

        # allow single item via bool method
        s = Series([True])
        self.assertTrue(s.bool())

        s = Series([False])
        self.assertFalse(s.bool())

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

        _metadata = Series._metadata
        _finalize = Series.__finalize__
        Series._metadata = ['name','filename']
        o.filename = 'foo'
        o2.filename = 'bar'

        def finalize(self, other, method=None, **kwargs):
            for name in self._metadata:
                if method == 'concat' and name == 'filename':
                    value = '+'.join([ getattr(o,name) for o in other.objs if getattr(o,name,None) ])
                    object.__setattr__(self, name, value)
                else:
                    object.__setattr__(self, name, getattr(other, name, None))

            return self

        Series.__finalize__ = finalize

        result = pd.concat([o, o2])
        self.assertEqual(result.filename,'foo+bar')
        self.assertIsNone(result.name)

        # reset
        Series._metadata = _metadata
        Series.__finalize__ = _finalize

    def test_interpolate(self):
        ts = Series(np.arange(len(self.ts), dtype=float), self.ts.index)

        ts_copy = ts.copy()
        ts_copy[5:10] = np.NaN

        linear_interp = ts_copy.interpolate(method='linear')
        self.assert_numpy_array_equal(linear_interp, ts)

        ord_ts = Series([d.toordinal() for d in self.ts.index],
                        index=self.ts.index).astype(float)

        ord_ts_copy = ord_ts.copy()
        ord_ts_copy[5:10] = np.NaN

        time_interp = ord_ts_copy.interpolate(method='time')
        self.assert_numpy_array_equal(time_interp, ord_ts)

        # try time interpolation on a non-TimeSeries
        # Only raises ValueError if there are NaNs.
        non_ts = self.series.copy()
        non_ts[0] = np.NaN
        self.assertRaises(ValueError, non_ts.interpolate, method='time')

    def test_interp_regression(self):
        tm._skip_if_no_scipy()
        _skip_if_no_pchip()

        ser = Series(np.sort(np.random.uniform(size=100)))

        # interpolate at new_index
        new_index = ser.index.union(Index([49.25, 49.5, 49.75, 50.25, 50.5, 50.75]))
        interp_s = ser.reindex(new_index).interpolate(method='pchip')
        # does not blow up, GH5977
        interp_s[49:51]

    def test_interpolate_corners(self):
        s = Series([np.nan, np.nan])
        assert_series_equal(s.interpolate(), s)

        s = Series([]).interpolate()
        assert_series_equal(s.interpolate(), s)

        tm._skip_if_no_scipy()
        s = Series([np.nan, np.nan])
        assert_series_equal(s.interpolate(method='polynomial', order=1), s)

        s = Series([]).interpolate()
        assert_series_equal(s.interpolate(method='polynomial', order=1), s)

    def test_interpolate_index_values(self):
        s = Series(np.nan, index=np.sort(np.random.rand(30)))
        s[::3] = np.random.randn(10)

        vals = s.index.values.astype(float)

        result = s.interpolate(method='index')

        expected = s.copy()
        bad = isnull(expected.values)
        good = ~bad
        expected = Series(
            np.interp(vals[bad], vals[good], s.values[good]), index=s.index[bad])

        assert_series_equal(result[bad], expected)

        # 'values' is synonymous with 'index' for the method kwarg
        other_result = s.interpolate(method='values')

        assert_series_equal(other_result, result)
        assert_series_equal(other_result[bad], expected)

    def test_interpolate_non_ts(self):
        s = Series([1, 3, np.nan, np.nan, np.nan, 11])
        with tm.assertRaises(ValueError):
            s.interpolate(method='time')

    # New interpolation tests
    def test_nan_interpolate(self):
        s = Series([0, 1, np.nan, 3])
        result = s.interpolate()
        expected = Series([0., 1., 2., 3.])
        assert_series_equal(result, expected)

        tm._skip_if_no_scipy()
        result = s.interpolate(method='polynomial', order=1)
        assert_series_equal(result, expected)

    def test_nan_irregular_index(self):
        s = Series([1, 2, np.nan, 4], index=[1, 3, 5, 9])
        result = s.interpolate()
        expected = Series([1., 2., 3., 4.], index=[1, 3, 5, 9])
        assert_series_equal(result, expected)

    def test_nan_str_index(self):
        s = Series([0, 1, 2, np.nan], index=list('abcd'))
        result = s.interpolate()
        expected = Series([0., 1., 2., 2.], index=list('abcd'))
        assert_series_equal(result, expected)

    def test_interp_quad(self):
        tm._skip_if_no_scipy()
        sq = Series([1, 4, np.nan, 16], index=[1, 2, 3, 4])
        result = sq.interpolate(method='quadratic')
        expected = Series([1., 4., 9., 16.], index=[1, 2, 3, 4])
        assert_series_equal(result, expected)

    def test_interp_scipy_basic(self):
        tm._skip_if_no_scipy()
        s = Series([1, 3, np.nan, 12, np.nan, 25])
        # slinear
        expected = Series([1., 3., 7.5, 12., 18.5, 25.])
        result = s.interpolate(method='slinear')
        assert_series_equal(result, expected)

        result = s.interpolate(method='slinear', downcast='infer')
        assert_series_equal(result, expected)
        # nearest
        expected = Series([1, 3, 3, 12, 12, 25])
        result = s.interpolate(method='nearest')
        assert_series_equal(result, expected.astype('float'))

        result = s.interpolate(method='nearest', downcast='infer')
        assert_series_equal(result, expected)
        # zero
        expected = Series([1, 3, 3, 12, 12, 25])
        result = s.interpolate(method='zero')
        assert_series_equal(result, expected.astype('float'))

        result = s.interpolate(method='zero', downcast='infer')
        assert_series_equal(result, expected)
        # quadratic
        expected = Series([1, 3., 6.769231, 12., 18.230769, 25.])
        result = s.interpolate(method='quadratic')
        assert_series_equal(result, expected)

        result = s.interpolate(method='quadratic', downcast='infer')
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

    def test_interp_limit_forward(self):
        s = Series([1, 3, np.nan, np.nan, np.nan, 11])

        # Provide 'forward' (the default) explicitly here.
        expected = Series([1., 3., 5., 7., np.nan, 11.])

        result = s.interpolate(
            method='linear', limit=2, limit_direction='forward')
        assert_series_equal(result, expected)

        result = s.interpolate(
            method='linear', limit=2, limit_direction='FORWARD')
        assert_series_equal(result, expected)

    def test_interp_limit_bad_direction(self):
        s = Series([1, 3, np.nan, np.nan, np.nan, 11])

        self.assertRaises(ValueError, s.interpolate,
                          method='linear', limit=2,
                          limit_direction='abc')

        # raises an error even if no limit is specified.
        self.assertRaises(ValueError, s.interpolate,
                          method='linear',
                          limit_direction='abc')

    def test_interp_limit_direction(self):
        # These tests are for issue #9218 -- fill NaNs in both directions.
        s = Series([1, 3, np.nan, np.nan, np.nan, 11])

        expected = Series([1., 3., np.nan, 7., 9., 11.])
        result = s.interpolate(
            method='linear', limit=2, limit_direction='backward')
        assert_series_equal(result, expected)

        expected = Series([1., 3., 5., np.nan, 9., 11.])
        result = s.interpolate(
            method='linear', limit=1, limit_direction='both')
        assert_series_equal(result, expected)

        # Check that this works on a longer series of nans.
        s = Series([1, 3, np.nan, np.nan, np.nan, 7, 9, np.nan, np.nan, 12, np.nan])

        expected = Series([1., 3., 4., 5., 6., 7., 9., 10., 11., 12., 12.])
        result = s.interpolate(
            method='linear', limit=2, limit_direction='both')
        assert_series_equal(result, expected)

        expected = Series([1., 3., 4., np.nan, 6., 7., 9., 10., 11., 12., 12.])
        result = s.interpolate(
            method='linear', limit=1, limit_direction='both')
        assert_series_equal(result, expected)

    def test_interp_limit_to_ends(self):
        # These test are for issue #10420 -- flow back to beginning.
        s = Series([np.nan, np.nan, 5, 7, 9, np.nan])

        expected = Series([5., 5., 5., 7., 9., np.nan])
        result = s.interpolate(
            method='linear', limit=2, limit_direction='backward')
        assert_series_equal(result, expected)

        expected = Series([5., 5., 5., 7., 9., 9.])
        result = s.interpolate(
            method='linear', limit=2, limit_direction='both')
        assert_series_equal(result, expected)

    def test_interp_limit_before_ends(self):
        # These test are for issue #11115 -- limit ends properly.
        s = Series([np.nan, np.nan, 5, 7, np.nan, np.nan])

        expected = Series([np.nan, np.nan, 5., 7., 7., np.nan])
        result = s.interpolate(
            method='linear', limit=1, limit_direction='forward')
        assert_series_equal(result, expected)

        expected = Series([np.nan, 5., 5., 7., np.nan, np.nan])
        result = s.interpolate(
            method='linear', limit=1, limit_direction='backward')
        assert_series_equal(result, expected)

        expected = Series([np.nan, 5., 5., 7., 7., np.nan])
        result = s.interpolate(
            method='linear', limit=1, limit_direction='both')
        assert_series_equal(result, expected)

    def test_interp_all_good(self):
        # scipy
        tm._skip_if_no_scipy()
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
        result = s.interpolate()
        assert_series_equal(result, expected)

        tm._skip_if_no_scipy()
        with tm.assertRaises(ValueError):
            s.interpolate(method='polynomial', order=1)

    def test_interp_nonmono_raise(self):
        tm._skip_if_no_scipy()
        s = Series([1, np.nan, 3], index=[0, 2, 1])
        with tm.assertRaises(ValueError):
            s.interpolate(method='krogh')

    def test_interp_datetime64(self):
        tm._skip_if_no_scipy()
        df = Series([1, np.nan, 3], index=date_range('1/1/2000', periods=3))
        result = df.interpolate(method='nearest')
        expected = Series([1., 1., 3.], index=date_range('1/1/2000', periods=3))
        assert_series_equal(result, expected)

    def test_interp_limit_no_nans(self):
        # GH 7173
        s = pd.Series([1., 2., 3.])
        result = s.interpolate(limit=1)
        expected = s
        assert_series_equal(result, expected)

    def test_describe(self):
        _ = self.series.describe()
        _ = self.ts.describe()

    def test_describe_objects(self):
        s = Series(['a', 'b', 'b', np.nan, np.nan, np.nan, 'c', 'd', 'a', 'a'])
        result = s.describe()
        expected = Series({'count': 7, 'unique': 4,
                           'top': 'a', 'freq': 3}, index=result.index)
        assert_series_equal(result, expected)

        dt = list(self.ts.index)
        dt.append(dt[0])
        ser = Series(dt)
        rs = ser.describe()
        min_date = min(dt)
        max_date = max(dt)
        xp = Series({'count': len(dt),
                     'unique': len(self.ts.index),
                     'first': min_date, 'last': max_date, 'freq': 2,
                     'top': min_date}, index=rs.index)
        assert_series_equal(rs, xp)

    def test_describe_empty(self):
        result = pd.Series().describe()

        self.assertEqual(result['count'], 0)
        self.assertTrue(result.drop('count').isnull().all())

        nanSeries = Series([np.nan])
        nanSeries.name = 'NaN'
        result = nanSeries.describe()
        self.assertEqual(result['count'], 0)
        self.assertTrue(result.drop('count').isnull().all())

    def test_describe_none(self):
        noneSeries = Series([None])
        noneSeries.name = 'None'
        expected = Series([0, 0], index=['count', 'unique'], name='None')
        assert_series_equal(noneSeries.describe(), expected)


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
        self.assertTrue(df.bool())

        df = DataFrame([[False]])
        self.assertFalse(df.bool())

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
        expected = DataFrame({'A': [1., 2., 3., 4.], 'B': [1., 4., 9., 9.],
                              'C': [1, 2, 3, 5], 'D': list('abcd')})
        result = df.interpolate()
        assert_frame_equal(result, expected)

        result = df.set_index('C').interpolate()
        expected = df.set_index('C')
        expected.loc[3,'A'] = 3
        expected.loc[5,'B'] = 9
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
        expected = Series([1., 2., 3., 4.], name='A')
        assert_series_equal(result, expected)

        result = df['A'].interpolate(downcast='infer')
        expected = Series([1, 2, 3, 4], name='A')
        assert_series_equal(result, expected)

    def test_interp_nan_idx(self):
        df = DataFrame({'A': [1, 2, np.nan, 4], 'B': [np.nan, 2, 3, 4]})
        df = df.set_index('A')
        with tm.assertRaises(NotImplementedError):
            df.interpolate(method='values')

    def test_interp_various(self):
        tm._skip_if_no_scipy()
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
        tm._skip_if_no_scipy()
        df = DataFrame({'A': [1, 2, np.nan, 4, 5, np.nan, 7],
                        'C': [1, 2, 3, 5, 8, 13, 21]})
        result = df.interpolate(method='barycentric')
        expected = df.copy()
        expected.ix[2,'A'] = 3
        expected.ix[5,'A'] = 6
        assert_frame_equal(result, expected)

        result = df.interpolate(method='barycentric', downcast='infer')
        assert_frame_equal(result, expected.astype(np.int64))

        result = df.interpolate(method='krogh')
        expectedk = df.copy()
        expectedk['A'] = expected['A']
        assert_frame_equal(result, expectedk)

        _skip_if_no_pchip()
        result = df.interpolate(method='pchip')
        expected.ix[2,'A'] = 3
        expected.ix[5,'A'] = 6.125
        assert_frame_equal(result, expected)

    def test_interp_rowwise(self):
        df = DataFrame({0: [1, 2, np.nan, 4],
                        1: [2, 3, 4, np.nan],
                        2: [np.nan, 4, 5, 6],
                        3: [4, np.nan, 6, 7],
                        4: [1, 2, 3, 4]})
        result = df.interpolate(axis=1)
        expected = df.copy()
        expected.loc[3,1] = 5
        expected.loc[0,2] = 3
        expected.loc[1,3] = 3
        expected[4] = expected[4].astype(np.float64)
        assert_frame_equal(result, expected)

        # scipy route
        tm._skip_if_no_scipy()
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

        tm._skip_if_no_scipy()
        result = df.interpolate(method='polynomial', order=1)
        assert_frame_equal(result, expected)

    def test_interp_raise_on_only_mixed(self):
        df = DataFrame({'A': [1, 2, np.nan, 4], 'B': ['a', 'b', 'c', 'd'],
                        'C': [np.nan, 2, 5, 7], 'D': [np.nan, np.nan, 9, 9],
                        'E': [1, 2, 3, 4]})
        with tm.assertRaises(TypeError):
            df.interpolate(axis=1)

    def test_interp_inplace(self):
        df = DataFrame({'a': [1., 2., np.nan, 4.]})
        expected = DataFrame({'a': [1., 2., 3., 4.]})
        result = df.copy()
        result['a'].interpolate(inplace=True)
        assert_frame_equal(result, expected)

        result = df.copy()
        result['a'].interpolate(inplace=True, downcast='infer')
        assert_frame_equal(result, expected.astype('int64'))

    def test_interp_inplace_row(self):
        # GH 10395
        result = DataFrame({'a': [1.,2.,3.,4.], 'b': [np.nan, 2., 3., 4.],
                            'c': [3, 2, 2, 2]})
        expected = result.interpolate(method='linear', axis=1, inplace=False)
        result.interpolate(method='linear', axis=1, inplace=True)
        assert_frame_equal(result, expected)

    def test_interp_ignore_all_good(self):
        # GH
        df = DataFrame({'A': [1, 2, np.nan, 4],
                        'B': [1, 2, 3, 4],
                        'C': [1., 2., np.nan, 4.],
                        'D': [1., 2., 3., 4.]})
        expected = DataFrame({'A': np.array([1, 2, 3, 4], dtype='float64'),
                              'B': np.array([1, 2, 3, 4], dtype='int64'),
                              'C': np.array([1., 2., 3, 4.], dtype='float64'),
                              'D': np.array([1., 2., 3., 4.], dtype='float64')})

        result = df.interpolate(downcast=None)
        assert_frame_equal(result, expected)

        # all good
        result = df[['B', 'D']].interpolate(downcast=None)
        assert_frame_equal(result, df[['B', 'D']])

    def test_describe(self):
        desc = tm.makeDataFrame().describe()
        desc = tm.makeMixedDataFrame().describe()
        desc = tm.makeTimeDataFrame().describe()

    def test_describe_percentiles_percent_or_raw(self):
        msg = 'percentiles should all be in the interval \\[0, 1\\]'

        df = tm.makeDataFrame()
        with tm.assertRaisesRegexp(ValueError, msg):
            df.describe(percentiles=[10, 50, 100])

        with tm.assertRaisesRegexp(ValueError, msg):
            df.describe(percentiles=[2])

        with tm.assertRaisesRegexp(ValueError, msg):
            df.describe(percentiles=[-2])

    def test_describe_percentiles_equivalence(self):
        df = tm.makeDataFrame()
        d1 = df.describe()
        d2 = df.describe(percentiles=[.25, .75])
        assert_frame_equal(d1, d2)

    def test_describe_percentiles_insert_median(self):
        df = tm.makeDataFrame()
        d1 = df.describe(percentiles=[.25, .75])
        d2 = df.describe(percentiles=[.25, .5, .75])
        assert_frame_equal(d1, d2)
        self.assertTrue('25%' in d1.index)
        self.assertTrue('75%' in d2.index)

        # none above
        d1 = df.describe(percentiles=[.25, .45])
        d2 = df.describe(percentiles=[.25, .45, .5])
        assert_frame_equal(d1, d2)
        self.assertTrue('25%' in d1.index)
        self.assertTrue('45%' in d2.index)

        # none below
        d1 = df.describe(percentiles=[.75, 1])
        d2 = df.describe(percentiles=[.5, .75, 1])
        assert_frame_equal(d1, d2)
        self.assertTrue('75%' in d1.index)
        self.assertTrue('100%' in d2.index)

        # edge
        d1 = df.describe(percentiles=[0, 1])
        d2 = df.describe(percentiles=[0, .5, 1])
        assert_frame_equal(d1, d2)
        self.assertTrue('0%' in d1.index)
        self.assertTrue('100%' in d2.index)

    def test_describe_no_numeric(self):
        df = DataFrame({'A': ['foo', 'foo', 'bar'] * 8,
                        'B': ['a', 'b', 'c', 'd'] * 6})
        desc = df.describe()
        expected = DataFrame(dict((k, v.describe())
                                  for k, v in compat.iteritems(df)),
                             columns=df.columns)
        assert_frame_equal(desc, expected)

        ts = tm.makeTimeSeries()
        df = DataFrame({'time': ts.index})
        desc = df.describe()
        self.assertEqual(desc.time['first'], min(ts.index))

    def test_describe_empty_int_columns(self):
        df = DataFrame([[0, 1], [1, 2]])
        desc = df[df[0] < 0].describe()  # works
        assert_series_equal(desc.xs('count'),
                            Series([0, 0], dtype=float, name='count'))
        self.assertTrue(isnull(desc.ix[1:]).all().all())

    def test_describe_objects(self):
        df = DataFrame({"C1": ['a', 'a', 'c'], "C2": ['d', 'd', 'f']})
        result = df.describe()
        expected = DataFrame({"C1": [3, 2, 'a', 2], "C2": [3, 2, 'd', 2]},
                             index=['count', 'unique', 'top', 'freq'])
        assert_frame_equal(result, expected)

        df = DataFrame({"C1": pd.date_range('2010-01-01', periods=4, freq='D')})
        df.loc[4] = pd.Timestamp('2010-01-04')
        result = df.describe()
        expected = DataFrame({"C1": [5, 4, pd.Timestamp('2010-01-04'), 2,
                                     pd.Timestamp('2010-01-01'),
                                     pd.Timestamp('2010-01-04')]},
                             index=['count', 'unique', 'top', 'freq',
                                    'first', 'last'])
        assert_frame_equal(result, expected)

        # mix time and str
        df['C2'] = ['a', 'a', 'b', 'c', 'a']
        result = df.describe()
        expected['C2'] = [5, 3, 'a', 3, np.nan, np.nan]
        assert_frame_equal(result, expected)

        # just str
        expected = DataFrame({'C2': [5, 3, 'a', 4]},
                             index=['count', 'unique', 'top', 'freq'])
        result = df[['C2']].describe()

        # mix of time, str, numeric
        df['C3'] = [2, 4, 6, 8, 2]
        result = df.describe()
        expected = DataFrame({"C3": [5., 4.4, 2.607681, 2., 2., 4., 6., 8.]},
                             index=['count', 'mean', 'std', 'min', '25%',
                                    '50%', '75%', 'max'])
        assert_frame_equal(result, expected)
        assert_frame_equal(df.describe(), df[['C3']].describe())

        assert_frame_equal(df[['C1', 'C3']].describe(), df[['C3']].describe())
        assert_frame_equal(df[['C2', 'C3']].describe(), df[['C3']].describe())

    def test_describe_typefiltering(self):
        df = DataFrame({'catA': ['foo', 'foo', 'bar'] * 8,
                        'catB': ['a', 'b', 'c', 'd'] * 6,
                        'numC': np.arange(24, dtype='int64'),
                        'numD': np.arange(24.) + .5,
                        'ts': tm.makeTimeSeries()[:24].index})

        descN = df.describe()
        expected_cols = ['numC', 'numD',]
        expected = DataFrame(dict((k, df[k].describe())
                                  for k in expected_cols),
                             columns=expected_cols)
        assert_frame_equal(descN, expected)

        desc = df.describe(include=['number'])
        assert_frame_equal(desc, descN)
        desc = df.describe(exclude=['object', 'datetime'])
        assert_frame_equal(desc, descN)
        desc = df.describe(include=['float'])
        assert_frame_equal(desc, descN.drop('numC',1))

        descC = df.describe(include=['O'])
        expected_cols = ['catA', 'catB']
        expected = DataFrame(dict((k, df[k].describe())
                                  for k in expected_cols),
                             columns=expected_cols)
        assert_frame_equal(descC, expected)

        descD = df.describe(include=['datetime'])
        assert_series_equal( descD.ts, df.ts.describe())

        desc = df.describe(include=['object','number', 'datetime'])
        assert_frame_equal(desc.loc[:,["numC","numD"]].dropna(), descN)
        assert_frame_equal(desc.loc[:,["catA","catB"]].dropna(), descC)
        descDs = descD.sort_index() # the index order change for mixed-types
        assert_frame_equal(desc.loc[:,"ts":].dropna().sort_index(), descDs)

        desc = df.loc[:,'catA':'catB'].describe(include='all')
        assert_frame_equal(desc, descC)
        desc = df.loc[:,'numC':'numD'].describe(include='all')
        assert_frame_equal(desc, descN)

        desc = df.describe(percentiles = [], include='all')
        cnt = Series(data=[4,4,6,6,6], index=['catA','catB','numC','numD','ts'])
        assert_series_equal( desc.count(), cnt)
        self.assertTrue('count' in desc.index)
        self.assertTrue('unique' in desc.index)
        self.assertTrue('50%' in desc.index)
        self.assertTrue('first' in desc.index)

        desc = df.drop("ts", 1).describe(percentiles = [], include='all')
        assert_series_equal( desc.count(), cnt.drop("ts"))
        self.assertTrue('first' not in desc.index)
        desc = df.drop(["numC","numD"], 1).describe(percentiles = [], include='all')
        assert_series_equal( desc.count(), cnt.drop(["numC","numD"]))
        self.assertTrue('50%' not in desc.index)

    def test_describe_typefiltering_category_bool(self):
        df = DataFrame({'A_cat': pd.Categorical(['foo', 'foo', 'bar'] * 8),
                        'B_str': ['a', 'b', 'c', 'd'] * 6,
                        'C_bool': [True] * 12 + [False] * 12,
                        'D_num': np.arange(24.) + .5,
                        'E_ts': tm.makeTimeSeries()[:24].index})

        # bool is considered numeric in describe, although not an np.number
        desc = df.describe()
        expected_cols = ['C_bool', 'D_num']
        expected = DataFrame(dict((k, df[k].describe())
                                  for k in expected_cols),
                             columns=expected_cols)
        assert_frame_equal(desc, expected)

        desc = df.describe(include=["category"])
        self.assertTrue(desc.columns.tolist() == ["A_cat"])

        # 'all' includes numpy-dtypes + category
        desc1 = df.describe(include="all")
        desc2 = df.describe(include=[np.generic, "category"])
        assert_frame_equal(desc1, desc2)

    def test_describe_timedelta(self):
        df = DataFrame({"td": pd.to_timedelta(np.arange(24)%20,"D")})
        self.assertTrue(df.describe().loc["mean"][0] == pd.to_timedelta("8d4h"))

    def test_describe_typefiltering_dupcol(self):
        df = DataFrame({'catA': ['foo', 'foo', 'bar'] * 8,
                        'catB': ['a', 'b', 'c', 'd'] * 6,
                        'numC': np.arange(24),
                        'numD': np.arange(24.) + .5,
                        'ts': tm.makeTimeSeries()[:24].index})
        s = df.describe(include='all').shape[1]
        df = pd.concat([df, df], axis=1)
        s2 = df.describe(include='all').shape[1]
        self.assertTrue(s2 == 2 * s)

    def test_describe_typefiltering_groupby(self):
        df = DataFrame({'catA': ['foo', 'foo', 'bar'] * 8,
                'catB': ['a', 'b', 'c', 'd'] * 6,
                'numC': np.arange(24),
                'numD': np.arange(24.) + .5,
                'ts': tm.makeTimeSeries()[:24].index})
        G = df.groupby('catA')
        self.assertTrue(G.describe(include=['number']).shape == (16, 2))
        self.assertTrue(G.describe(include=['number', 'object']).shape == (22, 3))
        self.assertTrue(G.describe(include='all').shape == (26, 4))

    def test_describe_multi_index_df_column_names(self):
        """ Test that column names persist after the describe operation."""

        df = pd.DataFrame({'A': ['foo', 'bar', 'foo', 'bar', 'foo', 'bar', 'foo', 'foo'],
                           'B': ['one', 'one', 'two', 'three', 'two', 'two', 'one', 'three'],
                           'C': np.random.randn(8),
                           'D': np.random.randn(8)})

        # GH 11517
        # test for hierarchical index
        hierarchical_index_df = df.groupby(['A', 'B']).mean().T
        self.assertTrue(hierarchical_index_df.columns.names == ['A', 'B'])
        self.assertTrue(hierarchical_index_df.describe().columns.names == ['A', 'B'])

        # test for non-hierarchical index
        non_hierarchical_index_df = df.groupby(['A']).mean().T
        self.assertTrue(non_hierarchical_index_df.columns.names == ['A'])
        self.assertTrue(non_hierarchical_index_df.describe().columns.names == ['A'])

    def test_no_order(self):
        tm._skip_if_no_scipy()
        s = Series([0, 1, np.nan, 3])
        with tm.assertRaises(ValueError):
            s.interpolate(method='polynomial')
        with tm.assertRaises(ValueError):
            s.interpolate(method='spline')

    def test_spline(self):
        tm._skip_if_no_scipy()
        s = Series([1, 2, np.nan, 4, 5, np.nan, 7])
        result = s.interpolate(method='spline', order=1)
        expected = Series([1., 2., 3., 4., 5., 6., 7.])
        assert_series_equal(result, expected)

    def test_spline_extrapolate(self):
        tm.skip_if_no_package('scipy', '0.15', 'setting ext on scipy.interpolate.UnivariateSpline')
        s = Series([1, 2, 3, 4, np.nan, 6, np.nan])
        result3 = s.interpolate(method='spline', order=1, ext=3)
        expected3 = Series([1., 2., 3., 4., 5., 6., 6.])
        assert_series_equal(result3, expected3)

        result1 = s.interpolate(method='spline', order=1, ext=0)
        expected1 = Series([1., 2., 3., 4., 5., 6., 7.])
        assert_series_equal(result1, expected1)

    def test_spline_smooth(self):
        tm._skip_if_no_scipy()
        s = Series([1, 2, np.nan, 4, 5.1, np.nan, 7])
        self.assertNotEqual(s.interpolate(method='spline', order=3, s=0)[5],
                            s.interpolate(method='spline', order=3)[5])

    def test_spline_interpolation(self):
        tm._skip_if_no_scipy()

        s = Series(np.arange(10)**2)
        s[np.random.randint(0,9,3)] = np.nan
        result1 = s.interpolate(method='spline', order=1)
        expected1 = s.interpolate(method='spline', order=1)
        assert_series_equal(result1, expected1)

    # GH #10633
    def test_spline_error(self):
        tm._skip_if_no_scipy()

        s = pd.Series(np.arange(10)**2)
        s[np.random.randint(0,9,3)] = np.nan
        with tm.assertRaises(ValueError):
            s.interpolate(method='spline')

        with tm.assertRaises(ValueError):
            s.interpolate(method='spline', order=0)

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

        # merging with override
        # GH 6923
        _metadata = DataFrame._metadata
        _finalize = DataFrame.__finalize__

        np.random.seed(10)
        df1 = DataFrame(np.random.randint(0, 4, (3, 2)), columns=['a', 'b'])
        df2 = DataFrame(np.random.randint(0, 4, (3, 2)), columns=['c', 'd'])
        DataFrame._metadata = ['filename']
        df1.filename = 'fname1.csv'
        df2.filename = 'fname2.csv'

        def finalize(self, other, method=None, **kwargs):

            for name in self._metadata:
                if method == 'merge':
                    left, right = other.left, other.right
                    value = getattr(left, name, '') + '|' + getattr(right, name, '')
                    object.__setattr__(self, name, value)
                else:
                    object.__setattr__(self, name, getattr(other, name, ''))

            return self

        DataFrame.__finalize__ = finalize
        result = df1.merge(df2, left_on=['a'], right_on=['c'], how='inner')
        self.assertEqual(result.filename,'fname1.csv|fname2.csv')

        # concat
        # GH 6927
        DataFrame._metadata = ['filename']
        df1 = DataFrame(np.random.randint(0, 4, (3, 2)), columns=list('ab'))
        df1.filename = 'foo'

        def finalize(self, other, method=None, **kwargs):
            for name in self._metadata:
                if method == 'concat':
                    value = '+'.join([ getattr(o,name) for o in other.objs if getattr(o,name,None) ])
                    object.__setattr__(self, name, value)
                else:
                    object.__setattr__(self, name, getattr(other, name, None))

            return self

        DataFrame.__finalize__ = finalize

        result = pd.concat([df1, df1])
        self.assertEqual(result.filename,'foo+foo')

        # reset
        DataFrame._metadata = _metadata
        DataFrame.__finalize__ = _finalize

    def test_tz_convert_and_localize(self):
        l0 = date_range('20140701', periods=5, freq='D')

        # TODO: l1 should be a PeriodIndex for testing
        #       after GH2106 is addressed
        with tm.assertRaises(NotImplementedError):
            period_range('20140701', periods=1).tz_convert('UTC')
        with tm.assertRaises(NotImplementedError):
            period_range('20140701', periods=1).tz_localize('UTC')
        # l1 = period_range('20140701', periods=5, freq='D')
        l1 = date_range('20140701', periods=5, freq='D')

        int_idx = Index(range(5))

        for fn in ['tz_localize', 'tz_convert']:

            if fn == 'tz_convert':
                l0 = l0.tz_localize('UTC')
                l1 = l1.tz_localize('UTC')

            for idx in [l0, l1]:

                l0_expected  = getattr(idx, fn)('US/Pacific')
                l1_expected  = getattr(idx, fn)('US/Pacific')

                df1 = DataFrame(np.ones(5), index=l0)
                df1 = getattr(df1, fn)('US/Pacific')
                self.assertTrue(df1.index.equals(l0_expected))

                # MultiIndex
                # GH7846
                df2 = DataFrame(np.ones(5),
                                MultiIndex.from_arrays([l0, l1]))

                df3 = getattr(df2, fn)('US/Pacific', level=0)
                self.assertFalse(df3.index.levels[0].equals(l0))
                self.assertTrue(df3.index.levels[0].equals(l0_expected))
                self.assertTrue(df3.index.levels[1].equals(l1))
                self.assertFalse(df3.index.levels[1].equals(l1_expected))

                df3 = getattr(df2, fn)('US/Pacific', level=1)
                self.assertTrue(df3.index.levels[0].equals(l0))
                self.assertFalse(df3.index.levels[0].equals(l0_expected))
                self.assertTrue(df3.index.levels[1].equals(l1_expected))
                self.assertFalse(df3.index.levels[1].equals(l1))

                df4 = DataFrame(np.ones(5),
                                MultiIndex.from_arrays([int_idx, l0]))

                df5 = getattr(df4, fn)('US/Pacific', level=1)
                self.assertTrue(df3.index.levels[0].equals(l0))
                self.assertFalse(df3.index.levels[0].equals(l0_expected))
                self.assertTrue(df3.index.levels[1].equals(l1_expected))
                self.assertFalse(df3.index.levels[1].equals(l1))

        # Bad Inputs
        for fn in ['tz_localize', 'tz_convert']:
            # Not DatetimeIndex / PeriodIndex
            with tm.assertRaisesRegexp(TypeError, 'DatetimeIndex'):
                df = DataFrame(index=int_idx)
                df = getattr(df, fn)('US/Pacific')

            # Not DatetimeIndex / PeriodIndex
            with tm.assertRaisesRegexp(TypeError, 'DatetimeIndex'):
                df = DataFrame(np.ones(5),
                            MultiIndex.from_arrays([int_idx, l0]))
                df = getattr(df, fn)('US/Pacific', level=0)

            # Invalid level
            with tm.assertRaisesRegexp(ValueError, 'not valid'):
                df = DataFrame(index=l0)
                df = getattr(df, fn)('US/Pacific', level=1)

    def test_set_attribute(self):
        # Test for consistent setattr behavior when an attribute and a column
        # have the same name (Issue #8994)
        df = DataFrame({'x':[1, 2, 3]})

        df.y = 2
        df['y'] = [2, 4, 6]
        df.y = 5

        assert_equal(df.y, 5)
        assert_series_equal(df['y'], Series([2, 4, 6], name='y'))

    def test_pct_change(self):
        # GH 11150
        pnl = DataFrame([np.arange(0, 40, 10), np.arange(0, 40, 10), np.arange(0, 40, 10)]).astype(np.float64)
        pnl.iat[1,0] = np.nan
        pnl.iat[1,1] = np.nan
        pnl.iat[2,3] = 60

        mask = pnl.isnull()

        for axis in range(2):
            expected = pnl.ffill(axis=axis)/pnl.ffill(axis=axis).shift(axis=axis) - 1
            expected[mask] = np.nan
            result = pnl.pct_change(axis=axis, fill_method='pad')

            self.assert_frame_equal(result, expected)

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

        # don't fail with 0 length dimensions GH11229 & GH8999
        empty_series=pd.Series([], name='five')
        empty_frame=pd.DataFrame([empty_series])
        empty_panel=pd.Panel({'six':empty_frame})

        [tm.assert_series_equal(empty_series, higher_dim.squeeze())
         for higher_dim in [empty_series, empty_frame, empty_panel]]


    def test_equals(self):
        s1 = pd.Series([1, 2, 3], index=[0, 2, 1])
        s2 = s1.copy()
        self.assertTrue(s1.equals(s2))

        s1[1] = 99
        self.assertFalse(s1.equals(s2))

        # NaNs compare as equal
        s1 = pd.Series([1, np.nan, 3, np.nan], index=[0, 2, 1, 3])
        s2 = s1.copy()
        self.assertTrue(s1.equals(s2))

        s2[0] = 9.9
        self.assertFalse(s1.equals(s2))

        idx = MultiIndex.from_tuples([(0, 'a'), (1, 'b'), (2, 'c')])
        s1 = Series([1, 2, np.nan], index=idx)
        s2 = s1.copy()
        self.assertTrue(s1.equals(s2))

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
        self.assertTrue(df1['text'].equals(df2['text']))
        self.assertTrue(df1['start'].equals(df2['start']))
        self.assertTrue(df1['end'].equals(df2['end']))
        self.assertTrue(df1['diff'].equals(df2['diff']))
        self.assertTrue(df1['bool'].equals(df2['bool']))
        self.assertTrue(df1.equals(df2))
        self.assertFalse(df1.equals(object))

        # different dtype
        different = df1.copy()
        different['floats'] = different['floats'].astype('float32')
        self.assertFalse(df1.equals(different))

        # different index
        different_index = -index
        different = df2.set_index(different_index)
        self.assertFalse(df1.equals(different))

        # different columns
        different = df2.copy()
        different.columns = df2.columns[::-1]
        self.assertFalse(df1.equals(different))

        # DatetimeIndex
        index = pd.date_range('2000-1-1', periods=10, freq='T')
        df1 = df1.set_index(index)
        df2 = df1.copy()
        self.assertTrue(df1.equals(df2))

        # MultiIndex
        df3 = df1.set_index(['text'], append=True)
        df2 = df1.set_index(['text'], append=True)
        self.assertTrue(df3.equals(df2))

        df2 = df1.set_index(['floats'], append=True)
        self.assertFalse(df3.equals(df2))

        # NaN in index
        df3 = df1.set_index(['floats'], append=True)
        df2 = df1.set_index(['floats'], append=True)
        self.assertTrue(df3.equals(df2))

        # GH 8437
        a = pd.Series([False, np.nan])
        b = pd.Series([False, np.nan])
        c = pd.Series(index=range(2))
        d = pd.Series(index=range(2))
        e = pd.Series(index=range(2))
        f = pd.Series(index=range(2))
        c[:-1] = d[:-1] = e[0] = f[0] = False
        self.assertTrue(a.equals(a))
        self.assertTrue(a.equals(b))
        self.assertTrue(a.equals(c))
        self.assertTrue(a.equals(d))
        self.assertFalse(a.equals(e))
        self.assertTrue(e.equals(f))

    def test_describe_raises(self):
        with tm.assertRaises(NotImplementedError):
            tm.makePanel().describe()

    def test_pipe(self):
        df = DataFrame({'A': [1, 2, 3]})
        f = lambda x, y: x ** y
        result = df.pipe(f, 2)
        expected = DataFrame({'A': [1, 4, 9]})
        self.assert_frame_equal(result, expected)

        result = df.A.pipe(f, 2)
        self.assert_series_equal(result, expected.A)

    def test_pipe_tuple(self):
        df = DataFrame({'A': [1, 2, 3]})
        f = lambda x, y: y
        result = df.pipe((f, 'y'), 0)
        self.assert_frame_equal(result, df)

        result = df.A.pipe((f, 'y'), 0)
        self.assert_series_equal(result, df.A)

    def test_pipe_tuple_error(self):
        df = DataFrame({"A": [1, 2, 3]})
        f = lambda x, y: y
        with tm.assertRaises(ValueError):
            result = df.pipe((f, 'y'), x=1, y=0)

        with tm.assertRaises(ValueError):
            result = df.A.pipe((f, 'y'), x=1, y=0)

    def test_pipe_panel(self):
        wp = Panel({'r1': DataFrame({"A": [1, 2, 3]})})
        f = lambda x, y: x + y
        result = wp.pipe(f, 2)
        expected = wp + 2
        assert_panel_equal(result, expected)

        result = wp.pipe((f, 'y'), x=1)
        expected = wp + 1
        assert_panel_equal(result, expected)

        with tm.assertRaises(ValueError):
            result = wp.pipe((f, 'y'), x=1, y=1)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
