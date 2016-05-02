# -*- coding: utf-8 -*-
# pylint: disable-msg=E1101,W0612

from operator import methodcaller
import nose
import numpy as np
from numpy import nan
import pandas as pd

from pandas import (Index, Series, DataFrame, Panel, isnull,
                    date_range, period_range, Panel4D)
from pandas.core.index import MultiIndex

import pandas.formats.printing as printing
import pandas.lib as lib

from pandas.compat import range, zip, PY3
from pandas import compat
from pandas.util.testing import (assertRaisesRegexp,
                                 assert_series_equal,
                                 assert_frame_equal,
                                 assert_panel_equal,
                                 assert_panel4d_equal,
                                 assert_almost_equal,
                                 assert_equal)

import pandas.util.testing as tm


# ----------------------------------------------------------------------
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

        if isinstance(shape, int):
            shape = tuple([shape] * self._ndim)
        if value is not None:
            if lib.isscalar(value):
                if value == 'empty':
                    arr = None

                    # remove the info axis
                    kwargs.pop(self._typ._info_axis_name, None)
                else:
                    arr = np.empty(shape, dtype=dtype)
                    arr.fill(value)
            else:
                fshape = np.prod(shape)
                arr = value.ravel()
                new_shape = fshape / arr.shape[0]
                if fshape % arr.shape[0] != 0:
                    raise Exception("invalid value passed in _construct")

                arr = np.repeat(arr, new_shape).reshape(shape)
        else:
            arr = np.random.randn(*shape)
        return self._typ(arr, dtype=dtype, **kwargs)

    def _compare(self, result, expected):
        self._comparator(result, expected)

    def test_rename(self):

        # single axis
        idx = list('ABCD')
        # relabeling values passed into self.rename
        args = [
            str.lower,
            {x: x.lower() for x in idx},
            Series({x: x.lower() for x in idx}),
        ]

        for axis in self._axes():
            kwargs = {axis: idx}
            obj = self._construct(4, **kwargs)

            for arg in args:
                # rename a single axis
                result = obj.rename(**{axis: arg})
                expected = obj.copy()
                setattr(expected, axis, list('abcd'))
                self._compare(result, expected)

        # multiple axes at once

    def test_rename_axis(self):
        idx = list('ABCD')
        # relabeling values passed into self.rename
        args = [
            str.lower,
            {x: x.lower() for x in idx},
            Series({x: x.lower() for x in idx}),
        ]

        for axis in self._axes():
            kwargs = {axis: idx}
            obj = self._construct(4, **kwargs)

            for arg in args:
                # rename a single axis
                result = obj.rename_axis(arg, axis=axis)
                expected = obj.copy()
                setattr(expected, axis, list('abcd'))
                self._compare(result, expected)
            # scalar values
            for arg in ['foo', None]:
                result = obj.rename_axis(arg, axis=axis)
                expected = obj.copy()
                getattr(expected, axis).name = arg
                self._compare(result, expected)

    def test_get_numeric_data(self):

        n = 4
        kwargs = {}
        for i in range(self._ndim):
            kwargs[self._typ._AXIS_NAMES[i]] = list(range(n))

        # get the numeric data
        o = self._construct(n, **kwargs)
        result = o._get_numeric_data()
        self._compare(result, o)

        # non-inclusion
        result = o._get_bool_data()
        expected = self._construct(n, value='empty', **kwargs)
        self._compare(result, expected)

        # get the bool data
        arr = np.array([True, True, False, True])
        o = self._construct(n, value=arr, **kwargs)
        result = o._get_numeric_data()
        self._compare(result, o)

        # _get_numeric_data is includes _get_bool_data, so can't test for
        # non-inclusion

    def test_get_default(self):

        # GH 7725
        d0 = "a", "b", "c", "d"
        d1 = np.arange(4, dtype='int64')
        others = "e", 10

        for data, index in ((d0, d1), (d1, d0)):
            s = Series(data, index=index)
            for i, d in zip(index, data):
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
        self.assertRaises(ValueError, lambda: bool(obj == 0))
        self.assertRaises(ValueError, lambda: bool(obj == 1))
        self.assertRaises(ValueError, lambda: bool(obj))

        obj = self._construct(shape=4, value=1)
        self.assertRaises(ValueError, lambda: bool(obj == 0))
        self.assertRaises(ValueError, lambda: bool(obj == 1))
        self.assertRaises(ValueError, lambda: bool(obj))

        obj = self._construct(shape=4, value=np.nan)
        self.assertRaises(ValueError, lambda: bool(obj == 0))
        self.assertRaises(ValueError, lambda: bool(obj == 1))
        self.assertRaises(ValueError, lambda: bool(obj))

        # empty
        obj = self._construct(shape=0)
        self.assertRaises(ValueError, lambda: bool(obj))

        # invalid behaviors

        obj1 = self._construct(shape=4, value=1)
        obj2 = self._construct(shape=4, value=1)

        def f():
            if obj1:
                printing.pprint_thing("this works and shouldn't")

        self.assertRaises(ValueError, f)
        self.assertRaises(ValueError, lambda: obj1 and obj2)
        self.assertRaises(ValueError, lambda: obj1 or obj2)
        self.assertRaises(ValueError, lambda: not obj1)

    def test_numpy_1_7_compat_numeric_methods(self):
        # GH 4435
        # numpy in 1.7 tries to pass addtional arguments to pandas functions

        o = self._construct(shape=4)
        for op in ['min', 'max', 'max', 'var', 'std', 'prod', 'sum', 'cumsum',
                   'cumprod', 'median', 'skew', 'kurt', 'compound', 'cummax',
                   'cummin', 'all', 'any']:
            f = getattr(np, op, None)
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

        self.assertRaises(NotImplementedError, f, [("A", "datetime64[h]"),
                                                   ("B", "str"),
                                                   ("C", "int32")])

        # these work (though results may be unexpected)
        f('int64')
        f('float64')
        f('M8[ns]')

    def check_metadata(self, x, y=None):
        for m in x._metadata:
            v = getattr(x, m, None)
            if y is None:
                self.assertIsNone(v)
            else:
                self.assertEqual(v, getattr(y, m, None))

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
        for op in ['__add__', '__sub__', '__truediv__', '__mul__']:
            result = getattr(o, op)(1)
            self.check_metadata(o, result)

        # ops with like
        for op in ['__add__', '__sub__', '__truediv__', '__mul__']:
            try:
                result = getattr(o, op)(o)
                self.check_metadata(o, result)
            except (ValueError, AttributeError):
                pass

        # simple boolean
        for op in ['__eq__', '__le__', '__ge__']:
            v1 = getattr(o, op)(o)
            self.check_metadata(o, v1)

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
            self.check_metadata(o, result)
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
        for op in ['__eq__', '__le__', '__ge__']:

            # this is a name matching op
            v1 = getattr(o, op)(o)

            v2 = getattr(o, op)(o2)
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
        for index in [tm.makeFloatIndex, tm.makeIntIndex, tm.makeStringIndex,
                      tm.makeUnicodeIndex, tm.makeDateIndex,
                      tm.makePeriodIndex]:
            axis = o._get_axis_name(0)
            setattr(o, axis, index(len(getattr(o, axis))))

            # Panel + dims
            try:
                o.head()
            except (NotImplementedError):
                raise nose.SkipTest('not implemented on {0}'.format(
                    o.__class__.__name__))

            self._compare(o.head(), o.iloc[:5])
            self._compare(o.tail(), o.iloc[-5:])

            # 0-len
            self._compare(o.head(0), o.iloc[0:0])
            self._compare(o.tail(0), o.iloc[0:0])

            # bounded
            self._compare(o.head(len(o) + 1), o)
            self._compare(o.tail(len(o) + 1), o)

            # neg index
            self._compare(o.head(-3), o.head(7))
            self._compare(o.tail(-3), o.tail(7))

    def test_sample(self):
        # Fixes issue: 2419

        o = self._construct(shape=10)

        ###
        # Check behavior of random_state argument
        ###

        # Check for stability when receives seed or random state -- run 10
        # times.
        for test in range(10):
            seed = np.random.randint(0, 100)
            self._compare(
                o.sample(n=4, random_state=seed), o.sample(n=4,
                                                           random_state=seed))
            self._compare(
                o.sample(frac=0.7, random_state=seed), o.sample(
                    frac=0.7, random_state=seed))

            self._compare(
                o.sample(n=4, random_state=np.random.RandomState(test)),
                o.sample(n=4, random_state=np.random.RandomState(test)))

            self._compare(
                o.sample(frac=0.7, random_state=np.random.RandomState(test)),
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
            o.sample(n=3.2)

        # Check lengths are right
        self.assertTrue(len(o.sample(n=4) == 4))
        self.assertTrue(len(o.sample(frac=0.34) == 3))
        self.assertTrue(len(o.sample(frac=0.36) == 4))

        ###
        # Check weights
        ###

        # Weight length must be right
        with tm.assertRaises(ValueError):
            o.sample(n=3, weights=[0, 1])

        with tm.assertRaises(ValueError):
            bad_weights = [0.5] * 11
            o.sample(n=3, weights=bad_weights)

        with tm.assertRaises(ValueError):
            bad_weight_series = Series([0, 0, 0.2])
            o.sample(n=4, weights=bad_weight_series)

        # Check won't accept negative weights
        with tm.assertRaises(ValueError):
            bad_weights = [-0.1] * 10
            o.sample(n=3, weights=bad_weights)

        # Check inf and -inf throw errors:
        with tm.assertRaises(ValueError):
            weights_with_inf = [0.1] * 10
            weights_with_inf[0] = np.inf
            o.sample(n=3, weights=weights_with_inf)

        with tm.assertRaises(ValueError):
            weights_with_ninf = [0.1] * 10
            weights_with_ninf[0] = -np.inf
            o.sample(n=3, weights=weights_with_ninf)

        # All zeros raises errors
        zero_weights = [0] * 10
        with tm.assertRaises(ValueError):
            o.sample(n=3, weights=zero_weights)

        # All missing weights
        nan_weights = [np.nan] * 10
        with tm.assertRaises(ValueError):
            o.sample(n=3, weights=nan_weights)

        # Check np.nan are replaced by zeros.
        weights_with_nan = [np.nan] * 10
        weights_with_nan[5] = 0.5
        self._compare(
            o.sample(n=1, axis=0, weights=weights_with_nan), o.iloc[5:6])

        # Check None are also replaced by zeros.
        weights_with_None = [None] * 10
        weights_with_None[5] = 0.5
        self._compare(
            o.sample(n=1, axis=0, weights=weights_with_None), o.iloc[5:6])

    def test_size_compat(self):
        # GH8846
        # size property should be defined

        o = self._construct(shape=10)
        self.assertTrue(o.size == np.prod(o.shape))
        self.assertTrue(o.size == 10 ** len(o.axes))

    def test_split_compat(self):
        # xref GH8846
        o = self._construct(shape=10)
        self.assertTrue(len(np.array_split(o, 5)) == 5)
        self.assertTrue(len(np.array_split(o, 2)) == 2)

    def test_unexpected_keyword(self):  # GH8597
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

    # See gh-12301
    def test_stat_unexpected_keyword(self):
        obj = self._construct(5)
        starwars = 'Star Wars'
        errmsg = 'unexpected keyword'

        with assertRaisesRegexp(TypeError, errmsg):
            obj.max(epic=starwars)  # stat_function
        with assertRaisesRegexp(TypeError, errmsg):
            obj.var(epic=starwars)  # stat_function_ddof
        with assertRaisesRegexp(TypeError, errmsg):
            obj.sum(epic=starwars)  # cum_function
        with assertRaisesRegexp(TypeError, errmsg):
            obj.any(epic=starwars)  # logical_function

    def test_api_compat(self):

        # GH 12021
        # compat for __name__, __qualname__

        obj = self._construct(5)
        for func in ['sum', 'cumsum', 'any', 'var']:
            f = getattr(obj, func)
            self.assertEqual(f.__name__, func)
            if PY3:
                self.assertTrue(f.__qualname__.endswith(func))

    def test_stat_non_defaults_args(self):
        obj = self._construct(5)
        out = np.array([0])
        errmsg = "the 'out' parameter is not supported"

        with assertRaisesRegexp(ValueError, errmsg):
            obj.max(out=out)  # stat_function
        with assertRaisesRegexp(ValueError, errmsg):
            obj.var(out=out)  # stat_function_ddof
        with assertRaisesRegexp(ValueError, errmsg):
            obj.sum(out=out)  # cum_function
        with assertRaisesRegexp(ValueError, errmsg):
            obj.any(out=out)  # logical_function

    def test_clip(self):
        lower = 1
        upper = 3
        col = np.arange(5)

        obj = self._construct(len(col), value=col)

        if isinstance(obj, Panel):
            msg = "clip is not supported yet for panels"
            tm.assertRaisesRegexp(NotImplementedError, msg,
                                  obj.clip, lower=lower,
                                  upper=upper)

        else:
            out = obj.clip(lower=lower, upper=upper)
            expected = self._construct(len(col), value=col
                                       .clip(lower, upper))
            self._compare(out, expected)

            bad_axis = 'foo'
            msg = ('No axis named {axis} '
                   'for object').format(axis=bad_axis)
            assertRaisesRegexp(ValueError, msg, obj.clip,
                               lower=lower, upper=upper,
                               axis=bad_axis)

    def test_numpy_clip(self):
        lower = 1
        upper = 3
        col = np.arange(5)

        obj = self._construct(len(col), value=col)

        if isinstance(obj, Panel):
            msg = "clip is not supported yet for panels"
            tm.assertRaisesRegexp(NotImplementedError, msg,
                                  np.clip, obj,
                                  lower, upper)
        else:
            out = np.clip(obj, lower, upper)
            expected = self._construct(len(col), value=col
                                       .clip(lower, upper))
            self._compare(out, expected)

            msg = "the 'out' parameter is not supported"
            tm.assertRaisesRegexp(ValueError, msg,
                                  np.clip, obj,
                                  lower, upper, out=col)


class TestSeries(tm.TestCase, Generic):
    _typ = Series
    _comparator = lambda self, x, y: assert_series_equal(x, y)

    def setUp(self):
        self.ts = tm.makeTimeSeries()  # Was at top level in test_series
        self.ts.name = 'ts'

        self.series = tm.makeStringSeries()
        self.series.name = 'series'

    def test_rename_mi(self):
        s = Series([11, 21, 31],
                   index=MultiIndex.from_tuples(
                       [("A", x) for x in ["a", "B", "c"]]))
        s.rename(str.lower)

    def test_set_axis_name(self):
        s = Series([1, 2, 3], index=['a', 'b', 'c'])
        funcs = ['rename_axis', '_set_axis_name']
        name = 'foo'
        for func in funcs:
            result = methodcaller(func, name)(s)
            self.assertTrue(s.index.name is None)
            self.assertEqual(result.index.name, name)

    def test_set_axis_name_mi(self):
        s = Series([11, 21, 31], index=MultiIndex.from_tuples(
            [("A", x) for x in ["a", "B", "c"]],
            names=['l1', 'l2'])
        )
        funcs = ['rename_axis', '_set_axis_name']
        for func in funcs:
            result = methodcaller(func, ['L1', 'L2'])(s)
            self.assertTrue(s.index.name is None)
            self.assertEqual(s.index.names, ['l1', 'l2'])
            self.assertTrue(result.index.name is None)
            self.assertTrue(result.index.names, ['L1', 'L2'])

    def test_set_axis_name_raises(self):
        s = pd.Series([1])
        with tm.assertRaises(ValueError):
            s._set_axis_name(name='a', axis=1)

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

        o = Series(date_range('20130101', periods=3))
        result = o._get_numeric_data()
        expected = Series([], dtype='M8[ns]', index=pd.Index([], dtype=object))
        self._compare(result, expected)

    def test_nonzero_single_element(self):

        # allow single item via bool method
        s = Series([True])
        self.assertTrue(s.bool())

        s = Series([False])
        self.assertFalse(s.bool())

        # single item nan to raise
        for s in [Series([np.nan]), Series([pd.NaT]), Series([True]),
                  Series([False])]:
            self.assertRaises(ValueError, lambda: bool(s))

        for s in [Series([np.nan]), Series([pd.NaT])]:
            self.assertRaises(ValueError, lambda: s.bool())

        # multiple bool are still an error
        for s in [Series([True, True]), Series([False, False])]:
            self.assertRaises(ValueError, lambda: bool(s))
            self.assertRaises(ValueError, lambda: s.bool())

        # single non-bool are an error
        for s in [Series([1]), Series([0]), Series(['a']), Series([0.0])]:
            self.assertRaises(ValueError, lambda: bool(s))
            self.assertRaises(ValueError, lambda: s.bool())

    def test_metadata_propagation_indiv(self):
        # check that the metadata matches up on the resulting ops

        o = Series(range(3), range(3))
        o.name = 'foo'
        o2 = Series(range(3), range(3))
        o2.name = 'bar'

        result = o.T
        self.check_metadata(o, result)

        # resample
        ts = Series(np.random.rand(1000),
                    index=date_range('20130101', periods=1000, freq='s'),
                    name='foo')
        result = ts.resample('1T').mean()
        self.check_metadata(ts, result)

        result = ts.resample('1T').min()
        self.check_metadata(ts, result)

        result = ts.resample('1T').apply(lambda x: x.sum())
        self.check_metadata(ts, result)

        _metadata = Series._metadata
        _finalize = Series.__finalize__
        Series._metadata = ['name', 'filename']
        o.filename = 'foo'
        o2.filename = 'bar'

        def finalize(self, other, method=None, **kwargs):
            for name in self._metadata:
                if method == 'concat' and name == 'filename':
                    value = '+'.join([getattr(
                        o, name) for o in other.objs if getattr(o, name, None)
                    ])
                    object.__setattr__(self, name, value)
                else:
                    object.__setattr__(self, name, getattr(other, name, None))

            return self

        Series.__finalize__ = finalize

        result = pd.concat([o, o2])
        self.assertEqual(result.filename, 'foo+bar')
        self.assertIsNone(result.name)

        # reset
        Series._metadata = _metadata
        Series.__finalize__ = _finalize

    def test_describe(self):
        self.series.describe()
        self.ts.describe()

    def test_describe_objects(self):
        s = Series(['a', 'b', 'b', np.nan, np.nan, np.nan, 'c', 'd', 'a', 'a'])
        result = s.describe()
        expected = Series({'count': 7, 'unique': 4,
                           'top': 'a', 'freq': 3, 'second': 'b',
                           'second_freq': 2}, index=result.index)
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

    def test_to_xarray(self):

        tm._skip_if_no_xarray()
        from xarray import DataArray

        s = Series([])
        s.index.name = 'foo'
        result = s.to_xarray()
        self.assertEqual(len(result), 0)
        self.assertEqual(len(result.coords), 1)
        assert_almost_equal(list(result.coords.keys()), ['foo'])
        self.assertIsInstance(result, DataArray)

        def testit(index, check_index_type=True):
            s = Series(range(6), index=index(6))
            s.index.name = 'foo'
            result = s.to_xarray()
            repr(result)
            self.assertEqual(len(result), 6)
            self.assertEqual(len(result.coords), 1)
            assert_almost_equal(list(result.coords.keys()), ['foo'])
            self.assertIsInstance(result, DataArray)

            # idempotency
            assert_series_equal(result.to_series(), s,
                                check_index_type=check_index_type)

        for index in [tm.makeFloatIndex, tm.makeIntIndex,
                      tm.makeStringIndex, tm.makeUnicodeIndex,
                      tm.makeDateIndex, tm.makePeriodIndex,
                      tm.makeTimedeltaIndex]:
            testit(index)

        # not idempotent
        testit(tm.makeCategoricalIndex, check_index_type=False)

        s = Series(range(6))
        s.index.name = 'foo'
        s.index = pd.MultiIndex.from_product([['a', 'b'], range(3)],
                                             names=['one', 'two'])
        result = s.to_xarray()
        self.assertEqual(len(result), 2)
        assert_almost_equal(list(result.coords.keys()), ['one', 'two'])
        self.assertIsInstance(result, DataArray)
        assert_series_equal(result.to_series(), s)


class TestDataFrame(tm.TestCase, Generic):
    _typ = DataFrame
    _comparator = lambda self, x, y: assert_frame_equal(x, y)

    def test_rename_mi(self):
        df = DataFrame([
            11, 21, 31
        ], index=MultiIndex.from_tuples([("A", x) for x in ["a", "B", "c"]]))
        df.rename(str.lower)

    def test_set_axis_name(self):
        df = pd.DataFrame([[1, 2], [3, 4]])
        funcs = ['_set_axis_name', 'rename_axis']
        for func in funcs:
            result = methodcaller(func, 'foo')(df)
            self.assertTrue(df.index.name is None)
            self.assertEqual(result.index.name, 'foo')

            result = methodcaller(func, 'cols', axis=1)(df)
            self.assertTrue(df.columns.name is None)
            self.assertEqual(result.columns.name, 'cols')

    def test_set_axis_name_mi(self):
        df = DataFrame(
            np.empty((3, 3)),
            index=MultiIndex.from_tuples([("A", x) for x in list('aBc')]),
            columns=MultiIndex.from_tuples([('C', x) for x in list('xyz')])
        )

        level_names = ['L1', 'L2']
        funcs = ['_set_axis_name', 'rename_axis']
        for func in funcs:
            result = methodcaller(func, level_names)(df)
            self.assertEqual(result.index.names, level_names)
            self.assertEqual(result.columns.names, [None, None])

            result = methodcaller(func, level_names, axis=1)(df)
            self.assertEqual(result.columns.names, ["L1", "L2"])
            self.assertEqual(result.index.names, [None, None])

    def test_nonzero_single_element(self):

        # allow single item via bool method
        df = DataFrame([[True]])
        self.assertTrue(df.bool())

        df = DataFrame([[False]])
        self.assertFalse(df.bool())

        df = DataFrame([[False, False]])
        self.assertRaises(ValueError, lambda: df.bool())
        self.assertRaises(ValueError, lambda: bool(df))

    def test_get_numeric_data_preserve_dtype(self):

        # get the numeric data
        o = DataFrame({'A': [1, '2', 3.]})
        result = o._get_numeric_data()
        expected = DataFrame(index=[0, 1, 2], dtype=object)
        self._compare(result, expected)

    def test_describe(self):
        tm.makeDataFrame().describe()
        tm.makeMixedDataFrame().describe()
        tm.makeTimeDataFrame().describe()

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

        df = DataFrame({"C1": pd.date_range('2010-01-01', periods=4, freq='D')
                        })
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
        expected_cols = ['numC', 'numD', ]
        expected = DataFrame(dict((k, df[k].describe())
                                  for k in expected_cols),
                             columns=expected_cols)
        assert_frame_equal(descN, expected)

        desc = df.describe(include=['number'])
        assert_frame_equal(desc, descN)
        desc = df.describe(exclude=['object', 'datetime'])
        assert_frame_equal(desc, descN)
        desc = df.describe(include=['float'])
        assert_frame_equal(desc, descN.drop('numC', 1))

        descC = df.describe(include=['O'])
        expected_cols = ['catA', 'catB']
        expected = DataFrame(dict((k, df[k].describe())
                                  for k in expected_cols),
                             columns=expected_cols)
        assert_frame_equal(descC, expected)

        descD = df.describe(include=['datetime'])
        assert_series_equal(descD.ts, df.ts.describe())

        desc = df.describe(include=['object', 'number', 'datetime'])
        assert_frame_equal(desc.loc[:, ["numC", "numD"]].dropna(), descN)
        assert_frame_equal(desc.loc[:, ["catA", "catB"]].dropna(), descC)
        descDs = descD.sort_index()  # the index order change for mixed-types
        assert_frame_equal(desc.loc[:, "ts":].dropna().sort_index(), descDs)

        desc = df.loc[:, 'catA':'catB'].describe(include='all')
        assert_frame_equal(desc, descC)
        desc = df.loc[:, 'numC':'numD'].describe(include='all')
        assert_frame_equal(desc, descN)

        desc = df.describe(percentiles=[], include='all')
        cnt = Series(data=[4, 4, 6, 6, 6],
                     index=['catA', 'catB', 'numC', 'numD', 'ts'])
        assert_series_equal(desc.count(), cnt)
        self.assertTrue('count' in desc.index)
        self.assertTrue('unique' in desc.index)
        self.assertTrue('50%' in desc.index)
        self.assertTrue('first' in desc.index)

        desc = df.drop("ts", 1).describe(percentiles=[], include='all')
        assert_series_equal(desc.count(), cnt.drop("ts"))
        self.assertTrue('first' not in desc.index)
        desc = df.drop(["numC", "numD"], 1).describe(percentiles=[],
                                                     include='all')
        assert_series_equal(desc.count(), cnt.drop(["numC", "numD"]))
        self.assertTrue('50%' not in desc.index)

    def test_describe_typefiltering_category_bool(self):
        df = DataFrame({'A_cat': pd.Categorical(['foo', 'foo', 'bar'] * 8),
                        'B_str': ['a', 'b', 'c', 'd'] * 6,
                        'C_bool': [True] * 12 + [False] * 12,
                        'D_num': np.arange(24.) + .5,
                        'E_ts': tm.makeTimeSeries()[:24].index})

        desc = df.describe()
        expected_cols = ['D_num']
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
        df = DataFrame({"td": pd.to_timedelta(np.arange(24) % 20, "D")})
        self.assertTrue(df.describe().loc["mean"][0] == pd.to_timedelta(
            "8d4h"))

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
        self.assertTrue(G.describe(include=['number', 'object']).shape == (22,
                                                                           3))
        self.assertTrue(G.describe(include='all').shape == (26, 4))

    def test_describe_multi_index_df_column_names(self):
        """ Test that column names persist after the describe operation."""

        df = pd.DataFrame(
            {'A': ['foo', 'bar', 'foo', 'bar', 'foo', 'bar', 'foo', 'foo'],
             'B': ['one', 'one', 'two', 'three', 'two', 'two', 'one', 'three'],
             'C': np.random.randn(8),
             'D': np.random.randn(8)})

        # GH 11517
        # test for hierarchical index
        hierarchical_index_df = df.groupby(['A', 'B']).mean().T
        self.assertTrue(hierarchical_index_df.columns.names == ['A', 'B'])
        self.assertTrue(hierarchical_index_df.describe().columns.names ==
                        ['A', 'B'])

        # test for non-hierarchical index
        non_hierarchical_index_df = df.groupby(['A']).mean().T
        self.assertTrue(non_hierarchical_index_df.columns.names == ['A'])
        self.assertTrue(non_hierarchical_index_df.describe().columns.names ==
                        ['A'])

    def test_metadata_propagation_indiv(self):

        # groupby
        df = DataFrame(
            {'A': ['foo', 'bar', 'foo', 'bar', 'foo', 'bar', 'foo', 'foo'],
             'B': ['one', 'one', 'two', 'three', 'two', 'two', 'one', 'three'],
             'C': np.random.randn(8),
             'D': np.random.randn(8)})
        result = df.groupby('A').sum()
        self.check_metadata(df, result)

        # resample
        df = DataFrame(np.random.randn(1000, 2),
                       index=date_range('20130101', periods=1000, freq='s'))
        result = df.resample('1T')
        self.check_metadata(df, result)

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
                    value = getattr(left, name, '') + '|' + getattr(right,
                                                                    name, '')
                    object.__setattr__(self, name, value)
                else:
                    object.__setattr__(self, name, getattr(other, name, ''))

            return self

        DataFrame.__finalize__ = finalize
        result = df1.merge(df2, left_on=['a'], right_on=['c'], how='inner')
        self.assertEqual(result.filename, 'fname1.csv|fname2.csv')

        # concat
        # GH 6927
        DataFrame._metadata = ['filename']
        df1 = DataFrame(np.random.randint(0, 4, (3, 2)), columns=list('ab'))
        df1.filename = 'foo'

        def finalize(self, other, method=None, **kwargs):
            for name in self._metadata:
                if method == 'concat':
                    value = '+'.join([getattr(
                        o, name) for o in other.objs if getattr(o, name, None)
                    ])
                    object.__setattr__(self, name, value)
                else:
                    object.__setattr__(self, name, getattr(other, name, None))

            return self

        DataFrame.__finalize__ = finalize

        result = pd.concat([df1, df1])
        self.assertEqual(result.filename, 'foo+foo')

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

                l0_expected = getattr(idx, fn)('US/Pacific')
                l1_expected = getattr(idx, fn)('US/Pacific')

                df1 = DataFrame(np.ones(5), index=l0)
                df1 = getattr(df1, fn)('US/Pacific')
                self.assertTrue(df1.index.equals(l0_expected))

                # MultiIndex
                # GH7846
                df2 = DataFrame(np.ones(5), MultiIndex.from_arrays([l0, l1]))

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

                # TODO: untested
                df5 = getattr(df4, fn)('US/Pacific', level=1)  # noqa

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
        df = DataFrame({'x': [1, 2, 3]})

        df.y = 2
        df['y'] = [2, 4, 6]
        df.y = 5

        assert_equal(df.y, 5)
        assert_series_equal(df['y'], Series([2, 4, 6], name='y'))

    def test_pct_change(self):
        # GH 11150
        pnl = DataFrame([np.arange(0, 40, 10), np.arange(0, 40, 10), np.arange(
            0, 40, 10)]).astype(np.float64)
        pnl.iat[1, 0] = np.nan
        pnl.iat[1, 1] = np.nan
        pnl.iat[2, 3] = 60

        mask = pnl.isnull()

        for axis in range(2):
            expected = pnl.ffill(axis=axis) / pnl.ffill(axis=axis).shift(
                axis=axis) - 1
            expected[mask] = np.nan
            result = pnl.pct_change(axis=axis, fill_method='pad')

            self.assert_frame_equal(result, expected)

    def test_to_xarray(self):

        tm._skip_if_no_xarray()
        from xarray import Dataset

        df = DataFrame({'a': list('abc'),
                        'b': list(range(1, 4)),
                        'c': np.arange(3, 6).astype('u1'),
                        'd': np.arange(4.0, 7.0, dtype='float64'),
                        'e': [True, False, True],
                        'f': pd.Categorical(list('abc')),
                        'g': pd.date_range('20130101', periods=3),
                        'h': pd.date_range('20130101',
                                           periods=3,
                                           tz='US/Eastern')}
                       )

        df.index.name = 'foo'
        result = df[0:0].to_xarray()
        self.assertEqual(result.dims['foo'], 0)
        self.assertIsInstance(result, Dataset)

        for index in [tm.makeFloatIndex, tm.makeIntIndex,
                      tm.makeStringIndex, tm.makeUnicodeIndex,
                      tm.makeDateIndex, tm.makePeriodIndex,
                      tm.makeCategoricalIndex, tm.makeTimedeltaIndex]:
            df.index = index(3)
            df.index.name = 'foo'
            df.columns.name = 'bar'
            result = df.to_xarray()
            self.assertEqual(result.dims['foo'], 3)
            self.assertEqual(len(result.coords), 1)
            self.assertEqual(len(result.data_vars), 8)
            assert_almost_equal(list(result.coords.keys()), ['foo'])
            self.assertIsInstance(result, Dataset)

            # idempotency
            # categoricals are not preserved
            # datetimes w/tz are not preserved
            # column names are lost
            expected = df.copy()
            expected['f'] = expected['f'].astype(object)
            expected['h'] = expected['h'].astype('datetime64[ns]')
            expected.columns.name = None
            assert_frame_equal(result.to_dataframe(),
                               expected,
                               check_index_type=False)

        # available in 0.7.1
        # MultiIndex
        df.index = pd.MultiIndex.from_product([['a'], range(3)],
                                              names=['one', 'two'])
        result = df.to_xarray()
        self.assertEqual(result.dims['one'], 1)
        self.assertEqual(result.dims['two'], 3)
        self.assertEqual(len(result.coords), 2)
        self.assertEqual(len(result.data_vars), 8)
        assert_almost_equal(list(result.coords.keys()), ['one', 'two'])
        self.assertIsInstance(result, Dataset)

        result = result.to_dataframe()
        expected = df.copy()
        expected['f'] = expected['f'].astype(object)
        expected['h'] = expected['h'].astype('datetime64[ns]')
        expected.columns.name = None
        assert_frame_equal(result,
                           expected,
                           check_index_type=False)


class TestPanel(tm.TestCase, Generic):
    _typ = Panel
    _comparator = lambda self, x, y: assert_panel_equal(x, y)

    def test_to_xarray(self):

        tm._skip_if_no_xarray()
        from xarray import DataArray

        p = tm.makePanel()

        result = p.to_xarray()
        self.assertIsInstance(result, DataArray)
        self.assertEqual(len(result.coords), 3)
        assert_almost_equal(list(result.coords.keys()),
                            ['items', 'major_axis', 'minor_axis'])
        self.assertEqual(len(result.dims), 3)

        # idempotency
        assert_panel_equal(result.to_pandas(), p)


class TestPanel4D(tm.TestCase, Generic):
    _typ = Panel4D
    _comparator = lambda self, x, y: assert_panel4d_equal(x, y)

    def test_sample(self):
        raise nose.SkipTest("sample on Panel4D")

    def test_to_xarray(self):

        tm._skip_if_no_xarray()
        from xarray import DataArray

        p = tm.makePanel4D()

        result = p.to_xarray()
        self.assertIsInstance(result, DataArray)
        self.assertEqual(len(result.coords), 4)
        assert_almost_equal(list(result.coords.keys()),
                            ['labels', 'items', 'major_axis', 'minor_axis'])
        self.assertEqual(len(result.dims), 4)

        # non-convertible
        self.assertRaises(ValueError, lambda: result.to_pandas())


class TestNDFrame(tm.TestCase):
    # tests that don't fit elsewhere

    def test_sample(sel):
        # Fixes issue: 2419
        # additional specific object based tests

        # A few dataframe test with degenerate weights.
        easy_weight_list = [0] * 10
        easy_weight_list[5] = 1

        df = pd.DataFrame({'col1': range(10, 20),
                           'col2': range(20, 30),
                           'colString': ['a'] * 10,
                           'easyweights': easy_weight_list})
        sample1 = df.sample(n=1, weights='easyweights')
        assert_frame_equal(sample1, df.iloc[5:6])

        # Ensure proper error if string given as weight for Series, panel, or
        # DataFrame with axis = 1.
        s = Series(range(10))
        with tm.assertRaises(ValueError):
            s.sample(n=3, weights='weight_column')

        panel = pd.Panel(items=[0, 1, 2], major_axis=[2, 3, 4],
                         minor_axis=[3, 4, 5])
        with tm.assertRaises(ValueError):
            panel.sample(n=1, weights='weight_column')

        with tm.assertRaises(ValueError):
            df.sample(n=1, weights='weight_column', axis=1)

        # Check weighting key error
        with tm.assertRaises(KeyError):
            df.sample(n=3, weights='not_a_real_column_name')

        # Check that re-normalizes weights that don't sum to one.
        weights_less_than_1 = [0] * 10
        weights_less_than_1[0] = 0.5
        tm.assert_frame_equal(
            df.sample(n=1, weights=weights_less_than_1), df.iloc[:1])

        ###
        # Test axis argument
        ###

        # Test axis argument
        df = pd.DataFrame({'col1': range(10), 'col2': ['a'] * 10})
        second_column_weight = [0, 1]
        assert_frame_equal(
            df.sample(n=1, axis=1, weights=second_column_weight), df[['col2']])

        # Different axis arg types
        assert_frame_equal(df.sample(n=1, axis='columns',
                                     weights=second_column_weight),
                           df[['col2']])

        weight = [0] * 10
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
            df.sample(n=1, axis=1, weights=[0.5] * 10)

        # Check weights with axis = 1
        easy_weight_list = [0] * 3
        easy_weight_list[2] = 1

        df = pd.DataFrame({'col1': range(10, 20),
                           'col2': range(20, 30),
                           'colString': ['a'] * 10})
        sample1 = df.sample(n=1, axis=1, weights=easy_weight_list)
        assert_frame_equal(sample1, df[['colString']])

        # Test default axes
        p = pd.Panel(items=['a', 'b', 'c'], major_axis=[2, 4, 6],
                     minor_axis=[1, 3, 5])
        assert_panel_equal(
            p.sample(n=3, random_state=42), p.sample(n=3, axis=1,
                                                     random_state=42))
        assert_frame_equal(
            df.sample(n=3, random_state=42), df.sample(n=3, axis=0,
                                                       random_state=42))

        # Test that function aligns weights with frame
        df = DataFrame(
            {'col1': [5, 6, 7],
             'col2': ['a', 'b', 'c'], }, index=[9, 5, 3])
        s = Series([1, 0, 0], index=[3, 5, 9])
        assert_frame_equal(df.loc[[3]], df.sample(1, weights=s))

        # Weights have index values to be dropped because not in
        # sampled DataFrame
        s2 = Series([0.001, 0, 10000], index=[3, 5, 10])
        assert_frame_equal(df.loc[[3]], df.sample(1, weights=s2))

        # Weights have empty values to be filed with zeros
        s3 = Series([0.01, 0], index=[3, 5])
        assert_frame_equal(df.loc[[3]], df.sample(1, weights=s3))

        # No overlap in weight and sampled DataFrame indices
        s4 = Series([1, 0], index=[1, 2])
        with tm.assertRaises(ValueError):
            df.sample(1, weights=s4)

    def test_squeeze(self):
        # noop
        for s in [tm.makeFloatSeries(), tm.makeStringSeries(),
                  tm.makeObjectSeries()]:
            tm.assert_series_equal(s.squeeze(), s)
        for df in [tm.makeTimeDataFrame()]:
            tm.assert_frame_equal(df.squeeze(), df)
        for p in [tm.makePanel()]:
            tm.assert_panel_equal(p.squeeze(), p)
        for p4d in [tm.makePanel4D()]:
            tm.assert_panel4d_equal(p4d.squeeze(), p4d)

        # squeezing
        df = tm.makeTimeDataFrame().reindex(columns=['A'])
        tm.assert_series_equal(df.squeeze(), df['A'])

        p = tm.makePanel().reindex(items=['ItemA'])
        tm.assert_frame_equal(p.squeeze(), p['ItemA'])

        p = tm.makePanel().reindex(items=['ItemA'], minor_axis=['A'])
        tm.assert_series_equal(p.squeeze(), p.ix['ItemA', :, 'A'])

        p4d = tm.makePanel4D().reindex(labels=['label1'])
        tm.assert_panel_equal(p4d.squeeze(), p4d['label1'])

        p4d = tm.makePanel4D().reindex(labels=['label1'], items=['ItemA'])
        tm.assert_frame_equal(p4d.squeeze(), p4d.ix['label1', 'ItemA'])

        # don't fail with 0 length dimensions GH11229 & GH8999
        empty_series = pd.Series([], name='five')
        empty_frame = pd.DataFrame([empty_series])
        empty_panel = pd.Panel({'six': empty_frame})

        [tm.assert_series_equal(empty_series, higher_dim.squeeze())
         for higher_dim in [empty_series, empty_frame, empty_panel]]

    def test_numpy_squeeze(self):
        s = tm.makeFloatSeries()
        tm.assert_series_equal(np.squeeze(s), s)

        df = tm.makeTimeDataFrame().reindex(columns=['A'])
        tm.assert_series_equal(np.squeeze(df), df['A'])

        msg = "the 'axis' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, msg,
                              np.squeeze, s, axis=0)

    def test_transpose(self):
        msg = ("transpose\(\) got multiple values for "
               "keyword argument 'axes'")
        for s in [tm.makeFloatSeries(), tm.makeStringSeries(),
                  tm.makeObjectSeries()]:
            # calls implementation in pandas/core/base.py
            tm.assert_series_equal(s.transpose(), s)
        for df in [tm.makeTimeDataFrame()]:
            tm.assert_frame_equal(df.transpose().transpose(), df)
        for p in [tm.makePanel()]:
            tm.assert_panel_equal(p.transpose(2, 0, 1)
                                  .transpose(1, 2, 0), p)
            tm.assertRaisesRegexp(TypeError, msg, p.transpose,
                                  2, 0, 1, axes=(2, 0, 1))
        for p4d in [tm.makePanel4D()]:
            tm.assert_panel4d_equal(p4d.transpose(2, 0, 3, 1)
                                    .transpose(1, 3, 0, 2), p4d)
            tm.assertRaisesRegexp(TypeError, msg, p4d.transpose,
                                  2, 0, 3, 1, axes=(2, 0, 3, 1))

    def test_numpy_transpose(self):
        msg = "the 'axes' parameter is not supported"

        s = tm.makeFloatSeries()
        tm.assert_series_equal(
            np.transpose(s), s)
        tm.assertRaisesRegexp(ValueError, msg,
                              np.transpose, s, axes=1)

        df = tm.makeTimeDataFrame()
        tm.assert_frame_equal(np.transpose(
            np.transpose(df)), df)
        tm.assertRaisesRegexp(ValueError, msg,
                              np.transpose, df, axes=1)

        p = tm.makePanel()
        tm.assert_panel_equal(np.transpose(
            np.transpose(p, axes=(2, 0, 1)),
            axes=(1, 2, 0)), p)

        p4d = tm.makePanel4D()
        tm.assert_panel4d_equal(np.transpose(
            np.transpose(p4d, axes=(2, 0, 3, 1)),
            axes=(1, 3, 0, 2)), p4d)

    def test_take(self):
        indices = [1, 5, -2, 6, 3, -1]
        for s in [tm.makeFloatSeries(), tm.makeStringSeries(),
                  tm.makeObjectSeries()]:
            out = s.take(indices)
            expected = Series(data=s.values.take(indices),
                              index=s.index.take(indices))
            tm.assert_series_equal(out, expected)
        for df in [tm.makeTimeDataFrame()]:
            out = df.take(indices)
            expected = DataFrame(data=df.values.take(indices, axis=0),
                                 index=df.index.take(indices),
                                 columns=df.columns)
            tm.assert_frame_equal(out, expected)

        indices = [-3, 2, 0, 1]
        for p in [tm.makePanel()]:
            out = p.take(indices)
            expected = Panel(data=p.values.take(indices, axis=0),
                             items=p.items.take(indices),
                             major_axis=p.major_axis,
                             minor_axis=p.minor_axis)
            tm.assert_panel_equal(out, expected)
        for p4d in [tm.makePanel4D()]:
            out = p4d.take(indices)
            expected = Panel4D(data=p4d.values.take(indices, axis=0),
                               labels=p4d.labels.take(indices),
                               major_axis=p4d.major_axis,
                               minor_axis=p4d.minor_axis,
                               items=p4d.items)
            tm.assert_panel4d_equal(out, expected)

    def test_take_invalid_kwargs(self):
        indices = [-3, 2, 0, 1]
        s = tm.makeFloatSeries()
        df = tm.makeTimeDataFrame()
        p = tm.makePanel()
        p4d = tm.makePanel4D()

        for obj in (s, df, p, p4d):
            msg = "take\(\) got an unexpected keyword argument 'foo'"
            tm.assertRaisesRegexp(TypeError, msg, obj.take,
                                  indices, foo=2)

            msg = "the 'out' parameter is not supported"
            tm.assertRaisesRegexp(ValueError, msg, obj.take,
                                  indices, out=indices)

            msg = "the 'mode' parameter is not supported"
            tm.assertRaisesRegexp(ValueError, msg, obj.take,
                                  indices, mode='clip')

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
        df1 = DataFrame(
            np.random.random(10, ), index=index, columns=['floats'])
        df1['text'] = 'the sky is so blue. we could use more chocolate.'.split(
        )
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
            df.pipe((f, 'y'), x=1, y=0)

        with tm.assertRaises(ValueError):
            df.A.pipe((f, 'y'), x=1, y=0)

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
