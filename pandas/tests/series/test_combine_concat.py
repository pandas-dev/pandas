# coding=utf-8
# pylint: disable-msg=E1101,W0612

from datetime import datetime

from numpy import nan
import numpy as np
import pandas as pd

from pandas import Series, DataFrame

from pandas import compat
from pandas.util.testing import assert_series_equal
import pandas.util.testing as tm

from .common import TestData


class TestSeriesCombine(TestData, tm.TestCase):

    _multiprocess_can_split_ = True

    def test_append(self):
        appendedSeries = self.series.append(self.objSeries)
        for idx, value in compat.iteritems(appendedSeries):
            if idx in self.series.index:
                self.assertEqual(value, self.series[idx])
            elif idx in self.objSeries.index:
                self.assertEqual(value, self.objSeries[idx])
            else:
                self.fail("orphaned index!")

        self.assertRaises(ValueError, self.ts.append, self.ts,
                          verify_integrity=True)

    def test_append_many(self):
        pieces = [self.ts[:5], self.ts[5:10], self.ts[10:]]

        result = pieces[0].append(pieces[1:])
        assert_series_equal(result, self.ts)

    def test_combine_first(self):
        values = tm.makeIntIndex(20).values.astype(float)
        series = Series(values, index=tm.makeIntIndex(20))

        series_copy = series * 2
        series_copy[::2] = np.NaN

        # nothing used from the input
        combined = series.combine_first(series_copy)

        self.assert_numpy_array_equal(combined, series)

        # Holes filled from input
        combined = series_copy.combine_first(series)
        self.assertTrue(np.isfinite(combined).all())

        self.assert_numpy_array_equal(combined[::2], series[::2])
        self.assert_numpy_array_equal(combined[1::2], series_copy[1::2])

        # mixed types
        index = tm.makeStringIndex(20)
        floats = Series(tm.randn(20), index=index)
        strings = Series(tm.makeStringIndex(10), index=index[::2])

        combined = strings.combine_first(floats)

        tm.assert_dict_equal(strings, combined, compare_keys=False)
        tm.assert_dict_equal(floats[1::2], combined, compare_keys=False)

        # corner case
        s = Series([1., 2, 3], index=[0, 1, 2])
        result = s.combine_first(Series([], index=[]))
        assert_series_equal(s, result)

    def test_update(self):
        s = Series([1.5, nan, 3., 4., nan])
        s2 = Series([nan, 3.5, nan, 5.])
        s.update(s2)

        expected = Series([1.5, 3.5, 3., 5., np.nan])
        assert_series_equal(s, expected)

        # GH 3217
        df = DataFrame([{"a": 1}, {"a": 3, "b": 2}])
        df['c'] = np.nan

        # this will fail as long as series is a sub-class of ndarray
        # df['c'].update(Series(['foo'],index=[0])) #####

    def test_concat_empty_series_dtypes_roundtrips(self):

        # round-tripping with self & like self
        dtypes = map(np.dtype, ['float64', 'int8', 'uint8', 'bool', 'm8[ns]',
                                'M8[ns]'])

        for dtype in dtypes:
            self.assertEqual(pd.concat([Series(dtype=dtype)]).dtype, dtype)
            self.assertEqual(pd.concat([Series(dtype=dtype),
                                        Series(dtype=dtype)]).dtype, dtype)

        def int_result_type(dtype, dtype2):
            typs = set([dtype.kind, dtype2.kind])
            if not len(typs - set(['i', 'u', 'b'])) and (dtype.kind == 'i' or
                                                         dtype2.kind == 'i'):
                return 'i'
            elif not len(typs - set(['u', 'b'])) and (dtype.kind == 'u' or
                                                      dtype2.kind == 'u'):
                return 'u'
            return None

        def float_result_type(dtype, dtype2):
            typs = set([dtype.kind, dtype2.kind])
            if not len(typs - set(['f', 'i', 'u'])) and (dtype.kind == 'f' or
                                                         dtype2.kind == 'f'):
                return 'f'
            return None

        def get_result_type(dtype, dtype2):
            result = float_result_type(dtype, dtype2)
            if result is not None:
                return result
            result = int_result_type(dtype, dtype2)
            if result is not None:
                return result
            return 'O'

        for dtype in dtypes:
            for dtype2 in dtypes:
                if dtype == dtype2:
                    continue

                expected = get_result_type(dtype, dtype2)
                result = pd.concat([Series(dtype=dtype), Series(dtype=dtype2)
                                    ]).dtype
                self.assertEqual(result.kind, expected)

    def test_concat_empty_series_dtypes(self):

        # bools
        self.assertEqual(pd.concat([Series(dtype=np.bool_),
                                    Series(dtype=np.int32)]).dtype, np.int32)
        self.assertEqual(pd.concat([Series(dtype=np.bool_),
                                    Series(dtype=np.float32)]).dtype,
                         np.object_)

        # datetimelike
        self.assertEqual(pd.concat([Series(dtype='m8[ns]'),
                                    Series(dtype=np.bool)]).dtype, np.object_)
        self.assertEqual(pd.concat([Series(dtype='m8[ns]'),
                                    Series(dtype=np.int64)]).dtype, np.object_)
        self.assertEqual(pd.concat([Series(dtype='M8[ns]'),
                                    Series(dtype=np.bool)]).dtype, np.object_)
        self.assertEqual(pd.concat([Series(dtype='M8[ns]'),
                                    Series(dtype=np.int64)]).dtype, np.object_)
        self.assertEqual(pd.concat([Series(dtype='M8[ns]'),
                                    Series(dtype=np.bool_),
                                    Series(dtype=np.int64)]).dtype, np.object_)

        # categorical
        self.assertEqual(pd.concat([Series(dtype='category'),
                                    Series(dtype='category')]).dtype,
                         'category')
        self.assertEqual(pd.concat([Series(dtype='category'),
                                    Series(dtype='float64')]).dtype,
                         np.object_)
        self.assertEqual(pd.concat([Series(dtype='category'),
                                    Series(dtype='object')]).dtype, 'category')

        # sparse
        result = pd.concat([Series(dtype='float64').to_sparse(), Series(
            dtype='float64').to_sparse()])
        self.assertEqual(result.dtype, np.float64)
        self.assertEqual(result.ftype, 'float64:sparse')

        result = pd.concat([Series(dtype='float64').to_sparse(), Series(
            dtype='float64')])
        self.assertEqual(result.dtype, np.float64)
        self.assertEqual(result.ftype, 'float64:sparse')

        result = pd.concat([Series(dtype='float64').to_sparse(), Series(
            dtype='object')])
        self.assertEqual(result.dtype, np.object_)
        self.assertEqual(result.ftype, 'object:dense')

    def test_combine_first_dt64(self):
        from pandas.tseries.tools import to_datetime
        s0 = to_datetime(Series(["2010", np.NaN]))
        s1 = to_datetime(Series([np.NaN, "2011"]))
        rs = s0.combine_first(s1)
        xp = to_datetime(Series(['2010', '2011']))
        assert_series_equal(rs, xp)

        s0 = to_datetime(Series(["2010", np.NaN]))
        s1 = Series([np.NaN, "2011"])
        rs = s0.combine_first(s1)
        xp = Series([datetime(2010, 1, 1), '2011'])
        assert_series_equal(rs, xp)
