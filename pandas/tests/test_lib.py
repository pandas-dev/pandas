# -*- coding: utf-8 -*-

from datetime import datetime, timedelta, date, time

import numpy as np
import pandas as pd
import pandas.lib as lib
import pandas.util.testing as tm

from pandas.compat import long, u, PY2


def _assert_same_values_and_dtype(res, exp):
    tm.assert_equal(res.dtype, exp.dtype)
    tm.assert_almost_equal(res, exp)


class TestMisc(tm.TestCase):

    def test_max_len_string_array(self):

        arr = a = np.array(['foo', 'b', np.nan], dtype='object')
        self.assertTrue(lib.max_len_string_array(arr), 3)

        # unicode
        arr = a.astype('U').astype(object)
        self.assertTrue(lib.max_len_string_array(arr), 3)

        # bytes for python3
        arr = a.astype('S').astype(object)
        self.assertTrue(lib.max_len_string_array(arr), 3)

        # raises
        tm.assertRaises(TypeError,
                        lambda: lib.max_len_string_array(arr.astype('U')))

    def test_infer_dtype_bytes(self):
        compare = 'string' if PY2 else 'bytes'

        # string array of bytes
        arr = np.array(list('abc'), dtype='S1')
        self.assertEqual(pd.lib.infer_dtype(arr), compare)

        # object array of bytes
        arr = arr.astype(object)
        self.assertEqual(pd.lib.infer_dtype(arr), compare)

    def test_maybe_indices_to_slice_left_edge(self):
        target = np.arange(100)

        # slice
        indices = np.array([], dtype=np.int64)
        maybe_slice = lib.maybe_indices_to_slice(indices, len(target))
        self.assertTrue(isinstance(maybe_slice, slice))
        self.assert_numpy_array_equal(target[indices], target[maybe_slice])

        for end in [1, 2, 5, 20, 99]:
            for step in [1, 2, 4]:
                indices = np.arange(0, end, step, dtype=np.int64)
                maybe_slice = lib.maybe_indices_to_slice(indices, len(target))
                self.assertTrue(isinstance(maybe_slice, slice))
                self.assert_numpy_array_equal(target[indices],
                                              target[maybe_slice])

                # reverse
                indices = indices[::-1]
                maybe_slice = lib.maybe_indices_to_slice(indices, len(target))
                self.assertTrue(isinstance(maybe_slice, slice))
                self.assert_numpy_array_equal(target[indices],
                                              target[maybe_slice])

        # not slice
        for case in [[2, 1, 2, 0], [2, 2, 1, 0], [0, 1, 2, 1], [-2, 0, 2],
                     [2, 0, -2]]:
            indices = np.array(case, dtype=np.int64)
            maybe_slice = lib.maybe_indices_to_slice(indices, len(target))
            self.assertFalse(isinstance(maybe_slice, slice))
            self.assert_numpy_array_equal(maybe_slice, indices)
            self.assert_numpy_array_equal(target[indices], target[maybe_slice])

    def test_maybe_indices_to_slice_right_edge(self):
        target = np.arange(100)

        # slice
        for start in [0, 2, 5, 20, 97, 98]:
            for step in [1, 2, 4]:
                indices = np.arange(start, 99, step, dtype=np.int64)
                maybe_slice = lib.maybe_indices_to_slice(indices, len(target))
                self.assertTrue(isinstance(maybe_slice, slice))
                self.assert_numpy_array_equal(target[indices],
                                              target[maybe_slice])

                # reverse
                indices = indices[::-1]
                maybe_slice = lib.maybe_indices_to_slice(indices, len(target))
                self.assertTrue(isinstance(maybe_slice, slice))
                self.assert_numpy_array_equal(target[indices],
                                              target[maybe_slice])

        # not slice
        indices = np.array([97, 98, 99, 100], dtype=np.int64)
        maybe_slice = lib.maybe_indices_to_slice(indices, len(target))
        self.assertFalse(isinstance(maybe_slice, slice))
        self.assert_numpy_array_equal(maybe_slice, indices)
        with self.assertRaises(IndexError):
            target[indices]
        with self.assertRaises(IndexError):
            target[maybe_slice]

        indices = np.array([100, 99, 98, 97], dtype=np.int64)
        maybe_slice = lib.maybe_indices_to_slice(indices, len(target))
        self.assertFalse(isinstance(maybe_slice, slice))
        self.assert_numpy_array_equal(maybe_slice, indices)
        with self.assertRaises(IndexError):
            target[indices]
        with self.assertRaises(IndexError):
            target[maybe_slice]

        for case in [[99, 97, 99, 96], [99, 99, 98, 97], [98, 98, 97, 96]]:
            indices = np.array(case, dtype=np.int64)
            maybe_slice = lib.maybe_indices_to_slice(indices, len(target))
            self.assertFalse(isinstance(maybe_slice, slice))
            self.assert_numpy_array_equal(maybe_slice, indices)
            self.assert_numpy_array_equal(target[indices], target[maybe_slice])

    def test_maybe_indices_to_slice_both_edges(self):
        target = np.arange(10)

        # slice
        for step in [1, 2, 4, 5, 8, 9]:
            indices = np.arange(0, 9, step, dtype=np.int64)
            maybe_slice = lib.maybe_indices_to_slice(indices, len(target))
            self.assertTrue(isinstance(maybe_slice, slice))
            self.assert_numpy_array_equal(target[indices], target[maybe_slice])

            # reverse
            indices = indices[::-1]
            maybe_slice = lib.maybe_indices_to_slice(indices, len(target))
            self.assertTrue(isinstance(maybe_slice, slice))
            self.assert_numpy_array_equal(target[indices], target[maybe_slice])

        # not slice
        for case in [[4, 2, 0, -2], [2, 2, 1, 0], [0, 1, 2, 1]]:
            indices = np.array(case, dtype=np.int64)
            maybe_slice = lib.maybe_indices_to_slice(indices, len(target))
            self.assertFalse(isinstance(maybe_slice, slice))
            self.assert_numpy_array_equal(maybe_slice, indices)
            self.assert_numpy_array_equal(target[indices], target[maybe_slice])

    def test_maybe_indices_to_slice_middle(self):
        target = np.arange(100)

        # slice
        for start, end in [(2, 10), (5, 25), (65, 97)]:
            for step in [1, 2, 4, 20]:
                indices = np.arange(start, end, step, dtype=np.int64)
                maybe_slice = lib.maybe_indices_to_slice(indices, len(target))
                self.assertTrue(isinstance(maybe_slice, slice))
                self.assert_numpy_array_equal(target[indices],
                                              target[maybe_slice])

                # reverse
                indices = indices[::-1]
                maybe_slice = lib.maybe_indices_to_slice(indices, len(target))
                self.assertTrue(isinstance(maybe_slice, slice))
                self.assert_numpy_array_equal(target[indices],
                                              target[maybe_slice])

        # not slice
        for case in [[14, 12, 10, 12], [12, 12, 11, 10], [10, 11, 12, 11]]:
            indices = np.array(case, dtype=np.int64)
            maybe_slice = lib.maybe_indices_to_slice(indices, len(target))
            self.assertFalse(isinstance(maybe_slice, slice))
            self.assert_numpy_array_equal(maybe_slice, indices)
            self.assert_numpy_array_equal(target[indices], target[maybe_slice])

    def test_isinf_scalar(self):
        # GH 11352
        self.assertTrue(lib.isposinf_scalar(float('inf')))
        self.assertTrue(lib.isposinf_scalar(np.inf))
        self.assertFalse(lib.isposinf_scalar(-np.inf))
        self.assertFalse(lib.isposinf_scalar(1))
        self.assertFalse(lib.isposinf_scalar('a'))

        self.assertTrue(lib.isneginf_scalar(float('-inf')))
        self.assertTrue(lib.isneginf_scalar(-np.inf))
        self.assertFalse(lib.isneginf_scalar(np.inf))
        self.assertFalse(lib.isneginf_scalar(1))
        self.assertFalse(lib.isneginf_scalar('a'))


class Testisscalar(tm.TestCase):

    def test_isscalar_builtin_scalars(self):
        self.assertTrue(lib.isscalar(None))
        self.assertTrue(lib.isscalar(True))
        self.assertTrue(lib.isscalar(False))
        self.assertTrue(lib.isscalar(0.))
        self.assertTrue(lib.isscalar(np.nan))
        self.assertTrue(lib.isscalar('foobar'))
        self.assertTrue(lib.isscalar(b'foobar'))
        self.assertTrue(lib.isscalar(u('efoobar')))
        self.assertTrue(lib.isscalar(datetime(2014, 1, 1)))
        self.assertTrue(lib.isscalar(date(2014, 1, 1)))
        self.assertTrue(lib.isscalar(time(12, 0)))
        self.assertTrue(lib.isscalar(timedelta(hours=1)))
        self.assertTrue(lib.isscalar(pd.NaT))

    def test_isscalar_builtin_nonscalars(self):
        self.assertFalse(lib.isscalar({}))
        self.assertFalse(lib.isscalar([]))
        self.assertFalse(lib.isscalar([1]))
        self.assertFalse(lib.isscalar(()))
        self.assertFalse(lib.isscalar((1, )))
        self.assertFalse(lib.isscalar(slice(None)))
        self.assertFalse(lib.isscalar(Ellipsis))

    def test_isscalar_numpy_array_scalars(self):
        self.assertTrue(lib.isscalar(np.int64(1)))
        self.assertTrue(lib.isscalar(np.float64(1.)))
        self.assertTrue(lib.isscalar(np.int32(1)))
        self.assertTrue(lib.isscalar(np.object_('foobar')))
        self.assertTrue(lib.isscalar(np.str_('foobar')))
        self.assertTrue(lib.isscalar(np.unicode_(u('foobar'))))
        self.assertTrue(lib.isscalar(np.bytes_(b'foobar')))
        self.assertTrue(lib.isscalar(np.datetime64('2014-01-01')))
        self.assertTrue(lib.isscalar(np.timedelta64(1, 'h')))

    def test_isscalar_numpy_zerodim_arrays(self):
        for zerodim in [np.array(1), np.array('foobar'),
                        np.array(np.datetime64('2014-01-01')),
                        np.array(np.timedelta64(1, 'h')),
                        np.array(np.datetime64('NaT'))]:
            self.assertFalse(lib.isscalar(zerodim))
            self.assertTrue(lib.isscalar(lib.item_from_zerodim(zerodim)))

    def test_isscalar_numpy_arrays(self):
        self.assertFalse(lib.isscalar(np.array([])))
        self.assertFalse(lib.isscalar(np.array([[]])))
        self.assertFalse(lib.isscalar(np.matrix('1; 2')))

    def test_isscalar_pandas_scalars(self):
        self.assertTrue(lib.isscalar(pd.Timestamp('2014-01-01')))
        self.assertTrue(lib.isscalar(pd.Timedelta(hours=1)))
        self.assertTrue(lib.isscalar(pd.Period('2014-01-01')))

    def test_lisscalar_pandas_containers(self):
        self.assertFalse(lib.isscalar(pd.Series()))
        self.assertFalse(lib.isscalar(pd.Series([1])))
        self.assertFalse(lib.isscalar(pd.DataFrame()))
        self.assertFalse(lib.isscalar(pd.DataFrame([[1]])))
        self.assertFalse(lib.isscalar(pd.Panel()))
        self.assertFalse(lib.isscalar(pd.Panel([[[1]]])))
        self.assertFalse(lib.isscalar(pd.Index([])))
        self.assertFalse(lib.isscalar(pd.Index([1])))


class TestParseSQL(tm.TestCase):

    def test_convert_sql_column_floats(self):
        arr = np.array([1.5, None, 3, 4.2], dtype=object)
        result = lib.convert_sql_column(arr)
        expected = np.array([1.5, np.nan, 3, 4.2], dtype='f8')
        _assert_same_values_and_dtype(result, expected)

    def test_convert_sql_column_strings(self):
        arr = np.array(['1.5', None, '3', '4.2'], dtype=object)
        result = lib.convert_sql_column(arr)
        expected = np.array(['1.5', np.nan, '3', '4.2'], dtype=object)
        _assert_same_values_and_dtype(result, expected)

    def test_convert_sql_column_unicode(self):
        arr = np.array([u('1.5'), None, u('3'), u('4.2')],
                       dtype=object)
        result = lib.convert_sql_column(arr)
        expected = np.array([u('1.5'), np.nan, u('3'), u('4.2')],
                            dtype=object)
        _assert_same_values_and_dtype(result, expected)

    def test_convert_sql_column_ints(self):
        arr = np.array([1, 2, 3, 4], dtype='O')
        arr2 = np.array([1, 2, 3, 4], dtype='i4').astype('O')
        result = lib.convert_sql_column(arr)
        result2 = lib.convert_sql_column(arr2)
        expected = np.array([1, 2, 3, 4], dtype='i8')
        _assert_same_values_and_dtype(result, expected)
        _assert_same_values_and_dtype(result2, expected)

        arr = np.array([1, 2, 3, None, 4], dtype='O')
        result = lib.convert_sql_column(arr)
        expected = np.array([1, 2, 3, np.nan, 4], dtype='f8')
        _assert_same_values_and_dtype(result, expected)

    def test_convert_sql_column_longs(self):
        arr = np.array([long(1), long(2), long(3), long(4)], dtype='O')
        result = lib.convert_sql_column(arr)
        expected = np.array([1, 2, 3, 4], dtype='i8')
        _assert_same_values_and_dtype(result, expected)

        arr = np.array([long(1), long(2), long(3), None, long(4)], dtype='O')
        result = lib.convert_sql_column(arr)
        expected = np.array([1, 2, 3, np.nan, 4], dtype='f8')
        _assert_same_values_and_dtype(result, expected)

    def test_convert_sql_column_bools(self):
        arr = np.array([True, False, True, False], dtype='O')
        result = lib.convert_sql_column(arr)
        expected = np.array([True, False, True, False], dtype=bool)
        _assert_same_values_and_dtype(result, expected)

        arr = np.array([True, False, None, False], dtype='O')
        result = lib.convert_sql_column(arr)
        expected = np.array([True, False, np.nan, False], dtype=object)
        _assert_same_values_and_dtype(result, expected)

    def test_convert_sql_column_decimals(self):
        from decimal import Decimal
        arr = np.array([Decimal('1.5'), None, Decimal('3'), Decimal('4.2')])
        result = lib.convert_sql_column(arr)
        expected = np.array([1.5, np.nan, 3, 4.2], dtype='f8')
        _assert_same_values_and_dtype(result, expected)

if __name__ == '__main__':
    import nose

    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
