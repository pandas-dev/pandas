# -*- coding: utf-8 -*-
from datetime import datetime, timedelta, date, time

import numpy as np

import pandas as pd
import pandas.lib as lib
import pandas.util.testing as tm
from pandas.compat import u, PY2


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
                self.assert_numpy_array_equal(target[indices], target[maybe_slice])

                # reverse
                indices = indices[::-1]
                maybe_slice = lib.maybe_indices_to_slice(indices, len(target))
                self.assertTrue(isinstance(maybe_slice, slice))
                self.assert_numpy_array_equal(target[indices], target[maybe_slice])

        # not slice
        for case in [[2, 1, 2, 0], [2, 2, 1, 0], [0, 1, 2, 1], [-2, 0, 2], [2, 0, -2]]:
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
                self.assert_numpy_array_equal(target[indices], target[maybe_slice])

                # reverse
                indices = indices[::-1]
                maybe_slice = lib.maybe_indices_to_slice(indices, len(target))
                self.assertTrue(isinstance(maybe_slice, slice))
                self.assert_numpy_array_equal(target[indices], target[maybe_slice])

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
                self.assert_numpy_array_equal(target[indices], target[maybe_slice])

                # reverse
                indices = indices[::-1]
                maybe_slice = lib.maybe_indices_to_slice(indices, len(target))
                self.assertTrue(isinstance(maybe_slice, slice))
                self.assert_numpy_array_equal(target[indices], target[maybe_slice])

        # not slice
        for case in [[14, 12, 10, 12], [12, 12, 11, 10], [10, 11, 12, 11]]:
            indices = np.array(case, dtype=np.int64)
            maybe_slice = lib.maybe_indices_to_slice(indices, len(target))
            self.assertFalse(isinstance(maybe_slice, slice))
            self.assert_numpy_array_equal(maybe_slice, indices)
            self.assert_numpy_array_equal(target[indices], target[maybe_slice])


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
        self.assertFalse(lib.isscalar((1,)))
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
        for zerodim in [np.array(1),
                        np.array('foobar'),
                        np.array(np.datetime64('2014-01-01')),
                        np.array(np.timedelta64(1, 'h'))]:
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


if __name__ == '__main__':
    import nose

    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)