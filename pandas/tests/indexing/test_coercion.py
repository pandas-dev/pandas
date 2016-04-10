# -*- coding: utf-8 -*-

import nose
import numpy as np

import pandas as pd
import pandas.util.testing as tm
import pandas.compat as compat


###############################################################
# Index / Series common tests which may trigger dtype coercions
###############################################################


class TestIndexCoercion(tm.TestCase):

    _multiprocess_can_split_ = True

    def test_setitem_index_numeric_coercion_int(self):
        # tests setitem with non-existing numeric key
        s = pd.Series([1, 2, 3, 4])
        self.assertEqual(s.index.dtype, np.int64)

        # int + int -> int
        temp = s.copy()
        temp[5] = 5
        tm.assert_series_equal(temp, pd.Series([1, 2, 3, 4, 5],
                                               index=[0, 1, 2, 3, 5]))
        self.assertEqual(temp.index.dtype, np.int64)

        # int + float -> float
        temp = s.copy()
        temp[1.1] = 5
        tm.assert_series_equal(temp, pd.Series([1, 2, 3, 4, 5],
                                               index=[0, 1, 2, 3, 1.1]))
        self.assertEqual(temp.index.dtype, np.float64)

    def test_setitem_index_numeric_coercion_float(self):
        # tests setitem with non-existing numeric key
        s = pd.Series([1, 2, 3, 4], index=[1.1, 2.1, 3.1, 4.1])
        self.assertEqual(s.index.dtype, np.float64)

        # float + int -> int
        temp = s.copy()
        # TODO_GH12747 The result must be float
        with tm.assertRaises(IndexError):
            temp[5] = 5

        # float + float -> float
        temp = s.copy()
        temp[5.1] = 5
        exp = pd.Series([1, 2, 3, 4, 5], index=[1.1, 2.1, 3.1, 4.1, 5.1])
        tm.assert_series_equal(temp, exp)
        self.assertEqual(temp.index.dtype, np.float64)

    def test_insert_numeric_coercion_int(self):
        idx = pd.Int64Index([1, 2, 3, 4])
        self.assertEqual(idx.dtype, np.int64)

        # int + int -> int
        res = idx.insert(1, 1)
        tm.assert_index_equal(res, pd.Index([1, 1, 2, 3, 4]))
        self.assertEqual(res.dtype, np.int64)

        # int + float -> float
        res = idx.insert(1, 1.1)
        tm.assert_index_equal(res, pd.Index([1, 1.1, 2, 3, 4]))
        self.assertEqual(res.dtype, np.float64)

        # int + bool -> int
        res = idx.insert(1, False)
        tm.assert_index_equal(res, pd.Index([1, 0, 2, 3, 4]))
        self.assertEqual(res.dtype, np.int64)

    def test_insert_numeric_coercion_float(self):
        idx = pd.Float64Index([1, 2, 3, 4])
        self.assertEqual(idx.dtype, np.float64)

        # float + int -> int
        res = idx.insert(1, 1)
        tm.assert_index_equal(res, pd.Index([1., 1., 2., 3., 4.]))
        self.assertEqual(res.dtype, np.float64)

        # float + float -> float
        res = idx.insert(1, 1.1)
        tm.assert_index_equal(res, pd.Index([1., 1.1, 2., 3., 4.]))
        self.assertEqual(res.dtype, np.float64)

        # float + bool -> float
        res = idx.insert(1, False)
        tm.assert_index_equal(res, pd.Index([1., 0., 2., 3., 4.]))
        self.assertEqual(res.dtype, np.float64)


class TestSeriesCoercion(tm.TestCase):

    _multiprocess_can_split_ = True

    def setUp(self):
        self.rep = {}
        self.rep['object'] = ['a', 'b']
        self.rep['int64'] = [4, 5]
        self.rep['float64'] = [1.1, 2.2]
        self.rep['complex128'] = [1 + 1j, 2 + 2j]
        self.rep['bool'] = [True, False]

    def test_setitem_numeric_coercion_int(self):
        s = pd.Series([1, 2, 3, 4])
        self.assertEqual(s.dtype, np.int64)

        # int + int -> int
        temp = s.copy()
        temp[1] = 1
        tm.assert_series_equal(temp, pd.Series([1, 1, 3, 4]))
        self.assertEqual(temp.dtype, np.int64)

        # int + float -> float
        # TODO_GH12747 The result must be float
        temp = s.copy()
        temp[1] = 1.1
        # tm.assert_series_equal(temp, pd.Series([1, 1.1, 3, 4]))
        # self.assertEqual(temp.dtype, np.float64)
        tm.assert_series_equal(temp, pd.Series([1, 1, 3, 4]))
        self.assertEqual(temp.dtype, np.int64)

        # int + complex -> complex
        temp = s.copy()
        temp[1] = 1 + 1j
        tm.assert_series_equal(temp, pd.Series([1, 1 + 1j, 3, 4]))
        self.assertEqual(temp.dtype, np.complex128)

        # int + bool -> int
        temp = s.copy()
        temp[1] = True
        tm.assert_series_equal(temp, pd.Series([1, 1, 3, 4]))
        self.assertEqual(temp.dtype, np.int64)

    def test_setitem_numeric_coercion_float(self):
        s = pd.Series([1.1, 2.2, 3.3, 4.4])
        self.assertEqual(s.dtype, np.float64)

        # float + int -> float
        temp = s.copy()
        temp[1] = 1
        tm.assert_series_equal(temp, pd.Series([1.1, 1.0, 3.3, 4.4]))
        self.assertEqual(temp.dtype, np.float64)

        # float + float -> float
        temp = s.copy()
        temp[1] = 1.1
        tm.assert_series_equal(temp, pd.Series([1.1, 1.1, 3.3, 4.4]))
        self.assertEqual(temp.dtype, np.float64)

        # float + complex -> complex
        temp = s.copy()
        temp[1] = 1 + 1j
        tm.assert_series_equal(temp, pd.Series([1.1, 1 + 1j, 3.3, 4.4]))
        self.assertEqual(temp.dtype, np.complex128)

        # float + bool -> float
        temp = s.copy()
        temp[1] = True
        tm.assert_series_equal(temp, pd.Series([1.1, 1.0, 3.3, 4.4]))
        self.assertEqual(temp.dtype, np.float64)

    def test_setitem_numeric_coercion_complex(self):
        s = pd.Series([1 + 1j, 2 + 2j, 3 + 3j, 4 + 4j])
        self.assertEqual(s.dtype, np.complex128)

        # complex + int -> complex
        temp = s.copy()
        temp[1] = 1
        tm.assert_series_equal(temp, pd.Series([1 + 1j, 1, 3 + 3j, 4 + 4j]))
        self.assertEqual(temp.dtype, np.complex128)

        # complex + float -> complex
        temp = s.copy()
        temp[1] = 1.1
        tm.assert_series_equal(temp, pd.Series([1 + 1j, 1.1, 3 + 3j, 4 + 4j]))
        self.assertEqual(temp.dtype, np.complex128)

        # complex + complex -> complex
        temp = s.copy()
        temp[1] = 1 + 1j
        tm.assert_series_equal(temp,
                               pd.Series([1 + 1j, 1 + 1j, 3 + 3j, 4 + 4j]))
        self.assertEqual(temp.dtype, np.complex128)

        # complex + bool -> complex
        temp = s.copy()
        temp[1] = True
        tm.assert_series_equal(temp, pd.Series([1 + 1j, 1, 3 + 3j, 4 + 4j]))
        self.assertEqual(temp.dtype, np.complex128)

    def test_setitem_numeric_coercion_bool(self):
        s = pd.Series([True, False, True, False])
        self.assertEqual(s.dtype, np.bool)

        # bool + int -> int
        # TODO_GH12747 The result must be int
        temp = s.copy()
        temp[1] = 1
        # tm.assert_series_equal(temp, pd.Series([1, 1, 1, 0]))
        # self.assertEqual(temp.dtype, np.int64)
        tm.assert_series_equal(temp, pd.Series([True, True, True, False]))
        self.assertEqual(temp.dtype, np.bool)

        # TODO_GH12747 The result must be int
        temp = s.copy()
        temp[1] = 3     # greater than bool
        # tm.assert_series_equal(temp, pd.Series([1, 3, 1, 0]))
        # self.assertEqual(temp.dtype, np.int64)
        tm.assert_series_equal(temp, pd.Series([True, True, True, False]))
        self.assertEqual(temp.dtype, np.bool)

        # bool + float -> float
        # TODO_GH12747 The result must be float
        temp = s.copy()
        temp[1] = 1.1
        # tm.assert_series_equal(temp, pd.Series([1., 1.1, 1., 0.]))
        # self.assertEqual(temp.dtype, np.float64)
        tm.assert_series_equal(temp, pd.Series([True, True, True, False]))
        self.assertEqual(temp.dtype, np.bool)

        # bool + complex -> complex (buggy, results in bool)
        # TODO_GH12747 The result must be complex
        temp = s.copy()
        temp[1] = 1 + 1j
        # tm.assert_series_equal(temp, pd.Series([1, 1 + 1j, 1, 0]))
        # self.assertEqual(temp.dtype, np.complex128)
        tm.assert_series_equal(temp, pd.Series([True, True, True, False]))
        self.assertEqual(temp.dtype, np.bool)

        # bool + bool -> int
        temp = s.copy()
        temp[1] = True
        tm.assert_series_equal(temp, pd.Series([True, True, True, False]))
        self.assertEqual(temp.dtype, np.bool)

    def test_where_numeric_coercion_int(self):
        s = pd.Series([1, 2, 3, 4])
        self.assertEqual(s.dtype, np.int64)
        cond = pd.Series([True, False, True, False])

        # int + int -> int
        res = s.where(cond, 1)
        tm.assert_series_equal(res, pd.Series([1, 1, 3, 1]))
        self.assertEqual(res.dtype, np.int64)
        res = s.where(cond, pd.Series([5, 6, 7, 8]))
        tm.assert_series_equal(res, pd.Series([1, 6, 3, 8]))
        self.assertEqual(res.dtype, np.int64)

        # int + float -> float
        res = s.where(cond, 1.1)
        tm.assert_series_equal(res, pd.Series([1, 1.1, 3, 1.1]))
        self.assertEqual(res.dtype, np.float64)
        res = s.where(cond, pd.Series([5.5, 6.6, 7.7, 8.8]))
        tm.assert_series_equal(res, pd.Series([1, 6.6, 3, 8.8]))
        self.assertEqual(res.dtype, np.float64)

        # int + complex -> complex
        res = s.where(cond, 1 + 1j)
        tm.assert_series_equal(res, pd.Series([1, 1 + 1j, 3, 1 + 1j]))
        self.assertEqual(res.dtype, np.complex128)
        res = s.where(cond, pd.Series([5 + 5j, 6 + 6j, 7 + 7j, 8 + 8j]))
        tm.assert_series_equal(res, pd.Series([1, 6 + 6j, 3, 8 + 8j]))
        self.assertEqual(res.dtype, np.complex128)

        # int + bool -> int
        res = s.where(cond, True)
        tm.assert_series_equal(res, pd.Series([1, 1, 3, 1]))
        self.assertEqual(res.dtype, np.int64)
        res = s.where(cond, pd.Series([True, False, True, True]))
        tm.assert_series_equal(res, pd.Series([1, 0, 3, 1]))
        self.assertEqual(res.dtype, np.int64)

    def test_where_numeric_coercion_float(self):
        s = pd.Series([1.1, 2.2, 3.3, 4.4])
        self.assertEqual(s.dtype, np.float64)
        cond = pd.Series([True, False, True, False])

        # float + int -> float
        res = s.where(cond, 1)
        tm.assert_series_equal(res, pd.Series([1.1, 1.0, 3.3, 1.0]))
        self.assertEqual(res.dtype, np.float64)
        res = s.where(cond, pd.Series([5, 6, 7, 8]))
        tm.assert_series_equal(res, pd.Series([1.1, 6.0, 3.3, 8.0]))
        self.assertEqual(res.dtype, np.float64)

        # float + float -> float
        res = s.where(cond, 1.1)
        tm.assert_series_equal(res, pd.Series([1.1, 1.1, 3.3, 1.1]))
        self.assertEqual(res.dtype, np.float64)
        res = s.where(cond, pd.Series([5.5, 6.6, 7.7, 8.8]))
        tm.assert_series_equal(res, pd.Series([1.1, 6.6, 3.3, 8.8]))
        self.assertEqual(res.dtype, np.float64)

        # float + complex -> complex
        res = s.where(cond, 1 + 1j)
        tm.assert_series_equal(res, pd.Series([1.1, 1 + 1j, 3.3, 1 + 1j]))
        self.assertEqual(res.dtype, np.complex128)
        res = s.where(cond, pd.Series([5 + 5j, 6 + 6j, 7 + 7j, 8 + 8j]))
        tm.assert_series_equal(res, pd.Series([1.1, 6 + 6j, 3.3, 8 + 8j]))
        self.assertEqual(res.dtype, np.complex128)

        # float + bool -> float
        res = s.where(cond, True)
        tm.assert_series_equal(res, pd.Series([1.1, 1.0, 3.3, 1.0]))
        self.assertEqual(res.dtype, np.float64)
        res = s.where(cond, pd.Series([True, False, True, True]))
        tm.assert_series_equal(res, pd.Series([1.1, 0.0, 3.3, 1.0]))
        self.assertEqual(res.dtype, np.float64)

    def test_where_numeric_coercion_complex(self):
        s = pd.Series([1 + 1j, 2 + 2j, 3 + 3j, 4 + 4j])
        self.assertEqual(s.dtype, np.complex128)
        cond = pd.Series([True, False, True, False])

        # complex + int -> float
        res = s.where(cond, 1)
        tm.assert_series_equal(res, pd.Series([1 + 1j, 1, 3 + 3j, 1]))
        self.assertEqual(res.dtype, np.complex128)
        res = s.where(cond, pd.Series([5, 6, 7, 8]))
        tm.assert_series_equal(res, pd.Series([1 + 1j, 6.0, 3 + 3j, 8.0]))
        self.assertEqual(res.dtype, np.complex128)

        # complex + float -> float
        res = s.where(cond, 1.1)
        tm.assert_series_equal(res, pd.Series([1 + 1j, 1.1, 3 + 3j, 1.1]))
        self.assertEqual(res.dtype, np.complex128)
        res = s.where(cond, pd.Series([5.5, 6.6, 7.7, 8.8]))
        tm.assert_series_equal(res, pd.Series([1 + 1j, 6.6, 3 + 3j, 8.8]))
        self.assertEqual(res.dtype, np.complex128)

        # complex + complex -> complex
        res = s.where(cond, 1 + 1j)
        tm.assert_series_equal(res,
                               pd.Series([1 + 1j, 1 + 1j, 3 + 3j, 1 + 1j]))
        self.assertEqual(res.dtype, np.complex128)
        res = s.where(cond, pd.Series([5 + 5j, 6 + 6j, 7 + 7j, 8 + 8j]))
        tm.assert_series_equal(res,
                               pd.Series([1 + 1j, 6 + 6j, 3 + 3j, 8 + 8j]))
        self.assertEqual(res.dtype, np.complex128)

        # complex + bool -> complex
        res = s.where(cond, True)
        tm.assert_series_equal(res, pd.Series([1 + 1j, 1, 3 + 3j, 1]))
        self.assertEqual(res.dtype, np.complex128)
        res = s.where(cond, pd.Series([True, False, True, True]))
        tm.assert_series_equal(res, pd.Series([1 + 1j, 0, 3 + 3j, 1]))
        self.assertEqual(res.dtype, np.complex128)

    def test_where_numeric_coercion_bool(self):
        s = pd.Series([True, False, True, False])
        self.assertEqual(s.dtype, np.bool)
        cond = pd.Series([True, False, True, False])

        # bool + int -> int
        res = s.where(cond, 1)
        tm.assert_series_equal(res, pd.Series([1, 1, 1, 1]))
        self.assertEqual(res.dtype, np.int64)
        res = s.where(cond, pd.Series([5, 6, 7, 8]))
        tm.assert_series_equal(res, pd.Series([1, 6, 1, 8]))
        self.assertEqual(res.dtype, np.int64)

        # bool + float -> float
        res = s.where(cond, 1.1)
        tm.assert_series_equal(res, pd.Series([1.0, 1.1, 1.0, 1.1]))
        self.assertEqual(res.dtype, np.float64)
        res = s.where(cond, pd.Series([5.5, 6.6, 7.7, 8.8]))
        tm.assert_series_equal(res, pd.Series([1.0, 6.6, 1.0, 8.8]))
        self.assertEqual(res.dtype, np.float64)

        # bool + complex -> complex
        res = s.where(cond, 1 + 1j)
        tm.assert_series_equal(res, pd.Series([1, 1 + 1j, 1, 1 + 1j]))
        self.assertEqual(res.dtype, np.complex128)
        res = s.where(cond, pd.Series([5 + 5j, 6 + 6j, 7 + 7j, 8 + 8j]))
        tm.assert_series_equal(res, pd.Series([1, 6 + 6j, 1, 8 + 8j]))
        self.assertEqual(res.dtype, np.complex128)

        # bool + bool -> bool
        res = s.where(cond, True)
        tm.assert_series_equal(res, pd.Series([True, True, True, True]))
        self.assertEqual(res.dtype, np.bool)
        res = s.where(cond, pd.Series([True, False, True, True]))
        tm.assert_series_equal(res, pd.Series([True, False, True, True]))
        self.assertEqual(res.dtype, np.bool)

    # not indexing, but place here for consisntency

    def test_fillna_numeric_coercion_int(self):
        # int can't hold NaN
        pass

    def test_fillna_numeric_coercion_float(self):
        s = pd.Series([1.1, np.nan, 3.3, 4.4])
        self.assertEqual(s.dtype, np.float64)

        # float + int -> float
        res = s.fillna(1)
        tm.assert_series_equal(res, pd.Series([1.1, 1.0, 3.3, 4.4]))
        self.assertEqual(res.dtype, np.float64)

        # float + float -> float
        res = s.fillna(1.1)
        tm.assert_series_equal(res, pd.Series([1.1, 1.1, 3.3, 4.4]))
        self.assertEqual(res.dtype, np.float64)

        # float + complex -> complex
        res = s.fillna(1 + 1j)
        tm.assert_series_equal(res, pd.Series([1.1, 1 + 1j, 3.3, 4.4]))
        self.assertEqual(res.dtype, np.complex128)

        # float + bool -> float
        res = s.fillna(True)
        tm.assert_series_equal(res, pd.Series([1.1, 1.0, 3.3, 4.4]))
        self.assertEqual(res.dtype, np.float64)

    def test_fillna_numeric_coercion_complex(self):
        s = pd.Series([1 + 1j, np.nan, 3 + 3j, 4 + 4j])
        self.assertEqual(s.dtype, np.complex128)

        # complex + int -> complex
        res = s.fillna(1)
        tm.assert_series_equal(res, pd.Series([1 + 1j, 1, 3 + 3j, 4 + 4j]))
        self.assertEqual(res.dtype, np.complex128)

        # complex + float -> complex
        res = s.fillna(1.1)
        tm.assert_series_equal(res, pd.Series([1 + 1j, 1.1, 3 + 3j, 4 + 4j]))
        self.assertEqual(res.dtype, np.complex128)

        # complex + complex -> complex
        res = s.fillna(1 + 1j)
        tm.assert_series_equal(res,
                               pd.Series([1 + 1j, 1 + 1j, 3 + 3j, 4 + 4j]))
        self.assertEqual(res.dtype, np.complex128)

        # complex + bool -> complex
        res = s.fillna(True)
        tm.assert_series_equal(res, pd.Series([1 + 1j, 1, 3 + 3j, 4 + 4j]))
        self.assertEqual(res.dtype, np.complex128)

    def test_fillna_numeric_coercion_bool(self):
        # bool can't hold NaN
        pass

    def _assert_replace_conversion(self, from_key, to_key, how):
        index = pd.Index([3, 4], name='xxx')
        s = pd.Series(self.rep[from_key], index=index, name='yyy')
        self.assertEqual(s.dtype, from_key)

        if how == 'dict':
            replacer = dict(zip(self.rep[from_key], self.rep[to_key]))
        elif how == 'series':
            replacer = pd.Series(self.rep[to_key], index=self.rep[from_key])
        else:
            raise ValueError

        result = s.replace(replacer)

        if ((from_key == 'float64' and
             to_key in ('bool', 'int64')) or

            (from_key == 'complex128' and
             to_key in ('bool', 'int64', 'float64')) or

            (from_key == 'int64' and
             to_key in ('bool')) or

            # TODO_GH12747 The result must be int?
           (from_key == 'bool' and to_key in ('int64'))):

            # Expected: do not downcast by replacement
            exp = pd.Series(self.rep[to_key], index=index,
                            name='yyy', dtype=from_key)

        else:
            exp = pd.Series(self.rep[to_key], index=index, name='yyy')
            self.assertEqual(exp.dtype, to_key)

        tm.assert_series_equal(result, exp)

    def test_replace_conversion_dict_from_object(self):
        from_key = 'object'
        for to_key in self.rep:
            self._assert_replace_conversion(from_key, to_key, how='dict')

    def test_replace_conversion_dict_from_int(self):
        from_key = 'int64'
        for to_key in self.rep:
            self._assert_replace_conversion(from_key, to_key, how='dict')

    def test_replace_conversion_dict_from_float(self):
        from_key = 'float64'
        for to_key in self.rep:
            self._assert_replace_conversion(from_key, to_key, how='dict')

    def test_replace_conversion_dict_from_complex(self):
        from_key = 'complex128'
        for to_key in self.rep:
            self._assert_replace_conversion(from_key, to_key, how='dict')

    def test_replace_conversion_dict_from_bool(self):
        from_key = 'bool'
        for to_key in self.rep:
            self._assert_replace_conversion(from_key, to_key, how='dict')

    # Series
    def test_replace_conversion_series_from_object(self):
        from_key = 'object'
        for to_key in self.rep:
            self._assert_replace_conversion(from_key, to_key, how='series')

    def test_replace_conversion_series_from_int(self):
        from_key = 'int64'
        for to_key in self.rep:
            self._assert_replace_conversion(from_key, to_key, how='series')

    def test_replace_conversion_series_from_float(self):
        from_key = 'float64'
        for to_key in self.rep:
            self._assert_replace_conversion(from_key, to_key, how='series')

    def test_replace_conversion_series_from_complex(self):
        from_key = 'complex128'
        for to_key in self.rep:
            self._assert_replace_conversion(from_key, to_key, how='series')

    def test_replace_conversion_series_from_bool(self):
        from_key = 'bool'
        for to_key in self.rep:

            if compat.PY3:
                # doesn't work in PY3, though ...dict_from_bool works fine
                raise nose.SkipTest("doesn't work as in PY3")

            self._assert_replace_conversion(from_key, to_key, how='series')
