# -*- coding: utf-8 -*-

import pytest
import numpy as np

import pandas as pd
import pandas.util.testing as tm
import pandas.compat as compat


###############################################################
# Index / Series common tests which may trigger dtype coercions
###############################################################


class CoercionBase(object):

    klasses = ['index', 'series']
    dtypes = ['object', 'int64', 'float64', 'complex128', 'bool',
              'datetime64', 'datetime64tz', 'timedelta64', 'period']

    @property
    def method(self):
        raise NotImplementedError(self)

    def _assert(self, left, right, dtype):
        # explicitly check dtype to avoid any unexpected result
        if isinstance(left, pd.Series):
            tm.assert_series_equal(left, right)
        elif isinstance(left, pd.Index):
            tm.assert_index_equal(left, right)
        else:
            raise NotImplementedError
        assert left.dtype == dtype
        assert right.dtype == dtype

    def test_has_comprehensive_tests(self):
        for klass in self.klasses:
            for dtype in self.dtypes:
                method_name = 'test_{0}_{1}_{2}'.format(self.method,
                                                        klass, dtype)
                if not hasattr(self, method_name):
                    msg = 'test method is not defined: {0}, {1}'
                    raise AssertionError(msg.format(type(self), method_name))


class TestSetitemCoercion(CoercionBase):

    method = 'setitem'

    def _assert_setitem_series_conversion(self, original_series, loc_value,
                                          expected_series, expected_dtype):
        """ test series value's coercion triggered by assignment """
        temp = original_series.copy()
        temp[1] = loc_value
        tm.assert_series_equal(temp, expected_series)
        # check dtype explicitly for sure
        assert temp.dtype == expected_dtype

        # .loc works different rule, temporary disable
        # temp = original_series.copy()
        # temp.loc[1] = loc_value
        # tm.assert_series_equal(temp, expected_series)

    def test_setitem_series_object(self):
        obj = pd.Series(list('abcd'))
        assert obj.dtype == np.object

        # object + int -> object
        exp = pd.Series(['a', 1, 'c', 'd'])
        self._assert_setitem_series_conversion(obj, 1, exp, np.object)

        # object + float -> object
        exp = pd.Series(['a', 1.1, 'c', 'd'])
        self._assert_setitem_series_conversion(obj, 1.1, exp, np.object)

        # object + complex -> object
        exp = pd.Series(['a', 1 + 1j, 'c', 'd'])
        self._assert_setitem_series_conversion(obj, 1 + 1j, exp, np.object)

        # object + bool -> object
        exp = pd.Series(['a', True, 'c', 'd'])
        self._assert_setitem_series_conversion(obj, True, exp, np.object)

    def test_setitem_series_int64(self):
        obj = pd.Series([1, 2, 3, 4])
        assert obj.dtype == np.int64

        # int + int -> int
        exp = pd.Series([1, 1, 3, 4])
        self._assert_setitem_series_conversion(obj, 1, exp, np.int64)

        # int + float -> float
        # TODO_GH12747 The result must be float
        # tm.assert_series_equal(temp, pd.Series([1, 1.1, 3, 4]))
        # assert temp.dtype == np.float64
        exp = pd.Series([1, 1, 3, 4])
        self._assert_setitem_series_conversion(obj, 1.1, exp, np.int64)

        # int + complex -> complex
        exp = pd.Series([1, 1 + 1j, 3, 4])
        self._assert_setitem_series_conversion(obj, 1 + 1j, exp, np.complex128)

        # int + bool -> object
        exp = pd.Series([1, True, 3, 4])
        self._assert_setitem_series_conversion(obj, True, exp, np.object)

    def test_setitem_series_int8(self):
        # integer dtype coercion (no change)
        obj = pd.Series([1, 2, 3, 4], dtype=np.int8)
        assert obj.dtype == np.int8

        exp = pd.Series([1, 1, 3, 4], dtype=np.int8)
        self._assert_setitem_series_conversion(obj, np.int32(1), exp, np.int8)

        # BUG: it must be Series([1, 1, 3, 4], dtype=np.int16)
        exp = pd.Series([1, 0, 3, 4], dtype=np.int8)
        self._assert_setitem_series_conversion(obj, np.int16(2**9), exp,
                                               np.int8)

    def test_setitem_series_float64(self):
        obj = pd.Series([1.1, 2.2, 3.3, 4.4])
        assert obj.dtype == np.float64

        # float + int -> float
        exp = pd.Series([1.1, 1.0, 3.3, 4.4])
        self._assert_setitem_series_conversion(obj, 1, exp, np.float64)

        # float + float -> float
        exp = pd.Series([1.1, 1.1, 3.3, 4.4])
        self._assert_setitem_series_conversion(obj, 1.1, exp, np.float64)

        # float + complex -> complex
        exp = pd.Series([1.1, 1 + 1j, 3.3, 4.4])
        self._assert_setitem_series_conversion(obj, 1 + 1j, exp,
                                               np.complex128)

        # float + bool -> object
        exp = pd.Series([1.1, True, 3.3, 4.4])
        self._assert_setitem_series_conversion(obj, True, exp, np.object)

    def test_setitem_series_complex128(self):
        obj = pd.Series([1 + 1j, 2 + 2j, 3 + 3j, 4 + 4j])
        assert obj.dtype == np.complex128

        # complex + int -> complex
        exp = pd.Series([1 + 1j, 1, 3 + 3j, 4 + 4j])
        self._assert_setitem_series_conversion(obj, 1, exp, np.complex128)

        # complex + float -> complex
        exp = pd.Series([1 + 1j, 1.1, 3 + 3j, 4 + 4j])
        self._assert_setitem_series_conversion(obj, 1.1, exp, np.complex128)

        # complex + complex -> complex
        exp = pd.Series([1 + 1j, 1 + 1j, 3 + 3j, 4 + 4j])
        self._assert_setitem_series_conversion(obj, 1 + 1j, exp, np.complex128)

        # complex + bool -> object
        exp = pd.Series([1 + 1j, True, 3 + 3j, 4 + 4j])
        self._assert_setitem_series_conversion(obj, True, exp, np.object)

    def test_setitem_series_bool(self):
        obj = pd.Series([True, False, True, False])
        assert obj.dtype == np.bool

        # bool + int -> int
        # TODO_GH12747 The result must be int
        # tm.assert_series_equal(temp, pd.Series([1, 1, 1, 0]))
        # assert temp.dtype == np.int64
        exp = pd.Series([True, True, True, False])
        self._assert_setitem_series_conversion(obj, 1, exp, np.bool)

        # TODO_GH12747 The result must be int
        # assigning int greater than bool
        # tm.assert_series_equal(temp, pd.Series([1, 3, 1, 0]))
        # assert temp.dtype == np.int64
        exp = pd.Series([True, True, True, False])
        self._assert_setitem_series_conversion(obj, 3, exp, np.bool)

        # bool + float -> float
        # TODO_GH12747 The result must be float
        # tm.assert_series_equal(temp, pd.Series([1., 1.1, 1., 0.]))
        # assert temp.dtype == np.float64
        exp = pd.Series([True, True, True, False])
        self._assert_setitem_series_conversion(obj, 1.1, exp, np.bool)

        # bool + complex -> complex (buggy, results in bool)
        # TODO_GH12747 The result must be complex
        # tm.assert_series_equal(temp, pd.Series([1, 1 + 1j, 1, 0]))
        # assert temp.dtype == np.complex128
        exp = pd.Series([True, True, True, False])
        self._assert_setitem_series_conversion(obj, 1 + 1j, exp, np.bool)

        # bool + bool -> bool
        exp = pd.Series([True, True, True, False])
        self._assert_setitem_series_conversion(obj, True, exp, np.bool)

    def test_setitem_series_datetime64(self):
        obj = pd.Series([pd.Timestamp('2011-01-01'),
                         pd.Timestamp('2011-01-02'),
                         pd.Timestamp('2011-01-03'),
                         pd.Timestamp('2011-01-04')])
        assert obj.dtype == 'datetime64[ns]'

        # datetime64 + datetime64 -> datetime64
        exp = pd.Series([pd.Timestamp('2011-01-01'),
                         pd.Timestamp('2012-01-01'),
                         pd.Timestamp('2011-01-03'),
                         pd.Timestamp('2011-01-04')])
        self._assert_setitem_series_conversion(obj, pd.Timestamp('2012-01-01'),
                                               exp, 'datetime64[ns]')

        # datetime64 + int -> object
        exp = pd.Series([pd.Timestamp('2011-01-01'),
                         1,
                         pd.Timestamp('2011-01-03'),
                         pd.Timestamp('2011-01-04')])
        self._assert_setitem_series_conversion(obj, 1, exp, 'object')

        # datetime64 + object -> object
        exp = pd.Series([pd.Timestamp('2011-01-01'),
                         'x',
                         pd.Timestamp('2011-01-03'),
                         pd.Timestamp('2011-01-04')])
        self._assert_setitem_series_conversion(obj, 'x', exp, np.object)

    def test_setitem_series_datetime64tz(self):
        tz = 'US/Eastern'
        obj = pd.Series([pd.Timestamp('2011-01-01', tz=tz),
                         pd.Timestamp('2011-01-02', tz=tz),
                         pd.Timestamp('2011-01-03', tz=tz),
                         pd.Timestamp('2011-01-04', tz=tz)])
        assert obj.dtype == 'datetime64[ns, US/Eastern]'

        # datetime64tz + datetime64tz -> datetime64tz
        exp = pd.Series([pd.Timestamp('2011-01-01', tz=tz),
                         pd.Timestamp('2012-01-01', tz=tz),
                         pd.Timestamp('2011-01-03', tz=tz),
                         pd.Timestamp('2011-01-04', tz=tz)])
        value = pd.Timestamp('2012-01-01', tz=tz)
        self._assert_setitem_series_conversion(obj, value, exp,
                                               'datetime64[ns, US/Eastern]')

        # datetime64tz + datetime64tz (different tz) -> object
        exp = pd.Series([pd.Timestamp('2011-01-01', tz=tz),
                         pd.Timestamp('2012-01-01', tz='US/Pacific'),
                         pd.Timestamp('2011-01-03', tz=tz),
                         pd.Timestamp('2011-01-04', tz=tz)])
        value = pd.Timestamp('2012-01-01', tz='US/Pacific')
        self._assert_setitem_series_conversion(obj, value, exp, np.object)

        # datetime64tz + datetime64 -> object
        exp = pd.Series([pd.Timestamp('2011-01-01', tz=tz),
                         pd.Timestamp('2012-01-01'),
                         pd.Timestamp('2011-01-03', tz=tz),
                         pd.Timestamp('2011-01-04', tz=tz)])
        value = pd.Timestamp('2012-01-01')
        self._assert_setitem_series_conversion(obj, value, exp, np.object)

        # datetime64 + int -> object
        exp = pd.Series([pd.Timestamp('2011-01-01', tz=tz),
                         1,
                         pd.Timestamp('2011-01-03', tz=tz),
                         pd.Timestamp('2011-01-04', tz=tz)])
        self._assert_setitem_series_conversion(obj, 1, exp, np.object)

        # ToDo: add more tests once the above issue has been fixed

    def test_setitem_series_timedelta64(self):
        obj = pd.Series([pd.Timedelta('1 day'),
                         pd.Timedelta('2 day'),
                         pd.Timedelta('3 day'),
                         pd.Timedelta('4 day')])
        assert obj.dtype == 'timedelta64[ns]'

        # timedelta64 + timedelta64 -> timedelta64
        exp = pd.Series([pd.Timedelta('1 day'),
                         pd.Timedelta('12 day'),
                         pd.Timedelta('3 day'),
                         pd.Timedelta('4 day')])
        self._assert_setitem_series_conversion(obj, pd.Timedelta('12 day'),
                                               exp, 'timedelta64[ns]')

        # timedelta64 + int -> object
        exp = pd.Series([pd.Timedelta('1 day'),
                         1,
                         pd.Timedelta('3 day'),
                         pd.Timedelta('4 day')])
        self._assert_setitem_series_conversion(obj, 1, exp, np.object)

        # timedelta64 + object -> object
        exp = pd.Series([pd.Timedelta('1 day'),
                         'x',
                         pd.Timedelta('3 day'),
                         pd.Timedelta('4 day')])
        self._assert_setitem_series_conversion(obj, 'x', exp, np.object)

    def test_setitem_series_period(self):
        pass

    def _assert_setitem_index_conversion(self, original_series, loc_key,
                                         expected_index, expected_dtype):
        """ test index's coercion triggered by assign key """
        temp = original_series.copy()
        temp[loc_key] = 5
        exp = pd.Series([1, 2, 3, 4, 5], index=expected_index)
        tm.assert_series_equal(temp, exp)
        # check dtype explicitly for sure
        assert temp.index.dtype == expected_dtype

        temp = original_series.copy()
        temp.loc[loc_key] = 5
        exp = pd.Series([1, 2, 3, 4, 5], index=expected_index)
        tm.assert_series_equal(temp, exp)
        # check dtype explicitly for sure
        assert temp.index.dtype == expected_dtype

    def test_setitem_index_object(self):
        obj = pd.Series([1, 2, 3, 4], index=list('abcd'))
        assert obj.index.dtype == np.object

        # object + object -> object
        exp_index = pd.Index(list('abcdx'))
        self._assert_setitem_index_conversion(obj, 'x', exp_index, np.object)

        # object + int -> IndexError, regarded as location
        temp = obj.copy()
        with pytest.raises(IndexError):
            temp[5] = 5

        # object + float -> object
        exp_index = pd.Index(['a', 'b', 'c', 'd', 1.1])
        self._assert_setitem_index_conversion(obj, 1.1, exp_index, np.object)

    def test_setitem_index_int64(self):
        # tests setitem with non-existing numeric key
        obj = pd.Series([1, 2, 3, 4])
        assert obj.index.dtype == np.int64

        # int + int -> int
        exp_index = pd.Index([0, 1, 2, 3, 5])
        self._assert_setitem_index_conversion(obj, 5, exp_index, np.int64)

        # int + float -> float
        exp_index = pd.Index([0, 1, 2, 3, 1.1])
        self._assert_setitem_index_conversion(obj, 1.1, exp_index, np.float64)

        # int + object -> object
        exp_index = pd.Index([0, 1, 2, 3, 'x'])
        self._assert_setitem_index_conversion(obj, 'x', exp_index, np.object)

    def test_setitem_index_float64(self):
        # tests setitem with non-existing numeric key
        obj = pd.Series([1, 2, 3, 4], index=[1.1, 2.1, 3.1, 4.1])
        assert obj.index.dtype == np.float64

        # float + int -> int
        temp = obj.copy()
        # TODO_GH12747 The result must be float
        with pytest.raises(IndexError):
            temp[5] = 5

        # float + float -> float
        exp_index = pd.Index([1.1, 2.1, 3.1, 4.1, 5.1])
        self._assert_setitem_index_conversion(obj, 5.1, exp_index, np.float64)

        # float + object -> object
        exp_index = pd.Index([1.1, 2.1, 3.1, 4.1, 'x'])
        self._assert_setitem_index_conversion(obj, 'x', exp_index, np.object)

    def test_setitem_index_complex128(self):
        pass

    def test_setitem_index_bool(self):
        pass

    def test_setitem_index_datetime64(self):
        pass

    def test_setitem_index_datetime64tz(self):
        pass

    def test_setitem_index_timedelta64(self):
        pass

    def test_setitem_index_period(self):
        pass


class TestInsertIndexCoercion(CoercionBase):

    klasses = ['index']
    method = 'insert'

    def _assert_insert_conversion(self, original, value,
                                  expected, expected_dtype):
        """ test coercion triggered by insert """
        target = original.copy()
        res = target.insert(1, value)
        tm.assert_index_equal(res, expected)
        assert res.dtype == expected_dtype

    @pytest.mark.parametrize("insert, coerced_val, coerced_dtype", [
        (1, 1, np.object),
        (1.1, 1.1, np.object),
        (False, False, np.object),
        ('x', 'x', np.object)])
    def test_insert_index_object(self, insert, coerced_val, coerced_dtype):
        obj = pd.Index(list('abcd'))
        assert obj.dtype == np.object

        exp = pd.Index(['a', coerced_val, 'b', 'c', 'd'])
        self._assert_insert_conversion(obj, insert, exp, coerced_dtype)

    @pytest.mark.parametrize("insert, coerced_val, coerced_dtype", [
        (1, 1, np.int64),
        (1.1, 1.1, np.float64),
        (False, 0, int),
        ('x', 'x', np.object)])
    def test_insert_index_int64(self, insert, coerced_val, coerced_dtype):
        obj = pd.Int64Index([1, 2, 3, 4])
        assert obj.dtype == np.int64

        exp = pd.Index([1, coerced_val, 2, 3, 4])
        self._assert_insert_conversion(obj, insert, exp, coerced_dtype)

    @pytest.mark.parametrize("insert, coerced_val, coerced_dtype", [
        (1, 1., np.float64),
        (1.1, 1.1, np.float64),
        (False, 0., np.float64),
        ('x', 'x', np.object)])
    def test_insert_index_float64(self, insert, coerced_val, coerced_dtype):
        obj = pd.Float64Index([1., 2., 3., 4.])
        assert obj.dtype == np.float64

        # float + int -> int
        exp = pd.Index([1., coerced_val, 2., 3., 4.])
        self._assert_insert_conversion(obj, insert, exp, coerced_dtype)

    def test_insert_index_complex128(self):
        pass

    def test_insert_index_bool(self):
        pass

    def test_insert_index_datetime64(self):
        obj = pd.DatetimeIndex(['2011-01-01', '2011-01-02', '2011-01-03',
                                '2011-01-04'])
        assert obj.dtype == 'datetime64[ns]'

        # datetime64 + datetime64 => datetime64
        exp = pd.DatetimeIndex(['2011-01-01', '2012-01-01', '2011-01-02',
                                '2011-01-03', '2011-01-04'])
        self._assert_insert_conversion(obj, pd.Timestamp('2012-01-01'),
                                       exp, 'datetime64[ns]')

        # ToDo: must coerce to object
        msg = "Passed item and index have different timezone"
        with tm.assert_raises_regex(ValueError, msg):
            obj.insert(1, pd.Timestamp('2012-01-01', tz='US/Eastern'))

        # ToDo: must coerce to object
        msg = "cannot insert DatetimeIndex with incompatible label"
        with tm.assert_raises_regex(TypeError, msg):
            obj.insert(1, 1)

    def test_insert_index_datetime64tz(self):
        obj = pd.DatetimeIndex(['2011-01-01', '2011-01-02', '2011-01-03',
                                '2011-01-04'], tz='US/Eastern')
        assert obj.dtype == 'datetime64[ns, US/Eastern]'

        # datetime64tz + datetime64tz => datetime64
        exp = pd.DatetimeIndex(['2011-01-01', '2012-01-01', '2011-01-02',
                                '2011-01-03', '2011-01-04'], tz='US/Eastern')
        val = pd.Timestamp('2012-01-01', tz='US/Eastern')
        self._assert_insert_conversion(obj, val, exp,
                                       'datetime64[ns, US/Eastern]')

        # ToDo: must coerce to object
        msg = "Passed item and index have different timezone"
        with tm.assert_raises_regex(ValueError, msg):
            obj.insert(1, pd.Timestamp('2012-01-01'))

        # ToDo: must coerce to object
        msg = "Passed item and index have different timezone"
        with tm.assert_raises_regex(ValueError, msg):
            obj.insert(1, pd.Timestamp('2012-01-01', tz='Asia/Tokyo'))

        # ToDo: must coerce to object
        msg = "cannot insert DatetimeIndex with incompatible label"
        with tm.assert_raises_regex(TypeError, msg):
            obj.insert(1, 1)

    def test_insert_index_timedelta64(self):
        obj = pd.TimedeltaIndex(['1 day', '2 day', '3 day', '4 day'])
        assert obj.dtype == 'timedelta64[ns]'

        # timedelta64 + timedelta64 => timedelta64
        exp = pd.TimedeltaIndex(['1 day', '10 day', '2 day', '3 day', '4 day'])
        self._assert_insert_conversion(obj, pd.Timedelta('10 day'),
                                       exp, 'timedelta64[ns]')

        # ToDo: must coerce to object
        msg = "cannot insert TimedeltaIndex with incompatible label"
        with tm.assert_raises_regex(TypeError, msg):
            obj.insert(1, pd.Timestamp('2012-01-01'))

        # ToDo: must coerce to object
        msg = "cannot insert TimedeltaIndex with incompatible label"
        with tm.assert_raises_regex(TypeError, msg):
            obj.insert(1, 1)

    @pytest.mark.parametrize("insert, coerced_val, coerced_dtype", [
        (pd.Period('2012-01', freq='M'), '2012-01', 'period[M]'),
        (pd.Timestamp('2012-01-01'), pd.Timestamp('2012-01-01'), np.object),
        (1, 1, np.object),
        ('x', 'x', np.object)])
    def test_insert_index_period(self, insert, coerced_val, coerced_dtype):
        obj = pd.PeriodIndex(['2011-01', '2011-02', '2011-03', '2011-04'],
                             freq='M')
        assert obj.dtype == 'period[M]'

        if isinstance(insert, pd.Period):
            index_type = pd.PeriodIndex
        else:
            index_type = pd.Index

        exp = index_type([pd.Period('2011-01', freq='M'),
                        coerced_val,
                        pd.Period('2011-02', freq='M'),
                        pd.Period('2011-03', freq='M'),
                        pd.Period('2011-04', freq='M')], freq='M')        
        self._assert_insert_conversion(obj, insert, exp, coerced_dtype)


class TestWhereCoercion(CoercionBase):

    method = 'where'

    def _assert_where_conversion(self, original, cond, values,
                                 expected, expected_dtype):
        """ test coercion triggered by where """
        target = original.copy()
        res = target.where(cond, values)
        self._assert(res, expected, expected_dtype)

    def _where_object_common(self, klass):
        obj = klass(list('abcd'))
        assert obj.dtype == np.object
        cond = klass([True, False, True, False])

        # object + int -> object
        exp = klass(['a', 1, 'c', 1])
        self._assert_where_conversion(obj, cond, 1, exp, np.object)

        values = klass([5, 6, 7, 8])
        exp = klass(['a', 6, 'c', 8])
        self._assert_where_conversion(obj, cond, values, exp, np.object)

        # object + float -> object
        exp = klass(['a', 1.1, 'c', 1.1])
        self._assert_where_conversion(obj, cond, 1.1, exp, np.object)

        values = klass([5.5, 6.6, 7.7, 8.8])
        exp = klass(['a', 6.6, 'c', 8.8])
        self._assert_where_conversion(obj, cond, values, exp, np.object)

        # object + complex -> object
        exp = klass(['a', 1 + 1j, 'c', 1 + 1j])
        self._assert_where_conversion(obj, cond, 1 + 1j, exp, np.object)

        values = klass([5 + 5j, 6 + 6j, 7 + 7j, 8 + 8j])
        exp = klass(['a', 6 + 6j, 'c', 8 + 8j])
        self._assert_where_conversion(obj, cond, values, exp, np.object)

        if klass is pd.Series:
            exp = klass(['a', 1, 'c', 1])
            self._assert_where_conversion(obj, cond, True, exp, np.object)

            values = klass([True, False, True, True])
            exp = klass(['a', 0, 'c', 1])
            self._assert_where_conversion(obj, cond, values, exp, np.object)
        elif klass is pd.Index:
            # object + bool -> object
            exp = klass(['a', True, 'c', True])
            self._assert_where_conversion(obj, cond, True, exp, np.object)

            values = klass([True, False, True, True])
            exp = klass(['a', False, 'c', True])
            self._assert_where_conversion(obj, cond, values, exp, np.object)
        else:
            NotImplementedError

    def test_where_series_object(self):
        self._where_object_common(pd.Series)

    def test_where_index_object(self):
        self._where_object_common(pd.Index)

    def _where_int64_common(self, klass):
        obj = klass([1, 2, 3, 4])
        assert obj.dtype == np.int64
        cond = klass([True, False, True, False])

        # int + int -> int
        exp = klass([1, 1, 3, 1])
        self._assert_where_conversion(obj, cond, 1, exp, np.int64)

        values = klass([5, 6, 7, 8])
        exp = klass([1, 6, 3, 8])
        self._assert_where_conversion(obj, cond, values, exp, np.int64)

        # int + float -> float
        exp = klass([1, 1.1, 3, 1.1])
        self._assert_where_conversion(obj, cond, 1.1, exp, np.float64)

        values = klass([5.5, 6.6, 7.7, 8.8])
        exp = klass([1, 6.6, 3, 8.8])
        self._assert_where_conversion(obj, cond, values, exp, np.float64)

        # int + complex -> complex
        if klass is pd.Series:
            exp = klass([1, 1 + 1j, 3, 1 + 1j])
            self._assert_where_conversion(obj, cond, 1 + 1j, exp,
                                          np.complex128)

            values = klass([5 + 5j, 6 + 6j, 7 + 7j, 8 + 8j])
            exp = klass([1, 6 + 6j, 3, 8 + 8j])
            self._assert_where_conversion(obj, cond, values, exp,
                                          np.complex128)

        # int + bool -> object
        exp = klass([1, True, 3, True])
        self._assert_where_conversion(obj, cond, True, exp, np.object)

        values = klass([True, False, True, True])
        exp = klass([1, False, 3, True])
        self._assert_where_conversion(obj, cond, values, exp, np.object)

    def test_where_series_int64(self):
        self._where_int64_common(pd.Series)

    def test_where_index_int64(self):
        self._where_int64_common(pd.Index)

    def _where_float64_common(self, klass):
        obj = klass([1.1, 2.2, 3.3, 4.4])
        assert obj.dtype == np.float64
        cond = klass([True, False, True, False])

        # float + int -> float
        exp = klass([1.1, 1.0, 3.3, 1.0])
        self._assert_where_conversion(obj, cond, 1, exp, np.float64)

        values = klass([5, 6, 7, 8])
        exp = klass([1.1, 6.0, 3.3, 8.0])
        self._assert_where_conversion(obj, cond, values, exp, np.float64)

        # float + float -> float
        exp = klass([1.1, 1.1, 3.3, 1.1])
        self._assert_where_conversion(obj, cond, 1.1, exp, np.float64)

        values = klass([5.5, 6.6, 7.7, 8.8])
        exp = klass([1.1, 6.6, 3.3, 8.8])
        self._assert_where_conversion(obj, cond, values, exp, np.float64)

        # float + complex -> complex
        if klass is pd.Series:
            exp = klass([1.1, 1 + 1j, 3.3, 1 + 1j])
            self._assert_where_conversion(obj, cond, 1 + 1j, exp,
                                          np.complex128)

            values = klass([5 + 5j, 6 + 6j, 7 + 7j, 8 + 8j])
            exp = klass([1.1, 6 + 6j, 3.3, 8 + 8j])
            self._assert_where_conversion(obj, cond, values, exp,
                                          np.complex128)

        # float + bool -> object
        exp = klass([1.1, True, 3.3, True])
        self._assert_where_conversion(obj, cond, True, exp, np.object)

        values = klass([True, False, True, True])
        exp = klass([1.1, False, 3.3, True])
        self._assert_where_conversion(obj, cond, values, exp, np.object)

    def test_where_series_float64(self):
        self._where_float64_common(pd.Series)

    def test_where_index_float64(self):
        self._where_float64_common(pd.Index)

    def test_where_series_complex128(self):
        obj = pd.Series([1 + 1j, 2 + 2j, 3 + 3j, 4 + 4j])
        assert obj.dtype == np.complex128
        cond = pd.Series([True, False, True, False])

        # complex + int -> complex
        exp = pd.Series([1 + 1j, 1, 3 + 3j, 1])
        self._assert_where_conversion(obj, cond, 1, exp, np.complex128)

        values = pd.Series([5, 6, 7, 8])
        exp = pd.Series([1 + 1j, 6.0, 3 + 3j, 8.0])
        self._assert_where_conversion(obj, cond, values, exp, np.complex128)

        # complex + float -> complex
        exp = pd.Series([1 + 1j, 1.1, 3 + 3j, 1.1])
        self._assert_where_conversion(obj, cond, 1.1, exp, np.complex128)

        values = pd.Series([5.5, 6.6, 7.7, 8.8])
        exp = pd.Series([1 + 1j, 6.6, 3 + 3j, 8.8])
        self._assert_where_conversion(obj, cond, values, exp, np.complex128)

        # complex + complex -> complex
        exp = pd.Series([1 + 1j, 1 + 1j, 3 + 3j, 1 + 1j])
        self._assert_where_conversion(obj, cond, 1 + 1j, exp, np.complex128)

        values = pd.Series([5 + 5j, 6 + 6j, 7 + 7j, 8 + 8j])
        exp = pd.Series([1 + 1j, 6 + 6j, 3 + 3j, 8 + 8j])
        self._assert_where_conversion(obj, cond, values, exp, np.complex128)

        # complex + bool -> object
        exp = pd.Series([1 + 1j, True, 3 + 3j, True])
        self._assert_where_conversion(obj, cond, True, exp, np.object)

        values = pd.Series([True, False, True, True])
        exp = pd.Series([1 + 1j, False, 3 + 3j, True])
        self._assert_where_conversion(obj, cond, values, exp, np.object)

    def test_where_index_complex128(self):
        pass

    def test_where_series_bool(self):

        obj = pd.Series([True, False, True, False])
        assert obj.dtype == np.bool
        cond = pd.Series([True, False, True, False])

        # bool + int -> object
        exp = pd.Series([True, 1, True, 1])
        self._assert_where_conversion(obj, cond, 1, exp, np.object)

        values = pd.Series([5, 6, 7, 8])
        exp = pd.Series([True, 6, True, 8])
        self._assert_where_conversion(obj, cond, values, exp, np.object)

        # bool + float -> object
        exp = pd.Series([True, 1.1, True, 1.1])
        self._assert_where_conversion(obj, cond, 1.1, exp, np.object)

        values = pd.Series([5.5, 6.6, 7.7, 8.8])
        exp = pd.Series([True, 6.6, True, 8.8])
        self._assert_where_conversion(obj, cond, values, exp, np.object)

        # bool + complex -> object
        exp = pd.Series([True, 1 + 1j, True, 1 + 1j])
        self._assert_where_conversion(obj, cond, 1 + 1j, exp, np.object)

        values = pd.Series([5 + 5j, 6 + 6j, 7 + 7j, 8 + 8j])
        exp = pd.Series([True, 6 + 6j, True, 8 + 8j])
        self._assert_where_conversion(obj, cond, values, exp, np.object)

        # bool + bool -> bool
        exp = pd.Series([True, True, True, True])
        self._assert_where_conversion(obj, cond, True, exp, np.bool)

        values = pd.Series([True, False, True, True])
        exp = pd.Series([True, False, True, True])
        self._assert_where_conversion(obj, cond, values, exp, np.bool)

    def test_where_index_bool(self):
        pass

    def test_where_series_datetime64(self):
        obj = pd.Series([pd.Timestamp('2011-01-01'),
                         pd.Timestamp('2011-01-02'),
                         pd.Timestamp('2011-01-03'),
                         pd.Timestamp('2011-01-04')])
        assert obj.dtype == 'datetime64[ns]'
        cond = pd.Series([True, False, True, False])

        # datetime64 + datetime64 -> datetime64
        exp = pd.Series([pd.Timestamp('2011-01-01'),
                         pd.Timestamp('2012-01-01'),
                         pd.Timestamp('2011-01-03'),
                         pd.Timestamp('2012-01-01')])
        self._assert_where_conversion(obj, cond, pd.Timestamp('2012-01-01'),
                                      exp, 'datetime64[ns]')

        values = pd.Series([pd.Timestamp('2012-01-01'),
                            pd.Timestamp('2012-01-02'),
                            pd.Timestamp('2012-01-03'),
                            pd.Timestamp('2012-01-04')])
        exp = pd.Series([pd.Timestamp('2011-01-01'),
                         pd.Timestamp('2012-01-02'),
                         pd.Timestamp('2011-01-03'),
                         pd.Timestamp('2012-01-04')])
        self._assert_where_conversion(obj, cond, values, exp, 'datetime64[ns]')

        # datetime64 + datetime64tz -> object
        exp = pd.Series([pd.Timestamp('2011-01-01'),
                         pd.Timestamp('2012-01-01', tz='US/Eastern'),
                         pd.Timestamp('2011-01-03'),
                         pd.Timestamp('2012-01-01', tz='US/Eastern')])
        self._assert_where_conversion(
            obj, cond,
            pd.Timestamp('2012-01-01', tz='US/Eastern'),
            exp, np.object)

        # ToDo: do not coerce to UTC, must be object
        values = pd.Series([pd.Timestamp('2012-01-01', tz='US/Eastern'),
                            pd.Timestamp('2012-01-02', tz='US/Eastern'),
                            pd.Timestamp('2012-01-03', tz='US/Eastern'),
                            pd.Timestamp('2012-01-04', tz='US/Eastern')])
        exp = pd.Series([pd.Timestamp('2011-01-01'),
                         pd.Timestamp('2012-01-02 05:00'),
                         pd.Timestamp('2011-01-03'),
                         pd.Timestamp('2012-01-04 05:00')])
        self._assert_where_conversion(obj, cond, values, exp, 'datetime64[ns]')

    def test_where_index_datetime64(self):
        obj = pd.Index([pd.Timestamp('2011-01-01'),
                        pd.Timestamp('2011-01-02'),
                        pd.Timestamp('2011-01-03'),
                        pd.Timestamp('2011-01-04')])
        assert obj.dtype == 'datetime64[ns]'
        cond = pd.Index([True, False, True, False])

        # datetime64 + datetime64 -> datetime64
        # must support scalar
        msg = "cannot coerce a Timestamp with a tz on a naive Block"
        with pytest.raises(TypeError):
            obj.where(cond, pd.Timestamp('2012-01-01'))

        values = pd.Index([pd.Timestamp('2012-01-01'),
                           pd.Timestamp('2012-01-02'),
                           pd.Timestamp('2012-01-03'),
                           pd.Timestamp('2012-01-04')])
        exp = pd.Index([pd.Timestamp('2011-01-01'),
                        pd.Timestamp('2012-01-02'),
                        pd.Timestamp('2011-01-03'),
                        pd.Timestamp('2012-01-04')])
        self._assert_where_conversion(obj, cond, values, exp, 'datetime64[ns]')

        # ToDo: coerce to object
        msg = ("Index\\(\\.\\.\\.\\) must be called with a collection "
               "of some kind")
        with tm.assert_raises_regex(TypeError, msg):
            obj.where(cond, pd.Timestamp('2012-01-01', tz='US/Eastern'))

        # ToDo: do not ignore timezone, must be object
        values = pd.Index([pd.Timestamp('2012-01-01', tz='US/Eastern'),
                           pd.Timestamp('2012-01-02', tz='US/Eastern'),
                           pd.Timestamp('2012-01-03', tz='US/Eastern'),
                           pd.Timestamp('2012-01-04', tz='US/Eastern')])
        exp = pd.Index([pd.Timestamp('2011-01-01'),
                        pd.Timestamp('2012-01-02'),
                        pd.Timestamp('2011-01-03'),
                        pd.Timestamp('2012-01-04')])
        self._assert_where_conversion(obj, cond, values, exp, 'datetime64[ns]')

    def test_where_series_datetime64tz(self):
        pass

    def test_where_series_timedelta64(self):
        pass

    def test_where_series_period(self):
        pass

    def test_where_index_datetime64tz(self):
        pass

    def test_where_index_timedelta64(self):
        pass

    def test_where_index_period(self):
        pass


class TestFillnaSeriesCoercion(CoercionBase):

    # not indexing, but place here for consisntency

    method = 'fillna'

    def test_has_comprehensive_tests(self):
        pass

    def _assert_fillna_conversion(self, original, value,
                                  expected, expected_dtype):
        """ test coercion triggered by fillna """
        target = original.copy()
        res = target.fillna(value)
        self._assert(res, expected, expected_dtype)

    @pytest.mark.parametrize("klass", [pd.Series, pd.Index])
    @pytest.mark.parametrize("fill_val, fill_dtype", [
        (1, np.object),
        (1.1, np.object),
        (1 + 1j, np.object),
        (True, np.object)])
    def test_fillna_object(self, klass, fill_val, fill_dtype):
        obj = klass(['a', np.nan, 'c', 'd'])
        assert obj.dtype == np.object

        exp = klass(['a', fill_val, 'c', 'd'])
        self._assert_fillna_conversion(obj, fill_val, exp, fill_dtype)

    @pytest.mark.parametrize("klass", [pd.Series, pd.Index])
    @pytest.mark.parametrize("fill_val,filled_val,fill_dtype", [
        (1, 1.0, np.float64),
        (1.1, 1.1, np.float64),
        (1 + 1j, 1 + 1j, np.complex128),
        (True, True, np.object)])
    def test_fillna_float64(self, klass, fill_val, filled_val, fill_dtype):
        obj = klass([1.1, np.nan, 3.3, 4.4])
        assert obj.dtype == np.float64

        exp = klass([1.1, filled_val, 3.3, 4.4])
        # float + complex -> we don't support a complex Index
        # complex for Series,
        # object for Index
        if fill_dtype == np.complex128 and klass == pd.Index:
            fill_dtype = np.object
        self._assert_fillna_conversion(obj, fill_val, exp, fill_dtype)

    @pytest.mark.parametrize("fill_val,fill_dtype", [
        (1, np.complex128),
        (1.1, np.complex128),
        (1 + 1j, np.complex128),
        (True, np.object)])
    def test_fillna_series_complex128(self, fill_val, fill_dtype):
        obj = pd.Series([1 + 1j, np.nan, 3 + 3j, 4 + 4j])
        assert obj.dtype == np.complex128

        exp = pd.Series([1 + 1j, fill_val, 3 + 3j, 4 + 4j])
        self._assert_fillna_conversion(obj, fill_val, exp, fill_dtype)

    @pytest.mark.parametrize("klass", [pd.Series, pd.Index])
    @pytest.mark.parametrize("fill_val,fill_dtype", [
        (pd.Timestamp('2012-01-01'), 'datetime64[ns]'),
        (pd.Timestamp('2012-01-01', tz='US/Eastern'), np.object),
        (1, np.object),
        ('x', np.object)])
    def test_fillna_datetime64(self, klass, fill_val, fill_dtype):
        obj = klass([pd.Timestamp('2011-01-01'),
                     pd.NaT,
                     pd.Timestamp('2011-01-03'),
                     pd.Timestamp('2011-01-04')])
        assert obj.dtype == 'datetime64[ns]'

        exp = klass([pd.Timestamp('2011-01-01'),
                     fill_val,
                     pd.Timestamp('2011-01-03'),
                     pd.Timestamp('2011-01-04')])
        self._assert_fillna_conversion(obj, fill_val, exp, fill_dtype)

    @pytest.mark.parametrize("klass", [pd.Series, pd.Index])
    @pytest.mark.parametrize("fill_val,fill_dtype", [
        (pd.Timestamp('2012-01-01', tz='US/Eastern'),
         'datetime64[ns, US/Eastern]'),
        (pd.Timestamp('2012-01-01'), np.object),
        (pd.Timestamp('2012-01-01', tz='Asia/Tokyo'), np.object),
        (1, np.object),
        ('x', np.object)])
    def test_fillna_datetime64tz(self, klass, fill_val, fill_dtype):
        tz = 'US/Eastern'

        obj = klass([pd.Timestamp('2011-01-01', tz=tz),
                     pd.NaT,
                     pd.Timestamp('2011-01-03', tz=tz),
                     pd.Timestamp('2011-01-04', tz=tz)])
        assert obj.dtype == 'datetime64[ns, US/Eastern]'

        exp = klass([pd.Timestamp('2011-01-01', tz=tz),
                     fill_val,
                     pd.Timestamp('2011-01-03', tz=tz),
                     pd.Timestamp('2011-01-04', tz=tz)])
        self._assert_fillna_conversion(obj, fill_val, exp, fill_dtype)


@pytest.mark.parametrize('how', ['dict', 'series'])
@pytest.mark.parametrize('to_key', [
    'object', 'int64', 'float64', 'complex128', 'bool', 'datetime64[ns]',
    'datetime64[ns, UTC]', 'datetime64[ns, US/Eastern]', 'timedelta64[ns]'])
@pytest.mark.parametrize('from_key', [
    'object', 'int64', 'float64', 'complex128', 'bool', 'datetime64[ns]',
    'datetime64[ns, UTC]', 'datetime64[ns, US/Eastern]', 'timedelta64[ns]'])
class TestReplaceSeriesCoercion(object):

    rep = {}
    rep['object'] = ['a', 'b']
    rep['int64'] = [4, 5]
    rep['float64'] = [1.1, 2.2]
    rep['complex128'] = [1 + 1j, 2 + 2j]
    rep['bool'] = [True, False]
    rep['datetime64[ns]'] = [pd.Timestamp('2011-01-01'),
                             pd.Timestamp('2011-01-03')]

    for tz in ['UTC', 'US/Eastern']:
        # to test tz => different tz replacement
        key = 'datetime64[ns, {0}]'.format(tz)
        rep[key] = [pd.Timestamp('2011-01-01', tz=tz),
                    pd.Timestamp('2011-01-03', tz=tz)]

    rep['timedelta64[ns]'] = [pd.Timedelta('1 day'),
                              pd.Timedelta('2 day')]

    def test_replace_series(self, how, to_key, from_key):
        if from_key == 'bool' and how == 'series' and compat.PY3:
            # doesn't work in PY3, though ...dict_from_bool works fine
            pytest.skip("doesn't work as in PY3")

        index = pd.Index([3, 4], name='xxx')
        obj = pd.Series(self.rep[from_key], index=index, name='yyy')
        assert obj.dtype == from_key

        if (from_key.startswith('datetime') and to_key.startswith('datetime')):
            # different tz, currently mask_missing raises SystemError
            return

        if how == 'dict':
            replacer = dict(zip(self.rep[from_key], self.rep[to_key]))
        elif how == 'series':
            replacer = pd.Series(self.rep[to_key], index=self.rep[from_key])
        else:
            raise ValueError

        result = obj.replace(replacer)

        if ((from_key == 'float64' and to_key in ('int64')) or
            (from_key == 'complex128' and
             to_key in ('int64', 'float64'))):

            if compat.is_platform_32bit() or compat.is_platform_windows():
                pytest.skip("32-bit platform buggy: {0} -> {1}".format
                            (from_key, to_key))

            # Expected: do not downcast by replacement
            exp = pd.Series(self.rep[to_key], index=index,
                            name='yyy', dtype=from_key)

        else:
            exp = pd.Series(self.rep[to_key], index=index, name='yyy')
            assert exp.dtype == to_key

        tm.assert_series_equal(result, exp)
