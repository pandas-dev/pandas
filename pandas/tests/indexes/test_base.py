# -*- coding: utf-8 -*-

import pytest

from datetime import datetime, timedelta

from collections import defaultdict

import pandas.util.testing as tm
from pandas.core.dtypes.generic import ABCIndex
from pandas.core.dtypes.common import is_unsigned_integer_dtype
from pandas.core.indexes.api import Index, MultiIndex
from pandas.tests.indexes.common import Base

from pandas.compat import (range, lrange, lzip, u,
                           text_type, zip, PY3, PY36, PYPY)
import operator
import numpy as np

from pandas import (period_range, date_range, Series,
                    DataFrame, Float64Index, Int64Index, UInt64Index,
                    CategoricalIndex, DatetimeIndex, TimedeltaIndex,
                    PeriodIndex, RangeIndex, isna)
from pandas.core.index import _get_combined_index, _ensure_index_from_sequences
from pandas.util.testing import assert_almost_equal
from pandas.compat.numpy import np_datetime64_compat

import pandas.core.config as cf

from pandas.core.indexes.datetimes import _to_m8

import pandas as pd
from pandas._libs.tslib import Timestamp


class TestIndex(Base):
    _holder = Index

    def setup_method(self, method):
        self.indices = dict(unicodeIndex=tm.makeUnicodeIndex(100),
                            strIndex=tm.makeStringIndex(100),
                            dateIndex=tm.makeDateIndex(100),
                            periodIndex=tm.makePeriodIndex(100),
                            tdIndex=tm.makeTimedeltaIndex(100),
                            intIndex=tm.makeIntIndex(100),
                            uintIndex=tm.makeUIntIndex(100),
                            rangeIndex=tm.makeRangeIndex(100),
                            floatIndex=tm.makeFloatIndex(100),
                            boolIndex=Index([True, False]),
                            catIndex=tm.makeCategoricalIndex(100),
                            empty=Index([]),
                            tuples=MultiIndex.from_tuples(lzip(
                                ['foo', 'bar', 'baz'], [1, 2, 3])),
                            repeats=Index([0, 0, 1, 1, 2, 2]))
        self.setup_indices()

    def create_index(self):
        return Index(list('abcde'))

    def generate_index_types(self, skip_index_keys=[]):
        """
        Return a generator of the various index types, leaving
        out the ones with a key in skip_index_keys
        """
        for key, index in self.indices.items():
            if key not in skip_index_keys:
                yield key, index

    def test_can_hold_identifiers(self):
        idx = self.create_index()
        key = idx[0]
        assert idx._can_hold_identifiers_and_holds_name(key) is True

    def test_new_axis(self):
        new_index = self.dateIndex[None, :]
        assert new_index.ndim == 2
        assert isinstance(new_index, np.ndarray)

    def test_copy_and_deepcopy(self, indices):
        super(TestIndex, self).test_copy_and_deepcopy(indices)

        new_copy2 = self.intIndex.copy(dtype=int)
        assert new_copy2.dtype.kind == 'i'

    @pytest.mark.parametrize("attr", ['strIndex', 'dateIndex'])
    def test_constructor_regular(self, attr):
        # regular instance creation
        index = getattr(self, attr)
        tm.assert_contains_all(index, index)

    def test_constructor_casting(self):
        # casting
        arr = np.array(self.strIndex)
        index = Index(arr)
        tm.assert_contains_all(arr, index)
        tm.assert_index_equal(self.strIndex, index)

    def test_constructor_copy(self):
        # copy
        arr = np.array(self.strIndex)
        index = Index(arr, copy=True, name='name')
        assert isinstance(index, Index)
        assert index.name == 'name'
        tm.assert_numpy_array_equal(arr, index.values)
        arr[0] = "SOMEBIGLONGSTRING"
        assert index[0] != "SOMEBIGLONGSTRING"

        # what to do here?
        # arr = np.array(5.)
        # pytest.raises(Exception, arr.view, Index)

    def test_constructor_corner(self):
        # corner case
        pytest.raises(TypeError, Index, 0)

    @pytest.mark.parametrize("index_vals", [
        [('A', 1), 'B'], ['B', ('A', 1)]])
    def test_construction_list_mixed_tuples(self, index_vals):
        # see gh-10697: if we are constructing from a mixed list of tuples,
        # make sure that we are independent of the sorting order.
        index = Index(index_vals)
        assert isinstance(index, Index)
        assert not isinstance(index, MultiIndex)

    @pytest.mark.parametrize('na_value', [None, np.nan])
    @pytest.mark.parametrize('vtype', [list, tuple, iter])
    def test_construction_list_tuples_nan(self, na_value, vtype):
        # GH 18505 : valid tuples containing NaN
        values = [(1, 'two'), (3., na_value)]
        result = Index(vtype(values))
        expected = MultiIndex.from_tuples(values)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("cast_as_obj", [True, False])
    @pytest.mark.parametrize("index", [
        pd.date_range('2015-01-01 10:00', freq='D', periods=3,
                      tz='US/Eastern'),  # DTI with tz
        pd.date_range('2015-01-01 10:00', freq='D', periods=3),  # DTI no tz
        pd.timedelta_range('1 days', freq='D', periods=3),  # td
        pd.period_range('2015-01-01', freq='D', periods=3)  # period
    ])
    def test_constructor_from_index_dtlike(self, cast_as_obj, index):
        if cast_as_obj:
            result = pd.Index(index.astype(object))
        else:
            result = pd.Index(index)

        tm.assert_index_equal(result, index)

        if isinstance(index, pd.DatetimeIndex) and hasattr(index, 'tz'):
            assert result.tz == index.tz

    @pytest.mark.parametrize("index,has_tz", [
        (pd.date_range('2015-01-01 10:00', freq='D', periods=3,
                       tz='US/Eastern'), True),  # datetimetz
        (pd.timedelta_range('1 days', freq='D', periods=3), False),  # td
        (pd.period_range('2015-01-01', freq='D', periods=3), False)  # period
    ])
    def test_constructor_from_series_dtlike(self, index, has_tz):
        result = pd.Index(pd.Series(index))
        tm.assert_index_equal(result, index)

        if has_tz:
            assert result.tz == index.tz

    @pytest.mark.parametrize("klass", [Index, DatetimeIndex])
    def test_constructor_from_series(self, klass):
        expected = DatetimeIndex([Timestamp('20110101'), Timestamp('20120101'),
                                  Timestamp('20130101')])
        s = Series([Timestamp('20110101'), Timestamp('20120101'),
                    Timestamp('20130101')])
        result = klass(s)
        tm.assert_index_equal(result, expected)

    def test_constructor_from_series_freq(self):
        # GH 6273
        # create from a series, passing a freq
        dts = ['1-1-1990', '2-1-1990', '3-1-1990', '4-1-1990', '5-1-1990']
        expected = DatetimeIndex(dts, freq='MS')

        s = Series(pd.to_datetime(dts))
        result = DatetimeIndex(s, freq='MS')

        tm.assert_index_equal(result, expected)

    def test_constructor_from_frame_series_freq(self):
        # GH 6273
        # create from a series, passing a freq
        dts = ['1-1-1990', '2-1-1990', '3-1-1990', '4-1-1990', '5-1-1990']
        expected = DatetimeIndex(dts, freq='MS')

        df = pd.DataFrame(np.random.rand(5, 3))
        df['date'] = dts
        result = DatetimeIndex(df['date'], freq='MS')

        assert df['date'].dtype == object
        expected.name = 'date'
        tm.assert_index_equal(result, expected)

        expected = pd.Series(dts, name='date')
        tm.assert_series_equal(df['date'], expected)

        # GH 6274
        # infer freq of same
        freq = pd.infer_freq(df['date'])
        assert freq == 'MS'

    @pytest.mark.parametrize("array", [
        np.arange(5), np.array(['a', 'b', 'c']), date_range(
            '2000-01-01', periods=3).values
    ])
    def test_constructor_ndarray_like(self, array):
        # GH 5460#issuecomment-44474502
        # it should be possible to convert any object that satisfies the numpy
        # ndarray interface directly into an Index
        class ArrayLike(object):
            def __init__(self, array):
                self.array = array

            def __array__(self, dtype=None):
                return self.array

        expected = pd.Index(array)
        result = pd.Index(ArrayLike(array))
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize('dtype', [
        int, 'int64', 'int32', 'int16', 'int8', 'uint64', 'uint32',
        'uint16', 'uint8'])
    def test_constructor_int_dtype_float(self, dtype):
        # GH 18400
        if is_unsigned_integer_dtype(dtype):
            index_type = UInt64Index
        else:
            index_type = Int64Index

        expected = index_type([0, 1, 2, 3])
        result = Index([0., 1., 2., 3.], dtype=dtype)
        tm.assert_index_equal(result, expected)

    def test_constructor_int_dtype_nan(self):
        # see gh-15187
        data = [np.nan]
        expected = Float64Index(data)
        result = Index(data, dtype='float')
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("dtype", ['int64', 'uint64'])
    def test_constructor_int_dtype_nan_raises(self, dtype):
        # see gh-15187
        data = [np.nan]
        msg = "cannot convert"
        with tm.assert_raises_regex(ValueError, msg):
            Index(data, dtype=dtype)

    @pytest.mark.parametrize("klass,dtype,na_val", [
        (pd.Float64Index, np.float64, np.nan),
        (pd.DatetimeIndex, 'datetime64[ns]', pd.NaT)
    ])
    def test_index_ctor_infer_nan_nat(self, klass, dtype, na_val):
        # GH 13467
        na_list = [na_val, na_val]
        expected = klass(na_list)
        assert expected.dtype == dtype

        result = Index(na_list)
        tm.assert_index_equal(result, expected)

        result = Index(np.array(na_list))
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("pos", [0, 1])
    @pytest.mark.parametrize("klass,dtype,ctor", [
        (pd.DatetimeIndex, 'datetime64[ns]', np.datetime64('nat')),
        (pd.TimedeltaIndex, 'timedelta64[ns]', np.timedelta64('nat'))
    ])
    def test_index_ctor_infer_nat_dt_like(self, pos, klass, dtype, ctor,
                                          nulls_fixture):
        expected = klass([pd.NaT, pd.NaT])
        assert expected.dtype == dtype
        data = [ctor]
        data.insert(pos, nulls_fixture)

        result = Index(data)
        tm.assert_index_equal(result, expected)

        result = Index(np.array(data, dtype=object))
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("swap_objs", [True, False])
    def test_index_ctor_nat_result(self, swap_objs):
        # mixed np.datetime64/timedelta64 nat results in object
        data = [np.datetime64('nat'), np.timedelta64('nat')]
        if swap_objs:
            data = data[::-1]

        expected = pd.Index(data, dtype=object)
        tm.assert_index_equal(Index(data), expected)
        tm.assert_index_equal(Index(np.array(data, dtype=object)), expected)

    def test_index_ctor_infer_periodindex(self):
        xp = period_range('2012-1-1', freq='M', periods=3)
        rs = Index(xp)
        tm.assert_index_equal(rs, xp)
        assert isinstance(rs, PeriodIndex)

    @pytest.mark.parametrize("vals,dtype", [
        ([1, 2, 3, 4, 5], 'int'), ([1.1, np.nan, 2.2, 3.0], 'float'),
        (['A', 'B', 'C', np.nan], 'obj')
    ])
    def test_constructor_simple_new(self, vals, dtype):
        index = Index(vals, name=dtype)
        result = index._simple_new(index, dtype)
        tm.assert_index_equal(result, index)

    @pytest.mark.parametrize("vals", [
        [1, 2, 3], np.array([1, 2, 3]), np.array([1, 2, 3], dtype=int),
        # below should coerce
        [1., 2., 3.], np.array([1., 2., 3.], dtype=float)
    ])
    def test_constructor_dtypes_to_int64(self, vals):
        index = Index(vals, dtype=int)
        assert isinstance(index, Int64Index)

    @pytest.mark.parametrize("vals", [
        [1, 2, 3], [1., 2., 3.], np.array([1., 2., 3.]),
        np.array([1, 2, 3], dtype=int), np.array([1., 2., 3.], dtype=float)
    ])
    def test_constructor_dtypes_to_float64(self, vals):
        index = Index(vals, dtype=float)
        assert isinstance(index, Float64Index)

    @pytest.mark.parametrize("cast_index", [True, False])
    @pytest.mark.parametrize("vals", [
        [True, False, True], np.array([True, False, True], dtype=bool)
    ])
    def test_constructor_dtypes_to_object(self, cast_index, vals):
        if cast_index:
            index = Index(vals, dtype=bool)
        else:
            index = Index(vals)

        assert isinstance(index, Index)
        assert index.dtype == object

    @pytest.mark.parametrize("vals", [
        [1, 2, 3], np.array([1, 2, 3], dtype=int),
        np.array([np_datetime64_compat('2011-01-01'),
                  np_datetime64_compat('2011-01-02')]),
        [datetime(2011, 1, 1), datetime(2011, 1, 2)]
    ])
    def test_constructor_dtypes_to_categorical(self, vals):
        index = Index(vals, dtype='category')
        assert isinstance(index, CategoricalIndex)

    @pytest.mark.parametrize("cast_index", [True, False])
    @pytest.mark.parametrize("vals", [
        Index(np.array([np_datetime64_compat('2011-01-01'),
                        np_datetime64_compat('2011-01-02')])),
        Index([datetime(2011, 1, 1), datetime(2011, 1, 2)])

    ])
    def test_constructor_dtypes_to_datetime(self, cast_index, vals):
        if cast_index:
            index = Index(vals, dtype=object)
            assert isinstance(index, Index)
            assert index.dtype == object
        else:
            index = Index(vals)
            assert isinstance(index, DatetimeIndex)

    @pytest.mark.parametrize("cast_index", [True, False])
    @pytest.mark.parametrize("vals", [
        np.array([np.timedelta64(1, 'D'), np.timedelta64(1, 'D')]),
        [timedelta(1), timedelta(1)]
    ])
    def test_constructor_dtypes_to_timedelta(self, cast_index, vals):
        if cast_index:
            index = Index(vals, dtype=object)
            assert isinstance(index, Index)
            assert index.dtype == object
        else:
            index = Index(vals)
            assert isinstance(index, TimedeltaIndex)

    @pytest.mark.parametrize("values", [
        # pass values without timezone, as DatetimeIndex localizes it
        pd.date_range('2011-01-01', periods=5).values,
        pd.date_range('2011-01-01', periods=5).asi8])
    @pytest.mark.parametrize("klass", [pd.Index, pd.DatetimeIndex])
    def test_constructor_dtypes_datetime(self, tz_naive_fixture, values,
                                         klass):
        index = pd.date_range('2011-01-01', periods=5, tz=tz_naive_fixture)
        dtype = index.dtype

        result = klass(values, tz=tz_naive_fixture)
        tm.assert_index_equal(result, index)

        result = klass(values, dtype=dtype)
        tm.assert_index_equal(result, index)

        result = klass(list(values), tz=tz_naive_fixture)
        tm.assert_index_equal(result, index)

        result = klass(list(values), dtype=dtype)
        tm.assert_index_equal(result, index)

    @pytest.mark.parametrize("attr", ['values', 'asi8'])
    @pytest.mark.parametrize("klass", [pd.Index, pd.TimedeltaIndex])
    def test_constructor_dtypes_timedelta(self, attr, klass):
        index = pd.timedelta_range('1 days', periods=5)
        dtype = index.dtype

        values = getattr(index, attr)

        result = klass(values, dtype=dtype)
        tm.assert_index_equal(result, index)

        result = klass(list(values), dtype=dtype)
        tm.assert_index_equal(result, index)

    def test_constructor_empty_gen(self):
        skip_index_keys = ["repeats", "periodIndex", "rangeIndex",
                           "tuples"]
        for key, index in self.generate_index_types(skip_index_keys):
            empty = index.__class__([])
            assert isinstance(empty, index.__class__)
            assert not len(empty)

    @pytest.mark.parametrize("empty,klass", [
        (PeriodIndex([], freq='B'), PeriodIndex),
        (RangeIndex(step=1), pd.RangeIndex),
        (MultiIndex(levels=[[1, 2], ['blue', 'red']],
                    labels=[[], []]), MultiIndex)
    ])
    def test_constructor_empty(self, empty, klass):
        assert isinstance(empty, klass)
        assert not len(empty)

    def test_constructor_nonhashable_name(self, indices):
        # GH 20527

        if isinstance(indices, MultiIndex):
            pytest.skip("multiindex handled in test_multi.py")

        name = ['0']
        message = "Index.name must be a hashable type"
        tm.assert_raises_regex(TypeError, message, name=name)

        # With .rename()
        renamed = [['1']]
        tm.assert_raises_regex(TypeError, message,
                               indices.rename, name=renamed)
        # With .set_names()
        tm.assert_raises_regex(TypeError, message,
                               indices.set_names, names=renamed)

    def test_view_with_args(self):

        restricted = ['unicodeIndex', 'strIndex', 'catIndex', 'boolIndex',
                      'empty']

        for i in restricted:
            ind = self.indices[i]

            # with arguments
            pytest.raises(TypeError, lambda: ind.view('i8'))

        # these are ok
        for i in list(set(self.indices.keys()) - set(restricted)):
            ind = self.indices[i]

            # with arguments
            ind.view('i8')

    def test_astype(self):
        casted = self.intIndex.astype('i8')

        # it works!
        casted.get_loc(5)

        # pass on name
        self.intIndex.name = 'foobar'
        casted = self.intIndex.astype('i8')
        assert casted.name == 'foobar'

    def test_equals_object(self):
        # same
        assert Index(['a', 'b', 'c']).equals(Index(['a', 'b', 'c']))

    @pytest.mark.parametrize("comp", [
        Index(['a', 'b']), Index(['a', 'b', 'd']), ['a', 'b', 'c']])
    def test_not_equals_object(self, comp):
        assert not Index(['a', 'b', 'c']).equals(comp)

    def test_insert(self):

        # GH 7256
        # validate neg/pos inserts
        result = Index(['b', 'c', 'd'])

        # test 0th element
        tm.assert_index_equal(Index(['a', 'b', 'c', 'd']),
                              result.insert(0, 'a'))

        # test Nth element that follows Python list behavior
        tm.assert_index_equal(Index(['b', 'c', 'e', 'd']),
                              result.insert(-1, 'e'))

        # test loc +/- neq (0, -1)
        tm.assert_index_equal(result.insert(1, 'z'), result.insert(-2, 'z'))

        # test empty
        null_index = Index([])
        tm.assert_index_equal(Index(['a']), null_index.insert(0, 'a'))

    def test_insert_missing(self, nulls_fixture):
        # GH 18295 (test missing)
        expected = Index(['a', np.nan, 'b', 'c'])
        result = Index(list('abc')).insert(1, nulls_fixture)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("pos,expected", [
        (0, Index(['b', 'c', 'd'], name='index')),
        (-1, Index(['a', 'b', 'c'], name='index'))
    ])
    def test_delete(self, pos, expected):
        index = Index(['a', 'b', 'c', 'd'], name='index')
        result = index.delete(pos)
        tm.assert_index_equal(result, expected)
        assert result.name == expected.name

    def test_delete_raises(self):
        index = Index(['a', 'b', 'c', 'd'], name='index')
        with pytest.raises((IndexError, ValueError)):
            # either depending on numpy version
            index.delete(5)

    def test_identical(self):

        # index
        i1 = Index(['a', 'b', 'c'])
        i2 = Index(['a', 'b', 'c'])

        assert i1.identical(i2)

        i1 = i1.rename('foo')
        assert i1.equals(i2)
        assert not i1.identical(i2)

        i2 = i2.rename('foo')
        assert i1.identical(i2)

        i3 = Index([('a', 'a'), ('a', 'b'), ('b', 'a')])
        i4 = Index([('a', 'a'), ('a', 'b'), ('b', 'a')], tupleize_cols=False)
        assert not i3.identical(i4)

    def test_is_(self):
        ind = Index(range(10))
        assert ind.is_(ind)
        assert ind.is_(ind.view().view().view().view())
        assert not ind.is_(Index(range(10)))
        assert not ind.is_(ind.copy())
        assert not ind.is_(ind.copy(deep=False))
        assert not ind.is_(ind[:])
        assert not ind.is_(np.array(range(10)))

        # quasi-implementation dependent
        assert ind.is_(ind.view())
        ind2 = ind.view()
        ind2.name = 'bob'
        assert ind.is_(ind2)
        assert ind2.is_(ind)
        # doesn't matter if Indices are *actually* views of underlying data,
        assert not ind.is_(Index(ind.values))
        arr = np.array(range(1, 11))
        ind1 = Index(arr, copy=False)
        ind2 = Index(arr, copy=False)
        assert not ind1.is_(ind2)

    def test_asof(self):
        d = self.dateIndex[0]
        assert self.dateIndex.asof(d) == d
        assert isna(self.dateIndex.asof(d - timedelta(1)))

        d = self.dateIndex[-1]
        assert self.dateIndex.asof(d + timedelta(1)) == d

        d = self.dateIndex[0].to_pydatetime()
        assert isinstance(self.dateIndex.asof(d), Timestamp)

    def test_asof_datetime_partial(self):
        index = pd.date_range('2010-01-01', periods=2, freq='m')
        expected = Timestamp('2010-02-28')
        result = index.asof('2010-02')
        assert result == expected
        assert not isinstance(result, Index)

    def test_nanosecond_index_access(self):
        s = Series([Timestamp('20130101')]).values.view('i8')[0]
        r = DatetimeIndex([s + 50 + i for i in range(100)])
        x = Series(np.random.randn(100), index=r)

        first_value = x.asof(x.index[0])

        # this does not yet work, as parsing strings is done via dateutil
        # assert first_value == x['2013-01-01 00:00:00.000000050+0000']

        expected_ts = np_datetime64_compat('2013-01-01 00:00:00.000000050+'
                                           '0000', 'ns')
        assert first_value == x[Timestamp(expected_ts)]

    @pytest.mark.parametrize("op", [
        operator.eq, operator.ne, operator.gt, operator.lt,
        operator.ge, operator.le
    ])
    def test_comparators(self, op):
        index = self.dateIndex
        element = index[len(index) // 2]
        element = _to_m8(element)

        arr = np.array(index)
        arr_result = op(arr, element)
        index_result = op(index, element)

        assert isinstance(index_result, np.ndarray)
        tm.assert_numpy_array_equal(arr_result, index_result)

    def test_booleanindex(self):
        boolIndex = np.repeat(True, len(self.strIndex)).astype(bool)
        boolIndex[5:30:2] = False

        subIndex = self.strIndex[boolIndex]

        for i, val in enumerate(subIndex):
            assert subIndex.get_loc(val) == i

        subIndex = self.strIndex[list(boolIndex)]
        for i, val in enumerate(subIndex):
            assert subIndex.get_loc(val) == i

    def test_fancy(self):
        sl = self.strIndex[[1, 2, 3]]
        for i in sl:
            assert i == sl[sl.get_loc(i)]

    @pytest.mark.parametrize("attr", [
        'strIndex', 'intIndex', 'floatIndex'])
    @pytest.mark.parametrize("dtype", [np.int_, np.bool_])
    def test_empty_fancy(self, attr, dtype):
        empty_arr = np.array([], dtype=dtype)
        index = getattr(self, attr)
        empty_index = index.__class__([])

        assert index[[]].identical(empty_index)
        assert index[empty_arr].identical(empty_index)

    @pytest.mark.parametrize("attr", [
        'strIndex', 'intIndex', 'floatIndex'])
    def test_empty_fancy_raises(self, attr):
        # pd.DatetimeIndex is excluded, because it overrides getitem and should
        # be tested separately.
        empty_farr = np.array([], dtype=np.float_)
        index = getattr(self, attr)
        empty_index = index.__class__([])

        assert index[[]].identical(empty_index)
        # np.ndarray only accepts ndarray of int & bool dtypes, so should Index
        pytest.raises(IndexError, index.__getitem__, empty_farr)

    @pytest.mark.parametrize("itm", [101, 'no_int'])
    def test_getitem_error(self, indices, itm):
        with pytest.raises(IndexError):
            indices[itm]

    def test_intersection(self):
        first = self.strIndex[:20]
        second = self.strIndex[:10]
        intersect = first.intersection(second)
        assert tm.equalContents(intersect, second)

        # Corner cases
        inter = first.intersection(first)
        assert inter is first

    @pytest.mark.parametrize("index2,keeps_name", [
        (Index([3, 4, 5, 6, 7], name="index"), True),  # preserve same name
        (Index([3, 4, 5, 6, 7], name="other"), False),  # drop diff names
        (Index([3, 4, 5, 6, 7]), False)])
    def test_intersection_name_preservation(self, index2, keeps_name):
        index1 = Index([1, 2, 3, 4, 5], name='index')
        expected = Index([3, 4, 5])
        result = index1.intersection(index2)

        if keeps_name:
            expected.name = 'index'

        assert result.name == expected.name
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("first_name,second_name,expected_name", [
        ('A', 'A', 'A'), ('A', 'B', None), (None, 'B', None)])
    def test_intersection_name_preservation2(self, first_name, second_name,
                                             expected_name):
        first = self.strIndex[5:20]
        second = self.strIndex[:10]
        first.name = first_name
        second.name = second_name
        intersect = first.intersection(second)
        assert intersect.name == expected_name

    @pytest.mark.parametrize("index2,keeps_name", [
        (Index([4, 7, 6, 5, 3], name='index'), True),
        (Index([4, 7, 6, 5, 3], name='other'), False)])
    def test_intersection_monotonic(self, index2, keeps_name):
        index1 = Index([5, 3, 2, 4, 1], name='index')
        expected = Index([5, 3, 4])

        if keeps_name:
            expected.name = "index"

        result = index1.intersection(index2)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("index2,expected_arr", [
        (Index(['B', 'D']), ['B']),
        (Index(['B', 'D', 'A']), ['A', 'B', 'A'])])
    def test_intersection_non_monotonic_non_unique(self, index2, expected_arr):
        # non-monotonic non-unique
        index1 = Index(['A', 'B', 'A', 'C'])
        expected = Index(expected_arr, dtype='object')
        result = index1.intersection(index2)
        tm.assert_index_equal(result, expected)

    def test_intersect_str_dates(self):
        dt_dates = [datetime(2012, 2, 9), datetime(2012, 2, 22)]

        i1 = Index(dt_dates, dtype=object)
        i2 = Index(['aa'], dtype=object)
        result = i2.intersection(i1)

        assert len(result) == 0

    def test_union(self):
        # TODO: Replace with fixturesult
        first = self.strIndex[5:20]
        second = self.strIndex[:10]
        everything = self.strIndex[:20]

        union = first.union(second)
        assert tm.equalContents(union, everything)

    @pytest.mark.parametrize("klass", [
        np.array, Series, list])
    def test_union_from_iterables(self, klass):
        # GH 10149
        # TODO: Replace with fixturesult
        first = self.strIndex[5:20]
        second = self.strIndex[:10]
        everything = self.strIndex[:20]

        case = klass(second.values)
        result = first.union(case)
        assert tm.equalContents(result, everything)

    def test_union_identity(self):
        # TODO: replace with fixturesult
        first = self.strIndex[5:20]

        union = first.union(first)
        assert union is first

        union = first.union([])
        assert union is first

        union = Index([]).union(first)
        assert union is first

    @pytest.mark.parametrize("first_list", [list('ab'), list()])
    @pytest.mark.parametrize("second_list", [list('ab'), list()])
    @pytest.mark.parametrize("first_name, second_name, expected_name", [
        ('A', 'B', None), (None, 'B', 'B'), ('A', None, 'A')])
    def test_union_name_preservation(self, first_list, second_list, first_name,
                                     second_name, expected_name):
        first = Index(first_list, name=first_name)
        second = Index(second_list, name=second_name)
        union = first.union(second)

        vals = sorted(set(first_list).union(second_list))
        expected = Index(vals, name=expected_name)
        tm.assert_index_equal(union, expected)

    def test_union_dt_as_obj(self):
        # TODO: Replace with fixturesult
        with tm.assert_produces_warning(RuntimeWarning):
            firstCat = self.strIndex.union(self.dateIndex)
        secondCat = self.strIndex.union(self.strIndex)

        if self.dateIndex.dtype == np.object_:
            appended = np.append(self.strIndex, self.dateIndex)
        else:
            appended = np.append(self.strIndex, self.dateIndex.astype('O'))

        assert tm.equalContents(firstCat, appended)
        assert tm.equalContents(secondCat, self.strIndex)
        tm.assert_contains_all(self.strIndex, firstCat)
        tm.assert_contains_all(self.strIndex, secondCat)
        tm.assert_contains_all(self.dateIndex, firstCat)

    def test_add(self):
        index = self.strIndex
        expected = Index(self.strIndex.values * 2)
        tm.assert_index_equal(index + index, expected)
        tm.assert_index_equal(index + index.tolist(), expected)
        tm.assert_index_equal(index.tolist() + index, expected)

        # test add and radd
        index = Index(list('abc'))
        expected = Index(['a1', 'b1', 'c1'])
        tm.assert_index_equal(index + '1', expected)
        expected = Index(['1a', '1b', '1c'])
        tm.assert_index_equal('1' + index, expected)

    def test_sub(self):
        index = self.strIndex
        pytest.raises(TypeError, lambda: index - 'a')
        pytest.raises(TypeError, lambda: index - index)
        pytest.raises(TypeError, lambda: index - index.tolist())
        pytest.raises(TypeError, lambda: index.tolist() - index)

    def test_map_identity_mapping(self):
        # GH 12766
        # TODO: replace with fixture
        for name, cur_index in self.indices.items():
            tm.assert_index_equal(cur_index, cur_index.map(lambda x: x))

    def test_map_with_tuples(self):
        # GH 12766

        # Test that returning a single tuple from an Index
        #   returns an Index.
        index = tm.makeIntIndex(3)
        result = tm.makeIntIndex(3).map(lambda x: (x,))
        expected = Index([(i,) for i in index])
        tm.assert_index_equal(result, expected)

        # Test that returning a tuple from a map of a single index
        #   returns a MultiIndex object.
        result = index.map(lambda x: (x, x == 1))
        expected = MultiIndex.from_tuples([(i, i == 1) for i in index])
        tm.assert_index_equal(result, expected)

    def test_map_with_tuples_mi(self):
        # Test that returning a single object from a MultiIndex
        #   returns an Index.
        first_level = ['foo', 'bar', 'baz']
        multi_index = MultiIndex.from_tuples(lzip(first_level, [1, 2, 3]))
        reduced_index = multi_index.map(lambda x: x[0])
        tm.assert_index_equal(reduced_index, Index(first_level))

    @pytest.mark.parametrize("attr", [
        'makeDateIndex', 'makePeriodIndex', 'makeTimedeltaIndex'])
    def test_map_tseries_indices_return_index(self, attr):
        index = getattr(tm, attr)(10)
        expected = Index([1] * 10)
        result = index.map(lambda x: 1)
        tm.assert_index_equal(expected, result)

    def test_map_tseries_indices_accsr_return_index(self):
        date_index = tm.makeDateIndex(24, freq='h', name='hourly')
        expected = Index(range(24), name='hourly')
        tm.assert_index_equal(expected, date_index.map(lambda x: x.hour))

    @pytest.mark.parametrize(
        "mapper",
        [
            lambda values, index: {i: e for e, i in zip(values, index)},
            lambda values, index: pd.Series(values, index)])
    def test_map_dictlike(self, mapper):
        # GH 12756
        expected = Index(['foo', 'bar', 'baz'])
        index = tm.makeIntIndex(3)
        result = index.map(mapper(expected.values, index))
        tm.assert_index_equal(result, expected)

        # TODO: replace with fixture
        for name in self.indices.keys():
            if name == 'catIndex':
                # Tested in test_categorical
                continue
            elif name == 'repeats':
                # Cannot map duplicated index
                continue

            index = self.indices[name]
            expected = Index(np.arange(len(index), 0, -1))

            # to match proper result coercion for uints
            if name == 'empty':
                expected = Index([])

            result = index.map(mapper(expected, index))
            tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("mapper", [
        Series(['foo', 2., 'baz'], index=[0, 2, -1]),
        {0: 'foo', 2: 2.0, -1: 'baz'}])
    def test_map_with_non_function_missing_values(self, mapper):
        # GH 12756
        expected = Index([2., np.nan, 'foo'])
        result = Index([2, 1, 0]).map(mapper)

        tm.assert_index_equal(expected, result)

    def test_map_na_exclusion(self):
        index = Index([1.5, np.nan, 3, np.nan, 5])

        result = index.map(lambda x: x * 2, na_action='ignore')
        expected = index * 2
        tm.assert_index_equal(result, expected)

    def test_map_defaultdict(self):
        index = Index([1, 2, 3])
        default_dict = defaultdict(lambda: 'blank')
        default_dict[1] = 'stuff'
        result = index.map(default_dict)
        expected = Index(['stuff', 'blank', 'blank'])
        tm.assert_index_equal(result, expected)

    def test_append_multiple(self):
        index = Index(['a', 'b', 'c', 'd', 'e', 'f'])

        foos = [index[:2], index[2:4], index[4:]]
        result = foos[0].append(foos[1:])
        tm.assert_index_equal(result, index)

        # empty
        result = index.append([])
        tm.assert_index_equal(result, index)

    @pytest.mark.parametrize("name,expected", [
        ('foo', 'foo'), ('bar', None)])
    def test_append_empty_preserve_name(self, name, expected):
        left = Index([], name='foo')
        right = Index([1, 2, 3], name=name)

        result = left.append(right)
        assert result.name == expected

    def test_add_string(self):
        # from bug report
        index = Index(['a', 'b', 'c'])
        index2 = index + 'foo'

        assert 'a' not in index2
        assert 'afoo' in index2

    def test_iadd_string(self):
        index = pd.Index(['a', 'b', 'c'])
        # doesn't fail test unless there is a check before `+=`
        assert 'a' in index

        index += '_x'
        assert 'a_x' in index

    @pytest.mark.parametrize("second_name,expected", [
        (None, None), ('name', 'name')])
    def test_difference_name_preservation(self, second_name, expected):
        # TODO: replace with fixturesult
        first = self.strIndex[5:20]
        second = self.strIndex[:10]
        answer = self.strIndex[10:20]

        first.name = 'name'
        second.name = second_name
        result = first.difference(second)

        assert tm.equalContents(result, answer)

        if expected is None:
            assert result.name is None
        else:
            assert result.name == expected

    def test_difference_empty_arg(self):
        first = self.strIndex[5:20]
        first.name == 'name'
        result = first.difference([])

        assert tm.equalContents(result, first)
        assert result.name == first.name

    def test_difference_identity(self):
        first = self.strIndex[5:20]
        first.name == 'name'
        result = first.difference(first)

        assert len(result) == 0
        assert result.name == first.name

    def test_symmetric_difference(self):
        # smoke
        index1 = Index([1, 2, 3, 4], name='index1')
        index2 = Index([2, 3, 4, 5])
        result = index1.symmetric_difference(index2)
        expected = Index([1, 5])
        assert tm.equalContents(result, expected)
        assert result.name is None

        # __xor__ syntax
        expected = index1 ^ index2
        assert tm.equalContents(result, expected)
        assert result.name is None

    def test_symmetric_difference_mi(self):
        index1 = MultiIndex.from_tuples(self.tuples)
        index2 = MultiIndex.from_tuples([('foo', 1), ('bar', 3)])
        result = index1.symmetric_difference(index2)
        expected = MultiIndex.from_tuples([('bar', 2), ('baz', 3), ('bar', 3)])
        assert tm.equalContents(result, expected)

    @pytest.mark.parametrize("index2,expected", [
        (Index([0, 1, np.nan]), Index([0.0, 2.0, 3.0])),
        (Index([0, 1]), Index([0.0, 2.0, 3.0, np.nan]))])
    def test_symmetric_difference_missing(self, index2, expected):
        # GH 13514 change: {nan} - {nan} == {}
        # (GH 6444, sorting of nans, is no longer an issue)
        index1 = Index([1, np.nan, 2, 3])

        result = index1.symmetric_difference(index2)
        tm.assert_index_equal(result, expected)

    def test_symmetric_difference_non_index(self):
        index1 = Index([1, 2, 3, 4], name='index1')
        index2 = np.array([2, 3, 4, 5])
        expected = Index([1, 5])
        result = index1.symmetric_difference(index2)
        assert tm.equalContents(result, expected)
        assert result.name == 'index1'

        result = index1.symmetric_difference(index2, result_name='new_name')
        assert tm.equalContents(result, expected)
        assert result.name == 'new_name'

    def test_difference_type(self):
        # GH 20040
        # If taking difference of a set and itself, it
        # needs to preserve the type of the index
        skip_index_keys = ['repeats']
        for key, index in self.generate_index_types(skip_index_keys):
            result = index.difference(index)
            expected = index.drop(index)
            tm.assert_index_equal(result, expected)

    def test_intersection_difference(self):
        # GH 20040
        # Test that the intersection of an index with an
        # empty index produces the same index as the difference
        # of an index with itself.  Test for all types
        skip_index_keys = ['repeats']
        for key, index in self.generate_index_types(skip_index_keys):
            inter = index.intersection(index.drop(index))
            diff = index.difference(index)
            tm.assert_index_equal(inter, diff)

    @pytest.mark.parametrize("attr,expected", [
        ('strIndex', False), ('boolIndex', False), ('catIndex', False),
        ('intIndex', True), ('dateIndex', False), ('floatIndex', True)])
    def test_is_numeric(self, attr, expected):
        assert getattr(self, attr).is_numeric() == expected

    @pytest.mark.parametrize("attr,expected", [
        ('strIndex', True), ('boolIndex', True), ('catIndex', False),
        ('intIndex', False), ('dateIndex', False), ('floatIndex', False)])
    def test_is_object(self, attr, expected):
        assert getattr(self, attr).is_object() == expected

    @pytest.mark.parametrize("attr,expected", [
        ('strIndex', False), ('boolIndex', False), ('catIndex', False),
        ('intIndex', False), ('dateIndex', True), ('floatIndex', False)])
    def test_is_all_dates(self, attr, expected):
        assert getattr(self, attr).is_all_dates == expected

    def test_summary(self):
        self._check_method_works(Index._summary)
        # GH3869
        ind = Index(['{other}%s', "~:{range}:0"], name='A')
        result = ind._summary()
        # shouldn't be formatted accidentally.
        assert '~:{range}:0' in result
        assert '{other}%s' in result

    # GH18217
    def test_summary_deprecated(self):
        ind = Index(['{other}%s', "~:{range}:0"], name='A')

        with tm.assert_produces_warning(FutureWarning):
            ind.summary()

    def test_format(self):
        self._check_method_works(Index.format)

        # GH 14626
        # windows has different precision on datetime.datetime.now (it doesn't
        # include us since the default for Timestamp shows these but Index
        # formatting does not we are skipping)
        now = datetime.now()
        if not str(now).endswith("000"):
            index = Index([now])
            formatted = index.format()
            expected = [str(index[0])]
            assert formatted == expected

        self.strIndex[:0].format()

    @pytest.mark.parametrize("vals", [
        [1, 2.0 + 3.0j, 4.], ['a', 'b', 'c']])
    def test_format_missing(self, vals, nulls_fixture):
        # 2845
        vals = list(vals)  # Copy for each iteration
        vals.append(nulls_fixture)
        index = Index(vals)

        formatted = index.format()
        expected = [str(index[0]), str(index[1]), str(index[2]), u('NaN')]

        assert formatted == expected
        assert index[3] is nulls_fixture

    def test_format_with_name_time_info(self):
        # bug I fixed 12/20/2011
        inc = timedelta(hours=4)
        dates = Index([dt + inc for dt in self.dateIndex], name='something')

        formatted = dates.format(name=True)
        assert formatted[0] == 'something'

    def test_format_datetime_with_time(self):
        t = Index([datetime(2012, 2, 7), datetime(2012, 2, 7, 23)])

        result = t.format()
        expected = ['2012-02-07 00:00:00', '2012-02-07 23:00:00']
        assert len(result) == 2
        assert result == expected

    @pytest.mark.parametrize("op", ['any', 'all'])
    def test_logical_compat(self, op):
        index = self.create_index()
        assert getattr(index, op)() == getattr(index.values, op)()

    def _check_method_works(self, method):
        # TODO: make this a dedicated test with parametrized methods
        method(self.empty)
        method(self.dateIndex)
        method(self.unicodeIndex)
        method(self.strIndex)
        method(self.intIndex)
        method(self.tuples)
        method(self.catIndex)

    def test_get_indexer(self):
        index1 = Index([1, 2, 3, 4, 5])
        index2 = Index([2, 4, 6])

        r1 = index1.get_indexer(index2)
        e1 = np.array([1, 3, -1], dtype=np.intp)
        assert_almost_equal(r1, e1)

    @pytest.mark.parametrize("reverse", [True, False])
    @pytest.mark.parametrize("expected,method", [
        (np.array([-1, 0, 0, 1, 1], dtype=np.intp), 'pad'),
        (np.array([-1, 0, 0, 1, 1], dtype=np.intp), 'ffill'),
        (np.array([0, 0, 1, 1, 2], dtype=np.intp), 'backfill'),
        (np.array([0, 0, 1, 1, 2], dtype=np.intp), 'bfill')])
    def test_get_indexer_methods(self, reverse, expected, method):
        index1 = Index([1, 2, 3, 4, 5])
        index2 = Index([2, 4, 6])

        if reverse:
            index1 = index1[::-1]
            expected = expected[::-1]

        result = index2.get_indexer(index1, method=method)
        assert_almost_equal(result, expected)

    def test_get_indexer_invalid(self):
        # GH10411
        index = Index(np.arange(10))

        with tm.assert_raises_regex(ValueError, 'tolerance argument'):
            index.get_indexer([1, 0], tolerance=1)

        with tm.assert_raises_regex(ValueError, 'limit argument'):
            index.get_indexer([1, 0], limit=1)

    @pytest.mark.parametrize(
        'method, tolerance, indexer, expected',
        [
            ('pad', None, [0, 5, 9], [0, 5, 9]),
            ('backfill', None, [0, 5, 9], [0, 5, 9]),
            ('nearest', None, [0, 5, 9], [0, 5, 9]),
            ('pad', 0, [0, 5, 9], [0, 5, 9]),
            ('backfill', 0, [0, 5, 9], [0, 5, 9]),
            ('nearest', 0, [0, 5, 9], [0, 5, 9]),

            ('pad', None, [0.2, 1.8, 8.5], [0, 1, 8]),
            ('backfill', None, [0.2, 1.8, 8.5], [1, 2, 9]),
            ('nearest', None, [0.2, 1.8, 8.5], [0, 2, 9]),
            ('pad', 1, [0.2, 1.8, 8.5], [0, 1, 8]),
            ('backfill', 1, [0.2, 1.8, 8.5], [1, 2, 9]),
            ('nearest', 1, [0.2, 1.8, 8.5], [0, 2, 9]),

            ('pad', 0.2, [0.2, 1.8, 8.5], [0, -1, -1]),
            ('backfill', 0.2, [0.2, 1.8, 8.5], [-1, 2, -1]),
            ('nearest', 0.2, [0.2, 1.8, 8.5], [0, 2, -1])])
    def test_get_indexer_nearest(self, method, tolerance, indexer, expected):
        index = Index(np.arange(10))

        actual = index.get_indexer(indexer, method=method, tolerance=tolerance)
        tm.assert_numpy_array_equal(actual, np.array(expected,
                                                     dtype=np.intp))

    @pytest.mark.parametrize('listtype', [list, tuple, Series, np.array])
    @pytest.mark.parametrize(
        'tolerance, expected',
        list(zip([[0.3, 0.3, 0.1], [0.2, 0.1, 0.1],
                  [0.1, 0.5, 0.5]],
                 [[0, 2, -1], [0, -1, -1],
                  [-1, 2, 9]])))
    def test_get_indexer_nearest_listlike_tolerance(self, tolerance,
                                                    expected, listtype):
        index = Index(np.arange(10))

        actual = index.get_indexer([0.2, 1.8, 8.5], method='nearest',
                                   tolerance=listtype(tolerance))
        tm.assert_numpy_array_equal(actual, np.array(expected,
                                                     dtype=np.intp))

    def test_get_indexer_nearest_error(self):
        index = Index(np.arange(10))
        with tm.assert_raises_regex(ValueError, 'limit argument'):
            index.get_indexer([1, 0], method='nearest', limit=1)

        with pytest.raises(ValueError, match='tolerance size must match'):
            index.get_indexer([1, 0], method='nearest',
                              tolerance=[1, 2, 3])

    @pytest.mark.parametrize("method,expected", [
        ('pad', [8, 7, 0]), ('backfill', [9, 8, 1]), ('nearest', [9, 7, 0])])
    def test_get_indexer_nearest_decreasing(self, method, expected):
        index = Index(np.arange(10))[::-1]

        actual = index.get_indexer([0, 5, 9], method=method)
        tm.assert_numpy_array_equal(actual, np.array([9, 4, 0], dtype=np.intp))

        actual = index.get_indexer([0.2, 1.8, 8.5], method=method)
        tm.assert_numpy_array_equal(actual, np.array(expected, dtype=np.intp))

    @pytest.mark.parametrize("method,expected", [
        ('pad', np.array([-1, 0, 1, 1], dtype=np.intp)),
        ('backfill', np.array([0, 0, 1, -1], dtype=np.intp))])
    def test_get_indexer_strings(self, method, expected):
        index = pd.Index(['b', 'c'])
        actual = index.get_indexer(['a', 'b', 'c', 'd'], method=method)

        tm.assert_numpy_array_equal(actual, expected)

    def test_get_indexer_strings_raises(self):
        index = pd.Index(['b', 'c'])

        with pytest.raises(TypeError):
            index.get_indexer(['a', 'b', 'c', 'd'], method='nearest')

        with pytest.raises(TypeError):
            index.get_indexer(['a', 'b', 'c', 'd'], method='pad', tolerance=2)

        with pytest.raises(TypeError):
            index.get_indexer(['a', 'b', 'c', 'd'], method='pad',
                              tolerance=[2, 2, 2, 2])

    def test_get_indexer_numeric_index_boolean_target(self):
        # GH 16877
        numeric_idx = pd.Index(range(4))
        result = numeric_idx.get_indexer([True, False, True])
        expected = np.array([-1, -1, -1], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    def test_get_loc(self):
        idx = pd.Index([0, 1, 2])
        all_methods = [None, 'pad', 'backfill', 'nearest']
        for method in all_methods:
            assert idx.get_loc(1, method=method) == 1
            if method is not None:
                assert idx.get_loc(1, method=method, tolerance=0) == 1
            with pytest.raises(TypeError):
                idx.get_loc([1, 2], method=method)

        for method, loc in [('pad', 1), ('backfill', 2), ('nearest', 1)]:
            assert idx.get_loc(1.1, method) == loc

        for method, loc in [('pad', 1), ('backfill', 2), ('nearest', 1)]:
            assert idx.get_loc(1.1, method, tolerance=1) == loc

        for method in ['pad', 'backfill', 'nearest']:
            with pytest.raises(KeyError):
                idx.get_loc(1.1, method, tolerance=0.05)

        with tm.assert_raises_regex(ValueError, 'must be numeric'):
            idx.get_loc(1.1, 'nearest', tolerance='invalid')
        with tm.assert_raises_regex(ValueError, 'tolerance .* valid if'):
            idx.get_loc(1.1, tolerance=1)
        with pytest.raises(ValueError, match='tolerance size must match'):
            idx.get_loc(1.1, 'nearest', tolerance=[1, 1])

        idx = pd.Index(['a', 'c'])
        with pytest.raises(TypeError):
            idx.get_loc('a', method='nearest')
        with pytest.raises(TypeError):
            idx.get_loc('a', method='pad', tolerance='invalid')

    def test_slice_locs(self):
        for dtype in [int, float]:
            idx = Index(np.array([0, 1, 2, 5, 6, 7, 9, 10], dtype=dtype))
            n = len(idx)

            assert idx.slice_locs(start=2) == (2, n)
            assert idx.slice_locs(start=3) == (3, n)
            assert idx.slice_locs(3, 8) == (3, 6)
            assert idx.slice_locs(5, 10) == (3, n)
            assert idx.slice_locs(end=8) == (0, 6)
            assert idx.slice_locs(end=9) == (0, 7)

            # reversed
            idx2 = idx[::-1]
            assert idx2.slice_locs(8, 2) == (2, 6)
            assert idx2.slice_locs(7, 3) == (2, 5)

        # float slicing
        idx = Index(np.array([0, 1, 2, 5, 6, 7, 9, 10], dtype=float))
        n = len(idx)
        assert idx.slice_locs(5.0, 10.0) == (3, n)
        assert idx.slice_locs(4.5, 10.5) == (3, 8)
        idx2 = idx[::-1]
        assert idx2.slice_locs(8.5, 1.5) == (2, 6)
        assert idx2.slice_locs(10.5, -1) == (0, n)

        # int slicing with floats
        # GH 4892, these are all TypeErrors
        idx = Index(np.array([0, 1, 2, 5, 6, 7, 9, 10], dtype=int))
        pytest.raises(TypeError,
                      lambda: idx.slice_locs(5.0, 10.0), (3, n))
        pytest.raises(TypeError,
                      lambda: idx.slice_locs(4.5, 10.5), (3, 8))
        idx2 = idx[::-1]
        pytest.raises(TypeError,
                      lambda: idx2.slice_locs(8.5, 1.5), (2, 6))
        pytest.raises(TypeError,
                      lambda: idx2.slice_locs(10.5, -1), (0, n))

    def test_slice_locs_dup(self):
        idx = Index(['a', 'a', 'b', 'c', 'd', 'd'])
        assert idx.slice_locs('a', 'd') == (0, 6)
        assert idx.slice_locs(end='d') == (0, 6)
        assert idx.slice_locs('a', 'c') == (0, 4)
        assert idx.slice_locs('b', 'd') == (2, 6)

        idx2 = idx[::-1]
        assert idx2.slice_locs('d', 'a') == (0, 6)
        assert idx2.slice_locs(end='a') == (0, 6)
        assert idx2.slice_locs('d', 'b') == (0, 4)
        assert idx2.slice_locs('c', 'a') == (2, 6)

        for dtype in [int, float]:
            idx = Index(np.array([10, 12, 12, 14], dtype=dtype))
            assert idx.slice_locs(12, 12) == (1, 3)
            assert idx.slice_locs(11, 13) == (1, 3)

            idx2 = idx[::-1]
            assert idx2.slice_locs(12, 12) == (1, 3)
            assert idx2.slice_locs(13, 11) == (1, 3)

    def test_slice_locs_na(self):
        idx = Index([np.nan, 1, 2])
        pytest.raises(KeyError, idx.slice_locs, start=1.5)
        pytest.raises(KeyError, idx.slice_locs, end=1.5)
        assert idx.slice_locs(1) == (1, 3)
        assert idx.slice_locs(np.nan) == (0, 3)

        idx = Index([0, np.nan, np.nan, 1, 2])
        assert idx.slice_locs(np.nan) == (1, 5)

    def test_slice_locs_negative_step(self):
        idx = Index(list('bcdxy'))

        SLC = pd.IndexSlice

        def check_slice(in_slice, expected):
            s_start, s_stop = idx.slice_locs(in_slice.start, in_slice.stop,
                                             in_slice.step)
            result = idx[s_start:s_stop:in_slice.step]
            expected = pd.Index(list(expected))
            tm.assert_index_equal(result, expected)

        for in_slice, expected in [
            (SLC[::-1], 'yxdcb'), (SLC['b':'y':-1], ''),
            (SLC['b'::-1], 'b'), (SLC[:'b':-1], 'yxdcb'),
            (SLC[:'y':-1], 'y'), (SLC['y'::-1], 'yxdcb'),
            (SLC['y'::-4], 'yb'),
            # absent labels
            (SLC[:'a':-1], 'yxdcb'), (SLC[:'a':-2], 'ydb'),
            (SLC['z'::-1], 'yxdcb'), (SLC['z'::-3], 'yc'),
            (SLC['m'::-1], 'dcb'), (SLC[:'m':-1], 'yx'),
            (SLC['a':'a':-1], ''), (SLC['z':'z':-1], ''),
            (SLC['m':'m':-1], '')
        ]:
            check_slice(in_slice, expected)

    def test_drop(self):
        n = len(self.strIndex)

        drop = self.strIndex[lrange(5, 10)]
        dropped = self.strIndex.drop(drop)
        expected = self.strIndex[lrange(5) + lrange(10, n)]
        tm.assert_index_equal(dropped, expected)

        pytest.raises(KeyError, self.strIndex.drop, ['foo', 'bar'])
        pytest.raises(KeyError, self.strIndex.drop, ['1', 'bar'])

        # errors='ignore'
        mixed = drop.tolist() + ['foo']
        dropped = self.strIndex.drop(mixed, errors='ignore')
        expected = self.strIndex[lrange(5) + lrange(10, n)]
        tm.assert_index_equal(dropped, expected)

        dropped = self.strIndex.drop(['foo', 'bar'], errors='ignore')
        expected = self.strIndex[lrange(n)]
        tm.assert_index_equal(dropped, expected)

        dropped = self.strIndex.drop(self.strIndex[0])
        expected = self.strIndex[1:]
        tm.assert_index_equal(dropped, expected)

        ser = Index([1, 2, 3])
        dropped = ser.drop(1)
        expected = Index([2, 3])
        tm.assert_index_equal(dropped, expected)

        # errors='ignore'
        pytest.raises(KeyError, ser.drop, [3, 4])

        dropped = ser.drop(4, errors='ignore')
        expected = Index([1, 2, 3])
        tm.assert_index_equal(dropped, expected)

        dropped = ser.drop([3, 4, 5], errors='ignore')
        expected = Index([1, 2])
        tm.assert_index_equal(dropped, expected)

    @pytest.mark.parametrize("values", [['a', 'b', ('c', 'd')],
                                        ['a', ('c', 'd'), 'b'],
                                        [('c', 'd'), 'a', 'b']])
    @pytest.mark.parametrize("to_drop", [[('c', 'd'), 'a'], ['a', ('c', 'd')]])
    def test_drop_tuple(self, values, to_drop):
        # GH 18304
        index = pd.Index(values)
        expected = pd.Index(['b'])

        result = index.drop(to_drop)
        tm.assert_index_equal(result, expected)

        removed = index.drop(to_drop[0])
        for drop_me in to_drop[1], [to_drop[1]]:
            result = removed.drop(drop_me)
            tm.assert_index_equal(result, expected)

        removed = index.drop(to_drop[1])
        for drop_me in to_drop[1], [to_drop[1]]:
            pytest.raises(KeyError, removed.drop, drop_me)

    def test_tuple_union_bug(self):
        import pandas
        import numpy as np

        aidx1 = np.array([(1, 'A'), (2, 'A'), (1, 'B'), (2, 'B')],
                         dtype=[('num', int), ('let', 'a1')])
        aidx2 = np.array([(1, 'A'), (2, 'A'), (1, 'B'),
                          (2, 'B'), (1, 'C'), (2, 'C')],
                         dtype=[('num', int), ('let', 'a1')])

        idx1 = pandas.Index(aidx1)
        idx2 = pandas.Index(aidx2)

        # intersection broken?
        int_idx = idx1.intersection(idx2)
        # needs to be 1d like idx1 and idx2
        expected = idx1[:4]  # pandas.Index(sorted(set(idx1) & set(idx2)))
        assert int_idx.ndim == 1
        tm.assert_index_equal(int_idx, expected)

        # union broken
        union_idx = idx1.union(idx2)
        expected = idx2
        assert union_idx.ndim == 1
        tm.assert_index_equal(union_idx, expected)

    def test_is_monotonic_incomparable(self):
        index = Index([5, datetime.now(), 7])
        assert not index.is_monotonic_increasing
        assert not index.is_monotonic_decreasing
        assert not index._is_strictly_monotonic_increasing
        assert not index._is_strictly_monotonic_decreasing

    def test_get_set_value(self):
        values = np.random.randn(100)
        date = self.dateIndex[67]

        assert_almost_equal(self.dateIndex.get_value(values, date), values[67])

        self.dateIndex.set_value(values, date, 10)
        assert values[67] == 10

    def test_isin(self):
        values = ['foo', 'bar', 'quux']

        idx = Index(['qux', 'baz', 'foo', 'bar'])
        result = idx.isin(values)
        expected = np.array([False, False, True, True])
        tm.assert_numpy_array_equal(result, expected)

        # set
        result = idx.isin(set(values))
        tm.assert_numpy_array_equal(result, expected)

        # empty, return dtype bool
        idx = Index([])
        result = idx.isin(values)
        assert len(result) == 0
        assert result.dtype == np.bool_

    @pytest.mark.skipif(PYPY, reason="np.nan is float('nan') on PyPy")
    def test_isin_nan_not_pypy(self):
        tm.assert_numpy_array_equal(Index(['a', np.nan]).isin([float('nan')]),
                                    np.array([False, False]))

    @pytest.mark.skipif(not PYPY, reason="np.nan is float('nan') on PyPy")
    def test_isin_nan_pypy(self):
        tm.assert_numpy_array_equal(Index(['a', np.nan]).isin([float('nan')]),
                                    np.array([False, True]))

    def test_isin_nan_common(self):
        tm.assert_numpy_array_equal(Index(['a', np.nan]).isin([np.nan]),
                                    np.array([False, True]))
        tm.assert_numpy_array_equal(Index(['a', pd.NaT]).isin([pd.NaT]),
                                    np.array([False, True]))
        tm.assert_numpy_array_equal(Index(['a', np.nan]).isin([pd.NaT]),
                                    np.array([False, False]))

        # Float64Index overrides isin, so must be checked separately
        tm.assert_numpy_array_equal(Float64Index([1.0, np.nan]).isin([np.nan]),
                                    np.array([False, True]))
        tm.assert_numpy_array_equal(
            Float64Index([1.0, np.nan]).isin([float('nan')]),
            np.array([False, True]))

        # we cannot compare NaT with NaN
        tm.assert_numpy_array_equal(Float64Index([1.0, np.nan]).isin([pd.NaT]),
                                    np.array([False, False]))

    def test_isin_level_kwarg(self):
        def check_idx(idx):
            values = idx.tolist()[-2:] + ['nonexisting']

            expected = np.array([False, False, True, True])
            tm.assert_numpy_array_equal(expected, idx.isin(values, level=0))
            tm.assert_numpy_array_equal(expected, idx.isin(values, level=-1))

            pytest.raises(IndexError, idx.isin, values, level=1)
            pytest.raises(IndexError, idx.isin, values, level=10)
            pytest.raises(IndexError, idx.isin, values, level=-2)

            pytest.raises(KeyError, idx.isin, values, level=1.0)
            pytest.raises(KeyError, idx.isin, values, level='foobar')

            idx.name = 'foobar'
            tm.assert_numpy_array_equal(expected,
                                        idx.isin(values, level='foobar'))

            pytest.raises(KeyError, idx.isin, values, level='xyzzy')
            pytest.raises(KeyError, idx.isin, values, level=np.nan)

        check_idx(Index(['qux', 'baz', 'foo', 'bar']))
        # Float64Index overrides isin, so must be checked separately
        check_idx(Float64Index([1.0, 2.0, 3.0, 4.0]))

    @pytest.mark.parametrize("empty", [[], Series(), np.array([])])
    def test_isin_empty(self, empty):
        # see gh-16991
        idx = Index(["a", "b"])
        expected = np.array([False, False])

        result = idx.isin(empty)
        tm.assert_numpy_array_equal(expected, result)

    def test_boolean_cmp(self):
        values = [1, 2, 3, 4]

        idx = Index(values)
        res = (idx == values)

        tm.assert_numpy_array_equal(res, np.array(
            [True, True, True, True], dtype=bool))

    def test_get_level_values(self):
        result = self.strIndex.get_level_values(0)
        tm.assert_index_equal(result, self.strIndex)

        # test for name (GH 17414)
        index_with_name = self.strIndex.copy()
        index_with_name.name = 'a'
        result = index_with_name.get_level_values('a')
        tm.assert_index_equal(result, index_with_name)

    def test_slice_keep_name(self):
        idx = Index(['a', 'b'], name='asdf')
        assert idx.name == idx[1:].name

    # instance attributes of the form self.<name>Index
    @pytest.mark.parametrize('index_kind',
                             ['unicode', 'str', 'date', 'int', 'float'])
    def test_join_self(self, join_type, index_kind):

        res = getattr(self, '{0}Index'.format(index_kind))

        joined = res.join(res, how=join_type)
        assert res is joined

    def test_str_attribute(self):
        # GH9068
        methods = ['strip', 'rstrip', 'lstrip']
        idx = Index([' jack', 'jill ', ' jesse ', 'frank'])
        for method in methods:
            expected = Index([getattr(str, method)(x) for x in idx.values])
            tm.assert_index_equal(
                getattr(Index.str, method)(idx.str), expected)

        # create a few instances that are not able to use .str accessor
        indices = [Index(range(5)), tm.makeDateIndex(10),
                   MultiIndex.from_tuples([('foo', '1'), ('bar', '3')]),
                   PeriodIndex(start='2000', end='2010', freq='A')]
        for idx in indices:
            with tm.assert_raises_regex(AttributeError,
                                        'only use .str accessor'):
                idx.str.repeat(2)

        idx = Index(['a b c', 'd e', 'f'])
        expected = Index([['a', 'b', 'c'], ['d', 'e'], ['f']])
        tm.assert_index_equal(idx.str.split(), expected)
        tm.assert_index_equal(idx.str.split(expand=False), expected)

        expected = MultiIndex.from_tuples([('a', 'b', 'c'), ('d', 'e', np.nan),
                                           ('f', np.nan, np.nan)])
        tm.assert_index_equal(idx.str.split(expand=True), expected)

        # test boolean case, should return np.array instead of boolean Index
        idx = Index(['a1', 'a2', 'b1', 'b2'])
        expected = np.array([True, True, False, False])
        tm.assert_numpy_array_equal(idx.str.startswith('a'), expected)
        assert isinstance(idx.str.startswith('a'), np.ndarray)
        s = Series(range(4), index=idx)
        expected = Series(range(2), index=['a1', 'a2'])
        tm.assert_series_equal(s[s.index.str.startswith('a')], expected)

    def test_tab_completion(self):
        # GH 9910
        idx = Index(list('abcd'))
        assert 'str' in dir(idx)

        idx = Index(range(4))
        assert 'str' not in dir(idx)

    def test_indexing_doesnt_change_class(self):
        idx = Index([1, 2, 3, 'a', 'b', 'c'])

        assert idx[1:3].identical(pd.Index([2, 3], dtype=np.object_))
        assert idx[[0, 1]].identical(pd.Index([1, 2], dtype=np.object_))

    def test_outer_join_sort(self):
        left_idx = Index(np.random.permutation(15))
        right_idx = tm.makeDateIndex(10)

        with tm.assert_produces_warning(RuntimeWarning):
            joined = left_idx.join(right_idx, how='outer')

        # right_idx in this case because DatetimeIndex has join precedence over
        # Int64Index
        with tm.assert_produces_warning(RuntimeWarning):
            expected = right_idx.astype(object).union(left_idx.astype(object))
        tm.assert_index_equal(joined, expected)

    def test_nan_first_take_datetime(self):
        idx = Index([pd.NaT, Timestamp('20130101'), Timestamp('20130102')])
        res = idx.take([-1, 0, 1])
        exp = Index([idx[-1], idx[0], idx[1]])
        tm.assert_index_equal(res, exp)

    def test_take_fill_value(self):
        # GH 12631
        idx = pd.Index(list('ABC'), name='xxx')
        result = idx.take(np.array([1, 0, -1]))
        expected = pd.Index(list('BAC'), name='xxx')
        tm.assert_index_equal(result, expected)

        # fill_value
        result = idx.take(np.array([1, 0, -1]), fill_value=True)
        expected = pd.Index(['B', 'A', np.nan], name='xxx')
        tm.assert_index_equal(result, expected)

        # allow_fill=False
        result = idx.take(np.array([1, 0, -1]), allow_fill=False,
                          fill_value=True)
        expected = pd.Index(['B', 'A', 'C'], name='xxx')
        tm.assert_index_equal(result, expected)

        msg = ('When allow_fill=True and fill_value is not None, '
               'all indices must be >= -1')
        with tm.assert_raises_regex(ValueError, msg):
            idx.take(np.array([1, 0, -2]), fill_value=True)
        with tm.assert_raises_regex(ValueError, msg):
            idx.take(np.array([1, 0, -5]), fill_value=True)

        with pytest.raises(IndexError):
            idx.take(np.array([1, -5]))

    def test_reindex_preserves_name_if_target_is_list_or_ndarray(self):
        # GH6552
        idx = pd.Index([0, 1, 2])

        dt_idx = pd.date_range('20130101', periods=3)

        idx.name = None
        assert idx.reindex([])[0].name is None
        assert idx.reindex(np.array([]))[0].name is None
        assert idx.reindex(idx.tolist())[0].name is None
        assert idx.reindex(idx.tolist()[:-1])[0].name is None
        assert idx.reindex(idx.values)[0].name is None
        assert idx.reindex(idx.values[:-1])[0].name is None

        # Must preserve name even if dtype changes.
        assert idx.reindex(dt_idx.values)[0].name is None
        assert idx.reindex(dt_idx.tolist())[0].name is None

        idx.name = 'foobar'
        assert idx.reindex([])[0].name == 'foobar'
        assert idx.reindex(np.array([]))[0].name == 'foobar'
        assert idx.reindex(idx.tolist())[0].name == 'foobar'
        assert idx.reindex(idx.tolist()[:-1])[0].name == 'foobar'
        assert idx.reindex(idx.values)[0].name == 'foobar'
        assert idx.reindex(idx.values[:-1])[0].name == 'foobar'

        # Must preserve name even if dtype changes.
        assert idx.reindex(dt_idx.values)[0].name == 'foobar'
        assert idx.reindex(dt_idx.tolist())[0].name == 'foobar'

    def test_reindex_preserves_type_if_target_is_empty_list_or_array(self):
        # GH7774
        idx = pd.Index(list('abc'))

        def get_reindex_type(target):
            return idx.reindex(target)[0].dtype.type

        assert get_reindex_type([]) == np.object_
        assert get_reindex_type(np.array([])) == np.object_
        assert get_reindex_type(np.array([], dtype=np.int64)) == np.object_

    def test_reindex_doesnt_preserve_type_if_target_is_empty_index(self):
        # GH7774
        idx = pd.Index(list('abc'))

        def get_reindex_type(target):
            return idx.reindex(target)[0].dtype.type

        assert get_reindex_type(pd.Int64Index([])) == np.int64
        assert get_reindex_type(pd.Float64Index([])) == np.float64
        assert get_reindex_type(pd.DatetimeIndex([])) == np.datetime64

        reindexed = idx.reindex(pd.MultiIndex(
            [pd.Int64Index([]), pd.Float64Index([])], [[], []]))[0]
        assert reindexed.levels[0].dtype.type == np.int64
        assert reindexed.levels[1].dtype.type == np.float64

    def test_groupby(self):
        idx = Index(range(5))
        groups = idx.groupby(np.array([1, 1, 2, 2, 2]))
        exp = {1: pd.Index([0, 1]), 2: pd.Index([2, 3, 4])}
        tm.assert_dict_equal(groups, exp)

    def test_equals_op_multiindex(self):
        # GH9785
        # test comparisons of multiindex
        from pandas.compat import StringIO
        df = pd.read_csv(StringIO('a,b,c\n1,2,3\n4,5,6'), index_col=[0, 1])
        tm.assert_numpy_array_equal(df.index == df.index,
                                    np.array([True, True]))

        mi1 = MultiIndex.from_tuples([(1, 2), (4, 5)])
        tm.assert_numpy_array_equal(df.index == mi1, np.array([True, True]))
        mi2 = MultiIndex.from_tuples([(1, 2), (4, 6)])
        tm.assert_numpy_array_equal(df.index == mi2, np.array([True, False]))
        mi3 = MultiIndex.from_tuples([(1, 2), (4, 5), (8, 9)])
        with tm.assert_raises_regex(ValueError, "Lengths must match"):
            df.index == mi3

        index_a = Index(['foo', 'bar', 'baz'])
        with tm.assert_raises_regex(ValueError, "Lengths must match"):
            df.index == index_a
        tm.assert_numpy_array_equal(index_a == mi3,
                                    np.array([False, False, False]))

    def test_conversion_preserves_name(self):
        # GH 10875
        i = pd.Index(['01:02:03', '01:02:04'], name='label')
        assert i.name == pd.to_datetime(i).name
        assert i.name == pd.to_timedelta(i).name

    def test_string_index_repr(self):
        # py3/py2 repr can differ because of "u" prefix
        # which also affects to displayed element size

        if PY3:
            coerce = lambda x: x
        else:
            coerce = unicode  # noqa

        # short
        idx = pd.Index(['a', 'bb', 'ccc'])
        if PY3:
            expected = u"""Index(['a', 'bb', 'ccc'], dtype='object')"""
            assert repr(idx) == expected
        else:
            expected = u"""Index([u'a', u'bb', u'ccc'], dtype='object')"""
            assert coerce(idx) == expected

        # multiple lines
        idx = pd.Index(['a', 'bb', 'ccc'] * 10)
        if PY3:
            expected = u"""\
Index(['a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc',
       'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc',
       'a', 'bb', 'ccc', 'a', 'bb', 'ccc'],
      dtype='object')"""

            assert repr(idx) == expected
        else:
            expected = u"""\
Index([u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a',
       u'bb', u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a', u'bb',
       u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc'],
      dtype='object')"""

            assert coerce(idx) == expected

        # truncated
        idx = pd.Index(['a', 'bb', 'ccc'] * 100)
        if PY3:
            expected = u"""\
Index(['a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a',
       ...
       'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc'],
      dtype='object', length=300)"""

            assert repr(idx) == expected
        else:
            expected = u"""\
Index([u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a',
       ...
       u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc'],
      dtype='object', length=300)"""

            assert coerce(idx) == expected

        # short
        idx = pd.Index([u'あ', u'いい', u'ううう'])
        if PY3:
            expected = u"""Index(['あ', 'いい', 'ううう'], dtype='object')"""
            assert repr(idx) == expected
        else:
            expected = u"""Index([u'あ', u'いい', u'ううう'], dtype='object')"""
            assert coerce(idx) == expected

        # multiple lines
        idx = pd.Index([u'あ', u'いい', u'ううう'] * 10)
        if PY3:
            expected = (u"Index(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', "
                        u"'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう',\n"
                        u"       'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', "
                        u"'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう',\n"
                        u"       'あ', 'いい', 'ううう', 'あ', 'いい', "
                        u"'ううう'],\n"
                        u"      dtype='object')")
            assert repr(idx) == expected
        else:
            expected = (u"Index([u'あ', u'いい', u'ううう', u'あ', u'いい', "
                        u"u'ううう', u'あ', u'いい', u'ううう', u'あ',\n"
                        u"       u'いい', u'ううう', u'あ', u'いい', u'ううう', "
                        u"u'あ', u'いい', u'ううう', u'あ', u'いい',\n"
                        u"       u'ううう', u'あ', u'いい', u'ううう', u'あ', "
                        u"u'いい', u'ううう', u'あ', u'いい', u'ううう'],\n"
                        u"      dtype='object')")
            assert coerce(idx) == expected

        # truncated
        idx = pd.Index([u'あ', u'いい', u'ううう'] * 100)
        if PY3:
            expected = (u"Index(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', "
                        u"'あ', 'いい', 'ううう', 'あ',\n"
                        u"       ...\n"
                        u"       'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', "
                        u"'ううう', 'あ', 'いい', 'ううう'],\n"
                        u"      dtype='object', length=300)")
            assert repr(idx) == expected
        else:
            expected = (u"Index([u'あ', u'いい', u'ううう', u'あ', u'いい', "
                        u"u'ううう', u'あ', u'いい', u'ううう', u'あ',\n"
                        u"       ...\n"
                        u"       u'ううう', u'あ', u'いい', u'ううう', u'あ', "
                        u"u'いい', u'ううう', u'あ', u'いい', u'ううう'],\n"
                        u"      dtype='object', length=300)")

            assert coerce(idx) == expected

        # Emable Unicode option -----------------------------------------
        with cf.option_context('display.unicode.east_asian_width', True):

            # short
            idx = pd.Index([u'あ', u'いい', u'ううう'])
            if PY3:
                expected = (u"Index(['あ', 'いい', 'ううう'], "
                            u"dtype='object')")
                assert repr(idx) == expected
            else:
                expected = (u"Index([u'あ', u'いい', u'ううう'], "
                            u"dtype='object')")
                assert coerce(idx) == expected

            # multiple lines
            idx = pd.Index([u'あ', u'いい', u'ううう'] * 10)
            if PY3:
                expected = (u"Index(['あ', 'いい', 'ううう', 'あ', 'いい', "
                            u"'ううう', 'あ', 'いい', 'ううう',\n"
                            u"       'あ', 'いい', 'ううう', 'あ', 'いい', "
                            u"'ううう', 'あ', 'いい', 'ううう',\n"
                            u"       'あ', 'いい', 'ううう', 'あ', 'いい', "
                            u"'ううう', 'あ', 'いい', 'ううう',\n"
                            u"       'あ', 'いい', 'ううう'],\n"
                            u"      dtype='object')""")

                assert repr(idx) == expected
            else:
                expected = (u"Index([u'あ', u'いい', u'ううう', u'あ', u'いい', "
                            u"u'ううう', u'あ', u'いい',\n"
                            u"       u'ううう', u'あ', u'いい', u'ううう', "
                            u"u'あ', u'いい', u'ううう', u'あ',\n"
                            u"       u'いい', u'ううう', u'あ', u'いい', "
                            u"u'ううう', u'あ', u'いい',\n"
                            u"       u'ううう', u'あ', u'いい', u'ううう', "
                            u"u'あ', u'いい', u'ううう'],\n"
                            u"      dtype='object')")

                assert coerce(idx) == expected

            # truncated
            idx = pd.Index([u'あ', u'いい', u'ううう'] * 100)
            if PY3:
                expected = (u"Index(['あ', 'いい', 'ううう', 'あ', 'いい', "
                            u"'ううう', 'あ', 'いい', 'ううう',\n"
                            u"       'あ',\n"
                            u"       ...\n"
                            u"       'ううう', 'あ', 'いい', 'ううう', 'あ', "
                            u"'いい', 'ううう', 'あ', 'いい',\n"
                            u"       'ううう'],\n"
                            u"      dtype='object', length=300)")

                assert repr(idx) == expected
            else:
                expected = (u"Index([u'あ', u'いい', u'ううう', u'あ', u'いい', "
                            u"u'ううう', u'あ', u'いい',\n"
                            u"       u'ううう', u'あ',\n"
                            u"       ...\n"
                            u"       u'ううう', u'あ', u'いい', u'ううう', "
                            u"u'あ', u'いい', u'ううう', u'あ',\n"
                            u"       u'いい', u'ううう'],\n"
                            u"      dtype='object', length=300)")

                assert coerce(idx) == expected

    @pytest.mark.parametrize('dtype', [np.int64, np.float64])
    @pytest.mark.parametrize('delta', [1, 0, -1])
    def test_addsub_arithmetic(self, dtype, delta):
        # GH 8142
        delta = dtype(delta)
        idx = pd.Index([10, 11, 12], dtype=dtype)
        result = idx + delta
        expected = pd.Index(idx.values + delta, dtype=dtype)
        tm.assert_index_equal(result, expected)

        # this subtraction used to fail
        result = idx - delta
        expected = pd.Index(idx.values - delta, dtype=dtype)
        tm.assert_index_equal(result, expected)

        tm.assert_index_equal(idx + idx, 2 * idx)
        tm.assert_index_equal(idx - idx, 0 * idx)
        assert not (idx - idx).empty

    def test_iadd_preserves_name(self):
        # GH#17067, GH#19723 __iadd__ and __isub__ should preserve index name
        ser = pd.Series([1, 2, 3])
        ser.index.name = 'foo'

        ser.index += 1
        assert ser.index.name == "foo"

        ser.index -= 1
        assert ser.index.name == "foo"

    def test_cached_properties_not_settable(self):
        idx = pd.Index([1, 2, 3])
        with tm.assert_raises_regex(AttributeError, "Can't set attribute"):
            idx.is_unique = False

    def test_get_duplicates_deprecated(self):
        idx = pd.Index([1, 2, 3])
        with tm.assert_produces_warning(FutureWarning):
            idx.get_duplicates()


class TestMixedIntIndex(Base):
    # Mostly the tests from common.py for which the results differ
    # in py2 and py3 because ints and strings are uncomparable in py3
    # (GH 13514)

    _holder = Index

    def setup_method(self, method):
        self.indices = dict(mixedIndex=Index([0, 'a', 1, 'b', 2, 'c']))
        self.setup_indices()

    def create_index(self):
        return self.mixedIndex

    def test_argsort(self):
        idx = self.create_index()
        if PY36:
            with tm.assert_raises_regex(TypeError, "'>|<' not supported"):
                result = idx.argsort()
        elif PY3:
            with tm.assert_raises_regex(TypeError, "unorderable types"):
                result = idx.argsort()
        else:
            result = idx.argsort()
            expected = np.array(idx).argsort()
            tm.assert_numpy_array_equal(result, expected, check_dtype=False)

    def test_numpy_argsort(self):
        idx = self.create_index()
        if PY36:
            with tm.assert_raises_regex(TypeError, "'>|<' not supported"):
                result = np.argsort(idx)
        elif PY3:
            with tm.assert_raises_regex(TypeError, "unorderable types"):
                result = np.argsort(idx)
        else:
            result = np.argsort(idx)
            expected = idx.argsort()
            tm.assert_numpy_array_equal(result, expected)

    def test_copy_name(self):
        # Check that "name" argument passed at initialization is honoured
        # GH12309
        idx = self.create_index()

        first = idx.__class__(idx, copy=True, name='mario')
        second = first.__class__(first, copy=False)

        # Even though "copy=False", we want a new object.
        assert first is not second
        # Not using tm.assert_index_equal() since names differ:
        assert idx.equals(first)

        assert first.name == 'mario'
        assert second.name == 'mario'

        s1 = Series(2, index=first)
        s2 = Series(3, index=second[:-1])

        warning_type = RuntimeWarning if PY3 else None
        with tm.assert_produces_warning(warning_type):
            # Python 3: Unorderable types
            s3 = s1 * s2

        assert s3.index.name == 'mario'

    def test_copy_name2(self):
        # Check that adding a "name" parameter to the copy is honored
        # GH14302
        idx = pd.Index([1, 2], name='MyName')
        idx1 = idx.copy()

        assert idx.equals(idx1)
        assert idx.name == 'MyName'
        assert idx1.name == 'MyName'

        idx2 = idx.copy(name='NewName')

        assert idx.equals(idx2)
        assert idx.name == 'MyName'
        assert idx2.name == 'NewName'

        idx3 = idx.copy(names=['NewName'])

        assert idx.equals(idx3)
        assert idx.name == 'MyName'
        assert idx.names == ['MyName']
        assert idx3.name == 'NewName'
        assert idx3.names == ['NewName']

    def test_union_base(self):
        idx = self.create_index()
        first = idx[3:]
        second = idx[:5]

        if PY3:
            with tm.assert_produces_warning(RuntimeWarning):
                # unorderable types
                result = first.union(second)
                expected = Index(['b', 2, 'c', 0, 'a', 1])
                tm.assert_index_equal(result, expected)
        else:
            result = first.union(second)
            expected = Index(['b', 2, 'c', 0, 'a', 1])
            tm.assert_index_equal(result, expected)

        # GH 10149
        cases = [klass(second.values)
                 for klass in [np.array, Series, list]]
        for case in cases:
            if PY3:
                with tm.assert_produces_warning(RuntimeWarning):
                    # unorderable types
                    result = first.union(case)
                    assert tm.equalContents(result, idx)
            else:
                result = first.union(case)
                assert tm.equalContents(result, idx)

    def test_intersection_base(self):
        # (same results for py2 and py3 but sortedness not tested elsewhere)
        idx = self.create_index()
        first = idx[:5]
        second = idx[:3]
        result = first.intersection(second)
        expected = Index([0, 'a', 1])
        tm.assert_index_equal(result, expected)

        # GH 10149
        cases = [klass(second.values)
                 for klass in [np.array, Series, list]]
        for case in cases:
            result = first.intersection(case)
            assert tm.equalContents(result, second)

    def test_difference_base(self):
        # (same results for py2 and py3 but sortedness not tested elsewhere)
        idx = self.create_index()
        first = idx[:4]
        second = idx[3:]

        result = first.difference(second)
        expected = Index([0, 1, 'a'])
        tm.assert_index_equal(result, expected)

    def test_symmetric_difference(self):
        # (same results for py2 and py3 but sortedness not tested elsewhere)
        idx = self.create_index()
        first = idx[:4]
        second = idx[3:]

        result = first.symmetric_difference(second)
        expected = Index([0, 1, 2, 'a', 'c'])
        tm.assert_index_equal(result, expected)

    def test_logical_compat(self):
        idx = self.create_index()
        assert idx.all() == idx.values.all()
        assert idx.any() == idx.values.any()

    def test_dropna(self):
        # GH 6194
        for dtype in [None, object, 'category']:
            idx = pd.Index([1, 2, 3], dtype=dtype)
            tm.assert_index_equal(idx.dropna(), idx)

            idx = pd.Index([1., 2., 3.], dtype=dtype)
            tm.assert_index_equal(idx.dropna(), idx)
            nanidx = pd.Index([1., 2., np.nan, 3.], dtype=dtype)
            tm.assert_index_equal(nanidx.dropna(), idx)

            idx = pd.Index(['A', 'B', 'C'], dtype=dtype)
            tm.assert_index_equal(idx.dropna(), idx)
            nanidx = pd.Index(['A', np.nan, 'B', 'C'], dtype=dtype)
            tm.assert_index_equal(nanidx.dropna(), idx)

            tm.assert_index_equal(nanidx.dropna(how='any'), idx)
            tm.assert_index_equal(nanidx.dropna(how='all'), idx)

        idx = pd.DatetimeIndex(['2011-01-01', '2011-01-02', '2011-01-03'])
        tm.assert_index_equal(idx.dropna(), idx)
        nanidx = pd.DatetimeIndex(['2011-01-01', '2011-01-02',
                                   '2011-01-03', pd.NaT])
        tm.assert_index_equal(nanidx.dropna(), idx)

        idx = pd.TimedeltaIndex(['1 days', '2 days', '3 days'])
        tm.assert_index_equal(idx.dropna(), idx)
        nanidx = pd.TimedeltaIndex([pd.NaT, '1 days', '2 days',
                                    '3 days', pd.NaT])
        tm.assert_index_equal(nanidx.dropna(), idx)

        idx = pd.PeriodIndex(['2012-02', '2012-04', '2012-05'], freq='M')
        tm.assert_index_equal(idx.dropna(), idx)
        nanidx = pd.PeriodIndex(['2012-02', '2012-04', 'NaT', '2012-05'],
                                freq='M')
        tm.assert_index_equal(nanidx.dropna(), idx)

        msg = "invalid how option: xxx"
        with tm.assert_raises_regex(ValueError, msg):
            pd.Index([1, 2, 3]).dropna(how='xxx')

    def test_get_combined_index(self):
        result = _get_combined_index([])
        tm.assert_index_equal(result, Index([]))

    def test_repeat(self):
        repeats = 2
        idx = pd.Index([1, 2, 3])
        expected = pd.Index([1, 1, 2, 2, 3, 3])

        result = idx.repeat(repeats)
        tm.assert_index_equal(result, expected)

        with tm.assert_produces_warning(FutureWarning):
            result = idx.repeat(n=repeats)
            tm.assert_index_equal(result, expected)

    def test_is_monotonic_na(self):
        examples = [pd.Index([np.nan]),
                    pd.Index([np.nan, 1]),
                    pd.Index([1, 2, np.nan]),
                    pd.Index(['a', 'b', np.nan]),
                    pd.to_datetime(['NaT']),
                    pd.to_datetime(['NaT', '2000-01-01']),
                    pd.to_datetime(['2000-01-01', 'NaT', '2000-01-02']),
                    pd.to_timedelta(['1 day', 'NaT']), ]
        for index in examples:
            assert not index.is_monotonic_increasing
            assert not index.is_monotonic_decreasing
            assert not index._is_strictly_monotonic_increasing
            assert not index._is_strictly_monotonic_decreasing

    def test_repr_summary(self):
        with cf.option_context('display.max_seq_items', 10):
            r = repr(pd.Index(np.arange(1000)))
            assert len(r) < 200
            assert "..." in r

    def test_int_name_format(self):
        index = Index(['a', 'b', 'c'], name=0)
        s = Series(lrange(3), index)
        df = DataFrame(lrange(3), index=index)
        repr(s)
        repr(df)

    def test_print_unicode_columns(self):
        df = pd.DataFrame({u("\u05d0"): [1, 2, 3],
                           "\u05d1": [4, 5, 6],
                           "c": [7, 8, 9]})
        repr(df.columns)  # should not raise UnicodeDecodeError

    def test_unicode_string_with_unicode(self):
        idx = Index(lrange(1000))

        if PY3:
            str(idx)
        else:
            text_type(idx)

    def test_bytestring_with_unicode(self):
        idx = Index(lrange(1000))
        if PY3:
            bytes(idx)
        else:
            str(idx)

    def test_intersect_str_dates(self):
        dt_dates = [datetime(2012, 2, 9), datetime(2012, 2, 22)]

        i1 = Index(dt_dates, dtype=object)
        i2 = Index(['aa'], dtype=object)
        res = i2.intersection(i1)

        assert len(res) == 0

    @pytest.mark.parametrize('op', [operator.eq, operator.ne,
                                    operator.gt, operator.ge,
                                    operator.lt, operator.le])
    def test_comparison_tzawareness_compat(self, op):
        # GH#18162
        dr = pd.date_range('2016-01-01', periods=6)
        dz = dr.tz_localize('US/Pacific')

        # Check that there isn't a problem aware-aware and naive-naive do not
        # raise
        naive_series = Series(dr)
        aware_series = Series(dz)
        with pytest.raises(TypeError):
            op(dz, naive_series)
        with pytest.raises(TypeError):
            op(dr, aware_series)

        # TODO: implement _assert_tzawareness_compat for the reverse
        # comparison with the Series on the left-hand side


class TestIndexUtils(object):

    @pytest.mark.parametrize('data, names, expected', [
        ([[1, 2, 3]], None, Index([1, 2, 3])),
        ([[1, 2, 3]], ['name'], Index([1, 2, 3], name='name')),
        ([['a', 'a'], ['c', 'd']], None,
         MultiIndex([['a'], ['c', 'd']], [[0, 0], [0, 1]])),
        ([['a', 'a'], ['c', 'd']], ['L1', 'L2'],
         MultiIndex([['a'], ['c', 'd']], [[0, 0], [0, 1]],
                    names=['L1', 'L2'])),
    ])
    def test_ensure_index_from_sequences(self, data, names, expected):
        result = _ensure_index_from_sequences(data, names)
        tm.assert_index_equal(result, expected)


@pytest.mark.parametrize('opname', ['eq', 'ne', 'le', 'lt', 'ge', 'gt',
                                    'add', 'radd', 'sub', 'rsub',
                                    'mul', 'rmul', 'truediv', 'rtruediv',
                                    'floordiv', 'rfloordiv',
                                    'pow', 'rpow', 'mod', 'divmod'])
def test_generated_op_names(opname, indices):
    index = indices
    if isinstance(index, ABCIndex) and opname == 'rsub':
        # pd.Index.__rsub__ does not exist; though the method does exist
        # for subclasses.  see GH#19723
        return
    opname = '__{name}__'.format(name=opname)
    method = getattr(index, opname)
    assert method.__name__ == opname


@pytest.mark.parametrize('idx_maker', tm.index_subclass_makers_generator())
def test_index_subclass_constructor_wrong_kwargs(idx_maker):
    # GH #19348
    with tm.assert_raises_regex(TypeError, 'unexpected keyword argument'):
        idx_maker(foo='bar')
