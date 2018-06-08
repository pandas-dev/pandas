# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import pandas.util.testing as tm
import pytest
from pandas import (CategoricalIndex, DatetimeIndex, Float64Index, Index,
                    Int64Index, IntervalIndex, MultiIndex, PeriodIndex,
                    RangeIndex, Series, TimedeltaIndex, UInt64Index, compat,
                    isna)
from pandas._libs.tslib import iNaT
from pandas.compat import PY3
from pandas.core.indexes.base import InvalidIndexError
from pandas.core.indexes.datetimelike import DatetimeIndexOpsMixin


def verify_pickle(indices):
    unpickled = tm.round_trip_pickle(indices)
    assert indices.equals(unpickled)


def test_pickle_compat_construction(_holder):
    # this is testing for pickle compat
    if _holder is None:
        return

    # need an object to create with
    pytest.raises(TypeError, _holder)


def test_to_series(_index):
    # assert that we are creating a copy of the index

    idx = _index
    s = idx.to_series()
    assert s.values is not idx.values
    assert s.index is not idx
    assert s.name == idx.name


def test_to_series_with_arguments(_index):
    # GH18699

    # index kwarg
    idx = _index
    s = idx.to_series(index=idx)

    assert s.values is not idx.values
    assert s.index is idx
    assert s.name == idx.name

    # name kwarg
    idx = _index
    s = idx.to_series(name='__test')

    assert s.values is not idx.values
    assert s.index is not idx
    assert s.name != idx.name


def test_shift(_index):

    # GH8083 test the base class for shift
    idx = _index
    pytest.raises(NotImplementedError, idx.shift, 1)
    pytest.raises(NotImplementedError, idx.shift, 1, 2)


def test_create_index_existing_name(_index):

    # GH11193, when an existing index is passed, and a new name is not
    # specified, the new index should inherit the previous object name
    expected = _index
    if not isinstance(expected, MultiIndex):
        expected.name = 'foo'
        result = pd.Index(expected)
        tm.assert_index_equal(result, expected)

        result = pd.Index(expected, name='bar')
        expected.name = 'bar'
        tm.assert_index_equal(result, expected)
    else:
        expected.names = ['foo', 'bar']
        result = pd.Index(expected)
        tm.assert_index_equal(
            result, Index(Index([('foo', 'one'), ('foo', 'two'),
                                 ('bar', 'one'), ('baz', 'two'),
                                 ('qux', 'one'), ('qux', 'two')],
                                dtype='object'),
                          names=['foo', 'bar']))

        result = pd.Index(expected, names=['A', 'B'])
        tm.assert_index_equal(
            result,
            Index(Index([('foo', 'one'), ('foo', 'two'), ('bar', 'one'),
                         ('baz', 'two'), ('qux', 'one'), ('qux', 'two')],
                        dtype='object'), names=['A', 'B']))


def test_numeric_compat(_index):

    idx = _index
    tm.assert_raises_regex(TypeError, "cannot perform __mul__",
                           lambda: idx * 1)
    tm.assert_raises_regex(TypeError, "cannot perform __rmul__",
                           lambda: 1 * idx)

    div_err = "cannot perform __truediv__" if PY3 \
        else "cannot perform __div__"
    tm.assert_raises_regex(TypeError, div_err, lambda: idx / 1)
    div_err = div_err.replace(' __', ' __r')
    tm.assert_raises_regex(TypeError, div_err, lambda: 1 / idx)
    tm.assert_raises_regex(TypeError, "cannot perform __floordiv__",
                           lambda: idx // 1)
    tm.assert_raises_regex(TypeError, "cannot perform __rfloordiv__",
                           lambda: 1 // idx)


def test_logical_compat(_index):
    idx = _index
    tm.assert_raises_regex(TypeError, 'cannot perform all',
                           lambda: idx.all())
    tm.assert_raises_regex(TypeError, 'cannot perform any',
                           lambda: idx.any())


def test_boolean_context_compat(_index):

    # boolean context compat
    idx = _index

    def f():
        if idx:
            pass

    tm.assert_raises_regex(ValueError, 'The truth value of a', f)


def test_reindex_base(_index):
    idx = _index
    expected = np.arange(idx.size, dtype=np.intp)

    actual = idx.get_indexer(idx)
    tm.assert_numpy_array_equal(expected, actual)

    with tm.assert_raises_regex(ValueError, 'Invalid fill method'):
        idx.get_indexer(idx, method='invalid')


def test_get_indexer_consistency(named_index):
    # See GH 16819
    for name, index in named_index.items():
        if isinstance(index, IntervalIndex):
            continue

        if index.is_unique or isinstance(index, CategoricalIndex):
            indexer = index.get_indexer(index[0:2])
            assert isinstance(indexer, np.ndarray)
            assert indexer.dtype == np.intp
        else:
            e = "Reindexing only valid with uniquely valued Index objects"
            with tm.assert_raises_regex(InvalidIndexError, e):
                indexer = index.get_indexer(index[0:2])

        indexer, _ = index.get_indexer_non_unique(index[0:2])
        assert isinstance(indexer, np.ndarray)
        assert indexer.dtype == np.intp


def test_ndarray_compat_properties(_index, _compat_props):
    idx = _index
    assert idx.T.equals(idx)
    assert idx.transpose().equals(idx)

    values = idx.values
    for prop in _compat_props:
        assert getattr(idx, prop) == getattr(values, prop)

    # test for validity
    idx.nbytes
    idx.values.nbytes


def test_dtype_str(indices):
    dtype = indices.dtype_str
    assert isinstance(dtype, compat.string_types)
    assert dtype == str(indices.dtype)


def test_repr_max_seq_item_setting(_index):
    # GH10182
    idx = _index
    idx = idx.repeat(50)
    with pd.option_context("display.max_seq_items", None):
        repr(idx)
        assert '...' not in str(idx)


def test_wrong_number_names(indices):
    def testit(ind):
        ind.names = ["apple", "banana", "carrot"]
    tm.assert_raises_regex(ValueError, "^Length", testit, indices)


def test_hash_error(indices):
    index = indices
    tm.assert_raises_regex(TypeError, "unhashable type: %r" %
                           type(index).__name__, hash, indices)


def test_copy_name(named_index):
    # gh-12309: Check that the "name" argument
    # passed at initialization is honored.

    for name, index in compat.iteritems(named_index):
        if isinstance(index, MultiIndex):
            continue

        first = index.__class__(index, copy=True, name='mario')
        second = first.__class__(first, copy=False)

        # Even though "copy=False", we want a new object.
        assert first is not second

        # Not using tm.assert_index_equal() since names differ.
        assert index.equals(first)

        assert first.name == 'mario'
        assert second.name == 'mario'

        s1 = Series(2, index=first)
        s2 = Series(3, index=second[:-1])

        if not isinstance(index, CategoricalIndex):
            # See gh-13365
            s3 = s1 * s2
            assert s3.index.name == 'mario'


def test_ensure_copied_data(named_index):
    # Check the "copy" argument of each Index.__new__ is honoured
    # GH12309
    for name, index in compat.iteritems(named_index):
        init_kwargs = {}
        if isinstance(index, PeriodIndex):
            # Needs "freq" specification:
            init_kwargs['freq'] = index.freq
        elif isinstance(index, (RangeIndex, MultiIndex, CategoricalIndex)):
            # RangeIndex cannot be initialized from data
            # MultiIndex and CategoricalIndex are tested separately
            continue

        index_type = index.__class__
        result = index_type(index.values, copy=True, **init_kwargs)
        tm.assert_index_equal(index, result)
        tm.assert_numpy_array_equal(index.values, result.values,
                                    check_same='copy')

        if isinstance(index, PeriodIndex):
            # .values an object array of Period, thus copied
            result = index_type(ordinal=index.asi8, copy=False,
                                **init_kwargs)
            tm.assert_numpy_array_equal(index._ndarray_values,
                                        result._ndarray_values,
                                        check_same='same')
        elif isinstance(index, IntervalIndex):
            # checked in test_interval.py
            pass
        else:
            result = index_type(index.values, copy=False, **init_kwargs)
            tm.assert_numpy_array_equal(index.values, result.values,
                                        check_same='same')
            tm.assert_numpy_array_equal(index._ndarray_values,
                                        result._ndarray_values,
                                        check_same='same')


def test_copy_and_deepcopy(indices):
    from copy import copy, deepcopy

    if isinstance(indices, MultiIndex):
        return
    for func in (copy, deepcopy):
        idx_copy = func(indices)
        assert idx_copy is not indices
        assert idx_copy.equals(indices)

    new_copy = indices.copy(deep=True, name="banana")
    assert new_copy.name == "banana"


def test_unique_na():
    idx = pd.Index([2, np.nan, 2, 1], name='my_index')
    expected = pd.Index([2, np.nan, 1], name='my_index')
    result = idx.unique()
    tm.assert_index_equal(result, expected)


def test_sort(indices):
    pytest.raises(TypeError, indices.sort)


def test_mutability(indices):
    if not len(indices):
        return
    pytest.raises(TypeError, indices.__setitem__, 0, indices[0])


def test_compat(indices):
    assert indices.tolist() == list(indices)


def test_memory_usage(named_index):
    for name, index in compat.iteritems(named_index):
        result = index.memory_usage()
        if len(index):
            index.get_loc(index[0])
            result2 = index.memory_usage()
            result3 = index.memory_usage(deep=True)

            # RangeIndex, IntervalIndex
            # don't have engines
            if not isinstance(index, (RangeIndex, IntervalIndex)):
                assert result2 > result

            if index.inferred_type == 'object':
                assert result3 > result2

        else:

            # we report 0 for no-length
            assert result == 0


def test_numpy_argsort(named_index):
    for k, ind in named_index.items():
        result = np.argsort(ind)
        expected = ind.argsort()
        tm.assert_numpy_array_equal(result, expected)

        # these are the only two types that perform
        # pandas compatibility input validation - the
        # rest already perform separate (or no) such
        # validation via their 'values' attribute as
        # defined in pandas.core.indexes/base.py - they
        # cannot be changed at the moment due to
        # backwards compatibility concerns
        if isinstance(type(ind), (CategoricalIndex, RangeIndex)):
            msg = "the 'axis' parameter is not supported"
            tm.assert_raises_regex(ValueError, msg,
                                   np.argsort, ind, axis=1)

            msg = "the 'kind' parameter is not supported"
            tm.assert_raises_regex(ValueError, msg, np.argsort,
                                   ind, kind='mergesort')

            msg = "the 'order' parameter is not supported"
            tm.assert_raises_regex(ValueError, msg, np.argsort,
                                   ind, order=('a', 'b'))


def test_pickle(indices):
    verify_pickle(indices)
    original_name, indices.name = indices.name, 'foo'
    verify_pickle(indices)
    indices.name = original_name


def test_take(named_index):
    indexer = [4, 3, 0, 2]
    for k, ind in named_index.items():

        # separate
        if k in ['boolIndex', 'tuples', 'empty']:
            continue

        result = ind.take(indexer)
        expected = ind[indexer]
        assert result.equals(expected)

        if not isinstance(ind,
                          (DatetimeIndex, PeriodIndex, TimedeltaIndex)):
            # GH 10791
            with pytest.raises(AttributeError):
                ind.freq


def test_take_invalid_kwargs(_index):
    idx = _index
    indices = [1, 2]

    msg = r"take\(\) got an unexpected keyword argument 'foo'"
    tm.assert_raises_regex(TypeError, msg, idx.take,
                           indices, foo=2)

    msg = "the 'out' parameter is not supported"
    tm.assert_raises_regex(ValueError, msg, idx.take,
                           indices, out=indices)

    msg = "the 'mode' parameter is not supported"
    tm.assert_raises_regex(ValueError, msg, idx.take,
                           indices, mode='clip')


def test_setops_errorcases(named_index):
    for name, idx in compat.iteritems(named_index):
        # # non-iterable input
        cases = [0.5, 'xxx']
        methods = [idx.intersection, idx.union, idx.difference,
                   idx.symmetric_difference]

        for method in methods:
            for case in cases:
                tm.assert_raises_regex(TypeError,
                                       "Input must be Index "
                                       "or array-like",
                                       method, case)


def test_intersection_base(named_index):
    for name, idx in compat.iteritems(named_index):
        first = idx[:5]
        second = idx[:3]
        intersect = first.intersection(second)

        if isinstance(idx, CategoricalIndex):
            pass
        else:
            assert tm.equalContents(intersect, second)

        # GH 10149
        cases = [klass(second.values)
                 for klass in [np.array, Series, list]]
        for case in cases:
            if isinstance(idx, PeriodIndex):
                msg = "can only call with other PeriodIndex-ed objects"
                with tm.assert_raises_regex(ValueError, msg):
                    result = first.intersection(case)
            elif isinstance(idx, CategoricalIndex):
                pass
            else:
                result = first.intersection(case)
                assert tm.equalContents(result, second)

        if isinstance(idx, MultiIndex):
            msg = "other must be a MultiIndex or a list of tuples"
            with tm.assert_raises_regex(TypeError, msg):
                result = first.intersection([1, 2, 3])


def test_union_base(named_index):
    for name, idx in compat.iteritems(named_index):
        first = idx[3:]
        second = idx[:5]
        everything = idx
        union = first.union(second)
        assert tm.equalContents(union, everything)

        # GH 10149
        cases = [klass(second.values)
                 for klass in [np.array, Series, list]]
        for case in cases:
            if isinstance(idx, PeriodIndex):
                msg = "can only call with other PeriodIndex-ed objects"
                with tm.assert_raises_regex(ValueError, msg):
                    result = first.union(case)
            elif isinstance(idx, CategoricalIndex):
                pass
            else:
                result = first.union(case)
                assert tm.equalContents(result, everything)

        if isinstance(idx, MultiIndex):
            msg = "other must be a MultiIndex or a list of tuples"
            with tm.assert_raises_regex(TypeError, msg):
                result = first.union([1, 2, 3])


def test_difference_base(named_index):
    for name, idx in compat.iteritems(named_index):
        first = idx[2:]
        second = idx[:4]
        answer = idx[4:]
        result = first.difference(second)

        if isinstance(idx, CategoricalIndex):
            pass
        else:
            assert tm.equalContents(result, answer)

        # GH 10149
        cases = [klass(second.values)
                 for klass in [np.array, Series, list]]
        for case in cases:
            if isinstance(idx, PeriodIndex):
                msg = "can only call with other PeriodIndex-ed objects"
                with tm.assert_raises_regex(ValueError, msg):
                    result = first.difference(case)
            elif isinstance(idx, CategoricalIndex):
                pass
            elif isinstance(idx, (DatetimeIndex, TimedeltaIndex)):
                assert result.__class__ == answer.__class__
                tm.assert_numpy_array_equal(result.sort_values().asi8,
                                            answer.sort_values().asi8)
            else:
                result = first.difference(case)
                assert tm.equalContents(result, answer)

        if isinstance(idx, MultiIndex):
            msg = "other must be a MultiIndex or a list of tuples"
            with tm.assert_raises_regex(TypeError, msg):
                result = first.difference([1, 2, 3])


def test_symmetric_difference(named_index):
    for name, idx in compat.iteritems(named_index):
        first = idx[1:]
        second = idx[:-1]
        if isinstance(idx, CategoricalIndex):
            pass
        else:
            answer = idx[[0, -1]]
            result = first.symmetric_difference(second)
            assert tm.equalContents(result, answer)

        # GH 10149
        cases = [klass(second.values)
                 for klass in [np.array, Series, list]]
        for case in cases:
            if isinstance(idx, PeriodIndex):
                msg = "can only call with other PeriodIndex-ed objects"
                with tm.assert_raises_regex(ValueError, msg):
                    result = first.symmetric_difference(case)
            elif isinstance(idx, CategoricalIndex):
                pass
            else:
                result = first.symmetric_difference(case)
                assert tm.equalContents(result, answer)

        if isinstance(idx, MultiIndex):
            msg = "other must be a MultiIndex or a list of tuples"
            with tm.assert_raises_regex(TypeError, msg):
                first.symmetric_difference([1, 2, 3])


def test_insert_base(named_index):

    for name, idx in compat.iteritems(named_index):
        result = idx[1:4]

        if not len(idx):
            continue

        # test 0th element
        assert idx[0:4].equals(result.insert(0, idx[0]))


def test_delete_base(named_index):

    for name, idx in compat.iteritems(named_index):

        if not len(idx):
            continue

        if isinstance(idx, RangeIndex):
            # tested in class
            continue

        expected = idx[1:]
        result = idx.delete(0)
        assert result.equals(expected)
        assert result.name == expected.name

        expected = idx[:-1]
        result = idx.delete(-1)
        assert result.equals(expected)
        assert result.name == expected.name

        with pytest.raises((IndexError, ValueError)):
            # either depending on numpy version
            result = idx.delete(len(idx))


def test_equals(named_index):

    for name, idx in compat.iteritems(named_index):
        assert idx.equals(idx)
        assert idx.equals(idx.copy())
        assert idx.equals(idx.astype(object))

        assert not idx.equals(list(idx))
        assert not idx.equals(np.array(idx))

        # Cannot pass in non-int64 dtype to RangeIndex
        if not isinstance(idx, RangeIndex):
            same_values = Index(idx, dtype=object)
            assert idx.equals(same_values)
            assert same_values.equals(idx)

        if idx.nlevels == 1:
            # do not test MultiIndex
            assert not idx.equals(pd.Series(idx))


def test_equals_op(_index):
    # GH9947, GH10637
    index_a = _index
    if isinstance(index_a, PeriodIndex):
        return

    n = len(index_a)
    index_b = index_a[0:-1]
    index_c = index_a[0:-1].append(index_a[-2:-1])
    index_d = index_a[0:1]
    with tm.assert_raises_regex(ValueError, "Lengths must match"):
        index_a == index_b
    expected1 = np.array([True] * n)
    expected2 = np.array([True] * (n - 1) + [False])
    tm.assert_numpy_array_equal(index_a == index_a, expected1)
    tm.assert_numpy_array_equal(index_a == index_c, expected2)

    # test comparisons with numpy arrays
    array_a = np.array(index_a)
    array_b = np.array(index_a[0:-1])
    array_c = np.array(index_a[0:-1].append(index_a[-2:-1]))
    array_d = np.array(index_a[0:1])
    with tm.assert_raises_regex(ValueError, "Lengths must match"):
        index_a == array_b
    tm.assert_numpy_array_equal(index_a == array_a, expected1)
    tm.assert_numpy_array_equal(index_a == array_c, expected2)

    # test comparisons with Series
    series_a = Series(array_a)
    series_b = Series(array_b)
    series_c = Series(array_c)
    series_d = Series(array_d)
    with tm.assert_raises_regex(ValueError, "Lengths must match"):
        index_a == series_b

    tm.assert_numpy_array_equal(index_a == series_a, expected1)
    tm.assert_numpy_array_equal(index_a == series_c, expected2)

    # cases where length is 1 for one of them
    with tm.assert_raises_regex(ValueError, "Lengths must match"):
        index_a == index_d
    with tm.assert_raises_regex(ValueError, "Lengths must match"):
        index_a == series_d
    with tm.assert_raises_regex(ValueError, "Lengths must match"):
        index_a == array_d
    msg = "Can only compare identically-labeled Series objects"
    with tm.assert_raises_regex(ValueError, msg):
        series_a == series_d
    with tm.assert_raises_regex(ValueError, "Lengths must match"):
        series_a == array_d

    # comparing with a scalar should broadcast; note that we are excluding
    # MultiIndex because in this case each item in the index is a tuple of
    # length 2, and therefore is considered an array of length 2 in the
    # comparison instead of a scalar
    if not isinstance(index_a, MultiIndex):
        expected3 = np.array([False] * (len(index_a) - 2) + [True, False])
        # assuming the 2nd to last item is unique in the data
        item = index_a[-2]
        tm.assert_numpy_array_equal(index_a == item, expected3)
        tm.assert_series_equal(series_a == item, Series(expected3))


def test_numpy_ufuncs(named_index):
    # test ufuncs of numpy 1.9.2. see:
    # http://docs.scipy.org/doc/numpy/reference/ufuncs.html

    # some functions are skipped because it may return different result
    # for unicode input depending on numpy version

    for name, idx in compat.iteritems(named_index):
        for func in [np.exp, np.exp2, np.expm1, np.log, np.log2, np.log10,
                     np.log1p, np.sqrt, np.sin, np.cos, np.tan, np.arcsin,
                     np.arccos, np.arctan, np.sinh, np.cosh, np.tanh,
                     np.arcsinh, np.arccosh, np.arctanh, np.deg2rad,
                     np.rad2deg]:
            if isinstance(idx, DatetimeIndexOpsMixin):
                # raise TypeError or ValueError (PeriodIndex)
                # PeriodIndex behavior should be changed in future version
                with pytest.raises(Exception):
                    with np.errstate(all='ignore'):
                        func(idx)
            elif isinstance(idx, (Float64Index, Int64Index, UInt64Index)):
                # coerces to float (e.g. np.sin)
                with np.errstate(all='ignore'):
                    result = func(idx)
                    exp = Index(func(idx.values), name=idx.name)

                tm.assert_index_equal(result, exp)
                assert isinstance(result, pd.Float64Index)
            else:
                # raise AttributeError or TypeError
                if len(idx) == 0:
                    continue
                else:
                    with pytest.raises(Exception):
                        with np.errstate(all='ignore'):
                            func(idx)

        for func in [np.isfinite, np.isinf, np.isnan, np.signbit]:
            if isinstance(idx, DatetimeIndexOpsMixin):
                # raise TypeError or ValueError (PeriodIndex)
                with pytest.raises(Exception):
                    func(idx)
            elif isinstance(idx, (Float64Index, Int64Index, UInt64Index)):
                # Results in bool array
                result = func(idx)
                assert isinstance(result, np.ndarray)
                assert not isinstance(result, Index)
            else:
                if len(idx) == 0:
                    continue
                else:
                    with pytest.raises(Exception):
                        func(idx)


def test_hasnans_isnans(named_index):
    # GH 11343, added tests for hasnans / isnans
    for name, index in named_index.items():
        if isinstance(index, MultiIndex):
            pass
        else:
            idx = index.copy()

            # cases in indices doesn't include NaN
            expected = np.array([False] * len(idx), dtype=bool)
            tm.assert_numpy_array_equal(idx._isnan, expected)
            assert not idx.hasnans

            idx = index.copy()
            values = idx.values

            if len(index) == 0:
                continue
            elif isinstance(index, DatetimeIndexOpsMixin):
                values[1] = iNaT
            elif isinstance(index, (Int64Index, UInt64Index)):
                continue
            else:
                values[1] = np.nan

            if isinstance(index, PeriodIndex):
                idx = index.__class__(values, freq=index.freq)
            else:
                idx = index.__class__(values)

            expected = np.array([False] * len(idx), dtype=bool)
            expected[1] = True
            tm.assert_numpy_array_equal(idx._isnan, expected)
            assert idx.hasnans


def test_fillna(named_index):
    # GH 11343
    for name, index in named_index.items():
        if len(index) == 0:
            pass
        elif isinstance(index, MultiIndex):
            idx = index.copy()
            msg = "isna is not defined for MultiIndex"
            with tm.assert_raises_regex(NotImplementedError, msg):
                idx.fillna(idx[0])
        else:
            idx = index.copy()
            result = idx.fillna(idx[0])
            tm.assert_index_equal(result, idx)
            assert result is not idx

            msg = "'value' must be a scalar, passed: "
            with tm.assert_raises_regex(TypeError, msg):
                idx.fillna([idx[0]])

            idx = index.copy()
            values = idx.values

            if isinstance(index, DatetimeIndexOpsMixin):
                values[1] = iNaT
            elif isinstance(index, (Int64Index, UInt64Index)):
                continue
            else:
                values[1] = np.nan

            if isinstance(index, PeriodIndex):
                idx = index.__class__(values, freq=index.freq)
            else:
                idx = index.__class__(values)

            expected = np.array([False] * len(idx), dtype=bool)
            expected[1] = True
            tm.assert_numpy_array_equal(idx._isnan, expected)
            assert idx.hasnans


def test_nulls(named_index):
    # this is really a smoke test for the methods
    # as these are adequately tested for function elsewhere

    for name, index in named_index.items():
        if len(index) == 0:
            tm.assert_numpy_array_equal(
                index.isna(), np.array([], dtype=bool))
        elif isinstance(index, MultiIndex):
            idx = index.copy()
            msg = "isna is not defined for MultiIndex"
            with tm.assert_raises_regex(NotImplementedError, msg):
                idx.isna()
        else:

            if not index.hasnans:
                tm.assert_numpy_array_equal(
                    index.isna(), np.zeros(len(index), dtype=bool))
                tm.assert_numpy_array_equal(
                    index.notna(), np.ones(len(index), dtype=bool))
            else:
                result = isna(index)
                tm.assert_numpy_array_equal(index.isna(), result)
                tm.assert_numpy_array_equal(index.notna(), ~result)


def test_empty(_index):
    # GH 15270
    index = _index
    assert not index.empty
    assert index[:0].empty


def test_join_self_unique(_index, join_type):
    index = _index
    if index.is_unique:
        joined = index.join(index, how=join_type)
        assert (index == joined).all()


def test_searchsorted_monotonic(indices):
    # GH17271
    # not implemented for tuple searches in MultiIndex
    # or Intervals searches in IntervalIndex
    if isinstance(indices, (MultiIndex, IntervalIndex)):
        return

    # nothing to test if the index is empty
    if indices.empty:
        return
    value = indices[0]

    # determine the expected results (handle dupes for 'right')
    expected_left, expected_right = 0, (indices == value).argmin()
    if expected_right == 0:
        # all values are the same, expected_right should be length
        expected_right = len(indices)

    # test _searchsorted_monotonic in all cases
    # test searchsorted only for increasing
    if indices.is_monotonic_increasing:
        ssm_left = indices._searchsorted_monotonic(value, side='left')
        assert expected_left == ssm_left

        ssm_right = indices._searchsorted_monotonic(value, side='right')
        assert expected_right == ssm_right

        ss_left = indices.searchsorted(value, side='left')
        assert expected_left == ss_left

        ss_right = indices.searchsorted(value, side='right')
        assert expected_right == ss_right

    elif indices.is_monotonic_decreasing:
        ssm_left = indices._searchsorted_monotonic(value, side='left')
        assert expected_left == ssm_left

        ssm_right = indices._searchsorted_monotonic(value, side='right')
        assert expected_right == ssm_right

    else:
        # non-monotonic should raise.
        with pytest.raises(ValueError):
            indices._searchsorted_monotonic(value, side='left')


def test_map(_index):
    # callable
    index = _index

    # we don't infer UInt64
    if isinstance(index, pd.UInt64Index):
        expected = index.astype('int64')
    else:
        expected = index

    result = index.map(lambda x: x)
    tm.assert_index_equal(result, expected)


@pytest.mark.parametrize(
    "mapper",
    [
        lambda values, index: {i: e for e, i in zip(values, index)},
        lambda values, index: pd.Series(values, index)])
def test_map_dictlike(_index, mapper):

    index = _index
    if isinstance(index, (pd.CategoricalIndex, pd.IntervalIndex)):
        pytest.skip("skipping tests for {}".format(type(index)))

    identity = mapper(index.values, index)

    # we don't infer to UInt64 for a dict
    if isinstance(index, pd.UInt64Index) and isinstance(identity, dict):
        expected = index.astype('int64')
    else:
        expected = index

    result = index.map(identity)
    tm.assert_index_equal(result, expected)

    # empty mappable
    expected = pd.Index([np.nan] * len(index))
    result = index.map(mapper(expected, index))
    tm.assert_index_equal(result, expected)


def test_putmask_with_wrong_mask(_index):
    # GH18368
    index = _index

    with pytest.raises(ValueError):
        index.putmask(np.ones(len(index) + 1, np.bool), 1)

    with pytest.raises(ValueError):
        index.putmask(np.ones(len(index) - 1, np.bool), 1)

    with pytest.raises(ValueError):
        index.putmask('foo', 1)
