import gc
from typing import Optional, Type

import numpy as np
import pytest

from pandas._libs.tslib import iNaT

from pandas.core.dtypes.common import is_datetime64tz_dtype
from pandas.core.dtypes.dtypes import CategoricalDtype

import pandas as pd
from pandas import (
    CategoricalIndex,
    DatetimeIndex,
    Index,
    Int64Index,
    IntervalIndex,
    MultiIndex,
    PeriodIndex,
    RangeIndex,
    Series,
    TimedeltaIndex,
    UInt64Index,
    isna,
)
import pandas._testing as tm
from pandas.core.indexes.base import InvalidIndexError
from pandas.core.indexes.datetimelike import DatetimeIndexOpsMixin


class Base:
    """ base class for index sub-class tests """

    _holder: Optional[Type[Index]] = None
    _compat_props = ["shape", "ndim", "size", "nbytes"]

    def test_pickle_compat_construction(self):
        # need an object to create with
        msg = (
            r"Index\(\.\.\.\) must be called with a collection of some "
            r"kind, None was passed|"
            r"__new__\(\) missing 1 required positional argument: 'data'|"
            r"__new__\(\) takes at least 2 arguments \(1 given\)"
        )
        with pytest.raises(TypeError, match=msg):
            self._holder()

    def test_to_series(self):
        # assert that we are creating a copy of the index

        idx = self.create_index()
        s = idx.to_series()
        assert s.values is not idx.values
        assert s.index is not idx
        assert s.name == idx.name

    def test_to_series_with_arguments(self):
        # GH18699

        # index kwarg
        idx = self.create_index()
        s = idx.to_series(index=idx)

        assert s.values is not idx.values
        assert s.index is idx
        assert s.name == idx.name

        # name kwarg
        idx = self.create_index()
        s = idx.to_series(name="__test")

        assert s.values is not idx.values
        assert s.index is not idx
        assert s.name != idx.name

    @pytest.mark.parametrize("name", [None, "new_name"])
    def test_to_frame(self, name):
        # see GH-15230, GH-22580
        idx = self.create_index()

        if name:
            idx_name = name
        else:
            idx_name = idx.name or 0

        df = idx.to_frame(name=idx_name)

        assert df.index is idx
        assert len(df.columns) == 1
        assert df.columns[0] == idx_name
        assert df[idx_name].values is not idx.values

        df = idx.to_frame(index=False, name=idx_name)
        assert df.index is not idx

    def test_shift(self):

        # GH8083 test the base class for shift
        idx = self.create_index()
        msg = f"Not supported for type {type(idx).__name__}"
        with pytest.raises(NotImplementedError, match=msg):
            idx.shift(1)
        with pytest.raises(NotImplementedError, match=msg):
            idx.shift(1, 2)

    def test_constructor_name_unhashable(self):
        # GH#29069 check that name is hashable
        # See also same-named test in tests.series.test_constructors
        idx = self.create_index()
        with pytest.raises(TypeError, match="Index.name must be a hashable type"):
            type(idx)(idx, name=[])

    def test_create_index_existing_name(self):

        # GH11193, when an existing index is passed, and a new name is not
        # specified, the new index should inherit the previous object name
        expected = self.create_index()
        if not isinstance(expected, MultiIndex):
            expected.name = "foo"
            result = pd.Index(expected)
            tm.assert_index_equal(result, expected)

            result = pd.Index(expected, name="bar")
            expected.name = "bar"
            tm.assert_index_equal(result, expected)
        else:
            expected.names = ["foo", "bar"]
            result = pd.Index(expected)
            tm.assert_index_equal(
                result,
                Index(
                    Index(
                        [
                            ("foo", "one"),
                            ("foo", "two"),
                            ("bar", "one"),
                            ("baz", "two"),
                            ("qux", "one"),
                            ("qux", "two"),
                        ],
                        dtype="object",
                    ),
                    names=["foo", "bar"],
                ),
            )

            result = pd.Index(expected, names=["A", "B"])
            tm.assert_index_equal(
                result,
                Index(
                    Index(
                        [
                            ("foo", "one"),
                            ("foo", "two"),
                            ("bar", "one"),
                            ("baz", "two"),
                            ("qux", "one"),
                            ("qux", "two"),
                        ],
                        dtype="object",
                    ),
                    names=["A", "B"],
                ),
            )

    def test_numeric_compat(self):

        idx = self.create_index()
        # Check that this doesn't cover MultiIndex case, if/when it does,
        #  we can remove multi.test_compat.test_numeric_compat
        assert not isinstance(idx, MultiIndex)

        with pytest.raises(TypeError, match="cannot perform __mul__"):
            idx * 1
        with pytest.raises(TypeError, match="cannot perform __rmul__"):
            1 * idx

        div_err = "cannot perform __truediv__"
        with pytest.raises(TypeError, match=div_err):
            idx / 1

        div_err = div_err.replace(" __", " __r")
        with pytest.raises(TypeError, match=div_err):
            1 / idx
        with pytest.raises(TypeError, match="cannot perform __floordiv__"):
            idx // 1
        with pytest.raises(TypeError, match="cannot perform __rfloordiv__"):
            1 // idx

    def test_logical_compat(self):
        idx = self.create_index()
        with pytest.raises(TypeError, match="cannot perform all"):
            idx.all()
        with pytest.raises(TypeError, match="cannot perform any"):
            idx.any()

    def test_boolean_context_compat(self):

        # boolean context compat
        idx = self.create_index()

        with pytest.raises(ValueError, match="The truth value of a"):
            if idx:
                pass

    def test_reindex_base(self):
        idx = self.create_index()
        expected = np.arange(idx.size, dtype=np.intp)

        actual = idx.get_indexer(idx)
        tm.assert_numpy_array_equal(expected, actual)

        with pytest.raises(ValueError, match="Invalid fill method"):
            idx.get_indexer(idx, method="invalid")

    def test_get_indexer_consistency(self, indices):
        # See GH 16819
        if isinstance(indices, IntervalIndex):
            return

        if indices.is_unique or isinstance(indices, CategoricalIndex):
            indexer = indices.get_indexer(indices[0:2])
            assert isinstance(indexer, np.ndarray)
            assert indexer.dtype == np.intp
        else:
            e = "Reindexing only valid with uniquely valued Index objects"
            with pytest.raises(InvalidIndexError, match=e):
                indices.get_indexer(indices[0:2])

        indexer, _ = indices.get_indexer_non_unique(indices[0:2])
        assert isinstance(indexer, np.ndarray)
        assert indexer.dtype == np.intp

    def test_ndarray_compat_properties(self):
        idx = self.create_index()
        assert idx.T.equals(idx)
        assert idx.transpose().equals(idx)

        values = idx.values
        for prop in self._compat_props:
            assert getattr(idx, prop) == getattr(values, prop)

        # test for validity
        idx.nbytes
        idx.values.nbytes

    def test_repr_roundtrip(self):

        idx = self.create_index()
        tm.assert_index_equal(eval(repr(idx)), idx)

    def test_str(self):

        # test the string repr
        idx = self.create_index()
        idx.name = "foo"
        assert "'foo'" in str(idx)
        assert type(idx).__name__ in str(idx)

    def test_repr_max_seq_item_setting(self):
        # GH10182
        idx = self.create_index()
        idx = idx.repeat(50)
        with pd.option_context("display.max_seq_items", None):
            repr(idx)
            assert "..." not in str(idx)

    def test_copy_name(self, indices):
        # gh-12309: Check that the "name" argument
        # passed at initialization is honored.
        if isinstance(indices, MultiIndex):
            return

        first = type(indices)(indices, copy=True, name="mario")
        second = type(first)(first, copy=False)

        # Even though "copy=False", we want a new object.
        assert first is not second

        # Not using tm.assert_index_equal() since names differ.
        assert indices.equals(first)

        assert first.name == "mario"
        assert second.name == "mario"

        s1 = Series(2, index=first)
        s2 = Series(3, index=second[:-1])

        if not isinstance(indices, CategoricalIndex):
            # See gh-13365
            s3 = s1 * s2
            assert s3.index.name == "mario"

    def test_ensure_copied_data(self, indices):
        # Check the "copy" argument of each Index.__new__ is honoured
        # GH12309
        init_kwargs = {}
        if isinstance(indices, PeriodIndex):
            # Needs "freq" specification:
            init_kwargs["freq"] = indices.freq
        elif isinstance(indices, (RangeIndex, MultiIndex, CategoricalIndex)):
            # RangeIndex cannot be initialized from data
            # MultiIndex and CategoricalIndex are tested separately
            return

        index_type = type(indices)
        result = index_type(indices.values, copy=True, **init_kwargs)
        if is_datetime64tz_dtype(indices.dtype):
            result = result.tz_localize("UTC").tz_convert(indices.tz)

        tm.assert_index_equal(indices, result)
        tm.assert_numpy_array_equal(
            indices._ndarray_values, result._ndarray_values, check_same="copy"
        )

        if isinstance(indices, PeriodIndex):
            # .values an object array of Period, thus copied
            result = index_type(ordinal=indices.asi8, copy=False, **init_kwargs)
            tm.assert_numpy_array_equal(
                indices._ndarray_values, result._ndarray_values, check_same="same"
            )
        elif isinstance(indices, IntervalIndex):
            # checked in test_interval.py
            pass
        else:
            result = index_type(indices.values, copy=False, **init_kwargs)
            tm.assert_numpy_array_equal(
                indices.values, result.values, check_same="same"
            )
            tm.assert_numpy_array_equal(
                indices._ndarray_values, result._ndarray_values, check_same="same"
            )

    def test_memory_usage(self, indices):
        indices._engine.clear_mapping()
        result = indices.memory_usage()
        if indices.empty:
            # we report 0 for no-length
            assert result == 0
            return

        # non-zero length
        indices.get_loc(indices[0])
        result2 = indices.memory_usage()
        result3 = indices.memory_usage(deep=True)

        # RangeIndex, IntervalIndex
        # don't have engines
        if not isinstance(indices, (RangeIndex, IntervalIndex)):
            assert result2 > result

        if indices.inferred_type == "object":
            assert result3 > result2

    def test_argsort(self, request, indices):
        # separately tested
        if isinstance(indices, CategoricalIndex):
            return

        result = indices.argsort()
        expected = np.array(indices).argsort()
        tm.assert_numpy_array_equal(result, expected, check_dtype=False)

    def test_numpy_argsort(self, indices):
        result = np.argsort(indices)
        expected = indices.argsort()
        tm.assert_numpy_array_equal(result, expected)

        # these are the only two types that perform
        # pandas compatibility input validation - the
        # rest already perform separate (or no) such
        # validation via their 'values' attribute as
        # defined in pandas.core.indexes/base.py - they
        # cannot be changed at the moment due to
        # backwards compatibility concerns
        if isinstance(type(indices), (CategoricalIndex, RangeIndex)):
            msg = "the 'axis' parameter is not supported"
            with pytest.raises(ValueError, match=msg):
                np.argsort(indices, axis=1)

            msg = "the 'kind' parameter is not supported"
            with pytest.raises(ValueError, match=msg):
                np.argsort(indices, kind="mergesort")

            msg = "the 'order' parameter is not supported"
            with pytest.raises(ValueError, match=msg):
                np.argsort(indices, order=("a", "b"))

    def test_take(self, indices):
        indexer = [4, 3, 0, 2]
        if len(indices) < 5:
            # not enough elements; ignore
            return

        result = indices.take(indexer)
        expected = indices[indexer]
        assert result.equals(expected)

        if not isinstance(indices, (DatetimeIndex, PeriodIndex, TimedeltaIndex)):
            # GH 10791
            with pytest.raises(AttributeError):
                indices.freq

    def test_take_invalid_kwargs(self):
        idx = self.create_index()
        indices = [1, 2]

        msg = r"take\(\) got an unexpected keyword argument 'foo'"
        with pytest.raises(TypeError, match=msg):
            idx.take(indices, foo=2)

        msg = "the 'out' parameter is not supported"
        with pytest.raises(ValueError, match=msg):
            idx.take(indices, out=indices)

        msg = "the 'mode' parameter is not supported"
        with pytest.raises(ValueError, match=msg):
            idx.take(indices, mode="clip")

    def test_repeat(self):
        rep = 2
        i = self.create_index()
        expected = pd.Index(i.values.repeat(rep), name=i.name)
        tm.assert_index_equal(i.repeat(rep), expected)

        i = self.create_index()
        rep = np.arange(len(i))
        expected = pd.Index(i.values.repeat(rep), name=i.name)
        tm.assert_index_equal(i.repeat(rep), expected)

    def test_numpy_repeat(self):
        rep = 2
        i = self.create_index()
        expected = i.repeat(rep)
        tm.assert_index_equal(np.repeat(i, rep), expected)

        msg = "the 'axis' parameter is not supported"
        with pytest.raises(ValueError, match=msg):
            np.repeat(i, rep, axis=0)

    @pytest.mark.parametrize("klass", [list, tuple, np.array, Series])
    def test_where(self, klass):
        i = self.create_index()

        cond = [True] * len(i)
        result = i.where(klass(cond))
        expected = i
        tm.assert_index_equal(result, expected)

        cond = [False] + [True] * len(i[1:])
        expected = pd.Index([i._na_value] + i[1:].tolist(), dtype=i.dtype)
        result = i.where(klass(cond))
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("case", [0.5, "xxx"])
    @pytest.mark.parametrize(
        "method", ["intersection", "union", "difference", "symmetric_difference"]
    )
    def test_set_ops_error_cases(self, case, method, indices):
        # non-iterable input
        msg = "Input must be Index or array-like"
        with pytest.raises(TypeError, match=msg):
            getattr(indices, method)(case)

    def test_intersection_base(self, indices):
        if isinstance(indices, CategoricalIndex):
            return

        first = indices[:5]
        second = indices[:3]
        intersect = first.intersection(second)
        assert tm.equalContents(intersect, second)

        if is_datetime64tz_dtype(indices.dtype):
            # The second.values below will drop tz, so the rest of this test
            #  is not applicable.
            return

        # GH 10149
        cases = [klass(second.values) for klass in [np.array, Series, list]]
        for case in cases:
            result = first.intersection(case)
            assert tm.equalContents(result, second)

        if isinstance(indices, MultiIndex):
            msg = "other must be a MultiIndex or a list of tuples"
            with pytest.raises(TypeError, match=msg):
                first.intersection([1, 2, 3])

    def test_union_base(self, indices):
        first = indices[3:]
        second = indices[:5]
        everything = indices
        union = first.union(second)
        assert tm.equalContents(union, everything)

        if is_datetime64tz_dtype(indices.dtype):
            # The second.values below will drop tz, so the rest of this test
            #  is not applicable.
            return

        # GH 10149
        cases = [klass(second.values) for klass in [np.array, Series, list]]
        for case in cases:
            if not isinstance(indices, CategoricalIndex):
                result = first.union(case)
                assert tm.equalContents(result, everything)

        if isinstance(indices, MultiIndex):
            msg = "other must be a MultiIndex or a list of tuples"
            with pytest.raises(TypeError, match=msg):
                first.union([1, 2, 3])

    def test_difference_base(self, sort, indices):
        first = indices[2:]
        second = indices[:4]
        if isinstance(indices, CategoricalIndex) or indices.is_boolean():
            answer = []
        else:
            answer = indices[4:]
        result = first.difference(second, sort)
        assert tm.equalContents(result, answer)

        # GH 10149
        cases = [klass(second.values) for klass in [np.array, Series, list]]
        for case in cases:
            if isinstance(indices, (DatetimeIndex, TimedeltaIndex)):
                assert type(result) == type(answer)
                tm.assert_numpy_array_equal(
                    result.sort_values().asi8, answer.sort_values().asi8
                )
            else:
                result = first.difference(case, sort)
                assert tm.equalContents(result, answer)

        if isinstance(indices, MultiIndex):
            msg = "other must be a MultiIndex or a list of tuples"
            with pytest.raises(TypeError, match=msg):
                first.difference([1, 2, 3], sort)

    def test_symmetric_difference(self, indices):
        if isinstance(indices, CategoricalIndex):
            return

        first = indices[1:]
        second = indices[:-1]
        answer = indices[[0, -1]]
        result = first.symmetric_difference(second)
        assert tm.equalContents(result, answer)

        # GH 10149
        cases = [klass(second.values) for klass in [np.array, Series, list]]
        for case in cases:
            result = first.symmetric_difference(case)
            assert tm.equalContents(result, answer)

        if isinstance(indices, MultiIndex):
            msg = "other must be a MultiIndex or a list of tuples"
            with pytest.raises(TypeError, match=msg):
                first.symmetric_difference([1, 2, 3])

    def test_insert_base(self, indices):
        result = indices[1:4]

        if not len(indices):
            return

        # test 0th element
        assert indices[0:4].equals(result.insert(0, indices[0]))

    def test_delete_base(self, indices):
        if not len(indices):
            return

        if isinstance(indices, RangeIndex):
            # tested in class
            return

        expected = indices[1:]
        result = indices.delete(0)
        assert result.equals(expected)
        assert result.name == expected.name

        expected = indices[:-1]
        result = indices.delete(-1)
        assert result.equals(expected)
        assert result.name == expected.name

        with pytest.raises((IndexError, ValueError)):
            # either depending on numpy version
            indices.delete(len(indices))

    def test_equals(self, indices):
        if isinstance(indices, IntervalIndex):
            # IntervalIndex tested separately
            return

        assert indices.equals(indices)
        assert indices.equals(indices.copy())
        assert indices.equals(indices.astype(object))

        assert not indices.equals(list(indices))
        assert not indices.equals(np.array(indices))

        # Cannot pass in non-int64 dtype to RangeIndex
        if not isinstance(indices, (RangeIndex, CategoricalIndex)):
            # TODO: CategoricalIndex can be re-allowed following GH#32167
            same_values = Index(indices, dtype=object)
            assert indices.equals(same_values)
            assert same_values.equals(indices)

        if indices.nlevels == 1:
            # do not test MultiIndex
            assert not indices.equals(Series(indices))

    def test_equals_op(self):
        # GH9947, GH10637
        index_a = self.create_index()
        if isinstance(index_a, PeriodIndex):
            pytest.skip("Skip check for PeriodIndex")

        n = len(index_a)
        index_b = index_a[0:-1]
        index_c = index_a[0:-1].append(index_a[-2:-1])
        index_d = index_a[0:1]

        msg = "Lengths must match|could not be broadcast"
        with pytest.raises(ValueError, match=msg):
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
        with pytest.raises(ValueError, match=msg):
            index_a == array_b
        tm.assert_numpy_array_equal(index_a == array_a, expected1)
        tm.assert_numpy_array_equal(index_a == array_c, expected2)

        # test comparisons with Series
        series_a = Series(array_a)
        series_b = Series(array_b)
        series_c = Series(array_c)
        series_d = Series(array_d)
        with pytest.raises(ValueError, match=msg):
            index_a == series_b

        tm.assert_numpy_array_equal(index_a == series_a, expected1)
        tm.assert_numpy_array_equal(index_a == series_c, expected2)

        # cases where length is 1 for one of them
        with pytest.raises(ValueError, match="Lengths must match"):
            index_a == index_d
        with pytest.raises(ValueError, match="Lengths must match"):
            index_a == series_d
        with pytest.raises(ValueError, match="Lengths must match"):
            index_a == array_d
        msg = "Can only compare identically-labeled Series objects"
        with pytest.raises(ValueError, match=msg):
            series_a == series_d
        with pytest.raises(ValueError, match="Lengths must match"):
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

    def test_hasnans_isnans(self, indices):
        # GH 11343, added tests for hasnans / isnans
        if isinstance(indices, MultiIndex):
            return

        # cases in indices doesn't include NaN
        idx = indices.copy(deep=True)
        expected = np.array([False] * len(idx), dtype=bool)
        tm.assert_numpy_array_equal(idx._isnan, expected)
        assert idx.hasnans is False

        idx = indices.copy(deep=True)
        values = np.asarray(idx.values)

        if len(indices) == 0:
            return
        elif isinstance(indices, DatetimeIndexOpsMixin):
            values[1] = iNaT
        elif isinstance(indices, (Int64Index, UInt64Index)):
            return
        else:
            values[1] = np.nan

        if isinstance(indices, PeriodIndex):
            idx = type(indices)(values, freq=indices.freq)
        else:
            idx = type(indices)(values)

            expected = np.array([False] * len(idx), dtype=bool)
            expected[1] = True
            tm.assert_numpy_array_equal(idx._isnan, expected)
            assert idx.hasnans is True

    def test_fillna(self, indices):
        # GH 11343
        if len(indices) == 0:
            pass
        elif isinstance(indices, MultiIndex):
            idx = indices.copy(deep=True)
            msg = "isna is not defined for MultiIndex"
            with pytest.raises(NotImplementedError, match=msg):
                idx.fillna(idx[0])
        else:
            idx = indices.copy(deep=True)
            result = idx.fillna(idx[0])
            tm.assert_index_equal(result, idx)
            assert result is not idx

            msg = "'value' must be a scalar, passed: "
            with pytest.raises(TypeError, match=msg):
                idx.fillna([idx[0]])

            idx = indices.copy(deep=True)
            values = np.asarray(idx.values)

            if isinstance(indices, DatetimeIndexOpsMixin):
                values[1] = iNaT
            elif isinstance(indices, (Int64Index, UInt64Index)):
                return
            else:
                values[1] = np.nan

            if isinstance(indices, PeriodIndex):
                idx = type(indices)(values, freq=indices.freq)
            else:
                idx = type(indices)(values)

            expected = np.array([False] * len(idx), dtype=bool)
            expected[1] = True
            tm.assert_numpy_array_equal(idx._isnan, expected)
            assert idx.hasnans is True

    def test_nulls(self, indices):
        # this is really a smoke test for the methods
        # as these are adequately tested for function elsewhere
        if len(indices) == 0:
            tm.assert_numpy_array_equal(indices.isna(), np.array([], dtype=bool))
        elif isinstance(indices, MultiIndex):
            idx = indices.copy()
            msg = "isna is not defined for MultiIndex"
            with pytest.raises(NotImplementedError, match=msg):
                idx.isna()
        elif not indices.hasnans:
            tm.assert_numpy_array_equal(
                indices.isna(), np.zeros(len(indices), dtype=bool)
            )
            tm.assert_numpy_array_equal(
                indices.notna(), np.ones(len(indices), dtype=bool)
            )
        else:
            result = isna(indices)
            tm.assert_numpy_array_equal(indices.isna(), result)
            tm.assert_numpy_array_equal(indices.notna(), ~result)

    def test_empty(self):
        # GH 15270
        index = self.create_index()
        assert not index.empty
        assert index[:0].empty

    def test_join_self_unique(self, join_type):
        index = self.create_index()
        if index.is_unique:
            joined = index.join(index, how=join_type)
            assert (index == joined).all()

    def test_map(self):
        # callable
        index = self.create_index()

        # we don't infer UInt64
        if isinstance(index, pd.UInt64Index):
            expected = index.astype("int64")
        else:
            expected = index

        result = index.map(lambda x: x)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize(
        "mapper",
        [
            lambda values, index: {i: e for e, i in zip(values, index)},
            lambda values, index: pd.Series(values, index),
        ],
    )
    def test_map_dictlike(self, mapper):

        index = self.create_index()
        if isinstance(index, (pd.CategoricalIndex, pd.IntervalIndex)):
            pytest.skip(f"skipping tests for {type(index)}")

        identity = mapper(index.values, index)

        # we don't infer to UInt64 for a dict
        if isinstance(index, pd.UInt64Index) and isinstance(identity, dict):
            expected = index.astype("int64")
        else:
            expected = index

        result = index.map(identity)
        tm.assert_index_equal(result, expected)

        # empty mappable
        expected = pd.Index([np.nan] * len(index))
        result = index.map(mapper(expected, index))
        tm.assert_index_equal(result, expected)

    def test_map_str(self):
        # GH 31202
        index = self.create_index()
        result = index.map(str)
        expected = Index([str(x) for x in index], dtype=object)
        tm.assert_index_equal(result, expected)

    def test_putmask_with_wrong_mask(self):
        # GH18368
        index = self.create_index()

        with pytest.raises(ValueError):
            index.putmask(np.ones(len(index) + 1, np.bool), 1)

        with pytest.raises(ValueError):
            index.putmask(np.ones(len(index) - 1, np.bool), 1)

        with pytest.raises(ValueError):
            index.putmask("foo", 1)

    @pytest.mark.parametrize("copy", [True, False])
    @pytest.mark.parametrize("name", [None, "foo"])
    @pytest.mark.parametrize("ordered", [True, False])
    def test_astype_category(self, copy, name, ordered):
        # GH 18630
        index = self.create_index()
        if name:
            index = index.rename(name)

        # standard categories
        dtype = CategoricalDtype(ordered=ordered)
        result = index.astype(dtype, copy=copy)
        expected = CategoricalIndex(index.values, name=name, ordered=ordered)
        tm.assert_index_equal(result, expected)

        # non-standard categories
        dtype = CategoricalDtype(index.unique().tolist()[:-1], ordered)
        result = index.astype(dtype, copy=copy)
        expected = CategoricalIndex(index.values, name=name, dtype=dtype)
        tm.assert_index_equal(result, expected)

        if ordered is False:
            # dtype='category' defaults to ordered=False, so only test once
            result = index.astype("category", copy=copy)
            expected = CategoricalIndex(index.values, name=name)
            tm.assert_index_equal(result, expected)

    def test_is_unique(self):
        # initialize a unique index
        index = self.create_index().drop_duplicates()
        assert index.is_unique is True

        # empty index should be unique
        index_empty = index[:0]
        assert index_empty.is_unique is True

        # test basic dupes
        index_dup = index.insert(0, index[0])
        assert index_dup.is_unique is False

        # single NA should be unique
        index_na = index.insert(0, np.nan)
        assert index_na.is_unique is True

        # multiple NA should not be unique
        index_na_dup = index_na.insert(0, np.nan)
        assert index_na_dup.is_unique is False

    def test_engine_reference_cycle(self):
        # GH27585
        index = self.create_index()
        nrefs_pre = len(gc.get_referrers(index))
        index._engine
        assert len(gc.get_referrers(index)) == nrefs_pre

    def test_getitem_2d_deprecated(self):
        # GH#30588
        idx = self.create_index()
        with tm.assert_produces_warning(DeprecationWarning, check_stacklevel=False):
            res = idx[:, None]

        assert isinstance(res, np.ndarray), type(res)

    def test_contains_requires_hashable_raises(self):
        idx = self.create_index()
        with pytest.raises(TypeError, match="unhashable type"):
            [] in idx

        with pytest.raises(TypeError):
            {} in idx._engine
