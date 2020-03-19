"""
Collection of tests asserting things that should be true for
any index subclass. Makes use of the `indices` fixture defined
in pandas/tests/indexes/conftest.py.
"""
import re

import numpy as np
import pytest

from pandas._libs.tslibs import iNaT

from pandas.core.dtypes.common import is_period_dtype, needs_i8_conversion

import pandas as pd
from pandas import CategoricalIndex, MultiIndex, RangeIndex
import pandas._testing as tm


class TestCommon:
    def test_droplevel(self, indices):
        # GH 21115
        if isinstance(indices, MultiIndex):
            # Tested separately in test_multi.py
            return

        assert indices.droplevel([]).equals(indices)

        for level in indices.name, [indices.name]:
            if isinstance(indices.name, tuple) and level is indices.name:
                # GH 21121 : droplevel with tuple name
                continue
            with pytest.raises(ValueError):
                indices.droplevel(level)

        for level in "wrong", ["wrong"]:
            with pytest.raises(
                KeyError,
                match=r"'Requested level \(wrong\) does not match index name \(None\)'",
            ):
                indices.droplevel(level)

    def test_constructor_non_hashable_name(self, indices):
        # GH 20527

        if isinstance(indices, MultiIndex):
            pytest.skip("multiindex handled in test_multi.py")

        message = "Index.name must be a hashable type"
        renamed = [["1"]]

        # With .rename()
        with pytest.raises(TypeError, match=message):
            indices.rename(name=renamed)

        # With .set_names()
        with pytest.raises(TypeError, match=message):
            indices.set_names(names=renamed)

    def test_constructor_unwraps_index(self, indices):
        if isinstance(indices, pd.MultiIndex):
            raise pytest.skip("MultiIndex has no ._data")
        a = indices
        b = type(a)(a)
        tm.assert_equal(a._data, b._data)

    @pytest.mark.parametrize("itm", [101, "no_int"])
    # FutureWarning from non-tuple sequence of nd indexing
    @pytest.mark.filterwarnings("ignore::FutureWarning")
    def test_getitem_error(self, indices, itm):
        with pytest.raises(IndexError):
            indices[itm]

    @pytest.mark.parametrize(
        "fname, sname, expected_name",
        [
            ("A", "A", "A"),
            ("A", "B", None),
            ("A", None, None),
            (None, "B", None),
            (None, None, None),
        ],
    )
    def test_corner_union(self, indices, fname, sname, expected_name):
        # GH 9943 9862
        # Test unions with various name combinations
        # Do not test MultiIndex or repeats

        if isinstance(indices, MultiIndex) or not indices.is_unique:
            pytest.skip("Not for MultiIndex or repeated indices")

        # Test copy.union(copy)
        first = indices.copy().set_names(fname)
        second = indices.copy().set_names(sname)
        union = first.union(second)
        expected = indices.copy().set_names(expected_name)
        tm.assert_index_equal(union, expected)

        # Test copy.union(empty)
        first = indices.copy().set_names(fname)
        second = indices.drop(indices).set_names(sname)
        union = first.union(second)
        expected = indices.copy().set_names(expected_name)
        tm.assert_index_equal(union, expected)

        # Test empty.union(copy)
        first = indices.drop(indices).set_names(fname)
        second = indices.copy().set_names(sname)
        union = first.union(second)
        expected = indices.copy().set_names(expected_name)
        tm.assert_index_equal(union, expected)

        # Test empty.union(empty)
        first = indices.drop(indices).set_names(fname)
        second = indices.drop(indices).set_names(sname)
        union = first.union(second)
        expected = indices.drop(indices).set_names(expected_name)
        tm.assert_index_equal(union, expected)

    def test_to_flat_index(self, indices):
        # 22866
        if isinstance(indices, MultiIndex):
            pytest.skip("Separate expectation for MultiIndex")

        result = indices.to_flat_index()
        tm.assert_index_equal(result, indices)

    def test_wrong_number_names(self, indices):
        with pytest.raises(ValueError, match="^Length"):
            indices.names = ["apple", "banana", "carrot"]

    def test_set_name_methods(self, indices):
        new_name = "This is the new name for this index"

        # don't tests a MultiIndex here (as its tested separated)
        if isinstance(indices, MultiIndex):
            pytest.skip("Skip check for MultiIndex")
        original_name = indices.name
        new_ind = indices.set_names([new_name])
        assert new_ind.name == new_name
        assert indices.name == original_name
        res = indices.rename(new_name, inplace=True)

        # should return None
        assert res is None
        assert indices.name == new_name
        assert indices.names == [new_name]
        # FIXME: dont leave commented-out
        # with pytest.raises(TypeError, match="list-like"):
        #    # should still fail even if it would be the right length
        #    ind.set_names("a")
        with pytest.raises(ValueError, match="Level must be None"):
            indices.set_names("a", level=0)

        # rename in place just leaves tuples and other containers alone
        name = ("A", "B")
        indices.rename(name, inplace=True)
        assert indices.name == name
        assert indices.names == [name]

    def test_copy_and_deepcopy(self, indices):
        from copy import copy, deepcopy

        if isinstance(indices, MultiIndex):
            pytest.skip("Skip check for MultiIndex")

        for func in (copy, deepcopy):
            idx_copy = func(indices)
            assert idx_copy is not indices
            assert idx_copy.equals(indices)

        new_copy = indices.copy(deep=True, name="banana")
        assert new_copy.name == "banana"

    def test_unique(self, indices):
        # don't test a MultiIndex here (as its tested separated)
        # don't test a CategoricalIndex because categories change (GH 18291)
        if isinstance(indices, (MultiIndex, CategoricalIndex)):
            pytest.skip("Skip check for MultiIndex/CategoricalIndex")

        # GH 17896
        expected = indices.drop_duplicates()
        for level in 0, indices.name, None:
            result = indices.unique(level=level)
            tm.assert_index_equal(result, expected)

        msg = "Too many levels: Index has only 1 level, not 4"
        with pytest.raises(IndexError, match=msg):
            indices.unique(level=3)

        msg = (
            fr"Requested level \(wrong\) does not match index name "
            fr"\({re.escape(indices.name.__repr__())}\)"
        )
        with pytest.raises(KeyError, match=msg):
            indices.unique(level="wrong")

    def test_get_unique_index(self, indices):
        # MultiIndex tested separately
        if not len(indices) or isinstance(indices, MultiIndex):
            pytest.skip("Skip check for empty Index and MultiIndex")

        idx = indices[[0] * 5]
        idx_unique = indices[[0]]

        # We test against `idx_unique`, so first we make sure it's unique
        # and doesn't contain nans.
        assert idx_unique.is_unique is True
        try:
            assert idx_unique.hasnans is False
        except NotImplementedError:
            pass

        for dropna in [False, True]:
            result = idx._get_unique_index(dropna=dropna)
            tm.assert_index_equal(result, idx_unique)

        # nans:
        if not indices._can_hold_na:
            pytest.skip("Skip na-check if index cannot hold na")

        if is_period_dtype(indices):
            vals = indices[[0] * 5]._data
            vals[0] = pd.NaT
        elif needs_i8_conversion(indices):
            vals = indices.asi8[[0] * 5]
            vals[0] = iNaT
        else:
            vals = indices.values[[0] * 5]
            vals[0] = np.nan

        vals_unique = vals[:2]
        idx_nan = indices._shallow_copy(vals)
        idx_unique_nan = indices._shallow_copy(vals_unique)
        assert idx_unique_nan.is_unique is True

        assert idx_nan.dtype == indices.dtype
        assert idx_unique_nan.dtype == indices.dtype

        for dropna, expected in zip([False, True], [idx_unique_nan, idx_unique]):
            for i in [idx_nan, idx_unique_nan]:
                result = i._get_unique_index(dropna=dropna)
                tm.assert_index_equal(result, expected)

    def test_mutability(self, indices):
        if not len(indices):
            pytest.skip("Skip check for empty Index")
        msg = "Index does not support mutable operations"
        with pytest.raises(TypeError, match=msg):
            indices[0] = indices[0]

    def test_view(self, indices):
        assert indices.view().name == indices.name

    def test_searchsorted_monotonic(self, indices):
        # GH17271
        # not implemented for tuple searches in MultiIndex
        # or Intervals searches in IntervalIndex
        if isinstance(indices, (MultiIndex, pd.IntervalIndex)):
            pytest.skip("Skip check for MultiIndex/IntervalIndex")

        # nothing to test if the index is empty
        if indices.empty:
            pytest.skip("Skip check for empty Index")
        value = indices[0]

        # determine the expected results (handle dupes for 'right')
        expected_left, expected_right = 0, (indices == value).argmin()
        if expected_right == 0:
            # all values are the same, expected_right should be length
            expected_right = len(indices)

        # test _searchsorted_monotonic in all cases
        # test searchsorted only for increasing
        if indices.is_monotonic_increasing:
            ssm_left = indices._searchsorted_monotonic(value, side="left")
            assert expected_left == ssm_left

            ssm_right = indices._searchsorted_monotonic(value, side="right")
            assert expected_right == ssm_right

            ss_left = indices.searchsorted(value, side="left")
            assert expected_left == ss_left

            ss_right = indices.searchsorted(value, side="right")
            assert expected_right == ss_right

        elif indices.is_monotonic_decreasing:
            ssm_left = indices._searchsorted_monotonic(value, side="left")
            assert expected_left == ssm_left

            ssm_right = indices._searchsorted_monotonic(value, side="right")
            assert expected_right == ssm_right
        else:
            # non-monotonic should raise.
            with pytest.raises(ValueError):
                indices._searchsorted_monotonic(value, side="left")

    def test_pickle(self, indices):
        original_name, indices.name = indices.name, "foo"
        unpickled = tm.round_trip_pickle(indices)
        assert indices.equals(unpickled)
        indices.name = original_name

    def test_drop_duplicates(self, indices, keep):
        if isinstance(indices, MultiIndex):
            pytest.skip("MultiIndex is tested separately")
        if isinstance(indices, RangeIndex):
            pytest.skip(
                "RangeIndex is tested in test_drop_duplicates_no_duplicates"
                " as it cannot hold duplicates"
            )
        if len(indices) == 0:
            pytest.skip(
                "empty index is tested in test_drop_duplicates_no_duplicates"
                " as it cannot hold duplicates"
            )

        # make unique index
        holder = type(indices)
        unique_values = list(set(indices))
        unique_idx = holder(unique_values)

        # make duplicated index
        n = len(unique_idx)
        duplicated_selection = np.random.choice(n, int(n * 1.5))
        idx = holder(unique_idx.values[duplicated_selection])

        # Series.duplicated is tested separately
        expected_duplicated = (
            pd.Series(duplicated_selection).duplicated(keep=keep).values
        )
        tm.assert_numpy_array_equal(idx.duplicated(keep=keep), expected_duplicated)

        # Series.drop_duplicates is tested separately
        expected_dropped = holder(pd.Series(idx).drop_duplicates(keep=keep))
        tm.assert_index_equal(idx.drop_duplicates(keep=keep), expected_dropped)

    def test_drop_duplicates_no_duplicates(self, indices):
        if isinstance(indices, MultiIndex):
            pytest.skip("MultiIndex is tested separately")

        # make unique index
        if isinstance(indices, RangeIndex):
            # RangeIndex cannot have duplicates
            unique_idx = indices
        else:
            holder = type(indices)
            unique_values = list(set(indices))
            unique_idx = holder(unique_values)

        # check on unique index
        expected_duplicated = np.array([False] * len(unique_idx), dtype="bool")
        tm.assert_numpy_array_equal(unique_idx.duplicated(), expected_duplicated)
        result_dropped = unique_idx.drop_duplicates()
        tm.assert_index_equal(result_dropped, unique_idx)
        # validate shallow copy
        assert result_dropped is not unique_idx

    def test_drop_duplicates_inplace(self, indices):
        msg = r"drop_duplicates\(\) got an unexpected keyword argument"
        with pytest.raises(TypeError, match=msg):
            indices.drop_duplicates(inplace=True)

    def test_has_duplicates(self, indices):
        holder = type(indices)
        if not len(indices) or isinstance(indices, (MultiIndex, RangeIndex)):
            # MultiIndex tested separately in:
            #   tests/indexes/multi/test_unique_and_duplicates.
            # RangeIndex is unique by definition.
            pytest.skip("Skip check for empty Index, MultiIndex, and RangeIndex")

        idx = holder([indices[0]] * 5)
        assert idx.is_unique is False
        assert idx.has_duplicates is True
