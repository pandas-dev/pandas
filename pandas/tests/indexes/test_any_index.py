"""
Tests that can be parametrized over _any_ Index object.
"""

import re

import numpy as np
import pytest

from pandas.errors import InvalidIndexError

import pandas as pd
import pandas._testing as tm


def test_boolean_context_compat(index):
    # GH#7897
    with pytest.raises(ValueError, match="The truth value of a"):
        if index:
            pass

    with pytest.raises(ValueError, match="The truth value of a"):
        bool(index)


def test_hash_error(index):
    with pytest.raises(TypeError, match=f"unhashable type: '{type(index).__name__}'"):
        hash(index)


def test_mutability(index):
    if not len(index):
        pytest.skip("Test doesn't make sense for empty index")
    msg = "Index does not support mutable operations"
    with pytest.raises(TypeError, match=msg):
        index[0] = index[0]


def test_map_identity_mapping(index, request, using_infer_string):
    # GH#12766
    if (
        not using_infer_string
        and isinstance(index.dtype, pd.StringDtype)
        and index.dtype.storage == "python"
    ):
        mark = pytest.mark.xfail(reason="Does not preserve dtype")
        request.applymarker(mark)

    result = index.map(lambda x: x)
    if index.dtype == object and (result.dtype in (bool, "string")):
        assert (index == result).all()
        # TODO: could work that into the 'exact="equiv"'?
        return  # FIXME: doesn't belong in this file anymore!
    tm.assert_index_equal(result, index, exact="equiv", check_freq=False)


@pytest.fixture
def two_value_index(index_flat):
    # CategoricalIndex differs, see categorical/test_map.py
    if not len(index_flat) or isinstance(index_flat, pd.CategoricalIndex):
        pytest.skip("Test doesn't apply to empty or categorical index")
    return index_flat.unique()[:2].rename("n")


def _xfail_if_tuples_become_arrow_lists(index, request):
    # ExtensionArray.map casts tuples to pyarrow lists since GH#64801
    if index.dtype == "string[pyarrow]":
        request.applymarker(pytest.mark.xfail(reason="tuples cast to pyarrow lists"))


def test_map_to_tuples_names_each_level(two_value_index, request):
    # GH#24800
    _xfail_if_tuples_become_arrow_lists(two_value_index, request)
    first, second = two_value_index.tolist()
    result = two_value_index.map(lambda val: (val, val))
    expected = pd.MultiIndex.from_tuples(
        [(first, first), (second, second)], names=["n", "n"]
    )
    tm.assert_index_equal(result, expected)


def test_map_to_tuples_of_different_lengths(two_value_index, request):
    # GH#50293
    _xfail_if_tuples_become_arrow_lists(two_value_index, request)
    first, second = two_value_index.tolist()
    result = two_value_index.map(lambda val: (val,) if val == first else (val, val))
    expected = pd.Index([(first,), (second, second)], name="n", tupleize_cols=False)
    tm.assert_index_equal(result, expected)


def test_map_to_tuple_or_non_tuple(two_value_index):
    # GH#24800
    first, _ = two_value_index.tolist()
    result = two_value_index.map(lambda val: (val, val) if val == first else "b")
    expected = pd.Index([(first, first), "b"], dtype=object, name="n")
    tm.assert_index_equal(result, expected)


@pytest.mark.parametrize("na_value", [None, np.nan])
def test_map_to_tuple_or_na(two_value_index, na_value, request):
    # GH#24800
    _xfail_if_tuples_become_arrow_lists(two_value_index, request)
    first, _ = two_value_index.tolist()
    result = two_value_index.map(lambda val: (val, val) if val == first else na_value)
    expected = pd.MultiIndex.from_tuples(
        [(first, first), (np.nan, np.nan)], names=["n", "n"]
    )
    tm.assert_index_equal(result, expected)


def test_map_dict_to_tuples_missing_key(two_value_index, request):
    # GH#24800
    _xfail_if_tuples_become_arrow_lists(two_value_index, request)
    first, _ = two_value_index.tolist()
    result = two_value_index.map({first: (first, first)})
    expected = pd.MultiIndex.from_tuples(
        [(first, first), (np.nan, np.nan)], names=["n", "n"]
    )
    tm.assert_index_equal(result, expected)


def test_wrong_number_names(index):
    names = index.nlevels * ["apple", "banana", "carrot"]
    with pytest.raises(ValueError, match="^Length"):
        index.names = names


def test_view_preserves_name(index):
    assert index.view().name == index.name


def test_ravel(index):
    # GH#19956 ravel returning ndarray is deprecated, in 2.0 returns a view on self
    res = index.ravel()
    tm.assert_index_equal(res, index)


class TestConversion:
    def test_to_series(self, index):
        # assert that we are creating a copy of the index

        ser = index.to_series()
        assert ser._values is not index._values
        assert ser.index is not index
        assert ser.name == index.name

    def test_to_series_with_arguments(self, index):
        # GH#18699

        # index kwarg
        ser = index.to_series(index=index)

        assert ser._values is not index._values
        assert ser.index is index
        assert ser.name == index.name

        # name kwarg
        ser = index.to_series(name="__test")

        assert ser._values is not index._values
        assert ser.index is not index
        assert ser.name != index.name

    def test_tolist_matches_list(self, index):
        assert index.tolist() == list(index)


class TestRoundTrips:
    def test_pickle_roundtrip(self, index, temp_file):
        result = tm.round_trip_pickle(index, temp_file)
        tm.assert_index_equal(result, index, exact=True)
        if result.nlevels > 1:
            # GH#8367 round-trip with timezone
            assert index.equal_levels(result)

    def test_pickle_preserves_name(self, index, temp_file):
        original_name, index.name = index.name, "foo"
        unpickled = tm.round_trip_pickle(index, temp_file)
        assert index.equals(unpickled)
        index.name = original_name


class TestIndexing:
    def test_getitem_0d_ndarray(self, index):
        # GH#55601
        if len(index) == 0:
            pytest.skip(reason="Test assumes non-empty index")
        key = np.array(0)
        result = index[key]

        assert result == index[0]

    def test_get_loc_listlike_raises_invalid_index_error(self, index):
        # and never TypeError
        key = np.array([0, 1], dtype=np.intp)

        with pytest.raises(InvalidIndexError, match=r"\[0 1\]"):
            index.get_loc(key)

        with pytest.raises(InvalidIndexError, match=r"\[False  True\]"):
            index.get_loc(key.astype(bool))

    def test_getitem_ellipsis(self, index):
        # GH#21282
        result = index[...]
        assert result.equals(index)
        assert result is not index

    def test_slice_keeps_name(self, index):
        assert index.name == index[1:].name

    @pytest.mark.parametrize("item", [101, "no_int", 2.5])
    def test_getitem_error(self, index, item):
        msg = "|".join(
            [
                r"index 101 is out of bounds for axis 0 with size [\d]+",
                re.escape(
                    "only integers, slices (`:`), ellipsis (`...`), "
                    "numpy.newaxis (`None`) and integer or boolean arrays "
                    "are valid indices"
                ),
                "index out of bounds",  # string[pyarrow]
            ]
        )
        with pytest.raises(IndexError, match=msg):
            index[item]


class TestRendering:
    def test_str(self, index):
        # test the string repr
        index.name = "foo"
        assert "'foo'" in str(index)
        assert type(index).__name__ in str(index)


class TestReductions:
    def test_argmax_axis_invalid(self, index):
        # GH#23081
        msg = r"`axis` must be less than the number of dimensions \(1\)"
        with pytest.raises(ValueError, match=msg):
            index.argmax(axis=1)
        with pytest.raises(ValueError, match=msg):
            index.argmin(axis=2)
        with pytest.raises(ValueError, match=msg):
            index.min(axis=-2)
        with pytest.raises(ValueError, match=msg):
            index.max(axis=-3)
