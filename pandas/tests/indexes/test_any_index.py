"""
Tests that can be parametrized over _any_ Index object.

TODO: consider using hypothesis for these.
"""
import pytest

import pandas._testing as tm


def test_boolean_context_compat(indices):
    with pytest.raises(ValueError, match="The truth value of a"):
        if indices:
            pass


def test_sort(indices):
    msg = "cannot sort an Index object in-place, use sort_values instead"
    with pytest.raises(TypeError, match=msg):
        indices.sort()


def test_hash_error(indices):
    index = indices
    with pytest.raises(TypeError, match=f"unhashable type: '{type(index).__name__}'"):
        hash(indices)


def test_mutability(indices):
    if not len(indices):
        return
    msg = "Index does not support mutable operations"
    with pytest.raises(TypeError, match=msg):
        indices[0] = indices[0]


def test_wrong_number_names(indices):
    names = indices.nlevels * ["apple", "banana", "carrot"]
    with pytest.raises(ValueError, match="^Length"):
        indices.names = names


class TestConversion:
    def test_to_series(self, indices):
        # assert that we are creating a copy of the index

        ser = indices.to_series()
        assert ser.values is not indices.values
        assert ser.index is not indices
        assert ser.name == indices.name

    def test_to_series_with_arguments(self, indices):
        # GH#18699

        # index kwarg
        ser = indices.to_series(index=indices)

        assert ser.values is not indices.values
        assert ser.index is indices
        assert ser.name == indices.name

        # name kwarg
        ser = indices.to_series(name="__test")

        assert ser.values is not indices.values
        assert ser.index is not indices
        assert ser.name != indices.name

    def test_tolist_matches_list(self, indices):
        assert indices.tolist() == list(indices)


class TestRoundTrips:
    def test_pickle_roundtrip(self, indices):
        result = tm.round_trip_pickle(indices)
        tm.assert_index_equal(result, indices)
        if result.nlevels > 1:
            # GH#8367 round-trip with timezone
            assert indices.equal_levels(result)


class TestIndexing:
    def test_slice_keeps_name(self, indices):
        assert indices.name == indices[1:].name


class TestRendering:
    def test_str(self, indices):
        # test the string repr
        indices.name = "foo"
        assert "'foo'" in str(indices)
        assert type(indices).__name__ in str(indices)
