"""
Tests that can be parametrized over _any_ Index object.

TODO: consider using hypothesis for these.
"""
import pytest


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
    with pytest.raises(ValueError, match="^Length"):
        indices.names = ["apple", "banana", "carrot"]


def test_tolist_matches_list(indices):
    assert indices.tolist() == list(indices)


def test_map_with_invalid_na_action_raises(indices):
    msg = "na_action must either be 'ignore' or None"
    with pytest.raises(ValueError, match=msg):
        indices.map(lambda x: x, na_action="____")
