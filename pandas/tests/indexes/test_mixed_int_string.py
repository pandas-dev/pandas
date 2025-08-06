import pytest

import pandas as pd


def test_mixed_int_string_index():
    idx = pd.Index([0, "a", 1, "b", 2, "c"])

    # Check if the index is of type Index
    assert len(idx) == 6
    assert idx[1] == "a"
    assert idx[-1] == "c"

    # Check if the index is sorted (it should not be)
    with pytest.raises(TypeError):
        idx.sort_values()

    # Check if the index is unique
    assert idx.is_unique

    # Check if the index contains a specific value
    assert idx.get_loc("a") == 1
    with pytest.raises(KeyError):
        idx.get_loc("z")
