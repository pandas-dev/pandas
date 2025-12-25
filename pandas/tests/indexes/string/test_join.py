import pandas as pd


def test_left_join_preserves_string_index_dtype_no_overlap():
    left = pd.Index(["1", "2", "3"], dtype="string")
    right = pd.Index([0, 1, 2])

    result = left.join(right, how="left")

    assert result.dtype == "string"
    assert result.equals(left)


def test_left_join_string_index_with_overlap_upcasts():
    left = pd.Index(["1", "2", "3"], dtype="string")
    right = pd.Index(["2", "4"])

    result = left.join(right, how="left")

    assert result.dtype == object


def test_str_cat_preserves_string_index_dtype_with_series():
    s = pd.Series(["a", "b", "c"], index=["1", "2", "3"])
    other = pd.Series(["A", "B", "C"])

    result = s.str.cat(other, sep=",", na_rep="-")

    assert result.index.dtype == "string"
