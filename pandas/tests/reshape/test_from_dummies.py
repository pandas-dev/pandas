import numpy as np
import pytest

from pandas import (
    DataFrame,
    Series,
)
from pandas.core.reshape.reshape import from_dummies


def test_from_dummies_to_series_basic():
    dummies = DataFrame({"a": [1, 0, 0, 1], "b": [0, 1, 0, 0], "c": [0, 0, 1, 0]})
    expected = Series(list("abca"))
    result = from_dummies(dummies, to_series=True)
    assert all(result == expected)


def test_from_dummies_to_series_dummy_na():
    dummies = DataFrame({"a": [1, 0, 0], "b": [0, 1, 0], "NaN": [0, 0, 1]})
    expected = Series(["a", "b", np.nan])
    result = from_dummies(dummies, to_series=True, dummy_na=True)
    assert all(result[:2] == expected[:2])
    assert result[2] is expected[2]


def test_from_dummies_to_series_contains_nan():
    dummies = DataFrame({"a": [1, 0, 0], "b": [0, 1, 0]})
    expected = Series(["a", "b", np.nan])
    result = from_dummies(dummies, to_series=True)
    assert all(result[:2] == expected[:2])
    assert result[2] is expected[2]


def test_from_dummies_to_series_dropped_first():
    dummies = DataFrame({"a": [1, 0, 0], "b": [0, 1, 0]})
    expected = Series(["a", "b", "c"])
    result = from_dummies(dummies, to_series=True, dropped_first="c")
    assert all(result == expected)


def test_from_dummies_to_series_wrong_dropped_first():
    dummies = DataFrame({"a": [1, 0, 1], "b": [0, 1, 1]})
    with pytest.raises(
        ValueError,
        match=r"Only one dropped first value possible in 1D dummy DataFrame.",
    ):
        from_dummies(dummies, to_series=True, dropped_first=["c", "d"])


def test_from_dummies_to_series_multi_assignment():
    dummies = DataFrame({"a": [1, 0, 1], "b": [0, 1, 1]})
    with pytest.raises(
        ValueError, match=r"Dummy DataFrame contains multi-assignment in row 2."
    ):
        from_dummies(dummies, to_series=True)


def test_from_dummies_to_series_unassigned_row():
    dummies = DataFrame({"a": [1, 0, 0], "b": [0, 1, 0]})
    with pytest.raises(
        ValueError, match=r"Dummy DataFrame contains no assignment in row 2."
    ):
        from_dummies(dummies, to_series=True, dummy_na=True)


def test_from_dummies_no_dummies():
    dummies = DataFrame(
        {"a": [1, 6, 3, 1], "b": [0, 1, 0, 2], "c": ["c1", "c2", "c3", "c4"]}
    )
    expected = dummies
    result = from_dummies(dummies)
    assert all(result == expected)


def test_from_dummies_to_df_basic():
    dummies = DataFrame(
        {
            "C": [1, 2, 3],
            "col1_a": [1, 0, 1],
            "col1_b": [0, 1, 0],
            "col2_a": [0, 1, 0],
            "col2_b": [1, 0, 0],
            "col2_c": [0, 0, 1],
        },
    )
    expected = DataFrame(
        {"C": [1, 2, 3], "col1": ["a", "b", "a"], "col2": ["b", "a", "c"]}
    )
    result = from_dummies(dummies)
    assert all(result == expected)


def test_from_dummies_to_df_variable_string():
    dummies = DataFrame(
        {
            "C": [1, 2, 3],
            "col1_a": [1, 0, 1],
            "col1_b": [0, 1, 0],
            "col2_a": [0, 1, 0],
            "col2_b": [1, 0, 0],
            "col2_c": [0, 0, 1],
        },
    )
    expected = DataFrame(
        {"C": [1, 2, 3], "varname0": ["a", "b", "a"], "varname1": ["b", "a", "c"]}
    )
    result = from_dummies(dummies, variables="varname")
    assert all(result == expected)


def test_from_dummies_to_df_variable_list():
    dummies = DataFrame(
        {
            "C": [1, 2, 3],
            "col1_a": [1, 0, 1],
            "col1_b": [0, 1, 0],
            "col2_a": [0, 1, 0],
            "col2_b": [1, 0, 0],
            "col2_c": [0, 0, 1],
        },
    )
    expected = DataFrame({"C": [1, 2, 3], "A": ["a", "b", "a"], "B": ["b", "a", "c"]})
    result = from_dummies(dummies, variables=["A", "B"])
    assert all(result == expected)


def test_from_dummies_to_df_variable_list_not_complete():
    dummies = DataFrame(
        {
            "C": [1, 2, 3],
            "col1_a": [1, 0, 1],
            "col1_b": [0, 1, 0],
            "col2_a": [0, 1, 0],
            "col2_b": [1, 0, 0],
            "col2_c": [0, 0, 1],
        },
    )
    with pytest.raises(
        ValueError,
        match=(
            r"Length of 'variables' \(1\) did not match "
            r"the length of the columns being encoded \(2\)."
        ),
    ):
        from_dummies(dummies, variables=["A"])


def test_from_dummies_to_df_variable_dict():
    dummies = DataFrame(
        {
            "C": [1, 2, 3],
            "col1_a": [1, 0, 1],
            "col1_b": [0, 1, 0],
            "col2_a": [0, 1, 0],
            "col2_b": [1, 0, 0],
            "col2_c": [0, 0, 1],
        },
    )
    expected = DataFrame({"C": [1, 2, 3], "A": ["b", "a", "c"], "B": ["a", "b", "a"]})
    result = from_dummies(dummies, variables={"col2": "A", "col1": "B"})
    assert all(result == expected)


def test_from_dummies_to_df_variable_dict_not_complete():
    dummies = DataFrame(
        {
            "C": [1, 2, 3],
            "col1_a": [1, 0, 1],
            "col1_b": [0, 1, 0],
            "col2_a": [0, 1, 0],
            "col2_b": [1, 0, 0],
            "col2_c": [0, 0, 1],
        },
    )
    with pytest.raises(
        ValueError,
        match=(
            r"Length of 'variables' \(1\) did not match "
            r"the length of the columns being encoded \(2\)."
        ),
    ):
        from_dummies(dummies, variables={"col1": "A"})


def test_from_dummies_to_df_prefix_sep_list():
    dummies = DataFrame(
        {
            "C": [1, 2, 3],
            "col1_a": [1, 0, 1],
            "col1_b": [0, 1, 0],
            "col2-a": [0, 1, 0],
            "col2-b": [1, 0, 0],
            "col2-c": [0, 0, 1],
        },
    )
    expected = DataFrame(
        {"C": [1, 2, 3], "col1": ["a", "b", "a"], "col2": ["b", "a", "c"]}
    )
    result = from_dummies(dummies, prefix_sep=["_", "-"])
    assert all(result == expected)


def test_from_dummies_to_df_prefix_sep_dict():
    dummies = DataFrame(
        {
            "C": [1, 2, 3],
            "col1_a-a": [1, 0, 1],
            "col1_b-a": [0, 1, 0],
            "col2-a_a": [0, 1, 0],
            "col2-b_b": [1, 0, 0],
            "col2-c_c": [0, 0, 1],
        },
    )
    expected = DataFrame(
        {"C": [1, 2, 3], "col1": ["a-a", "b-b", "a-a"], "col2": ["b_b", "a_a", "c_c"]}
    )
    result = from_dummies(
        dummies,
        prefix_sep={
            "col1": "_",
            "col2": "-",
        },
    )
    assert all(result == expected)


def test_from_dummies_to_df_dummy_na():
    dummies = DataFrame(
        {
            "C": [1, 2, 3],
            "col1_a": [1, 0, 0],
            "col1_b": [0, 1, 0],
            "col1_NaN": [0, 0, 1],
            "col2_a": [0, 1, 0],
            "col2_b": [0, 0, 0],
            "col2_c": [0, 0, 1],
            "col2_NAN": [1, 0, 0],
        },
    )
    expected = DataFrame(
        {"C": [1, 2, 3], "col1": ["a", "b", np.nan], "col2": [np.nan, "a", "c"]}
    )
    result = from_dummies(dummies, dummy_na=True)
    assert all(result["C"] == expected["C"])
    assert all(result["col1"][:2] == expected["col1"][:2])
    assert all(result["col2"][1:] == expected["col2"][1:])
    assert result["col1"][2] is expected["col1"][2]
    assert result["col2"][1] is expected["col2"][1]


def test_from_dummies_to_df_contains_nan():
    dummies = DataFrame(
        {
            "C": [1, 2, 3],
            "col1_a": [1, 0, 0],
            "col1_b": [0, 1, 0],
            "col2_a": [0, 1, 0],
            "col2_b": [0, 0, 0],
            "col2_c": [0, 0, 1],
        },
    )
    expected = DataFrame(
        {"C": [1, 2, 3], "col1": ["a", "b", np.nan], "col2": [np.nan, "a", "c"]}
    )
    result = from_dummies(dummies)
    assert all(result["C"] == expected["C"])
    assert all(result["col1"][:2] == expected["col1"][:2])
    assert all(result["col2"][1:] == expected["col2"][1:])
    assert result["col1"][2] is expected["col1"][2]
    assert result["col2"][1] is expected["col2"][1]


def test_from_dummies_to_df_columns():
    dummies = DataFrame(
        {
            "C": [1, 2, 3],
            "col1_a": [1, 0, 1],
            "col1_b": [0, 1, 0],
            "col2_a": [0, 1, 0],
            "col2_b": [1, 0, 0],
            "col2_c": [0, 0, 1],
        },
    )
    expected = DataFrame({"col1": ["a", "b", "a"], "col2": ["b", "a", "c"]})
    result = from_dummies(
        dummies, columns=["col1_a", "col1_b", "col2_a", "col2_b", "col2_c"]
    )
    assert all(result == expected)


def test_from_dummies_to_df_dropped_first_str():
    dummies = DataFrame(
        {
            "C": [1, 2, 3],
            "col1_a": [1, 0, 0],
            "col1_b": [0, 1, 0],
            "col2_a": [0, 1, 0],
            "col2_b": [0, 0, 0],
            "col2_c": [0, 0, 1],
        },
    )
    expected = DataFrame(
        {"C": [1, 2, 3], "col1": ["a", "b", "x"], "col2": ["x", "a", "c"]}
    )
    result = from_dummies(dummies, dropped_first="x")
    assert all(result == expected)


def test_from_dummies_to_df_dropped_first_list():
    dummies = DataFrame(
        {
            "C": [1, 2, 3],
            "col1_a": [1, 0, 0],
            "col1_b": [0, 1, 0],
            "col2_a": [0, 1, 0],
            "col2_b": [0, 0, 0],
            "col2_c": [0, 0, 1],
        },
    )
    expected = DataFrame(
        {"C": [1, 2, 3], "col1": ["a", "b", "x"], "col2": ["y", "a", "c"]}
    )
    result = from_dummies(dummies, dropped_first=["x", "y"])
    assert all(result == expected)


def test_from_dummies_to_df_dropped_first_list_not_complete():
    dummies = DataFrame(
        {
            "C": [1, 2, 3],
            "col1_a": [1, 0, 0],
            "col1_b": [0, 1, 0],
            "col2_a": [0, 1, 0],
            "col2_b": [0, 0, 0],
            "col2_c": [0, 0, 1],
        },
    )
    with pytest.raises(
        ValueError,
        match=(
            r"Length of 'dropped_first' \(1\) did not match "
            r"the length of the columns being encoded \(2\)."
        ),
    ):
        from_dummies(dummies, dropped_first=["x"])


def test_from_dummies_to_df_dropped_first_dict():
    dummies = DataFrame(
        {
            "C": [1, 2, 3],
            "col1_a": [1, 0, 0],
            "col1_b": [0, 1, 0],
            "col2_a": [0, 1, 0],
            "col2_b": [0, 0, 0],
            "col2_c": [0, 0, 1],
        },
    )
    expected = DataFrame(
        {"C": [1, 2, 3], "col1": ["a", "b", "y"], "col2": ["x", "a", "c"]}
    )
    result = from_dummies(dummies, dropped_first={"col2": "x", "col1": "y"})
    assert all(result == expected)


def test_from_dummies_to_df_dropped_first_dict_not_complete():
    dummies = DataFrame(
        {
            "C": [1, 2, 3],
            "col1_a": [1, 0, 0],
            "col1_b": [0, 1, 0],
            "col2_a": [0, 1, 0],
            "col2_b": [0, 0, 0],
            "col2_c": [0, 0, 1],
        },
    )
    with pytest.raises(
        ValueError,
        match=(
            r"Length of 'dropped_first' \(1\) did not match "
            r"the length of the columns being encoded \(2\)."
        ),
    ):
        from_dummies(dummies, dropped_first={"col1": "x"})


def test_from_dummies_to_df_wrong_column_type():
    dummies = DataFrame(
        {
            "C": [1, 2, 3],
            "col1_a": [1, 0, 0],
            "col1_b": [0, 1, 0],
            "col2_a": [0, 1, 0],
            "col2_b": [1, 0, 0],
            "col2_c": [0, 0, 1],
        },
    )
    with pytest.raises(
        TypeError,
        match=r"Input must be a list-like for parameter 'columns'",
    ):
        from_dummies(dummies, columns="col1_a")


def test_from_dummies_to_df_double_assignment():
    dummies = DataFrame(
        {
            "C": [1, 2, 3],
            "col1_a": [1, 0, 0],
            "col1_b": [1, 1, 0],
            "col2_a": [0, 1, 0],
            "col2_b": [1, 0, 0],
            "col2_c": [0, 0, 1],
        },
    )
    with pytest.raises(
        ValueError,
        match=(
            r"Dummy DataFrame contains multi-assignment\(s\) for prefix: "
            r"'col1' in row 0."
        ),
    ):
        from_dummies(dummies)


def test_from_dummies_to_df_no_assignment():
    dummies = DataFrame(
        {
            "C": [1, 2, 3],
            "col1_a": [0, 0, 0],
            "col1_b": [0, 1, 0],
            "col2_a": [0, 1, 0],
            "col2_b": [1, 0, 0],
            "col2_c": [0, 0, 1],
        },
    )
    with pytest.raises(
        ValueError,
        match=r"Dummy DataFrame contains no assignment for prefix: 'col1' in row 0.",
    ):
        from_dummies(dummies, dummy_na=True)
