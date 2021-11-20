import numpy as np
import pytest

from pandas import DataFrame
import pandas._testing as tm
from pandas.core.reshape.reshape import from_dummies


@pytest.fixture
def dummies_basic():
    return DataFrame(
        {
            "col1_a": [1, 0, 1],
            "col1_b": [0, 1, 0],
            "col2_a": [0, 1, 0],
            "col2_b": [1, 0, 0],
            "col2_c": [0, 0, 1],
        },
    )


@pytest.fixture
def dummies_with_unassigned():
    return DataFrame(
        {
            "col1_a": [1, 0, 0],
            "col1_b": [0, 1, 0],
            "col2_a": [0, 1, 0],
            "col2_b": [0, 0, 0],
            "col2_c": [0, 0, 1],
        },
    )


def test_from_dummies_no_prefix_string_cats_basic():
    dummies = DataFrame({"a": [1, 0, 0, 1], "b": [0, 1, 0, 0], "c": [0, 0, 1, 0]})
    expected = DataFrame({"": ["a", "b", "c", "a"]})
    result = from_dummies(dummies)
    tm.assert_frame_equal(result, expected)


def test_from_dummies_no_prefix_int_cats_basic():
    dummies = DataFrame(
        {1: [1, 0, 0, 0], 25: [0, 1, 0, 0], 2: [0, 0, 1, 0], 5: [0, 0, 0, 1]}
    )
    expected = DataFrame({"": [1, 25, 2, 5]}, dtype="object")
    result = from_dummies(dummies)
    tm.assert_frame_equal(result, expected)


def test_from_dummies_no_prefix_float_cats_basic():
    dummies = DataFrame(
        {1.0: [1, 0, 0, 0], 25.0: [0, 1, 0, 0], 2.5: [0, 0, 1, 0], 5.84: [0, 0, 0, 1]}
    )
    expected = DataFrame({"": [1.0, 25.0, 2.5, 5.84]}, dtype="object")
    result = from_dummies(dummies)
    tm.assert_frame_equal(result, expected)


def test_from_dummies_no_prefix_mixed_cats_basic():
    dummies = DataFrame(
        {
            1.23: [1, 0, 0, 0, 0],
            "c": [0, 1, 0, 0, 0],
            2: [0, 0, 1, 0, 0],
            False: [0, 0, 0, 1, 0],
            None: [0, 0, 0, 0, 1],
        }
    )
    expected = DataFrame({"": [1.23, "c", 2, False, None]}, dtype="object")
    result = from_dummies(dummies)
    tm.assert_frame_equal(result, expected)


def test_from_dummies_no_prefix_string_cats_contains_get_dummies_NaN_column():
    dummies = DataFrame({"a": [1, 0, 0], "b": [0, 1, 0], "NaN": [0, 0, 1]})
    expected = DataFrame({"": ["a", "b", "NaN"]})
    result = from_dummies(dummies)
    tm.assert_frame_equal(result, expected)


def test_from_dummies_no_prefix_contains_unassigned():
    dummies = DataFrame({"a": [1, 0, 0], "b": [0, 1, 0]})
    with pytest.raises(
        ValueError,
        match=(
            r"Dummy DataFrame contains unassigned value\(s\); "
            r"First instance in row: 2"
        ),
    ):
        from_dummies(dummies)


def test_from_dummies_no_prefix_string_cats_dropped_first():
    dummies = DataFrame({"a": [1, 0, 0], "b": [0, 1, 0]})
    expected = DataFrame({"": ["a", "b", "c"]})
    result = from_dummies(dummies, dropped_first="c")
    tm.assert_frame_equal(result, expected)


def test_from_dummies_no_prefix_wrong_dropped_first_type():
    dummies = DataFrame({"a": [1, 0, 1], "b": [0, 1, 1]})
    with pytest.raises(
        TypeError,
        match=(
            r"Expected 'dropped_first' to be of type 'Hashable' or 'dict'; "
            r"Received 'dropped_first' of type: list"
        ),
    ):
        from_dummies(dummies, dropped_first=["c", "d"])


def test_from_dummies_no_prefix_multi_assignment():
    dummies = DataFrame({"a": [1, 0, 1], "b": [0, 1, 1]})
    with pytest.raises(
        ValueError,
        match=(
            r"Dummy DataFrame contains multi-assignment\(s\); "
            r"First instance in row: 2"
        ),
    ):
        from_dummies(dummies)


def test_from_dummies_no_prefix_contains_nan():
    dummies = DataFrame({"a": [1, 0, 0], "b": [0, 1, np.nan]})
    with pytest.raises(
        ValueError, match=r"Dummy DataFrame contains NA value in column: 'b'"
    ):
        from_dummies(dummies)


def test_from_dummies_no_dummies():
    dummies = DataFrame(
        {"a": [1, 6, 3, 1], "b": [0, 1, 0, 2], "c": ["c1", "c2", "c3", "c4"]}
    )
    with pytest.raises(
        TypeError,
        match=r"Passed DataFrame contains non-dummy data",
    ):
        from_dummies(dummies, sep="_")


def test_from_dummies_with_prefix_basic(dummies_basic):
    expected = DataFrame({"col1": ["a", "b", "a"], "col2": ["b", "a", "c"]})
    result = from_dummies(dummies_basic, sep="_")
    tm.assert_frame_equal(result, expected)


def test_from_dummies_with_prefix_prefix_multiple_seperators():
    dummies = DataFrame(
        {
            "col1_a": [1, 0, 1],
            "col1_b": [0, 1, 0],
            "col2-a": [0, 1, 0],
            "col2-b": [1, 0, 1],
        },
    )
    with pytest.raises(
        ValueError,
        match=(
            r"Separator not specified for all columns; "
            r"First instance column: col2-a"
        ),
    ):
        from_dummies(dummies, sep="_")


def test_from_dummies_with_prefix_prefix_sep_dict():
    dummies = DataFrame(
        {
            "col1_a-a": [1, 0, 1],
            "col1_b-b": [0, 1, 0],
            "col_2-a": [0, 1, 0],
            "col_2-b": [1, 0, 0],
            "col_2-c": [0, 0, 1],
        },
    )
    expected = DataFrame({"col1": ["a-a", "b-b", "a-a"], "col_2": ["b", "a", "c"]})
    result = from_dummies(
        dummies,
        sep={
            "col1": "_",
            "col_2": "-",
        },
    )
    tm.assert_frame_equal(result, expected)


def test_from_dummies_with_prefix_prefix_partial_sep_dict():
    dummies = DataFrame(
        {
            "col1_a-a": [1, 0, 1],
            "col1_b-b": [0, 1, 0],
            "col2-a_a": [0, 1, 0],
            "col2-b_b": [1, 0, 1],
        },
    )
    with pytest.raises(
        ValueError,
        match=(
            r"Separator not specified for all columns; "
            r"First instance column: col2-a_a"
        ),
    ):
        from_dummies(dummies, sep={"col1": "_"})


def test_from_dummies_with_prefix_sep_wrong_type(dummies_basic):

    with pytest.raises(
        TypeError,
        match=(
            r"Expected 'sep' to be of type 'str' or 'dict'; "
            r"Received 'sep' of type: list"
        ),
    ):
        from_dummies(dummies_basic, sep=["_"])


def test_from_dummies_with_prefix_contains_get_dummies_NaN_column():
    dummies = DataFrame(
        {
            "col1_a": [1, 0, 0],
            "col1_b": [0, 1, 0],
            "col1_NaN": [0, 0, 1],
            "col2_a": [0, 1, 0],
            "col2_b": [0, 0, 0],
            "col2_c": [0, 0, 1],
            "col2_NaN": [1, 0, 0],
        },
    )
    expected = DataFrame({"col1": ["a", "b", "NaN"], "col2": ["NaN", "a", "c"]})
    result = from_dummies(dummies, sep="_")
    tm.assert_frame_equal(result, expected)


def test_from_dummies_with_prefix_contains_unassigned(dummies_with_unassigned):
    with pytest.raises(
        ValueError,
        match=(
            r"Dummy DataFrame contains unassigned value\(s\); "
            r"First instance in row: 2"
        ),
    ):
        from_dummies(dummies_with_unassigned, sep="_")


def test_from_dummies_with_prefix_subset(dummies_basic):
    expected = DataFrame({"col1": ["a", "b", "a"], "col2": ["b", "a", "c"]})
    result = from_dummies(
        dummies_basic,
        subset=["col1_a", "col1_b", "col2_a", "col2_b", "col2_c"],
        sep="_",
    )
    tm.assert_frame_equal(result, expected)


def test_from_dummies_with_prefix_dropped_first_str(dummies_with_unassigned):
    expected = DataFrame({"col1": ["a", "b", "x"], "col2": ["x", "a", "c"]})
    result = from_dummies(dummies_with_unassigned, sep="_", dropped_first="x")
    tm.assert_frame_equal(result, expected)


def test_from_dummies_with_prefix_dropped_first_wrong_type(dummies_with_unassigned):
    with pytest.raises(
        TypeError,
        match=(
            r"Expected 'dropped_first' to be of type 'Hashable' or 'dict'; "
            r"Received 'dropped_first' of type: list"
        ),
    ):
        from_dummies(dummies_with_unassigned, sep="_", dropped_first=["x", "y"])


def test_from_dummies_with_prefix_dropped_first_dict(dummies_with_unassigned):
    expected = DataFrame({"col1": ["a", "b", "y"], "col2": ["x", "a", "c"]})
    result = from_dummies(
        dummies_with_unassigned,
        sep="_",
        dropped_first={"col2": "x", "col1": "y"},
    )
    tm.assert_frame_equal(result, expected)


def test_from_dummies_with_prefix_dropped_first_int_and_float(dummies_with_unassigned):
    expected = DataFrame({"col1": ["a", "b", 2.5], "col2": [1, "a", "c"]})
    result = from_dummies(
        dummies_with_unassigned,
        sep="_",
        dropped_first={"col2": 1, "col1": 2.5},
    )
    tm.assert_frame_equal(result, expected)


def test_from_dummies_with_prefix_dropped_first_bool_and_none(dummies_with_unassigned):
    expected = DataFrame({"col1": ["a", "b", False], "col2": [None, "a", "c"]})
    result = from_dummies(
        dummies_with_unassigned,
        sep="_",
        dropped_first={"col2": None, "col1": False},
    )
    tm.assert_frame_equal(result, expected)


def test_from_dummies_with_prefix_dropped_first_dict_not_complete(
    dummies_with_unassigned,
):
    with pytest.raises(
        ValueError,
        match=(
            r"Length of 'dropped_first' \(1\) did not match "
            r"the length of the columns being encoded \(2\)"
        ),
    ):
        from_dummies(dummies_with_unassigned, sep="_", dropped_first={"col1": "x"})


def test_from_dummies_with_prefix_wrong_subset_type(dummies_basic):
    with pytest.raises(
        TypeError,
        match=(
            r"Expected 'subset' to be of type 'Index', or 'list'; "
            r"Received 'subset' of type: str"
        ),
    ):
        from_dummies(
            dummies_basic,
            subset="col1_a",
            sep="_",
        )


def test_from_dummies_with_prefix_contains_nan(dummies_basic):
    dummies_basic["col2_c"][2] = np.nan
    with pytest.raises(
        ValueError, match=r"Dummy DataFrame contains NA value in column: 'col2_c'"
    ):
        from_dummies(dummies_basic, sep="_")


def test_from_dummies_with_prefix_double_assignment():
    dummies = DataFrame(
        {
            "col1_a": [1, 0, 1],
            "col1_b": [1, 1, 0],
            "col2_a": [0, 1, 0],
            "col2_b": [1, 0, 0],
            "col2_c": [0, 0, 1],
        },
    )
    with pytest.raises(
        ValueError,
        match=(
            r"Dummy DataFrame contains multi-assignment\(s\); "
            r"First instance in row: 0"
        ),
    ):
        from_dummies(dummies, sep="_")


def test_from_dummies_collate_prefix_sep_and_dropped_first_dict():
    dummies = DataFrame(
        {
            "col1_a-a": [1, 0, 0],
            "col1_b-b": [0, 1, 0],
            "col2-a_a": [0, 1, 0],
            "col2-b_b": [0, 0, 0],
            "col2-c_c": [0, 0, 1],
        },
    )
    expected = DataFrame({"col1": ["a-a", "b-b", "x"], "col2": ["y", "a_a", "c_c"]})
    result = from_dummies(
        dummies,
        sep={
            "col1": "_",
            "col2": "-",
        },
        dropped_first={"col1": "x", "col2": "y"},
    )
    tm.assert_frame_equal(result, expected)
