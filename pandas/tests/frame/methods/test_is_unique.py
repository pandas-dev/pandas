import pytest

from pandas import DataFrame, Series
import pandas._testing as tm


@pytest.mark.parametrize(
    "frame, expected",
    [
        # single column
        [DataFrame(), Series(dtype=bool)],
        [DataFrame({"a": ["x"]}), Series({"a": True})],
        [DataFrame({"a": ["x", "y"]}), Series({"a": True})],
        [DataFrame({"a": ["x", "x"]}), Series({"a": False})],
        [DataFrame({"a": ["x", "y", "y"]}), Series({"a": False})],
        # multiple columns
        [DataFrame(columns=["a", "b"]), Series({"a": True, "b": True})],
        [DataFrame({"a": ["x"], "b": ["y"]}), Series({"a": True, "b": True})],
        [
            DataFrame({"a": ["x", "y"], "b": ["x", "x"]}),
            Series({"a": True, "b": False}),
        ],
        # multiple columns, same column name
        [DataFrame(columns=["a", "a"]), Series([True, True], index=["a", "a"])],
        [
            DataFrame([["x", "y"]], columns=["a", "a"]),
            Series([True, True], index=["a", "a"]),
        ],
        [
            DataFrame([["x", "y"], ["y", "y"]], columns=["a", "a"]),
            Series([True, False], index=["a", "a"]),
        ],
    ],
)
def test_is_unique(frame, expected):
    # GH37565
    result = frame.is_unique()
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "frame, subset, expected",
    [
        [DataFrame(columns=["a", "b"]), ["a"], Series({"a": True})],
        [DataFrame({"a": ["x"], "b": ["y"]}), "a", Series({"a": True})],
        [DataFrame({"a": ["x"], "b": ["y"]}), ["a"], Series({"a": True})],
        [
            DataFrame({"a": ["x", "y"], "b": ["x", "x"]}),
            ["a", "b"],
            Series({"a": True, "b": False}),
        ],
    ],
)
def test_is_unique_subsetting(frame, subset, expected):
    # GH37565
    result = frame.is_unique(subset=subset)
    tm.assert_series_equal(result, expected)
