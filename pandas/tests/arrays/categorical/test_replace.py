import re

import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize(
    "to_replace,value,expected",
    [
        # one-to-one
        (4, 1, [1, 2, 3]),
        (3, 1, [1, 2, 1]),
        # many-to-one
        ((5, 6), 2, [1, 2, 3]),
        ((3, 2), 1, [1, 1, 1]),
    ],
)
def test_replace_categorical_series(to_replace, value, expected):
    # GH 31720
    ser = pd.Series([1, 2, 3], dtype="category")
    result = ser.replace(to_replace, value)
    expected = pd.Series(pd.Categorical(expected, categories=[1, 2, 3]))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "to_replace,value",
    [
        # one-to-one
        (3, 5),
        # many-to-one
        ((3, 2), 5),
    ],
)
def test_replace_categorical_series_new_category_raises(to_replace, value):
    # GH 31720
    ser = pd.Series([1, 2, 3], dtype="category")
    with pytest.raises(
        TypeError, match="Cannot setitem on a Categorical with a new category"
    ):
        ser.replace(to_replace, value)


def test_replace_maintain_ordering():
    # GH51016
    dtype = pd.CategoricalDtype([0, 1, 2], ordered=True)
    ser = pd.Series([0, 1, 2], dtype=dtype)
    result = ser.replace(0, 2)
    expected = pd.Series([2, 1, 2], dtype=dtype)
    tm.assert_series_equal(expected, result, check_category_order=True)


def test_replace_categorical_ea_dtype():
    # GH49404
    cat = pd.Categorical(pd.array(["a", "b", "c"], dtype="string"))
    result = pd.Series(cat).replace(["a", "b"], ["c", "c"])._values
    expected = pd.Categorical(
        pd.array(["c"] * 3, dtype="string"),
        categories=pd.array(["a", "b", "c"], dtype="string"),
    )
    tm.assert_categorical_equal(result, expected)


def test_replace_categorical_ea_dtype_different_cats_raises():
    # GH49404
    cat = pd.Categorical(pd.array(["a", "b"], dtype="string"))
    with pytest.raises(
        TypeError, match="Cannot setitem on a Categorical with a new category"
    ):
        pd.Series(cat).replace(["a", "b"], ["c", pd.NA])


@pytest.mark.parametrize(
    "kwargs",
    [
        {"regex": {"^a": "b"}},
        {"to_replace": "^a", "value": "b", "regex": True},
        {"to_replace": ["^a"], "value": ["b"], "regex": True},
        {"to_replace": re.compile("^a"), "value": "b"},
    ],
)
def test_replace_regex_existing_category(kwargs):
    # GH#38447 the dict, list and compiled spellings used to discard the replacement
    ser = pd.Series(pd.Categorical(["a", "b", "c"]))
    result = ser.replace(**kwargs)
    expected = pd.Series(pd.Categorical(["b", "b", "c"], categories=["a", "b", "c"]))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"regex": {"^a": "b"}},
        {"to_replace": "^a", "value": "b", "regex": True},
        {"to_replace": re.compile("^a"), "value": "b"},
    ],
)
def test_replace_regex_existing_category_arrow_strings(kwargs):
    # GH#69026 the gate resolves through _regex_target_dtype, so a Categorical of
    #  arrow strings stops being a silent no-op like every other Categorical
    pa = pytest.importorskip("pyarrow")
    cats = pd.array(["a", "b", "c"], dtype=pd.ArrowDtype(pa.string()))
    ser = pd.Series(pd.Categorical(cats))

    result = ser.replace(**kwargs)

    expected = pd.Series(pd.Categorical(cats.take([1, 1, 2]), categories=cats))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"regex": {"^a": "z"}},
        {"to_replace": "^a", "value": "z", "regex": True},
        {"to_replace": ["^a"], "value": ["z"], "regex": True},
        {"to_replace": re.compile("^a"), "value": "z"},
    ],
)
def test_replace_regex_new_category_raises(kwargs):
    # GH#38447 matches the non-regex spelling, which refuses to widen
    ser = pd.Series(pd.Categorical(["a", "b", "c"]))
    with pytest.raises(
        TypeError, match="Cannot setitem on a Categorical with a new category"
    ):
        ser.replace(**kwargs)


def test_replace_regex_non_string_categories():
    # GH#38447 the pattern must not be matched against a non-string category's repr
    ser = pd.Series(pd.Categorical([1, 2, 3]))
    expected = ser.copy()
    tm.assert_series_equal(ser.replace(regex={"^1": "9"}), expected)
    tm.assert_series_equal(ser.replace(re.compile("^1"), "9"), expected)
    tm.assert_series_equal(ser, expected)


def test_replace_regex_frame():
    # GH#38447
    df = pd.DataFrame({"A": pd.Categorical(["a", "b", "c"])})
    result = df.replace(regex={"^a": "b"})
    expected = pd.DataFrame(
        {"A": pd.Categorical(["b", "b", "c"], categories=["a", "b", "c"])}
    )
    tm.assert_frame_equal(result, expected)
