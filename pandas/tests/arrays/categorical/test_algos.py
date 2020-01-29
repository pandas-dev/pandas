import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize("ordered", [True, False])
@pytest.mark.parametrize("categories", [["b", "a", "c"], ["a", "b", "c", "d"]])
def test_factorize(categories, ordered):
    cat = pd.Categorical(
        ["b", "b", "a", "c", None], categories=categories, ordered=ordered
    )
    codes, uniques = pd.factorize(cat)
    expected_codes = np.array([0, 0, 1, 2, -1], dtype=np.intp)
    expected_uniques = pd.Categorical(
        ["b", "a", "c"], categories=categories, ordered=ordered
    )

    tm.assert_numpy_array_equal(codes, expected_codes)
    tm.assert_categorical_equal(uniques, expected_uniques)


def test_factorized_sort():
    cat = pd.Categorical(["b", "b", None, "a"])
    codes, uniques = pd.factorize(cat, sort=True)
    expected_codes = np.array([1, 1, -1, 0], dtype=np.intp)
    expected_uniques = pd.Categorical(["a", "b"])

    tm.assert_numpy_array_equal(codes, expected_codes)
    tm.assert_categorical_equal(uniques, expected_uniques)


def test_factorized_sort_ordered():
    cat = pd.Categorical(
        ["b", "b", None, "a"], categories=["c", "b", "a"], ordered=True
    )

    codes, uniques = pd.factorize(cat, sort=True)
    expected_codes = np.array([0, 0, -1, 1], dtype=np.intp)
    expected_uniques = pd.Categorical(
        ["b", "a"], categories=["c", "b", "a"], ordered=True
    )

    tm.assert_numpy_array_equal(codes, expected_codes)
    tm.assert_categorical_equal(uniques, expected_uniques)


def test_isin_cats():
    # GH2003
    cat = pd.Categorical(["a", "b", np.nan])

    result = cat.isin(["a", np.nan])
    expected = np.array([True, False, True], dtype=bool)
    tm.assert_numpy_array_equal(expected, result)

    result = cat.isin(["a", "c"])
    expected = np.array([True, False, False], dtype=bool)
    tm.assert_numpy_array_equal(expected, result)


@pytest.mark.parametrize(
    "to_replace, value, result, expected_error_msg",
    [
        ("b", "c", ["a", "c"], "Categorical.categories are different"),
        ("c", "d", ["a", "b"], None),
        ("b", None, ["a", None], "Categorical.categories length are different"),
    ],
)
def test_replace(to_replace, value, result, expected_error_msg):
    # GH 26988
    cat = pd.Categorical(["a", "b"])
    expected = pd.Categorical(result)
    result = cat.replace(to_replace, value)
    tm.assert_categorical_equal(result, expected)
    if to_replace == "b":  # the "c" test is supposed to be unchanged
        with pytest.raises(AssertionError, match=expected_error_msg):
            # ensure non-inplace call does not affect original
            tm.assert_categorical_equal(cat, expected)
    cat.replace(to_replace, value, inplace=True)
    tm.assert_categorical_equal(cat, expected)


@pytest.mark.parametrize("empty", [[], pd.Series(dtype=object), np.array([])])
def test_isin_empty(empty):
    s = pd.Categorical(["a", "b"])
    expected = np.array([False, False], dtype=bool)

    result = s.isin(empty)
    tm.assert_numpy_array_equal(expected, result)


def test_diff():
    s = pd.Series([1, 2, 3], dtype="category")
    with tm.assert_produces_warning(FutureWarning):
        result = s.diff()
    expected = pd.Series([np.nan, 1, 1])
    tm.assert_series_equal(result, expected)

    expected = expected.to_frame(name="A")
    df = s.to_frame(name="A")
    with tm.assert_produces_warning(FutureWarning):
        result = df.diff()

    tm.assert_frame_equal(result, expected)


class TestTake:
    # https://github.com/pandas-dev/pandas/issues/20664

    def test_take_default_allow_fill(self):
        cat = pd.Categorical(["a", "b"])
        with tm.assert_produces_warning(None):
            result = cat.take([0, -1])

        assert result.equals(cat)

    def test_take_positive_no_warning(self):
        cat = pd.Categorical(["a", "b"])
        with tm.assert_produces_warning(None):
            cat.take([0, 0])

    def test_take_bounds(self, allow_fill):
        # https://github.com/pandas-dev/pandas/issues/20664
        cat = pd.Categorical(["a", "b", "a"])
        if allow_fill:
            msg = "indices are out-of-bounds"
        else:
            msg = "index 4 is out of bounds for( axis 0 with)? size 3"
        with pytest.raises(IndexError, match=msg):
            cat.take([4, 5], allow_fill=allow_fill)

    def test_take_empty(self, allow_fill):
        # https://github.com/pandas-dev/pandas/issues/20664
        cat = pd.Categorical([], categories=["a", "b"])
        if allow_fill:
            msg = "indices are out-of-bounds"
        else:
            msg = "cannot do a non-empty take from an empty axes"
        with pytest.raises(IndexError, match=msg):
            cat.take([0], allow_fill=allow_fill)

    def test_positional_take(self, ordered_fixture):
        cat = pd.Categorical(
            ["a", "a", "b", "b"], categories=["b", "a"], ordered=ordered_fixture
        )
        result = cat.take([0, 1, 2], allow_fill=False)
        expected = pd.Categorical(
            ["a", "a", "b"], categories=cat.categories, ordered=ordered_fixture
        )
        tm.assert_categorical_equal(result, expected)

    def test_positional_take_unobserved(self, ordered_fixture):
        cat = pd.Categorical(
            ["a", "b"], categories=["a", "b", "c"], ordered=ordered_fixture
        )
        result = cat.take([1, 0], allow_fill=False)
        expected = pd.Categorical(
            ["b", "a"], categories=cat.categories, ordered=ordered_fixture
        )
        tm.assert_categorical_equal(result, expected)

    def test_take_allow_fill(self):
        # https://github.com/pandas-dev/pandas/issues/23296
        cat = pd.Categorical(["a", "a", "b"])
        result = cat.take([0, -1, -1], allow_fill=True)
        expected = pd.Categorical(["a", np.nan, np.nan], categories=["a", "b"])
        tm.assert_categorical_equal(result, expected)

    def test_take_fill_with_negative_one(self):
        # -1 was a category
        cat = pd.Categorical([-1, 0, 1])
        result = cat.take([0, -1, 1], allow_fill=True, fill_value=-1)
        expected = pd.Categorical([-1, -1, 0], categories=[-1, 0, 1])
        tm.assert_categorical_equal(result, expected)

    def test_take_fill_value(self):
        # https://github.com/pandas-dev/pandas/issues/23296
        cat = pd.Categorical(["a", "b", "c"])
        result = cat.take([0, 1, -1], fill_value="a", allow_fill=True)
        expected = pd.Categorical(["a", "b", "a"], categories=["a", "b", "c"])
        tm.assert_categorical_equal(result, expected)

    def test_take_fill_value_new_raises(self):
        # https://github.com/pandas-dev/pandas/issues/23296
        cat = pd.Categorical(["a", "b", "c"])
        xpr = r"'fill_value' \('d'\) is not in this Categorical's categories."
        with pytest.raises(TypeError, match=xpr):
            cat.take([0, 1, -1], fill_value="d", allow_fill=True)

    def test_take_nd_deprecated(self):
        cat = pd.Categorical(["a", "b", "c"])
        with tm.assert_produces_warning(FutureWarning):
            cat.take_nd([0, 1])

        ci = pd.Index(cat)
        with tm.assert_produces_warning(FutureWarning):
            ci.take_nd([0, 1])
