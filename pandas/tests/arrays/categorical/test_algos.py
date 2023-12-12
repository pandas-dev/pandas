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


@pytest.mark.parametrize("value", [[""], [None, ""], [pd.NaT, ""]])
def test_isin_cats_corner_cases(value):
    # GH36550
    cat = pd.Categorical([""])
    result = cat.isin(value)
    expected = np.array([True], dtype=bool)
    tm.assert_numpy_array_equal(expected, result)


@pytest.mark.parametrize("empty", [[], pd.Series(dtype=object), np.array([])])
def test_isin_empty(empty):
    s = pd.Categorical(["a", "b"])
    expected = np.array([False, False], dtype=bool)

    result = s.isin(empty)
    tm.assert_numpy_array_equal(expected, result)


def test_diff():
    ser = pd.Series([1, 2, 3], dtype="category")

    msg = "Convert to a suitable dtype"
    with pytest.raises(TypeError, match=msg):
        ser.diff()

    df = ser.to_frame(name="A")
    with pytest.raises(TypeError, match=msg):
        df.diff()


def test_factorize_with_original_factorization():
    # Original factorization setup
    original_values = np.array(['apple', 'banana', 'apple', 'cherry'])
    original_codes, original_uniques = pd.factorize(original_values)
    original_factorization = (original_uniques, original_codes)

    # New values including some from original and some new
    new_values = np.array(['banana', 'apple', 'durian', 'apple', 'elderberry'])

    # Expected outcome
    expected_new_codes = np.array([1, 0, 3, 0, 4])  # 'banana' and 'apple' should have original codes
    expected_new_uniques = np.array(['apple', 'banana', 'cherry', 'durian', 'elderberry'])

    # Run factorize with original_factorization
    new_codes, new_uniques = pd.factorize(new_values, original_factorization=original_factorization)

    # Assertions
    np.testing.assert_array_equal(new_codes, expected_new_codes)
    np.testing.assert_array_equal(new_uniques, expected_new_uniques)


def test_factorize_with_empty_original_factorization():
    values = np.array(['x', 'y', 'z'])
    original_factorization = (np.array([]), np.array([]))

    expected_codes = np.array([0, 1, 2])
    expected_uniques = np.array(['x', 'y', 'z'])

    codes, uniques = pd.factorize(values, original_factorization=original_factorization)

    np.testing.assert_array_equal(codes, expected_codes)
    np.testing.assert_array_equal(uniques, expected_uniques)


def test_factorize_with_non_overlapping_original_factorization():
    original_values = np.array(['a', 'b', 'c'])
    new_values = np.array(['x', 'y', 'z'])
    original_codes, original_uniques = pd.factorize(original_values)
    original_factorization = (original_uniques, original_codes)

    expected_new_codes = np.array([3, 4, 5])  # completely new codes
    expected_new_uniques = np.array(['a', 'b', 'c', 'x', 'y', 'z'])

    new_codes, new_uniques = pd.factorize(new_values, original_factorization=original_factorization)

    np.testing.assert_array_equal(new_codes, expected_new_codes)
    np.testing.assert_array_equal(new_uniques, expected_new_uniques)

def test_factorize_with_nan_in_original_factorization():
    original_values = np.array(['a', 'b', np.nan, 'd'])
    new_values = np.array(['a', 'b', 'c', np.nan])
    original_codes, original_uniques = pd.factorize(original_values)
    original_factorization = (original_uniques, original_codes)

    expected_new_codes = np.array([0, 1, 4, 2])  # 'c' is new, NaN should match original
    expected_new_uniques = np.array(['a', 'b', np.nan, 'd', 'c'])

    new_codes, new_uniques = pd.factorize(new_values, original_factorization=original_factorization)

    np.testing.assert_array_equal(new_codes, expected_new_codes)
    np.testing.assert_array_equal(new_uniques, expected_new_uniques)

def test_factorize_with_different_data_types():
    original_values = np.array([1, 2, 3])  # integers
    new_values = np.array(['1', '4', '2'])  # strings representing numbers
    original_codes, original_uniques = pd.factorize(original_values)
    original_factorization = (original_uniques, original_codes)

    expected_new_codes = np.array([3, 4, 1])  # '1' and '2' as strings are new
    expected_new_uniques = np.array([1, 2, 3, '1', '4'])

    new_codes, new_uniques = pd.factorize(new_values, original_factorization=original_factorization)
