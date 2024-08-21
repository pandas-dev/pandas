import numpy as np

from pandas import (
    DataFrame,
    Index,
    MultiIndex,
    Series,
    _testing as tm,
)


def test_get_dummies(any_string_dtype):
    s = Series(["a|b", "a|c", np.nan], dtype=any_string_dtype)
    result = s.str.get_dummies("|")
    expected = DataFrame([[1, 1, 0], [1, 0, 1], [0, 0, 0]], columns=list("abc"))
    tm.assert_frame_equal(result, expected)

    s = Series(["a;b", "a", 7], dtype=any_string_dtype)
    result = s.str.get_dummies(";")
    expected = DataFrame([[0, 1, 1], [0, 1, 0], [1, 0, 0]], columns=list("7ab"))
    tm.assert_frame_equal(result, expected)


def test_get_dummies_index():
    # GH9980, GH8028
    idx = Index(["a|b", "a|c", "b|c"])
    result = idx.str.get_dummies("|")

    expected = MultiIndex.from_tuples(
        [(1, 1, 0), (1, 0, 1), (0, 1, 1)], names=("a", "b", "c")
    )
    tm.assert_index_equal(result, expected)


def test_get_dummies_with_name_dummy(any_string_dtype):
    # GH 12180
    # Dummies named 'name' should work as expected
    s = Series(["a", "b,name", "b"], dtype=any_string_dtype)
    result = s.str.get_dummies(",")
    expected = DataFrame([[1, 0, 0], [0, 1, 1], [0, 1, 0]], columns=["a", "b", "name"])
    tm.assert_frame_equal(result, expected)


def test_get_dummies_with_name_dummy_index():
    # GH 12180
    # Dummies named 'name' should work as expected
    idx = Index(["a|b", "name|c", "b|name"])
    result = idx.str.get_dummies("|")

    expected = MultiIndex.from_tuples(
        [(1, 1, 0, 0), (0, 0, 1, 1), (0, 1, 0, 1)], names=("a", "b", "c", "name")
    )
    tm.assert_index_equal(result, expected)


def test_get_dummies_with_prefix(any_string_dtype):
    s = Series(["a|b", "a|c", np.nan], dtype=any_string_dtype)
    result = s.str.get_dummies(sep="|", prefix="prefix")
    expected = DataFrame(
        [[1, 1, 0], [1, 0, 1], [0, 0, 0]],
        columns=["prefix_a", "prefix_b", "prefix_c"],
    )
    tm.assert_frame_equal(result, expected)


def test_get_dummies_with_prefix_sep(any_string_dtype):
    s = Series(["a|b", "a|c", np.nan], dtype=any_string_dtype)
    result = s.str.get_dummies(sep="|", prefix=None, prefix_sep="__")
    expected = DataFrame([[1, 1, 0], [1, 0, 1], [0, 0, 0]], columns=["a", "b", "c"])
    tm.assert_frame_equal(result, expected)

    result = s.str.get_dummies(sep="|", prefix="col", prefix_sep="__")
    expected = DataFrame(
        [[1, 1, 0], [1, 0, 1], [0, 0, 0]],
        columns=["col__a", "col__b", "col__c"],
    )
    tm.assert_frame_equal(result, expected)


def test_get_dummies_with_dummy_na(any_string_dtype):
    s = Series(["a|b", "a|c", np.nan], dtype=any_string_dtype)
    result = s.str.get_dummies(sep="|", dummy_na=True)
    expected = DataFrame(
        [[1, 1, 0, 0], [1, 0, 1, 0], [0, 0, 0, 1]],
        columns=["a", "b", "c", np.nan],
    )
    tm.assert_frame_equal(result, expected)


def test_get_dummies_with_dtype(any_string_dtype):
    s = Series(["a|b", "a|c", np.nan], dtype=any_string_dtype)
    result = s.str.get_dummies(sep="|", dtype=bool)
    expected = DataFrame(
        [[True, True, False], [True, False, True], [False, False, False]],
        columns=["a", "b", "c"],
    )
    tm.assert_frame_equal(result, expected)
    assert (result.dtypes == bool).all()


def test_get_dummies_with_prefix_dict(any_string_dtype):
    s = Series(["a|b", "a|c", np.nan], dtype=any_string_dtype)
    prefix = {"a": "alpha", "b": "beta", "c": "gamma"}
    result = s.str.get_dummies(sep="|", prefix=prefix)
    expected = DataFrame(
        [[1, 1, 0], [1, 0, 1], [0, 0, 0]],
        columns=["alpha_a", "beta_b", "gamma_c"],
    )
    tm.assert_frame_equal(result, expected)
