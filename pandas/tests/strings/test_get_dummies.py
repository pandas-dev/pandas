import numpy as np

from pandas import (
    DataFrame,
    Index,
    MultiIndex,
    Series,
    _testing as tm,
    generate_fake_dataframe
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


def test_generate_fake_dataframe():
    # Assert that a dataframe with all different dtypes can be created
    df = generate_fake_dataframe(10, "icfd")
    assert df.shape == (10 , 4)


def test_generate_fake_dataframe_of_integers_has_all_zeroes():
    # Assert that a dataframe of integers between (0,1) has all zeroes
    df = generate_fake_dataframe(10, "iiii", intervals={"i" : (0 , 1)})
    assert (df.values == np.zeros((10 , 4))).all()


def test_generate_fake_dataframe_set_dates_interval():
    # Assert that a dataframe can be created with dates formatted as %Y–%m–%d
    df = generate_fake_dataframe(1, "d", intervals={"d" : ("2020-09-01", "2020-09-02")})
    assert isinstance(df.values[0][0], np.datetime64)


def test_generate_fake_dataframe_with_none_in_intervals():
    # Assert that a dataframe can be created with an interval list with Nones in it.
    df = generate_fake_dataframe(1, "iii", intervals=[None, (4 , 6), (6 , 10)])
    assert df.shape == (1 , 3)
