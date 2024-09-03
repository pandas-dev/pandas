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


def test_get_dummies_int8_dtype():
    s = Series(["1|2", "1|3", np.nan], dtype="string")
    result = s.str.get_dummies("|", dtype=np.int8)
    expected = DataFrame(
        [[1, 1, 0], [1, 0, 1], [0, 0, 0]], columns=list("123"), dtype=np.int8
    )
    tm.assert_frame_equal(result, expected)
    assert (result.dtypes == np.int8).all()


def test_get_dummies_uint8_dtype():
    s = Series(["a|b", "a|c", np.nan], dtype="string")
    result = s.str.get_dummies("|", dtype=np.uint8)
    expected = DataFrame(
        [[1, 1, 0], [1, 0, 1], [0, 0, 0]], columns=list("abc"), dtype=np.uint8
    )
    tm.assert_frame_equal(result, expected)
    assert (result.dtypes == np.uint8).all()


def test_get_dummies_int16_dtype():
    s = Series(["a|b", "a|c", np.nan], dtype="string")
    result = s.str.get_dummies("|", dtype=np.int16)
    expected = DataFrame(
        [[1, 1, 0], [1, 0, 1], [0, 0, 0]], columns=list("abc"), dtype=np.int16
    )
    tm.assert_frame_equal(result, expected)
    assert (result.dtypes == np.int16).all()


def test_get_dummies_uint16_dtype():
    s = Series(["a|b", "a|c", np.nan], dtype="string")
    result = s.str.get_dummies("|", dtype=np.uint16)
    expected = DataFrame(
        [[1, 1, 0], [1, 0, 1], [0, 0, 0]], columns=list("abc"), dtype=np.uint16
    )
    tm.assert_frame_equal(result, expected)
    assert (result.dtypes == np.uint16).all()


def test_get_dummies_int32_dtype():
    s = Series(["x|y", "x|z", np.nan], dtype="string")
    result = s.str.get_dummies("|", dtype=np.int32)
    expected = DataFrame(
        [[1, 1, 0], [1, 0, 1], [0, 0, 0]], columns=list("xyz"), dtype=np.int32
    )
    tm.assert_frame_equal(result, expected)
    assert (result.dtypes == np.int32).all()


def test_get_dummies_uint32_dtype():
    s = Series(["x|y", "x|z", np.nan], dtype="string")
    result = s.str.get_dummies("|", dtype=np.uint32)
    expected = DataFrame(
        [[1, 1, 0], [1, 0, 1], [0, 0, 0]], columns=list("xyz"), dtype=np.uint32
    )
    tm.assert_frame_equal(result, expected)
    assert (result.dtypes == np.uint32).all()


def test_get_dummies_int64_dtype():
    s = Series(["foo|bar", "foo|baz", np.nan], dtype="string")
    result = s.str.get_dummies("|", dtype=np.int64)
    expected = DataFrame(
        [[1, 0, 1], [0, 1, 1], [0, 0, 0]], columns=["bar", "baz", "foo"], dtype=np.int64
    )
    tm.assert_frame_equal(result, expected)
    assert (result.dtypes == np.int64).all()


def test_get_dummies_uint64_dtype():
    s = Series(["foo|bar", "foo|baz", np.nan], dtype="string")
    result = s.str.get_dummies("|", dtype=np.uint64)
    expected = DataFrame(
        [[1, 0, 1], [0, 1, 1], [0, 0, 0]],
        columns=["bar", "baz", "foo"],
        dtype=np.uint64,
    )
    tm.assert_frame_equal(result, expected)
    assert (result.dtypes == np.uint64).all()


def test_get_dummies_bool_dtype():
    s = Series(["a|b", "a|c", np.nan], dtype="string")
    result = s.str.get_dummies("|", dtype=bool)
    expected = DataFrame(
        [[True, True, False], [True, False, True], [False, False, False]],
        columns=["a", "b", "c"],
    )
    tm.assert_frame_equal(result, expected)
    assert (result.dtypes == bool).all()
