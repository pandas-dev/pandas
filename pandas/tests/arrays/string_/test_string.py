import operator

import numpy as np
import pytest

import pandas as pd
import pandas.util.testing as tm


def test_none_to_nan():
    a = pd.arrays.TextArray._from_sequence(["a", None, "b"])
    assert a[1] is not None
    assert np.isnan(a[1])


def test_setitem_validates():
    a = pd.arrays.TextArray._from_sequence(["a", "b"])
    with pytest.raises(ValueError, match="10"):
        a[0] = 10

    with pytest.raises(ValueError, match="strings"):
        a[:] = np.array([1, 2])


@pytest.mark.parametrize(
    "input, method",
    [
        (["a", "b", "c"], operator.methodcaller("capitalize")),
        (["a", "b", "c"], operator.methodcaller("capitalize")),
        (["a b", "a bc. de"], operator.methodcaller("capitalize")),
    ],
)
def test_string_methods(input, method):
    a = pd.Series(input, dtype="text")
    b = pd.Series(input, dtype="object")
    result = method(a.str)
    expected = method(b.str)

    assert result.dtype.name == "text"
    tm.assert_series_equal(result.astype(object), expected)


def test_astype_roundtrip():
    s = pd.Series(pd.date_range("2000", periods=12))
    s[0] = None

    result = s.astype("text").astype("datetime64[ns]")
    tm.assert_series_equal(result, s)


def test_add():
    a = pd.Series(["a", "b", "c", None, None], dtype="text")
    b = pd.Series(["x", "y", None, "z", None], dtype="text")

    result = a + b
    expected = pd.Series(["ax", "by", None, None, None], dtype="text")
    tm.assert_series_equal(result, expected)

    result = a.add(b)
    tm.assert_series_equal(result, expected)

    result = a.radd(b)
    expected = pd.Series(["xa", "yb", None, None, None], dtype="text")
    tm.assert_series_equal(result, expected)

    result = a.add(b, fill_value="-")
    expected = pd.Series(["ax", "by", "c-", "-z", None], dtype="text")
    tm.assert_series_equal(result, expected)


def test_add_sequence():
    a = pd.array(["a", "b", None, None], dtype="text")
    other = ["x", None, "y", None]

    result = a + other
    expected = pd.array(["ax", None, None, None], dtype="text")
    tm.assert_extension_array_equal(result, expected)

    result = other + a
    expected = pd.array(["xa", None, None, None], dtype="text")
    tm.assert_extension_array_equal(result, expected)


def test_mul():
    a = pd.array(["a", "b", None], dtype="text")
    result = a * 2
    expected = pd.array(["aa", "bb", None], dtype="text")
    tm.assert_extension_array_equal(result, expected)

    result = 2 * a
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.xfail(reason="GH-28527")
def test_add_strings():
    array = pd.array(["a", "b", "c", "d"], dtype="text")
    df = pd.DataFrame([["t", "u", "v", "w"]])
    assert array.__add__(df) is NotImplemented

    result = array + df
    expected = pd.DataFrame([["at", "bu", "cv", "dw"]]).astype("text")
    tm.assert_frame_equal(result, expected)

    result = df + array
    expected = pd.DataFrame([["ta", "ub", "vc", "wd"]]).astype("text")
    tm.assert_frame_equal(result, expected)


@pytest.mark.xfail(reason="GH-28527")
def test_add_frame():
    array = pd.array(["a", "b", np.nan, np.nan], dtype="text")
    df = pd.DataFrame([["x", np.nan, "y", np.nan]])

    assert array.__add__(df) is NotImplemented

    result = array + df
    expected = pd.DataFrame([["ax", np.nan, np.nan, np.nan]]).astype("text")
    tm.assert_frame_equal(result, expected)

    result = df + array
    expected = pd.DataFrame([["xa", np.nan, np.nan, np.nan]]).astype("text")
    tm.assert_frame_equal(result, expected)


def test_constructor_raises():
    with pytest.raises(ValueError, match="object-dtype ndarray"):
        pd.arrays.TextArray(np.array(["a", "b"], dtype="S1"))

    with pytest.raises(ValueError, match="object-dtype ndarray"):
        pd.arrays.TextArray(np.array([]))


@pytest.mark.parametrize("skipna", [True, False])
def test_reduce(skipna):
    arr = pd.Series(["a", "b", "c"], dtype="text")
    result = arr.sum(skipna=skipna)
    assert result == "abc"


@pytest.mark.parametrize("skipna", [True, False])
def test_reduce_missing(skipna):
    arr = pd.Series([None, "a", None, "b", "c", None], dtype="text")
    result = arr.sum(skipna=skipna)
    if skipna:
        assert result == "abc"
    else:
        assert pd.isna(result)
