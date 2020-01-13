import operator

import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd
import pandas._testing as tm


def test_repr():
    df = pd.DataFrame({"A": pd.array(["a", pd.NA, "b"], dtype="string")})
    expected = "      A\n0     a\n1  <NA>\n2     b"
    assert repr(df) == expected

    expected = "0       a\n1    <NA>\n2       b\nName: A, dtype: string"
    assert repr(df.A) == expected

    expected = "<StringArray>\n['a', <NA>, 'b']\nLength: 3, dtype: string"
    assert repr(df.A.array) == expected


def test_none_to_nan():
    a = pd.arrays.StringArray._from_sequence(["a", None, "b"])
    assert a[1] is not None
    assert a[1] is pd.NA


def test_setitem_validates():
    a = pd.arrays.StringArray._from_sequence(["a", "b"])
    with pytest.raises(ValueError, match="10"):
        a[0] = 10

    with pytest.raises(ValueError, match="strings"):
        a[:] = np.array([1, 2])


def test_setitem_with_scalar_string():
    # is_float_dtype considers some strings, like 'd', to be floats
    # which can cause issues.
    arr = pd.array(["a", "c"], dtype="string")
    arr[0] = "d"
    expected = pd.array(["d", "c"], dtype="string")
    tm.assert_extension_array_equal(arr, expected)


@pytest.mark.parametrize(
    "input, method",
    [
        (["a", "b", "c"], operator.methodcaller("capitalize")),
        (["a", "b", "c"], operator.methodcaller("capitalize")),
        (["a b", "a bc. de"], operator.methodcaller("capitalize")),
    ],
)
def test_string_methods(input, method):
    a = pd.Series(input, dtype="string")
    b = pd.Series(input, dtype="object")
    result = method(a.str)
    expected = method(b.str)

    assert result.dtype.name == "string"
    tm.assert_series_equal(result.astype(object), expected)


def test_astype_roundtrip():
    s = pd.Series(pd.date_range("2000", periods=12))
    s[0] = None

    result = s.astype("string").astype("datetime64[ns]")
    tm.assert_series_equal(result, s)


def test_add():
    a = pd.Series(["a", "b", "c", None, None], dtype="string")
    b = pd.Series(["x", "y", None, "z", None], dtype="string")

    result = a + b
    expected = pd.Series(["ax", "by", None, None, None], dtype="string")
    tm.assert_series_equal(result, expected)

    result = a.add(b)
    tm.assert_series_equal(result, expected)

    result = a.radd(b)
    expected = pd.Series(["xa", "yb", None, None, None], dtype="string")
    tm.assert_series_equal(result, expected)

    result = a.add(b, fill_value="-")
    expected = pd.Series(["ax", "by", "c-", "-z", None], dtype="string")
    tm.assert_series_equal(result, expected)


def test_add_2d():
    a = pd.array(["a", "b", "c"], dtype="string")
    b = np.array([["a", "b", "c"]], dtype=object)
    with pytest.raises(ValueError, match="3 != 1"):
        a + b

    s = pd.Series(a)
    with pytest.raises(ValueError, match="3 != 1"):
        s + b


def test_add_sequence():
    a = pd.array(["a", "b", None, None], dtype="string")
    other = ["x", None, "y", None]

    result = a + other
    expected = pd.array(["ax", None, None, None], dtype="string")
    tm.assert_extension_array_equal(result, expected)

    result = other + a
    expected = pd.array(["xa", None, None, None], dtype="string")
    tm.assert_extension_array_equal(result, expected)


def test_mul():
    a = pd.array(["a", "b", None], dtype="string")
    result = a * 2
    expected = pd.array(["aa", "bb", None], dtype="string")
    tm.assert_extension_array_equal(result, expected)

    result = 2 * a
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.xfail(reason="GH-28527")
def test_add_strings():
    array = pd.array(["a", "b", "c", "d"], dtype="string")
    df = pd.DataFrame([["t", "u", "v", "w"]])
    assert array.__add__(df) is NotImplemented

    result = array + df
    expected = pd.DataFrame([["at", "bu", "cv", "dw"]]).astype("string")
    tm.assert_frame_equal(result, expected)

    result = df + array
    expected = pd.DataFrame([["ta", "ub", "vc", "wd"]]).astype("string")
    tm.assert_frame_equal(result, expected)


@pytest.mark.xfail(reason="GH-28527")
def test_add_frame():
    array = pd.array(["a", "b", np.nan, np.nan], dtype="string")
    df = pd.DataFrame([["x", np.nan, "y", np.nan]])

    assert array.__add__(df) is NotImplemented

    result = array + df
    expected = pd.DataFrame([["ax", np.nan, np.nan, np.nan]]).astype("string")
    tm.assert_frame_equal(result, expected)

    result = df + array
    expected = pd.DataFrame([["xa", np.nan, np.nan, np.nan]]).astype("string")
    tm.assert_frame_equal(result, expected)


def test_comparison_methods_scalar(all_compare_operators):
    op_name = all_compare_operators

    a = pd.array(["a", None, "c"], dtype="string")
    other = "a"
    result = getattr(a, op_name)(other)
    expected = np.array([getattr(item, op_name)(other) for item in a], dtype=object)
    expected = pd.array(expected, dtype="boolean")
    tm.assert_extension_array_equal(result, expected)

    result = getattr(a, op_name)(pd.NA)
    expected = pd.array([None, None, None], dtype="boolean")
    tm.assert_extension_array_equal(result, expected)


def test_comparison_methods_array(all_compare_operators):
    op_name = all_compare_operators

    a = pd.array(["a", None, "c"], dtype="string")
    other = [None, None, "c"]
    result = getattr(a, op_name)(other)
    expected = np.empty_like(a, dtype="object")
    expected[-1] = getattr(other[-1], op_name)(a[-1])
    expected = pd.array(expected, dtype="boolean")
    tm.assert_extension_array_equal(result, expected)

    result = getattr(a, op_name)(pd.NA)
    expected = pd.array([None, None, None], dtype="boolean")
    tm.assert_extension_array_equal(result, expected)


def test_constructor_raises():
    with pytest.raises(ValueError, match="sequence of strings"):
        pd.arrays.StringArray(np.array(["a", "b"], dtype="S1"))

    with pytest.raises(ValueError, match="sequence of strings"):
        pd.arrays.StringArray(np.array([]))


@pytest.mark.parametrize("skipna", [True, False])
@pytest.mark.xfail(reason="Not implemented StringArray.sum")
def test_reduce(skipna):
    arr = pd.Series(["a", "b", "c"], dtype="string")
    result = arr.sum(skipna=skipna)
    assert result == "abc"


@pytest.mark.parametrize("skipna", [True, False])
@pytest.mark.xfail(reason="Not implemented StringArray.sum")
def test_reduce_missing(skipna):
    arr = pd.Series([None, "a", None, "b", "c", None], dtype="string")
    result = arr.sum(skipna=skipna)
    if skipna:
        assert result == "abc"
    else:
        assert pd.isna(result)


@td.skip_if_no("pyarrow", min_version="0.15.0")
def test_arrow_array():
    # protocol added in 0.15.0
    import pyarrow as pa

    data = pd.array(["a", "b", "c"], dtype="string")
    arr = pa.array(data)
    expected = pa.array(list(data), type=pa.string(), from_pandas=True)
    assert arr.equals(expected)


@td.skip_if_no("pyarrow", min_version="0.15.1.dev")
def test_arrow_roundtrip():
    # roundtrip possible from arrow 1.0.0
    import pyarrow as pa

    data = pd.array(["a", "b", None], dtype="string")
    df = pd.DataFrame({"a": data})
    table = pa.table(df)
    assert table.field("a").type == "string"
    result = table.to_pandas()
    assert isinstance(result["a"].dtype, pd.StringDtype)
    tm.assert_frame_equal(result, df)
    # ensure the missing value is represented by NA and not np.nan or None
    assert result.loc[2, "a"] is pd.NA


def test_value_counts_na():
    arr = pd.array(["a", "b", "a", pd.NA], dtype="string")
    result = arr.value_counts(dropna=False)
    expected = pd.Series([2, 1, 1], index=["a", "b", pd.NA], dtype="Int64")
    tm.assert_series_equal(result, expected)

    result = arr.value_counts(dropna=True)
    expected = pd.Series([2, 1], index=["a", "b"], dtype="Int64")
    tm.assert_series_equal(result, expected)
