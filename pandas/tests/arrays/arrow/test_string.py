from __future__ import annotations

import re
import sys
import unicodedata

import numpy as np
import pytest

from pandas.core.dtypes.dtypes import (
    ArrowDtype,
)

import pandas as pd
import pandas._testing as tm

pa = pytest.importorskip("pyarrow")

from pandas.core.arrays.arrow.array import ArrowExtensionArray


def test_arrow_string_multiplication():
    # GH 56537
    binary = pd.Series(["abc", "defg"], dtype=ArrowDtype(pa.string()))
    repeat = pd.Series([2, -2], dtype="int64[pyarrow]")
    result = binary * repeat
    expected = pd.Series(["abcabc", ""], dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)
    reflected_result = repeat * binary
    tm.assert_series_equal(result, reflected_result)


def test_arrow_string_multiplication_scalar_repeat():
    binary = pd.Series(["abc", "defg"], dtype=ArrowDtype(pa.string()))
    result = binary * 2
    expected = pd.Series(["abcabc", "defgdefg"], dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)
    reflected_result = 2 * binary
    tm.assert_series_equal(reflected_result, expected)


def test_arrow_string_addition_mixed_string_types():
    # https://github.com/pandas-dev/pandas/issues/65220
    left = pd.Series(["a", None], dtype=ArrowDtype(pa.string()))
    right = pd.Series(["b", "c"], dtype=ArrowDtype(pa.large_string()))

    result = left + right
    expected = pd.Series(["ab", None], dtype=ArrowDtype(pa.large_string()))
    tm.assert_series_equal(result, expected)

    reflected_result = right + left
    expected_reflected = pd.Series(["ba", None], dtype=ArrowDtype(pa.large_string()))
    tm.assert_series_equal(reflected_result, expected_reflected)


@pytest.mark.parametrize("string_type", [pa.string(), pa.large_string()])
def test_arrow_string_addition_mixed_with_binary_raises(string_type):
    left = pd.Series(["a", None], dtype=ArrowDtype(string_type))
    right = pd.Series([b"b", b"c"], dtype=ArrowDtype(pa.binary()))

    msg = (
        f"operation 'add' not supported for dtype '{left.dtype}' "
        f"with dtype '{right.dtype}'"
    )
    with pytest.raises(TypeError, match=re.escape(msg)):
        left + right

    reflected_msg = (
        f"operation 'add' not supported for dtype '{right.dtype}' "
        f"with dtype '{left.dtype}'"
    )
    with pytest.raises(TypeError, match=re.escape(reflected_msg)):
        right + left


@pytest.mark.parametrize("pat", ["abc", "a[a-z]{2}"])
def test_str_count(pat):
    ser = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.count(pat)
    expected = pd.Series([1, None], dtype=ArrowDtype(pa.int32()))
    tm.assert_series_equal(result, expected)


def test_str_count_flags():
    # GH#66348
    ser = pd.Series(["abc", "éxy", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.count(r"\w", flags=re.ASCII)
    expected = pd.Series([3, 2, None], dtype=ArrowDtype(pa.int32()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "side, str_func", [["left", "rjust"], ["right", "ljust"], ["both", "center"]]
)
def test_str_pad(side, str_func):
    ser = pd.Series(["a", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.pad(width=3, side=side, fillchar="x")
    expected = pd.Series(
        [getattr("a", str_func)(3, "x"), None], dtype=ArrowDtype(pa.string())
    )
    tm.assert_series_equal(result, expected)


def test_str_pad_invalid_side():
    ser = pd.Series(["a", None], dtype=ArrowDtype(pa.string()))
    with pytest.raises(ValueError, match="Invalid side: foo"):
        ser.str.pad(3, "foo", "x")


@pytest.mark.parametrize(
    "pat, case, na, regex, exp",
    [
        ["ab", False, None, False, [True, None]],
        ["Ab", True, None, False, [False, None]],
        ["ab", False, True, False, [True, True]],
        ["a[a-z]{1}", False, None, True, [True, None]],
        ["A[a-z]{1}", True, None, True, [False, None]],
    ],
)
def test_str_contains(pat, case, na, regex, exp):
    ser = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.contains(pat, case=case, na=na, regex=regex)
    expected = pd.Series(exp, dtype=ArrowDtype(pa.bool_()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("method", ["contains", "match", "fullmatch"])
def test_str_contains_match_flags(method):
    # GH#66348 re.ASCII and re.UNICODE are mutually exclusive, so the flag has
    #  to reach `re` rather than being silently dropped
    ser = pd.Series(["abc", "éxy", None], dtype=ArrowDtype(pa.string()))
    result = getattr(ser.str, method)(r"^\w+$", flags=re.ASCII)
    expected = pd.Series([True, False, None], dtype=ArrowDtype(pa.bool_()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("method", ["contains", "match", "fullmatch"])
def test_str_contains_match_compiled_pattern(method):
    # GH#66348
    ser = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    result = getattr(ser.str, method)(re.compile("ABC", re.IGNORECASE))
    expected = pd.Series([True, None], dtype=ArrowDtype(pa.bool_()))
    tm.assert_series_equal(result, expected)


def test_str_match_multiline_flag():
    # GH#66348 the anchors added for match/fullmatch must not become line
    #  anchors under re.MULTILINE
    ser = pd.Series(["a\nb", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.fullmatch("b", flags=re.MULTILINE)
    expected = pd.Series([False, None], dtype=ArrowDtype(pa.bool_()))
    tm.assert_series_equal(result, expected)


def test_str_contains_unicode_escape():
    # GH 63901, GH#63683 patterns are Python regular expressions
    ser = pd.Series(["a", "\u0e01", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.contains(r"[\u0e00-\u0e7f]")
    expected = pd.Series([False, True, None], dtype=ArrowDtype(pa.bool_()))
    tm.assert_series_equal(result, expected)

    with pytest.raises(re.error, match="incomplete escape"):
        ser.str.contains(r"[\x{0e00}-\x{0e7f}]")


@pytest.mark.parametrize(
    "side, pat, na, exp",
    [
        ["startswith", "ab", None, [True, None, False]],
        ["startswith", "b", False, [False, False, False]],
        ["endswith", "b", True, [False, True, False]],
        ["endswith", "bc", None, [True, None, False]],
        ["startswith", ("a", "e", "g"), None, [True, None, True]],
        ["endswith", ("a", "c", "g"), None, [True, None, True]],
        ["startswith", (), None, [False, None, False]],
        ["endswith", (), None, [False, None, False]],
    ],
)
def test_str_start_ends_with(side, pat, na, exp):
    ser = pd.Series(["abc", None, "efg"], dtype=ArrowDtype(pa.string()))
    result = getattr(ser.str, side)(pat, na=na)
    expected = pd.Series(exp, dtype=ArrowDtype(pa.bool_()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("side", ("startswith", "endswith"))
def test_str_starts_ends_with_all_nulls_empty_tuple(side):
    ser = pd.Series([None, None], dtype=ArrowDtype(pa.string()))
    result = getattr(ser.str, side)(())

    # bool datatype preserved for all nulls.
    expected = pd.Series([None, None], dtype=ArrowDtype(pa.bool_()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "kwargs, exp",
    [
        [{"pat": re.compile("b")}, ["axc", None]],
        [{"repl": lambda match: match.group().upper()}, ["aBc", None]],
        [{"pat": "B", "case": False}, ["axc", None]],
        [{"pat": "B", "flags": re.IGNORECASE}, ["axc", None]],
        [{"pat": "(?P<mid>b)", "repl": r"[\g<mid>]"}, ["a[b]c", None]],
    ],
)
def test_str_replace_re_fallback(kwargs, exp):
    # GH#66348 these are not expressible with pyarrow's kernel, so they are
    #  evaluated with `re` instead of raising
    ser = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.replace(**{"pat": "b", "repl": "x", "regex": True, **kwargs})
    expected = pd.Series(exp, dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "pat, repl, n, regex, exp",
    [
        ["a", "x", -1, False, ["xbxc", None]],
        ["a", "x", 1, False, ["xbac", None]],
        ["[a-b]", "x", -1, True, ["xxxc", None]],
    ],
)
def test_str_replace(pat, repl, n, regex, exp):
    ser = pd.Series(["abac", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.replace(pat, repl, n=n, regex=regex)
    expected = pd.Series(exp, dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


def test_str_replace_re2_unicode_property_raises():
    # GH#63683 patterns are Python regular expressions
    ser = pd.Series(["Jan", "Feb", None], dtype=ArrowDtype(pa.string()))
    with pytest.raises(re.error, match="bad escape"):
        ser.str.replace(r"\p{Lu}", "U", regex=True)


def test_str_replace_negative_n():
    # GH 56404
    ser = pd.Series(["abc", "aaaaaa"], dtype=ArrowDtype(pa.string()))
    actual = ser.str.replace("a", "", -3, True)
    expected = pd.Series(["bc", ""], dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(expected, actual)

    # Same bug for pyarrow-backed StringArray GH#59628
    ser2 = ser.astype(pd.StringDtype(storage="pyarrow"))
    actual2 = ser2.str.replace("a", "", -3, True)
    expected2 = expected.astype(ser2.dtype)
    tm.assert_series_equal(expected2, actual2)

    ser3 = ser.astype(pd.StringDtype(storage="pyarrow", na_value=np.nan))
    actual3 = ser3.str.replace("a", "", -3, True)
    expected3 = expected.astype(ser3.dtype)
    tm.assert_series_equal(expected3, actual3)


def test_str_replace_empty_pattern():
    # https://github.com/pandas-dev/pandas/issues/64941
    ser = pd.Series(["abcd"], dtype=ArrowDtype(pa.string()))

    result = ser.str.replace("", "")
    expected = pd.Series(["abcd"], dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)

    result = ser.str.replace("", "X")
    expected = pd.Series(["XaXbXcXdX"], dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


def test_str_repeat_unsupported():
    ser = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    with pytest.raises(NotImplementedError, match="repeat is not"):
        ser.str.repeat([1, 2])


def test_str_repeat():
    ser = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.repeat(2)
    expected = pd.Series(["abcabc", None], dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "pat, case, na, exp",
    [
        ["ab", False, None, [True, None]],
        ["Ab", True, None, [False, None]],
        ["bc", True, None, [False, None]],
        ["ab", False, True, [True, True]],
        ["a[a-z]{1}", False, None, [True, None]],
        ["A[a-z]{1}", True, None, [False, None]],
    ],
)
def test_str_match(pat, case, na, exp):
    ser = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.match(pat, case=case, na=na)
    expected = pd.Series(exp, dtype=ArrowDtype(pa.bool_()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "pat, case, na, exp",
    # Note: keep cases in sync with
    # pandas/tests/strings/test_find_replace.py::test_str_fullmatch_extra_cases
    [
        ["abc", False, None, [True, False, False, None]],
        ["Abc", True, None, [False, False, False, None]],
        ["bc", True, None, [False, False, False, None]],
        ["ab", False, None, [False, False, False, None]],
        ["a[a-z]{2}", False, None, [True, False, False, None]],
        ["A[a-z]{1}", True, None, [False, False, False, None]],
        # GH Issue: #56652
        ["abc$", False, None, [True, False, False, None]],
        ["abc\\$", False, None, [False, True, False, None]],
        ["Abc$", True, None, [False, False, False, None]],
        ["Abc\\$", True, None, [False, False, False, None]],
        # https://github.com/pandas-dev/pandas/issues/61072
        ["(abc)|(abx)", True, None, [True, False, False, None]],
        ["((abc)|(abx))", True, None, [True, False, False, None]],
    ],
)
def test_str_fullmatch(pat, case, na, exp):
    ser = pd.Series(["abc", "abc$", "$abc", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.fullmatch(pat, case=case, na=na)
    expected = pd.Series(exp, dtype=ArrowDtype(pa.bool_()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "sub, start, end, exp",
    [
        ["ab", 0, None, [0, None]],
        ["bc", 1, 3, [1, None]],
        ["ab", 1, 3, [-1, None]],
        ["ab", -3, -3, [-1, None]],
    ],
)
def test_str_find(sub, start, end, exp):
    ser = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.find(sub, start=start, end=end)
    expected = pd.Series(exp, dtype=ArrowDtype(pa.int64()))
    tm.assert_series_equal(result, expected)


def test_str_find_negative_start():
    # GH 56411
    ser = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.find(sub="b", start=-1000, end=3)
    expected = pd.Series([1, None], dtype=ArrowDtype(pa.int64()))
    tm.assert_series_equal(result, expected)


def test_str_find_no_end():
    ser = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.find("ab", start=1)
    expected = pd.Series([-1, None], dtype="int64[pyarrow]")
    tm.assert_series_equal(result, expected)


def test_str_find_negative_start_negative_end():
    # GH 56791
    ser = pd.Series(["abcdefg", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.find(sub="d", start=-6, end=-3)
    expected = pd.Series([3, None], dtype=ArrowDtype(pa.int64()))
    tm.assert_series_equal(result, expected)


def test_str_find_large_start():
    # GH 56791
    ser = pd.Series(["abcdefg", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.find(sub="d", start=16)
    expected = pd.Series([-1, None], dtype=ArrowDtype(pa.int64()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("start", [-15, -3, 0, 1, 15, None])
@pytest.mark.parametrize("end", [-15, -1, 0, 3, 15, None])
@pytest.mark.parametrize("sub", ["", "az", "abce", "a", "caa"])
def test_str_find_e2e(start, end, sub):
    s = pd.Series(
        ["abcaadef", "abc", "abcdeddefgj8292", "ab", "a", ""],
        dtype=ArrowDtype(pa.string()),
    )
    object_series = s.astype(pd.StringDtype(storage="python"))
    result = s.str.find(sub, start, end)
    expected = object_series.str.find(sub, start, end).astype(result.dtype)
    tm.assert_series_equal(result, expected)

    arrow_str_series = s.astype(pd.StringDtype(storage="pyarrow"))
    result2 = arrow_str_series.str.find(sub, start, end).astype(result.dtype)
    tm.assert_series_equal(result2, expected)


def test_str_find_negative_start_negative_end_no_match():
    # GH 56791
    ser = pd.Series(["abcdefg", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.find(sub="d", start=-3, end=-6)
    expected = pd.Series([-1, None], dtype=ArrowDtype(pa.int64()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "idx, expected_values",
    [
        [0, ["a", "d", None]],
        [1, ["b", "e", None]],
        [-1, ["c", "e", None]],
        [2, ["c", None, None]],
        [-3, ["a", None, None]],
        [4, [None, None, None]],
    ],
)
def test_str_get(idx, expected_values):
    ser = pd.Series(["abc", "de", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.get(idx)
    expected = pd.Series(expected_values, dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "idx, expected_values",
    [
        [0, ["a", "d", None]],
        [1, ["b", "e", None]],
        [-1, ["c", "e", None]],
        [2, ["c", None, None]],
        [-3, ["a", None, None]],
        [4, [None, None, None]],
    ],
)
def test_str_getitem(idx, expected_values):
    # GH 65112
    ser = pd.Series(["abc", "de", None], dtype=ArrowDtype(pa.string()))
    result = ser.str[idx]
    expected = pd.Series(expected_values, dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


@pytest.mark.xfail(
    reason="TODO: StringMethods._validate should support Arrow list types",
    raises=AttributeError,
)
def test_str_join():
    ser = pd.Series(ArrowExtensionArray(pa.array([list("abc"), list("123"), None])))
    result = ser.str.join("=")
    expected = pd.Series(["a=b=c", "1=2=3", None], dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


def test_str_join_string_type():
    ser = pd.Series(ArrowExtensionArray(pa.array(["abc", "123", None])))
    result = ser.str.join("=")
    expected = pd.Series(["a=b=c", "1=2=3", None], dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "start, stop, step, expected_values",
    [
        [None, 2, None, ["ab", None]],
        [None, 2, 1, ["ab", None]],
        [1, 3, 1, ["bc", None]],
        (None, None, -1, ["dcba", None]),
        (1, None, 2, ["bd", None]),
        (None, None, None, ["abcd", None]),
    ],
)
def test_str_slice(start, stop, step, expected_values):
    ser = pd.Series(["abcd", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.slice(start, stop, step)
    expected = pd.Series(expected_values, dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "start, stop, step, expected_values",
    [
        [None, 2, None, ["ab", None]],
        [None, 2, 1, ["ab", None]],
        [1, 3, 1, ["bc", None]],
        (None, None, -1, ["dcba", None]),
        (1, None, 2, ["bd", None]),
        (None, None, None, ["abcd", None]),
    ],
)
def test_str_getitem_range(start, stop, step, expected_values):
    # GH 65112
    ser = pd.Series(["abcd", None], dtype=ArrowDtype(pa.string()))
    result = ser.str[slice(start, stop, step)]
    expected = pd.Series(expected_values, dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "start, stop, repl, exp",
    [
        [1, 2, "x", ["axcd", None]],
        [None, 2, "x", ["xcd", None]],
        [None, 2, None, ["cd", None]],
    ],
)
def test_str_slice_replace(start, stop, repl, exp):
    ser = pd.Series(["abcd", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.slice_replace(start, stop, repl)
    expected = pd.Series(exp, dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "value, method, exp",
    [
        ["a1c", "isalnum", True],
        ["!|,", "isalnum", False],
        ["aaa", "isalpha", True],
        ["!!!", "isalpha", False],
        ["٠", "isdecimal", True],  # noqa: RUF001
        ["~!", "isdecimal", False],
        ["2", "isdigit", True],
        ["~", "isdigit", False],
        ["aaa", "islower", True],
        ["aaA", "islower", False],
        ["123", "isnumeric", True],
        ["11I", "isnumeric", False],
        [" ", "isspace", True],
        ["", "isspace", False],
        ["The That", "istitle", True],
        ["the That", "istitle", False],
        ["AAA", "isupper", True],
        ["AAc", "isupper", False],
    ],
)
def test_str_is_functions(value, method, exp):
    ser = pd.Series([value, None], dtype=ArrowDtype(pa.string()))
    result = getattr(ser.str, method)()
    expected = pd.Series([exp, None], dtype=ArrowDtype(pa.bool_()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "method, exp",
    [
        ["capitalize", "Abc def"],
        ["title", "Abc Def"],
        ["swapcase", "AbC Def"],
        ["lower", "abc def"],
        ["upper", "ABC DEF"],
        ["casefold", "abc def"],
    ],
)
def test_str_transform_functions(method, exp):
    ser = pd.Series(["aBc dEF", None], dtype=ArrowDtype(pa.string()))
    result = getattr(ser.str, method)()
    expected = pd.Series([exp, None], dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


def test_str_len():
    ser = pd.Series(["abcd", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.len()
    expected = pd.Series([4, None], dtype=ArrowDtype(pa.int32()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "method, to_strip, val",
    [
        ["strip", None, " abc "],
        ["strip", "x", "xabcx"],
        ["lstrip", None, " abc"],
        ["lstrip", "x", "xabc"],
        ["rstrip", None, "abc "],
        ["rstrip", "x", "abcx"],
    ],
)
def test_str_strip(method, to_strip, val):
    ser = pd.Series([val, None], dtype=ArrowDtype(pa.string()))
    result = getattr(ser.str, method)(to_strip=to_strip)
    expected = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("val", ["abc123", "abc"])
def test_str_removesuffix(val):
    ser = pd.Series([val, None], dtype=ArrowDtype(pa.string()))
    result = ser.str.removesuffix("123")
    expected = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("val", ["123abc", "abc"])
def test_str_removeprefix(val):
    ser = pd.Series([val, None], dtype=ArrowDtype(pa.string()))
    result = ser.str.removeprefix("123")
    expected = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("errors", ["ignore", "strict"])
@pytest.mark.parametrize(
    "encoding, exp",
    [
        ("utf8", {"little": b"abc", "big": "abc"}),
        (
            "utf32",
            {
                "little": b"\xff\xfe\x00\x00a\x00\x00\x00b\x00\x00\x00c\x00\x00\x00",
                "big": b"\x00\x00\xfe\xff\x00\x00\x00a\x00\x00\x00b\x00\x00\x00c",
            },
        ),
    ],
    ids=["utf8", "utf32"],
)
def test_str_encode(errors, encoding, exp):
    ser = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.encode(encoding, errors)
    expected = pd.Series([exp[sys.byteorder], None], dtype=ArrowDtype(pa.binary()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("flags", [0, 2])
def test_str_findall(flags):
    ser = pd.Series(["abc", "efg", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.findall("b", flags=flags)
    expected = pd.Series([["b"], [], None], dtype=ArrowDtype(pa.list_(pa.string())))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("method", ["index", "rindex"])
@pytest.mark.parametrize(
    "start, end",
    [
        [0, None],
        [1, 4],
    ],
)
def test_str_r_index(method, start, end):
    ser = pd.Series(["abcba", None], dtype=ArrowDtype(pa.string()))
    result = getattr(ser.str, method)("c", start, end)
    expected = pd.Series([2, None], dtype=ArrowDtype(pa.int64()))
    tm.assert_series_equal(result, expected)

    with pytest.raises(ValueError, match="substring not found"):
        getattr(ser.str, method)("foo", start, end)


@pytest.mark.parametrize("form", ["NFC", "NFD", "NFKC", "NFKD"])
def test_str_normalize(form):
    # GH#64359 the composing forms (NFC/NFKC) must apply canonical composition,
    #  not just decomposition: "e" + U+0301 must compose to U+00E9, and a
    #  pre-composed character must round-trip unchanged. The compatibility
    #  forms additionally fold e.g. the U+FB01 "fi" ligature.
    data = ["abc", "e\u0301", "\u00e9", "\u212b", "\ufb01", None]
    ser = pd.Series(data, dtype=ArrowDtype(pa.string()))
    result = ser.str.normalize(form)
    expected = pd.Series(
        [unicodedata.normalize(form, val) if val is not None else None for val in data],
        dtype=ArrowDtype(pa.string()),
    )
    tm.assert_series_equal(result, expected)
    if form in ("NFC", "NFKC"):
        # decomposed and pre-composed inputs both yield the single U+00E9
        assert result[1] == "\u00e9"
        assert result[2] == "\u00e9"


@pytest.mark.parametrize(
    "start, end",
    [
        [0, None],
        [1, 4],
    ],
)
def test_str_rfind(start, end):
    ser = pd.Series(["abcba", "foo", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.rfind("c", start, end)
    expected = pd.Series([2, -1, None], dtype=ArrowDtype(pa.int64()))
    tm.assert_series_equal(result, expected)


def test_str_translate():
    ser = pd.Series(["abcba", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.translate({97: "b"})
    expected = pd.Series(["bbcbb", None], dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


def test_str_wrap():
    ser = pd.Series(["abcba", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.wrap(3)
    expected = pd.Series(["abc\nba", None], dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


def test_get_dummies():
    ser = pd.Series(["a|b", None, "a|c"], dtype=ArrowDtype(pa.string()))
    result = ser.str.get_dummies()
    expected = pd.DataFrame(
        [[True, True, False], [False, False, False], [True, False, True]],
        dtype=ArrowDtype(pa.bool_()),
        columns=["a", "b", "c"],
    )
    tm.assert_frame_equal(result, expected)


def test_str_partition():
    ser = pd.Series(["abcba", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.partition("b")
    expected = pd.DataFrame(
        [["a", "b", "cba"], [None, None, None]],
        dtype=ArrowDtype(pa.string()),
        columns=pd.RangeIndex(3),
    )
    tm.assert_frame_equal(result, expected, check_column_type=True)

    result = ser.str.partition("b", expand=False)
    expected = pd.Series(ArrowExtensionArray(pa.array([["a", "b", "cba"], None])))
    tm.assert_series_equal(result, expected)

    result = ser.str.rpartition("b")
    expected = pd.DataFrame(
        [["abc", "b", "a"], [None, None, None]],
        dtype=ArrowDtype(pa.string()),
        columns=pd.RangeIndex(3),
    )
    tm.assert_frame_equal(result, expected, check_column_type=True)

    result = ser.str.rpartition("b", expand=False)
    expected = pd.Series(ArrowExtensionArray(pa.array([["abc", "b", "a"], None])))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("pa_type", [pa.string(), pa.large_string()])
def test_str_partition_chunked(pa_type):
    # GH#63602 each chunk is partitioned on its own, so the result keeps the
    #  input's chunking instead of concatenating into one oversized array
    arr = ArrowExtensionArray(
        pa.chunked_array(
            [pa.array(["abcba"], type=pa_type), pa.array(["a", None], type=pa_type)]
        )
    )
    result = pd.Series(arr).str.partition("b")
    expected = pd.DataFrame(
        [["a", "b", "cba"], ["a", "", ""], [None, None, None]],
        dtype=ArrowDtype(pa.string()),
        columns=pd.RangeIndex(3),
    )
    tm.assert_frame_equal(result, expected, check_column_type=True)

    assert [len(chunk) for chunk in arr._str_partition("b", True)._pa_array.chunks] == [
        1,
        2,
    ]


@pytest.mark.parametrize("method", ["partition", "split"])
@pytest.mark.parametrize("pa_type", [pa.string(), pa.large_string()])
@pytest.mark.parametrize("data", [[], [None, None]], ids=["empty", "all-na"])
def test_str_expand_no_width(data, pa_type, method):
    # GH#63602 no non-null row to take a width from; used to raise
    ser = pd.Series(data, dtype=ArrowDtype(pa_type))
    result = getattr(ser.str, method)("b", expand=True)
    expected = pd.DataFrame(
        {} if len(data) == 0 else {0: data},
        dtype=ArrowDtype(pa_type),
        index=pd.RangeIndex(len(data)),
        columns=pd.RangeIndex(0 if len(data) == 0 else 1),
    )
    tm.assert_frame_equal(result, expected, check_column_type=True)


# (method, args, kwargs, result pa type as a function of the input's own type)
_ELEMENTWISE_STR_FALLBACKS = [
    ("casefold", (), {}, lambda pa_type: pa_type),
    ("normalize", ("NFC",), {}, lambda pa_type: pa_type),
    ("translate", ({97: "b"},), {}, lambda pa_type: pa_type),
    ("wrap", (3,), {}, lambda pa_type: pa_type),
    ("join", ("-",), {}, lambda pa_type: pa_type),
    (
        "encode",
        ("utf-8",),
        {},
        lambda pa_type: (
            pa.large_binary() if pa.types.is_large_string(pa_type) else pa.binary()
        ),
    ),
    ("partition", ("b",), {"expand": False}, lambda pa_type: pa.list_(pa_type)),
    ("rpartition", ("b",), {"expand": False}, lambda pa_type: pa.list_(pa_type)),
    ("findall", ("b",), {}, lambda pa_type: pa.list_(pa_type)),
    ("index", ("b",), {}, lambda pa_type: pa.int64()),
    ("rindex", ("b",), {}, lambda pa_type: pa.int64()),
    ("rfind", ("b",), {}, lambda pa_type: pa.int64()),
]

_elementwise_str_fallback_params = pytest.mark.parametrize(
    "method, args, kwargs, result_pa_type",
    _ELEMENTWISE_STR_FALLBACKS,
    ids=[entry[0] for entry in _ELEMENTWISE_STR_FALLBACKS],
)


@_elementwise_str_fallback_params
@pytest.mark.parametrize(
    "chunks",
    [[], [[None, None]], [[None], ["abcba"]]],
    ids=["no-chunks", "all-na", "na-chunk-first"],
)
def test_str_elementwise_fallback_degenerate_chunks(
    method, args, kwargs, result_pa_type, chunks
):
    # GH#66706 the elementwise fallbacks rebuild the result with an explicit
    #  type, so chunkings that give pyarrow nothing to infer from -- no chunks
    #  at all, or a leading all-null chunk -- no longer raise or come back as
    #  null[pyarrow]
    arr = pa.chunked_array(
        [pa.array(chunk, type=pa.string()) for chunk in chunks], type=pa.string()
    )
    result = getattr(pd.Series(ArrowExtensionArray(arr)).str, method)(*args, **kwargs)
    assert result.dtype == ArrowDtype(result_pa_type(pa.string()))

    data = [val for chunk in chunks for val in chunk]
    unchunked = pd.Series(data, dtype=ArrowDtype(pa.string()))
    expected = getattr(unchunked.str, method)(*args, **kwargs)
    tm.assert_series_equal(result, expected)


@_elementwise_str_fallback_params
def test_str_elementwise_fallback_keeps_large_string(
    method, args, kwargs, result_pa_type
):
    # GH#66221 a large_string input must not silently come back as string; the
    #  string-valued results keep large_string, encode gives the matching
    #  large_binary, and the list-valued results nest large_string just like
    #  the native pc.split_pattern kernels do
    ser = pd.Series(["abcba", None], dtype=ArrowDtype(pa.large_string()))
    result = getattr(ser.str, method)(*args, **kwargs)
    assert result.dtype == ArrowDtype(result_pa_type(pa.large_string()))


@pytest.mark.parametrize("method", ["rsplit", "split"])
def test_str_split_pat_none(method):
    # GH 56271
    ser = pd.Series(["a1 cbc\nb", None], dtype=ArrowDtype(pa.string()))
    result = getattr(ser.str, method)()
    expected = pd.Series(ArrowExtensionArray(pa.array([["a1", "cbc", "b"], None])))
    tm.assert_series_equal(result, expected)


def test_str_split():
    # GH 52401
    ser = pd.Series(["a1cbcb", "a2cbcb", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.split("c")
    expected = pd.Series(
        ArrowExtensionArray(pa.array([["a1", "b", "b"], ["a2", "b", "b"], None]))
    )
    tm.assert_series_equal(result, expected)

    result = ser.str.split("c", n=1)
    expected = pd.Series(
        ArrowExtensionArray(pa.array([["a1", "bcb"], ["a2", "bcb"], None]))
    )
    tm.assert_series_equal(result, expected)

    result = ser.str.split("[1-2]", regex=True)
    expected = pd.Series(
        ArrowExtensionArray(pa.array([["a", "cbcb"], ["a", "cbcb"], None]))
    )
    tm.assert_series_equal(result, expected)

    result = ser.str.split("[1-2]", regex=True, expand=True)
    expected = pd.DataFrame(
        {
            0: ArrowExtensionArray(pa.array(["a", "a", None])),
            1: ArrowExtensionArray(pa.array(["cbcb", "cbcb", None])),
        }
    )
    tm.assert_frame_equal(result, expected)

    result = ser.str.split("1", expand=True)
    expected = pd.DataFrame(
        {
            0: ArrowExtensionArray(pa.array(["a", "a2cbcb", None])),
            1: ArrowExtensionArray(pa.array(["cbcb", None, None])),
        }
    )
    tm.assert_frame_equal(result, expected)


def test_str_rsplit():
    # GH 52401
    ser = pd.Series(["a1cbcb", "a2cbcb", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.rsplit("c")
    expected = pd.Series(
        ArrowExtensionArray(pa.array([["a1", "b", "b"], ["a2", "b", "b"], None]))
    )
    tm.assert_series_equal(result, expected)

    result = ser.str.rsplit("c", n=1)
    expected = pd.Series(
        ArrowExtensionArray(pa.array([["a1cb", "b"], ["a2cb", "b"], None]))
    )
    tm.assert_series_equal(result, expected)

    result = ser.str.rsplit("c", n=1, expand=True)
    expected = pd.DataFrame(
        {
            0: ArrowExtensionArray(pa.array(["a1cb", "a2cb", None])),
            1: ArrowExtensionArray(pa.array(["b", "b", None])),
        }
    )
    tm.assert_frame_equal(result, expected)

    result = ser.str.rsplit("1", expand=True)
    expected = pd.DataFrame(
        {
            0: ArrowExtensionArray(pa.array(["a", "a2cbcb", None])),
            1: ArrowExtensionArray(pa.array(["cbcb", None, None])),
        }
    )
    tm.assert_frame_equal(result, expected)


def test_str_extract_non_symbolic():
    # GH#63683
    ser = pd.Series(["a1", "b2", "c3"], dtype=ArrowDtype(pa.string()))
    result = ser.str.extract(r"[ab](\d)")
    expected = pd.DataFrame({0: ["1", "2", None]}, dtype=ArrowDtype(pa.string()))
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("expand", [True, False])
def test_str_extract(expand):
    ser = pd.Series(["a1", "b2", "c3"], dtype=ArrowDtype(pa.string()))
    result = ser.str.extract(r"(?P<letter>[ab])(?P<digit>\d)", expand=expand)
    expected = pd.DataFrame(
        {
            "letter": ArrowExtensionArray(pa.array(["a", "b", None])),
            "digit": ArrowExtensionArray(pa.array(["1", "2", None])),
        }
    )
    tm.assert_frame_equal(result, expected)


def test_str_extract_expand():
    ser = pd.Series(["a1", "b2", "c3"], dtype=ArrowDtype(pa.string()))
    result = ser.str.extract(r"[ab](?P<digit>\d)", expand=True)
    expected = pd.DataFrame(
        {
            "digit": ArrowExtensionArray(pa.array(["1", "2", None])),
        }
    )
    tm.assert_frame_equal(result, expected)

    result = ser.str.extract(r"[ab](?P<digit>\d)", expand=False)
    expected = pd.Series(ArrowExtensionArray(pa.array(["1", "2", None])), name="digit")
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("expand", [True, False])
def test_str_extract_flags(expand):
    # GH#66348
    ser = pd.Series(["a1", "A2", "b3"], dtype=ArrowDtype(pa.string()))
    result = ser.str.extract(r"a(?P<digit>\d)", flags=re.IGNORECASE, expand=expand)
    expected = pd.Series(ArrowExtensionArray(pa.array(["1", "2", None])), name="digit")
    if expand:
        expected = expected.to_frame()
    tm.assert_equal(result, expected)
