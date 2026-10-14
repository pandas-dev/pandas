"""
The pyarrow-backed dtypes evaluate regular expressions with RE2 (via pyarrow)
where they can, and must give the same results as Python's ``re`` does for
object dtype. GH#63683
"""

import re

import pytest

import pandas as pd
import pandas._testing as tm

pa = pytest.importorskip("pyarrow")


@pytest.fixture(
    params=[
        pd.StringDtype("pyarrow", na_value=float("nan")),
        pd.StringDtype("pyarrow"),
        pd.ArrowDtype(pa.string()),
    ],
    ids=["str[pyarrow]", "string[pyarrow]", "string[pyarrow-arrowdtype]"],
)
def arrow_string_dtype(request):
    return request.param


def assert_same_as_object(data, dtype, method, *args, **kwargs):
    expected_ser = pd.Series(data, dtype=object)
    ser = pd.Series(data, dtype=dtype)
    try:
        expected = getattr(expected_ser.str, method)(*args, **kwargs)
    except Exception as err:
        with pytest.raises(type(err)):
            getattr(ser.str, method)(*args, **kwargs)
        return
    result = getattr(ser.str, method)(*args, **kwargs)
    # compare values only: the result dtypes and missing value markers differ
    tm.assert_equal(
        result.astype(object).mask(result.isna(), None),
        expected.astype(object).mask(expected.isna(), None),
    )


@pytest.mark.parametrize(
    "data, method, args, kwargs",
    [
        # \w \W \d \D \s \S \b \B match Unicode characters
        (["é", "東", "٣"], "contains", (r"\w",), {}),
        (["é", "東", "٣"], "fullmatch", (r"\W",), {}),
        (["٣", "\uff11"], "contains", (r"\d",), {}),
        (["٣", "\uff11"], "fullmatch", (r"\D",), {}),
        (["\xa0", "\x0b", "\u3000", "\x1c"], "contains", (r"\s",), {}),
        (["\xa0", "\x0b", "\u3000", "\x1c"], "fullmatch", (r"\S",), {}),
        (["é", "ab cd", "é b"], "count", (r"\b",), {}),
        (["é", "é b"], "count", (r"\B",), {}),
        (["é", "a b", ""], "contains", (r"\B",), {}),
        (["é"], "match", (r"[\w]",), {}),
        (["é"], "match", (r"[^\W\d]",), {}),
        (["café crème!"], "replace", (r"[^\w\s]", ""), {"regex": True}),
        (["é"], "contains", (re.compile(r"\w"),), {}),
        (["é", "a"], "contains", (r"(?a)\w",), {}),
        # $ also matches before a trailing newline
        (["a\n", "a"], "contains", ("a$",), {}),
        (["a\n", "a"], "match", ("a$",), {}),
        (["a\n", "a"], "fullmatch", ("a$",), {}),
        (["a\n", "a"], "replace", ("$", "X"), {"regex": True}),
        (["a\n", "a\nb"], "count", ("(?m)$",), {}),
        # case-insensitive matching folds as re does
        (["İ", "\u0131", "I", "K", "\u017f"], "contains", ("i",), {"case": False}),
        (["İ", "\u0131"], "match", ("i",), {"case": False}),
        (["İ", "\u0131", "\u212a"], "count", ("(?i)[ik]",), {}),
        (["İ"], "contains", ("(?i)[a-z]",), {}),
        (["É"], "contains", (r"\w",), {"case": False}),
        (["\u212a", "k"], "contains", ("(?ai)k",), {}),
        (["Straße", "STRASSE"], "replace", ("ß", "ss"), {"case": False}),
        # {,n} quantifier
        (["aab", "b"], "contains", ("a{,2}b",), {}),
        # no POSIX bracket classes
        (["a", "[:a]"], "contains", ("[[:alpha:]]",), {}),
        (["a", "]"], "contains", ("[[:^alpha:]]",), {}),
        # replacement templates
        (["ab"], "replace", ("(a)(b)", r"\0"), {"regex": True}),
        (["abcdefghij"], "replace", ("(a)(b)(c)(d)(e)(f)(g)(h)(i)(j)", r"\10"), {}),
        (["ab"], "replace", ("(a)(b)", r"\2\1\g<1>\\"), {"regex": True}),
        (["ab"], "replace", ("(?P<x>a)", r"[\g<x>]"), {"regex": True}),
        (["a"], "replace", ("a", r"\n\t"), {"regex": True}),
        (["ab"], "replace", ("a", r"\g<1>"), {"regex": False}),
        (["ab"], "replace", ("(a)|(c)", r"[\2]"), {"regex": True}),
        # zero-width matches
        (["ab"], "count", ("^",), {}),
        (["ab"], "count", ("$",), {}),
        (["a\nb"], "count", ("(?m)^",), {}),
        (["東京"], "count", ("x*",), {}),
        (["baaa"], "replace", ("a*", "-"), {"regex": True}),
        (["baaa"], "replace", ("a?", "-"), {"regex": True}),
        (["abc"], "replace", ("x*", "-"), {"n": 2, "regex": True}),
        (["baaa"], "replace", ("a?", "-"), {"n": 2, "regex": True}),
        (["ab ab ab"], "replace", (r"\bab", "X"), {"n": 2, "regex": True}),
        (["ab ab"], "count", (r"\bab",), {}),
        # Python syntax that RE2 does not have
        (["é"], "contains", (r"é",), {}),
        (["é"], "contains", (r"[à-ÿ]",), {}),
        (["é"], "contains", (r"\U000000e9",), {}),
        (["é"], "contains", (r"\N{LATIN SMALL LETTER E WITH ACUTE}",), {}),
        (["é"], "contains", (r"\é",), {}),
        (["\x08"], "contains", (r"[\b]",), {}),
        (["\x01"], "contains", (r"[\1]",), {}),
        (["aab"], "contains", ("a*+b",), {}),
        (["aab"], "contains", ("(?>a+)b",), {}),
        (["a" * 1001, "a"], "contains", ("a{1001}",), {}),
        (["a"], "contains", ("(?:a{1000}){2}",), {}),
        (["a"], "contains", (r"\w{3000}",), {}),
        (["c", "ab"], "contains", ("(a)?(?(1)b|c)",), {}),
        (["a"], "contains", ("(?#note)a",), {}),
        (["ab"], "contains", ("(?x) a b # comment",), {}),
        (["a"], "contains", ("(?u)a",), {}),
        (["a\nb"], "contains", ("(?s)a.b",), {}),
        (["a\nb"], "contains", ("a.b",), {}),
        # RE2 syntax that Python does not have
        (["é"], "contains", (r"\pL",), {}),
        (["\u03b1"], "contains", (r"\p{Greek}",), {}),
        (["é"], "contains", (r"\x{e9}",), {}),
        (["a.*"], "contains", (r"\Q.*\E",), {}),
        (["a"], "contains", (r"a\z",), {}),
        (["a"], "contains", ("(?<n>a)",), {}),
        (["aB"], "contains", ("a(?i)b",), {}),
        (["aa"], "contains", ("(?U)a+",), {}),
        (["-"], "contains", (r"[\d-z]",), {}),
        (["ab"], "contains", ("(?P<n>a)(?P<n>b)",), {}),
        # invalid patterns and templates
        (["a"], "contains", ("(",), {}),
        (["a"], "contains", ("a**",), {}),
        (["a"], "contains", (r"\9",), {}),
        (["a"], "replace", ("a", r"\q"), {"regex": True}),
        (["a"], "replace", ("(a)", r"\9"), {"regex": True}),
        # split and extract
        (["café crème"], "split", (r"\w+",), {"regex": True}),
        (["café crème"], "split", (r"(é)",), {"regex": True}),
        (["a\n"], "split", ("a$",), {"regex": True}),
        (["baaa"], "split", ("(a*)",), {"regex": True}),
        (["a b  c"], "split", (r"\s+",), {"regex": True, "n": 1}),
        (["café"], "extract", (r"(\w+)",), {"expand": False}),
        (["٣ 4"], "extract", (r"(?P<n>\d)",), {"expand": False}),
        (["ab", "b", "x"], "extract", (r"(a)?(b)",), {"expand": True}),
        (["ab", "x"], "extract", (r"(?P<one>a)(b)",), {"expand": True}),
        (["aab"], "extract", (r"(a*)+b",), {"expand": False}),
    ],
)
def test_regex_same_as_object(arrow_string_dtype, data, method, args, kwargs):
    if method in ("split", "extract") and isinstance(
        arrow_string_dtype, pd.StringDtype
    ):
        pytest.skip("evaluated with re for StringDtype")
    assert_same_as_object(data, arrow_string_dtype, method, *args, **kwargs)
