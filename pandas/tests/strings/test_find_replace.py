from datetime import datetime
import re

import numpy as np
import pytest

import pandas as pd
from pandas import (
    Index,
    Series,
    _testing as tm,
)


def test_contains():
    values = np.array(
        ["foo", np.nan, "fooommm__foo", "mmm_", "foommm[_]+bar"], dtype=np.object_
    )
    values = Series(values)
    pat = "mmm[_]+"

    result = values.str.contains(pat)
    expected = Series(np.array([False, np.nan, True, True, False], dtype=np.object_))
    tm.assert_series_equal(result, expected)

    result = values.str.contains(pat, regex=False)
    expected = Series(np.array([False, np.nan, False, False, True], dtype=np.object_))
    tm.assert_series_equal(result, expected)

    values = Series(np.array(["foo", "xyz", "fooommm__foo", "mmm_"], dtype=object))
    result = values.str.contains(pat)
    expected = Series(np.array([False, False, True, True]))
    assert result.dtype == np.bool_
    tm.assert_series_equal(result, expected)

    # case insensitive using regex
    values = Series(np.array(["Foo", "xYz", "fOOomMm__fOo", "MMM_"], dtype=object))
    result = values.str.contains("FOO|mmm", case=False)
    expected = Series(np.array([True, False, True, True]))
    tm.assert_series_equal(result, expected)

    # case insensitive without regex
    result = Series(values).str.contains("foo", regex=False, case=False)
    expected = Series(np.array([True, False, True, False]))
    tm.assert_series_equal(result, expected)

    # mixed
    mixed = Series(
        np.array(
            ["a", np.nan, "b", True, datetime.today(), "foo", None, 1, 2.0],
            dtype=object,
        )
    )
    rs = mixed.str.contains("o")
    xp = Series(
        np.array(
            [False, np.nan, False, np.nan, np.nan, True, np.nan, np.nan, np.nan],
            dtype=np.object_,
        )
    )
    tm.assert_series_equal(rs, xp)

    rs = mixed.str.contains("o")
    xp = Series([False, np.nan, False, np.nan, np.nan, True, np.nan, np.nan, np.nan])
    assert isinstance(rs, Series)
    tm.assert_series_equal(rs, xp)

    # unicode
    values = Series(np.array(["foo", np.nan, "fooommm__foo", "mmm_"], dtype=np.object_))
    pat = "mmm[_]+"

    result = values.str.contains(pat)
    expected = Series(np.array([False, np.nan, True, True], dtype=np.object_))
    tm.assert_series_equal(result, expected)

    result = values.str.contains(pat, na=False)
    expected = Series(np.array([False, False, True, True]))
    tm.assert_series_equal(result, expected)

    values = Series(np.array(["foo", "xyz", "fooommm__foo", "mmm_"], dtype=np.object_))
    result = values.str.contains(pat)
    expected = Series(np.array([False, False, True, True]))
    assert result.dtype == np.bool_
    tm.assert_series_equal(result, expected)


def test_contains_for_object_category():
    # gh 22158

    # na for category
    values = Series(["a", "b", "c", "a", np.nan], dtype="category")
    result = values.str.contains("a", na=True)
    expected = Series([True, False, False, True, True])
    tm.assert_series_equal(result, expected)

    result = values.str.contains("a", na=False)
    expected = Series([True, False, False, True, False])
    tm.assert_series_equal(result, expected)

    # na for objects
    values = Series(["a", "b", "c", "a", np.nan])
    result = values.str.contains("a", na=True)
    expected = Series([True, False, False, True, True])
    tm.assert_series_equal(result, expected)

    result = values.str.contains("a", na=False)
    expected = Series([True, False, False, True, False])
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("dtype", [None, "category"])
@pytest.mark.parametrize("null_value", [None, np.nan, pd.NA])
@pytest.mark.parametrize("na", [True, False])
def test_startswith(dtype, null_value, na):
    # add category dtype parametrizations for GH-36241
    values = Series(
        ["om", null_value, "foo_nom", "nom", "bar_foo", null_value, "foo"],
        dtype=dtype,
    )

    result = values.str.startswith("foo")
    exp = Series([False, np.nan, True, False, False, np.nan, True])
    tm.assert_series_equal(result, exp)

    result = values.str.startswith("foo", na=na)
    exp = Series([False, na, True, False, False, na, True])
    tm.assert_series_equal(result, exp)

    # mixed
    mixed = np.array(
        ["a", np.nan, "b", True, datetime.today(), "foo", None, 1, 2.0],
        dtype=np.object_,
    )
    rs = Series(mixed).str.startswith("f")
    xp = Series([False, np.nan, False, np.nan, np.nan, True, np.nan, np.nan, np.nan])
    tm.assert_series_equal(rs, xp)


@pytest.mark.parametrize("dtype", [None, "category"])
@pytest.mark.parametrize("null_value", [None, np.nan, pd.NA])
@pytest.mark.parametrize("na", [True, False])
def test_endswith(dtype, null_value, na):
    # add category dtype parametrizations for GH-36241
    values = Series(
        ["om", null_value, "foo_nom", "nom", "bar_foo", null_value, "foo"],
        dtype=dtype,
    )

    result = values.str.endswith("foo")
    exp = Series([False, np.nan, False, False, True, np.nan, True])
    tm.assert_series_equal(result, exp)

    result = values.str.endswith("foo", na=na)
    exp = Series([False, na, False, False, True, na, True])
    tm.assert_series_equal(result, exp)

    # mixed
    mixed = np.array(
        ["a", np.nan, "b", True, datetime.today(), "foo", None, 1, 2.0],
        dtype=object,
    )
    rs = Series(mixed).str.endswith("f")
    xp = Series([False, np.nan, False, np.nan, np.nan, False, np.nan, np.nan, np.nan])
    tm.assert_series_equal(rs, xp)


def test_replace():
    values = Series(["fooBAD__barBAD", np.nan])

    result = values.str.replace("BAD[_]*", "", regex=True)
    exp = Series(["foobar", np.nan])
    tm.assert_series_equal(result, exp)

    result = values.str.replace("BAD[_]*", "", n=1, regex=True)
    exp = Series(["foobarBAD", np.nan])
    tm.assert_series_equal(result, exp)

    # mixed
    mixed = Series(
        ["aBAD", np.nan, "bBAD", True, datetime.today(), "fooBAD", None, 1, 2.0]
    )

    rs = Series(mixed).str.replace("BAD[_]*", "", regex=True)
    xp = Series(["a", np.nan, "b", np.nan, np.nan, "foo", np.nan, np.nan, np.nan])
    assert isinstance(rs, Series)
    tm.assert_almost_equal(rs, xp)

    # flags + unicode
    values = Series([b"abcd,\xc3\xa0".decode("utf-8")])
    exp = Series([b"abcd, \xc3\xa0".decode("utf-8")])
    result = values.str.replace(r"(?<=\w),(?=\w)", ", ", flags=re.UNICODE, regex=True)
    tm.assert_series_equal(result, exp)

    # GH 13438
    msg = "repl must be a string or callable"
    for klass in (Series, Index):
        for repl in (None, 3, {"a": "b"}):
            for data in (["a", "b", None], ["a", "b", "c", "ad"]):
                values = klass(data)
                with pytest.raises(TypeError, match=msg):
                    values.str.replace("a", repl)


def test_replace_callable():
    # GH 15055
    values = Series(["fooBAD__barBAD", np.nan])

    # test with callable
    repl = lambda m: m.group(0).swapcase()
    result = values.str.replace("[a-z][A-Z]{2}", repl, n=2, regex=True)
    exp = Series(["foObaD__baRbaD", np.nan])
    tm.assert_series_equal(result, exp)

    # test with wrong number of arguments, raising an error
    p_err = (
        r"((takes)|(missing)) (?(2)from \d+ to )?\d+ "
        r"(?(3)required )positional arguments?"
    )

    repl = lambda: None
    with pytest.raises(TypeError, match=p_err):
        values.str.replace("a", repl)

    repl = lambda m, x: None
    with pytest.raises(TypeError, match=p_err):
        values.str.replace("a", repl)

    repl = lambda m, x, y=None: None
    with pytest.raises(TypeError, match=p_err):
        values.str.replace("a", repl)

    # test regex named groups
    values = Series(["Foo Bar Baz", np.nan])
    pat = r"(?P<first>\w+) (?P<middle>\w+) (?P<last>\w+)"
    repl = lambda m: m.group("middle").swapcase()
    result = values.str.replace(pat, repl, regex=True)
    exp = Series(["bAR", np.nan])
    tm.assert_series_equal(result, exp)


def test_replace_compiled_regex():
    # GH 15446
    values = Series(["fooBAD__barBAD", np.nan])

    # test with compiled regex
    pat = re.compile(r"BAD_*")
    result = values.str.replace(pat, "", regex=True)
    exp = Series(["foobar", np.nan])
    tm.assert_series_equal(result, exp)

    result = values.str.replace(pat, "", n=1, regex=True)
    exp = Series(["foobarBAD", np.nan])
    tm.assert_series_equal(result, exp)

    # mixed
    mixed = Series(
        ["aBAD", np.nan, "bBAD", True, datetime.today(), "fooBAD", None, 1, 2.0]
    )

    rs = Series(mixed).str.replace(pat, "", regex=True)
    xp = Series(["a", np.nan, "b", np.nan, np.nan, "foo", np.nan, np.nan, np.nan])
    assert isinstance(rs, Series)
    tm.assert_almost_equal(rs, xp)

    # flags + unicode
    values = Series([b"abcd,\xc3\xa0".decode("utf-8")])
    exp = Series([b"abcd, \xc3\xa0".decode("utf-8")])
    pat = re.compile(r"(?<=\w),(?=\w)", flags=re.UNICODE)
    result = values.str.replace(pat, ", ")
    tm.assert_series_equal(result, exp)

    # case and flags provided to str.replace will have no effect
    # and will produce warnings
    values = Series(["fooBAD__barBAD__bad", np.nan])
    pat = re.compile(r"BAD_*")

    with pytest.raises(ValueError, match="case and flags cannot be"):
        result = values.str.replace(pat, "", flags=re.IGNORECASE)

    with pytest.raises(ValueError, match="case and flags cannot be"):
        result = values.str.replace(pat, "", case=False)

    with pytest.raises(ValueError, match="case and flags cannot be"):
        result = values.str.replace(pat, "", case=True)

    # test with callable
    values = Series(["fooBAD__barBAD", np.nan])
    repl = lambda m: m.group(0).swapcase()
    pat = re.compile("[a-z][A-Z]{2}")
    result = values.str.replace(pat, repl, n=2)
    exp = Series(["foObaD__baRbaD", np.nan])
    tm.assert_series_equal(result, exp)


def test_replace_literal():
    # GH16808 literal replace (regex=False vs regex=True)
    values = Series(["f.o", "foo", np.nan])
    exp = Series(["bao", "bao", np.nan])
    result = values.str.replace("f.", "ba", regex=True)
    tm.assert_series_equal(result, exp)

    exp = Series(["bao", "foo", np.nan])
    result = values.str.replace("f.", "ba", regex=False)
    tm.assert_series_equal(result, exp)

    # Cannot do a literal replace if given a callable repl or compiled
    # pattern
    callable_repl = lambda m: m.group(0).swapcase()
    compiled_pat = re.compile("[a-z][A-Z]{2}")

    msg = "Cannot use a callable replacement when regex=False"
    with pytest.raises(ValueError, match=msg):
        values.str.replace("abc", callable_repl, regex=False)

    msg = "Cannot use a compiled regex as replacement pattern with regex=False"
    with pytest.raises(ValueError, match=msg):
        values.str.replace(compiled_pat, "", regex=False)


def test_match():
    # New match behavior introduced in 0.13
    values = Series(["fooBAD__barBAD", np.nan, "foo"])
    result = values.str.match(".*(BAD[_]+).*(BAD)")
    exp = Series([True, np.nan, False])
    tm.assert_series_equal(result, exp)

    values = Series(["fooBAD__barBAD", "BAD_BADleroybrown", np.nan, "foo"])
    result = values.str.match(".*BAD[_]+.*BAD")
    exp = Series([True, True, np.nan, False])
    tm.assert_series_equal(result, exp)

    # mixed
    mixed = Series(
        [
            "aBAD_BAD",
            np.nan,
            "BAD_b_BAD",
            True,
            datetime.today(),
            "foo",
            None,
            1,
            2.0,
        ]
    )
    rs = Series(mixed).str.match(".*(BAD[_]+).*(BAD)")
    xp = Series([True, np.nan, True, np.nan, np.nan, False, np.nan, np.nan, np.nan])
    assert isinstance(rs, Series)
    tm.assert_series_equal(rs, xp)

    # na GH #6609
    res = Series(["a", 0, np.nan]).str.match("a", na=False)
    exp = Series([True, False, False])
    tm.assert_series_equal(exp, res)
    res = Series(["a", 0, np.nan]).str.match("a")
    exp = Series([True, np.nan, np.nan])
    tm.assert_series_equal(exp, res)

    values = Series(["ab", "AB", "abc", "ABC"])
    result = values.str.match("ab", case=False)
    expected = Series([True, True, True, True])
    tm.assert_series_equal(result, expected)


def test_fullmatch():
    # GH 32806
    values = Series(["fooBAD__barBAD", "BAD_BADleroybrown", np.nan, "foo"])
    result = values.str.fullmatch(".*BAD[_]+.*BAD")
    exp = Series([True, False, np.nan, False])
    tm.assert_series_equal(result, exp)

    # Make sure that the new string arrays work
    string_values = Series(
        ["fooBAD__barBAD", "BAD_BADleroybrown", np.nan, "foo"], dtype="string"
    )
    result = string_values.str.fullmatch(".*BAD[_]+.*BAD")
    # Result is nullable boolean with StringDtype
    string_exp = Series([True, False, np.nan, False], dtype="boolean")
    tm.assert_series_equal(result, string_exp)

    values = Series(["ab", "AB", "abc", "ABC"])
    result = values.str.fullmatch("ab", case=False)
    expected = Series([True, True, False, False])
    tm.assert_series_equal(result, expected)


def test_findall():
    values = Series(["fooBAD__barBAD", np.nan, "foo", "BAD"])

    result = values.str.findall("BAD[_]*")
    exp = Series([["BAD__", "BAD"], np.nan, [], ["BAD"]])
    tm.assert_almost_equal(result, exp)

    # mixed
    mixed = Series(
        [
            "fooBAD__barBAD",
            np.nan,
            "foo",
            True,
            datetime.today(),
            "BAD",
            None,
            1,
            2.0,
        ]
    )

    rs = Series(mixed).str.findall("BAD[_]*")
    xp = Series(
        [
            ["BAD__", "BAD"],
            np.nan,
            [],
            np.nan,
            np.nan,
            ["BAD"],
            np.nan,
            np.nan,
            np.nan,
        ]
    )

    assert isinstance(rs, Series)
    tm.assert_almost_equal(rs, xp)


def test_find():
    values = Series(["ABCDEFG", "BCDEFEF", "DEFGHIJEF", "EFGHEF", "XXXX"])
    result = values.str.find("EF")
    tm.assert_series_equal(result, Series([4, 3, 1, 0, -1]))
    expected = np.array([v.find("EF") for v in values.values], dtype=np.int64)
    tm.assert_numpy_array_equal(result.values, expected)

    result = values.str.rfind("EF")
    tm.assert_series_equal(result, Series([4, 5, 7, 4, -1]))
    expected = np.array([v.rfind("EF") for v in values.values], dtype=np.int64)
    tm.assert_numpy_array_equal(result.values, expected)

    result = values.str.find("EF", 3)
    tm.assert_series_equal(result, Series([4, 3, 7, 4, -1]))
    expected = np.array([v.find("EF", 3) for v in values.values], dtype=np.int64)
    tm.assert_numpy_array_equal(result.values, expected)

    result = values.str.rfind("EF", 3)
    tm.assert_series_equal(result, Series([4, 5, 7, 4, -1]))
    expected = np.array([v.rfind("EF", 3) for v in values.values], dtype=np.int64)
    tm.assert_numpy_array_equal(result.values, expected)

    result = values.str.find("EF", 3, 6)
    tm.assert_series_equal(result, Series([4, 3, -1, 4, -1]))
    expected = np.array([v.find("EF", 3, 6) for v in values.values], dtype=np.int64)
    tm.assert_numpy_array_equal(result.values, expected)

    result = values.str.rfind("EF", 3, 6)
    tm.assert_series_equal(result, Series([4, 3, -1, 4, -1]))
    expected = np.array([v.rfind("EF", 3, 6) for v in values.values], dtype=np.int64)
    tm.assert_numpy_array_equal(result.values, expected)

    with pytest.raises(TypeError, match="expected a string object, not int"):
        result = values.str.find(0)

    with pytest.raises(TypeError, match="expected a string object, not int"):
        result = values.str.rfind(0)


def test_find_nan():
    values = Series(["ABCDEFG", np.nan, "DEFGHIJEF", np.nan, "XXXX"])
    result = values.str.find("EF")
    tm.assert_series_equal(result, Series([4, np.nan, 1, np.nan, -1]))

    result = values.str.rfind("EF")
    tm.assert_series_equal(result, Series([4, np.nan, 7, np.nan, -1]))

    result = values.str.find("EF", 3)
    tm.assert_series_equal(result, Series([4, np.nan, 7, np.nan, -1]))

    result = values.str.rfind("EF", 3)
    tm.assert_series_equal(result, Series([4, np.nan, 7, np.nan, -1]))

    result = values.str.find("EF", 3, 6)
    tm.assert_series_equal(result, Series([4, np.nan, -1, np.nan, -1]))

    result = values.str.rfind("EF", 3, 6)
    tm.assert_series_equal(result, Series([4, np.nan, -1, np.nan, -1]))


def test_translate():
    def _check(result, expected):
        if isinstance(result, Series):
            tm.assert_series_equal(result, expected)
        else:
            tm.assert_index_equal(result, expected)

    for klass in [Series, Index]:
        s = klass(["abcdefg", "abcc", "cdddfg", "cdefggg"])
        table = str.maketrans("abc", "cde")
        result = s.str.translate(table)
        expected = klass(["cdedefg", "cdee", "edddfg", "edefggg"])
        _check(result, expected)

    # Series with non-string values
    s = Series(["a", "b", "c", 1.2])
    expected = Series(["c", "d", "e", np.nan])
    result = s.str.translate(table)
    tm.assert_series_equal(result, expected)


def test_contains_moar():
    # PR #1179
    s = Series(["A", "B", "C", "Aaba", "Baca", "", np.nan, "CABA", "dog", "cat"])

    result = s.str.contains("a")
    expected = Series(
        [False, False, False, True, True, False, np.nan, False, False, True]
    )
    tm.assert_series_equal(result, expected)

    result = s.str.contains("a", case=False)
    expected = Series(
        [True, False, False, True, True, False, np.nan, True, False, True]
    )
    tm.assert_series_equal(result, expected)

    result = s.str.contains("Aa")
    expected = Series(
        [False, False, False, True, False, False, np.nan, False, False, False]
    )
    tm.assert_series_equal(result, expected)

    result = s.str.contains("ba")
    expected = Series(
        [False, False, False, True, False, False, np.nan, False, False, False]
    )
    tm.assert_series_equal(result, expected)

    result = s.str.contains("ba", case=False)
    expected = Series(
        [False, False, False, True, True, False, np.nan, True, False, False]
    )
    tm.assert_series_equal(result, expected)


def test_contains_nan():
    # PR #14171
    s = Series([np.nan, np.nan, np.nan], dtype=np.object_)

    result = s.str.contains("foo", na=False)
    expected = Series([False, False, False], dtype=np.bool_)
    tm.assert_series_equal(result, expected)

    result = s.str.contains("foo", na=True)
    expected = Series([True, True, True], dtype=np.bool_)
    tm.assert_series_equal(result, expected)

    result = s.str.contains("foo", na="foo")
    expected = Series(["foo", "foo", "foo"], dtype=np.object_)
    tm.assert_series_equal(result, expected)

    result = s.str.contains("foo")
    expected = Series([np.nan, np.nan, np.nan], dtype=np.object_)
    tm.assert_series_equal(result, expected)


def test_replace_moar():
    # PR #1179
    s = Series(["A", "B", "C", "Aaba", "Baca", "", np.nan, "CABA", "dog", "cat"])

    result = s.str.replace("A", "YYY")
    expected = Series(
        ["YYY", "B", "C", "YYYaba", "Baca", "", np.nan, "CYYYBYYY", "dog", "cat"]
    )
    tm.assert_series_equal(result, expected)

    result = s.str.replace("A", "YYY", case=False)
    expected = Series(
        [
            "YYY",
            "B",
            "C",
            "YYYYYYbYYY",
            "BYYYcYYY",
            "",
            np.nan,
            "CYYYBYYY",
            "dog",
            "cYYYt",
        ]
    )
    tm.assert_series_equal(result, expected)

    result = s.str.replace("^.a|dog", "XX-XX ", case=False, regex=True)
    expected = Series(
        [
            "A",
            "B",
            "C",
            "XX-XX ba",
            "XX-XX ca",
            "",
            np.nan,
            "XX-XX BA",
            "XX-XX ",
            "XX-XX t",
        ]
    )
    tm.assert_series_equal(result, expected)


def test_match_findall_flags():
    data = {
        "Dave": "dave@google.com",
        "Steve": "steve@gmail.com",
        "Rob": "rob@gmail.com",
        "Wes": np.nan,
    }
    data = Series(data)

    pat = r"([A-Z0-9._%+-]+)@([A-Z0-9.-]+)\.([A-Z]{2,4})"

    result = data.str.extract(pat, flags=re.IGNORECASE, expand=True)
    assert result.iloc[0].tolist() == ["dave", "google", "com"]

    result = data.str.match(pat, flags=re.IGNORECASE)
    assert result[0]

    result = data.str.fullmatch(pat, flags=re.IGNORECASE)
    assert result[0]

    result = data.str.findall(pat, flags=re.IGNORECASE)
    assert result[0][0] == ("dave", "google", "com")

    result = data.str.count(pat, flags=re.IGNORECASE)
    assert result[0] == 1

    with tm.assert_produces_warning(UserWarning):
        result = data.str.contains(pat, flags=re.IGNORECASE)
    assert result[0]
