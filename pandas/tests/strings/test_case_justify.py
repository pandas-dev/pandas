from datetime import datetime

import numpy as np
import pytest

from pandas import (
    Series,
    _testing as tm,
)


def test_title():
    values = Series(["FOO", "BAR", np.nan, "Blah", "blurg"])

    result = values.str.title()
    exp = Series(["Foo", "Bar", np.nan, "Blah", "Blurg"])
    tm.assert_series_equal(result, exp)

    # mixed
    mixed = Series(["FOO", np.nan, "bar", True, datetime.today(), "blah", None, 1, 2.0])
    mixed = mixed.str.title()
    exp = Series(["Foo", np.nan, "Bar", np.nan, np.nan, "Blah", np.nan, np.nan, np.nan])
    tm.assert_almost_equal(mixed, exp)


def test_lower_upper():
    values = Series(["om", np.nan, "nom", "nom"])

    result = values.str.upper()
    exp = Series(["OM", np.nan, "NOM", "NOM"])
    tm.assert_series_equal(result, exp)

    result = result.str.lower()
    tm.assert_series_equal(result, values)

    # mixed
    mixed = Series(["a", np.nan, "b", True, datetime.today(), "foo", None, 1, 2.0])
    mixed = mixed.str.upper()
    rs = Series(mixed).str.lower()
    xp = Series(["a", np.nan, "b", np.nan, np.nan, "foo", np.nan, np.nan, np.nan])
    assert isinstance(rs, Series)
    tm.assert_series_equal(rs, xp)


def test_capitalize():
    values = Series(["FOO", "BAR", np.nan, "Blah", "blurg"])
    result = values.str.capitalize()
    exp = Series(["Foo", "Bar", np.nan, "Blah", "Blurg"])
    tm.assert_series_equal(result, exp)

    # mixed
    mixed = Series(["FOO", np.nan, "bar", True, datetime.today(), "blah", None, 1, 2.0])
    mixed = mixed.str.capitalize()
    exp = Series(["Foo", np.nan, "Bar", np.nan, np.nan, "Blah", np.nan, np.nan, np.nan])
    tm.assert_almost_equal(mixed, exp)


def test_swapcase():
    values = Series(["FOO", "BAR", np.nan, "Blah", "blurg"])
    result = values.str.swapcase()
    exp = Series(["foo", "bar", np.nan, "bLAH", "BLURG"])
    tm.assert_series_equal(result, exp)

    # mixed
    mixed = Series(["FOO", np.nan, "bar", True, datetime.today(), "Blah", None, 1, 2.0])
    mixed = mixed.str.swapcase()
    exp = Series(["foo", np.nan, "BAR", np.nan, np.nan, "bLAH", np.nan, np.nan, np.nan])
    tm.assert_almost_equal(mixed, exp)


def test_casemethods():
    values = ["aaa", "bbb", "CCC", "Dddd", "eEEE"]
    s = Series(values)
    assert s.str.lower().tolist() == [v.lower() for v in values]
    assert s.str.upper().tolist() == [v.upper() for v in values]
    assert s.str.title().tolist() == [v.title() for v in values]
    assert s.str.capitalize().tolist() == [v.capitalize() for v in values]
    assert s.str.swapcase().tolist() == [v.swapcase() for v in values]


def test_pad():
    values = Series(["a", "b", np.nan, "c", np.nan, "eeeeee"])

    result = values.str.pad(5, side="left")
    exp = Series(["    a", "    b", np.nan, "    c", np.nan, "eeeeee"])
    tm.assert_almost_equal(result, exp)

    result = values.str.pad(5, side="right")
    exp = Series(["a    ", "b    ", np.nan, "c    ", np.nan, "eeeeee"])
    tm.assert_almost_equal(result, exp)

    result = values.str.pad(5, side="both")
    exp = Series(["  a  ", "  b  ", np.nan, "  c  ", np.nan, "eeeeee"])
    tm.assert_almost_equal(result, exp)

    # mixed
    mixed = Series(["a", np.nan, "b", True, datetime.today(), "ee", None, 1, 2.0])

    rs = Series(mixed).str.pad(5, side="left")
    xp = Series(
        ["    a", np.nan, "    b", np.nan, np.nan, "   ee", np.nan, np.nan, np.nan]
    )

    assert isinstance(rs, Series)
    tm.assert_almost_equal(rs, xp)

    mixed = Series(["a", np.nan, "b", True, datetime.today(), "ee", None, 1, 2.0])

    rs = Series(mixed).str.pad(5, side="right")
    xp = Series(
        ["a    ", np.nan, "b    ", np.nan, np.nan, "ee   ", np.nan, np.nan, np.nan]
    )

    assert isinstance(rs, Series)
    tm.assert_almost_equal(rs, xp)

    mixed = Series(["a", np.nan, "b", True, datetime.today(), "ee", None, 1, 2.0])

    rs = Series(mixed).str.pad(5, side="both")
    xp = Series(
        ["  a  ", np.nan, "  b  ", np.nan, np.nan, "  ee ", np.nan, np.nan, np.nan]
    )

    assert isinstance(rs, Series)
    tm.assert_almost_equal(rs, xp)


def test_pad_fillchar():

    values = Series(["a", "b", np.nan, "c", np.nan, "eeeeee"])

    result = values.str.pad(5, side="left", fillchar="X")
    exp = Series(["XXXXa", "XXXXb", np.nan, "XXXXc", np.nan, "eeeeee"])
    tm.assert_almost_equal(result, exp)

    result = values.str.pad(5, side="right", fillchar="X")
    exp = Series(["aXXXX", "bXXXX", np.nan, "cXXXX", np.nan, "eeeeee"])
    tm.assert_almost_equal(result, exp)

    result = values.str.pad(5, side="both", fillchar="X")
    exp = Series(["XXaXX", "XXbXX", np.nan, "XXcXX", np.nan, "eeeeee"])
    tm.assert_almost_equal(result, exp)

    msg = "fillchar must be a character, not str"
    with pytest.raises(TypeError, match=msg):
        result = values.str.pad(5, fillchar="XY")

    msg = "fillchar must be a character, not int"
    with pytest.raises(TypeError, match=msg):
        result = values.str.pad(5, fillchar=5)


@pytest.mark.parametrize("f", ["center", "ljust", "rjust", "zfill", "pad"])
def test_pad_width(f):
    # see gh-13598
    s = Series(["1", "22", "a", "bb"])
    msg = "width must be of integer type, not*"

    with pytest.raises(TypeError, match=msg):
        getattr(s.str, f)("f")


def test_center_ljust_rjust():
    values = Series(["a", "b", np.nan, "c", np.nan, "eeeeee"])

    result = values.str.center(5)
    exp = Series(["  a  ", "  b  ", np.nan, "  c  ", np.nan, "eeeeee"])
    tm.assert_almost_equal(result, exp)

    result = values.str.ljust(5)
    exp = Series(["a    ", "b    ", np.nan, "c    ", np.nan, "eeeeee"])
    tm.assert_almost_equal(result, exp)

    result = values.str.rjust(5)
    exp = Series(["    a", "    b", np.nan, "    c", np.nan, "eeeeee"])
    tm.assert_almost_equal(result, exp)

    # mixed
    mixed = Series(["a", np.nan, "b", True, datetime.today(), "c", "eee", None, 1, 2.0])

    rs = Series(mixed).str.center(5)
    xp = Series(
        [
            "  a  ",
            np.nan,
            "  b  ",
            np.nan,
            np.nan,
            "  c  ",
            " eee ",
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    assert isinstance(rs, Series)
    tm.assert_almost_equal(rs, xp)

    rs = Series(mixed).str.ljust(5)
    xp = Series(
        [
            "a    ",
            np.nan,
            "b    ",
            np.nan,
            np.nan,
            "c    ",
            "eee  ",
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    assert isinstance(rs, Series)
    tm.assert_almost_equal(rs, xp)

    rs = Series(mixed).str.rjust(5)
    xp = Series(
        [
            "    a",
            np.nan,
            "    b",
            np.nan,
            np.nan,
            "    c",
            "  eee",
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    assert isinstance(rs, Series)
    tm.assert_almost_equal(rs, xp)


def test_center_ljust_rjust_fillchar():
    values = Series(["a", "bb", "cccc", "ddddd", "eeeeee"])

    result = values.str.center(5, fillchar="X")
    expected = Series(["XXaXX", "XXbbX", "Xcccc", "ddddd", "eeeeee"])
    tm.assert_series_equal(result, expected)
    expected = np.array([v.center(5, "X") for v in values.values], dtype=np.object_)
    tm.assert_numpy_array_equal(result.values, expected)

    result = values.str.ljust(5, fillchar="X")
    expected = Series(["aXXXX", "bbXXX", "ccccX", "ddddd", "eeeeee"])
    tm.assert_series_equal(result, expected)
    expected = np.array([v.ljust(5, "X") for v in values.values], dtype=np.object_)
    tm.assert_numpy_array_equal(result.values, expected)

    result = values.str.rjust(5, fillchar="X")
    expected = Series(["XXXXa", "XXXbb", "Xcccc", "ddddd", "eeeeee"])
    tm.assert_series_equal(result, expected)
    expected = np.array([v.rjust(5, "X") for v in values.values], dtype=np.object_)
    tm.assert_numpy_array_equal(result.values, expected)

    # If fillchar is not a charatter, normal str raises TypeError
    # 'aaa'.ljust(5, 'XY')
    # TypeError: must be char, not str
    template = "fillchar must be a character, not {dtype}"

    with pytest.raises(TypeError, match=template.format(dtype="str")):
        values.str.center(5, fillchar="XY")

    with pytest.raises(TypeError, match=template.format(dtype="str")):
        values.str.ljust(5, fillchar="XY")

    with pytest.raises(TypeError, match=template.format(dtype="str")):
        values.str.rjust(5, fillchar="XY")

    with pytest.raises(TypeError, match=template.format(dtype="int")):
        values.str.center(5, fillchar=1)

    with pytest.raises(TypeError, match=template.format(dtype="int")):
        values.str.ljust(5, fillchar=1)

    with pytest.raises(TypeError, match=template.format(dtype="int")):
        values.str.rjust(5, fillchar=1)


def test_zfill():
    values = Series(["1", "22", "aaa", "333", "45678"])

    result = values.str.zfill(5)
    expected = Series(["00001", "00022", "00aaa", "00333", "45678"])
    tm.assert_series_equal(result, expected)
    expected = np.array([v.zfill(5) for v in values.values], dtype=np.object_)
    tm.assert_numpy_array_equal(result.values, expected)

    result = values.str.zfill(3)
    expected = Series(["001", "022", "aaa", "333", "45678"])
    tm.assert_series_equal(result, expected)
    expected = np.array([v.zfill(3) for v in values.values], dtype=np.object_)
    tm.assert_numpy_array_equal(result.values, expected)

    values = Series(["1", np.nan, "aaa", np.nan, "45678"])
    result = values.str.zfill(5)
    expected = Series(["00001", np.nan, "00aaa", np.nan, "45678"])
    tm.assert_series_equal(result, expected)


def test_wrap():
    # test values are: two words less than width, two words equal to width,
    # two words greater than width, one word less than width, one word
    # equal to width, one word greater than width, multiple tokens with
    # trailing whitespace equal to width
    values = Series(
        [
            "hello world",
            "hello world!",
            "hello world!!",
            "abcdefabcde",
            "abcdefabcdef",
            "abcdefabcdefa",
            "ab ab ab ab ",
            "ab ab ab ab a",
            "\t",
        ]
    )

    # expected values
    xp = Series(
        [
            "hello world",
            "hello world!",
            "hello\nworld!!",
            "abcdefabcde",
            "abcdefabcdef",
            "abcdefabcdef\na",
            "ab ab ab ab",
            "ab ab ab ab\na",
            "",
        ]
    )

    rs = values.str.wrap(12, break_long_words=True)
    tm.assert_series_equal(rs, xp)

    # test with pre and post whitespace (non-unicode), NaN, and non-ascii
    # Unicode
    values = Series(["  pre  ", np.nan, "\xac\u20ac\U00008000 abadcafe"])
    xp = Series(["  pre", np.nan, "\xac\u20ac\U00008000 ab\nadcafe"])
    rs = values.str.wrap(6)
    tm.assert_series_equal(rs, xp)
