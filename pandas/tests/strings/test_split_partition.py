from datetime import datetime
import re

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.tests.strings import (
    _convert_na_value,
    is_object_or_nan_string_dtype,
)


@pytest.mark.parametrize("method", ["split", "rsplit"])
def test_split(any_string_dtype, method):
    values = pd.Series(["a_b_c", "c_d_e", np.nan, "f_g_h"], dtype=any_string_dtype)

    result = getattr(values.str, method)("_")
    exp = pd.Series([["a", "b", "c"], ["c", "d", "e"], np.nan, ["f", "g", "h"]])
    exp = _convert_na_value(values, exp)
    tm.assert_series_equal(result, exp)


@pytest.mark.parametrize("method", ["split", "rsplit"])
def test_split_more_than_one_char(any_string_dtype, method):
    # more than one char
    values = pd.Series(
        ["a__b__c", "c__d__e", np.nan, "f__g__h"], dtype=any_string_dtype
    )
    result = getattr(values.str, method)("__")
    exp = pd.Series([["a", "b", "c"], ["c", "d", "e"], np.nan, ["f", "g", "h"]])
    exp = _convert_na_value(values, exp)
    tm.assert_series_equal(result, exp)

    result = getattr(values.str, method)("__", expand=False)
    tm.assert_series_equal(result, exp)


def test_split_more_regex_split(any_string_dtype):
    # regex split
    values = pd.Series(["a,b_c", "c_d,e", np.nan, "f,g,h"], dtype=any_string_dtype)
    result = values.str.split("[,_]")
    exp = pd.Series([["a", "b", "c"], ["c", "d", "e"], np.nan, ["f", "g", "h"]])
    exp = _convert_na_value(values, exp)
    tm.assert_series_equal(result, exp)


def test_split_regex(any_string_dtype):
    # GH 43563
    # explicit regex = True split
    values = pd.Series("xxxjpgzzz.jpg", dtype=any_string_dtype)
    result = values.str.split(r"\.jpg", regex=True)
    exp = pd.Series([["xxxjpgzzz", ""]])
    tm.assert_series_equal(result, exp)


def test_split_regex_explicit(any_string_dtype):
    # explicit regex = True split with compiled regex
    regex_pat = re.compile(r".jpg")
    values = pd.Series("xxxjpgzzz.jpg", dtype=any_string_dtype)
    result = values.str.split(regex_pat)
    exp = pd.Series([["xx", "zzz", ""]])
    tm.assert_series_equal(result, exp)

    # explicit regex = False split
    result = values.str.split(r"\.jpg", regex=False)
    exp = pd.Series([["xxxjpgzzz.jpg"]])
    tm.assert_series_equal(result, exp)

    # non explicit regex split, pattern length == 1
    result = values.str.split(r".")
    exp = pd.Series([["xxxjpgzzz", "jpg"]])
    tm.assert_series_equal(result, exp)

    # non explicit regex split, pattern length != 1
    result = values.str.split(r".jpg")
    exp = pd.Series([["xx", "zzz", ""]])
    tm.assert_series_equal(result, exp)

    # regex=False with pattern compiled regex raises error
    with pytest.raises(
        ValueError,
        match="Cannot use a compiled regex as replacement pattern with regex=False",
    ):
        values.str.split(regex_pat, regex=False)


@pytest.mark.parametrize("expand", [None, False])
@pytest.mark.parametrize("method", ["split", "rsplit"])
def test_split_object_mixed(expand, method):
    mixed = pd.Series(
        ["a_b_c", np.nan, "d_e_f", True, datetime(2011, 1, 1), None, 1, 2.0]
    )
    result = getattr(mixed.str, method)("_", expand=expand)
    exp = pd.Series(
        [
            ["a", "b", "c"],
            np.nan,
            ["d", "e", "f"],
            np.nan,
            np.nan,
            None,
            np.nan,
            np.nan,
        ]
    )
    assert isinstance(result, pd.Series)
    tm.assert_almost_equal(result, exp)


@pytest.mark.parametrize("method", ["split", "rsplit"])
@pytest.mark.parametrize("n", [None, 0])
def test_split_n(any_string_dtype, method, n):
    s = pd.Series(["a b", pd.NA, "b c"], dtype=any_string_dtype)
    expected = pd.Series([["a", "b"], pd.NA, ["b", "c"]])
    result = getattr(s.str, method)(" ", n=n)
    expected = _convert_na_value(s, expected)
    tm.assert_series_equal(result, expected)


def test_rsplit(any_string_dtype):
    values = pd.Series(["a,b_c", "c_d,e", np.nan, "f,g,h"], dtype=any_string_dtype)
    result = values.str.rsplit("[,_]")
    exp = pd.Series([["a,b_c"], ["c_d,e"], np.nan, ["f,g,h"]])
    exp = _convert_na_value(values, exp)
    tm.assert_series_equal(result, exp)


def test_rsplit_with_regex(any_string_dtype):
    # regex split is not supported by rsplit
    values = pd.Series(["a,b_c", "c_d,e", np.nan, "f,g,h"], dtype=any_string_dtype)
    with pytest.raises(TypeError, match="expected a string object, not Pattern"):
        values.str.rsplit(re.compile("[,_]"))


def test_rsplit_max_number(any_string_dtype):
    # setting max number of splits, make sure it's from reverse
    values = pd.Series(["a_b_c", "c_d_e", np.nan, "f_g_h"], dtype=any_string_dtype)
    result = values.str.rsplit("_", n=1)
    exp = pd.Series([["a_b", "c"], ["c_d", "e"], np.nan, ["f_g", "h"]])
    exp = _convert_na_value(values, exp)
    tm.assert_series_equal(result, exp)


def test_split_blank_string(any_string_dtype):
    # expand blank split GH 20067
    values = pd.Series([""], name="test", dtype=any_string_dtype)
    result = values.str.split(expand=True)
    exp = pd.DataFrame([[]], dtype=any_string_dtype)  # NOTE: this is NOT an empty df
    tm.assert_frame_equal(result, exp)


def test_split_blank_string_with_non_empty(any_string_dtype):
    values = pd.Series(["a b c", "a b", "", " "], name="test", dtype=any_string_dtype)
    result = values.str.split(expand=True)
    exp = pd.DataFrame(
        [
            ["a", "b", "c"],
            ["a", "b", None],
            [None, None, None],
            [None, None, None],
        ],
        dtype=any_string_dtype,
    )
    tm.assert_frame_equal(result, exp)


@pytest.mark.parametrize("method", ["split", "rsplit"])
def test_split_noargs(any_string_dtype, method):
    # #1859
    s = pd.Series(["Wes McKinney", "Travis  Oliphant"], dtype=any_string_dtype)
    result = getattr(s.str, method)()
    expected = ["Travis", "Oliphant"]
    assert result[1] == expected


@pytest.mark.parametrize(
    "data, pat",
    [
        (["bd asdf jfg", "kjasdflqw asdfnfk"], None),
        (["bd asdf jfg", "kjasdflqw asdfnfk"], "asdf"),
        (["bd_asdf_jfg", "kjasdflqw_asdfnfk"], "_"),
    ],
)
@pytest.mark.parametrize("n", [-1, 0])
def test_split_maxsplit(data, pat, any_string_dtype, n):
    # re.split 0, str.split -1
    s = pd.Series(data, dtype=any_string_dtype)

    result = s.str.split(pat=pat, n=n)
    xp = s.str.split(pat=pat)
    tm.assert_series_equal(result, xp)


@pytest.mark.parametrize(
    "data, pat, expected_val",
    [
        (
            ["split once", "split once too!"],
            None,
            "once too!",
        ),
        (
            ["split_once", "split_once_too!"],
            "_",
            "once_too!",
        ),
    ],
)
def test_split_no_pat_with_nonzero_n(data, pat, expected_val, any_string_dtype):
    s = pd.Series(data, dtype=any_string_dtype)
    result = s.str.split(pat=pat, n=1)
    expected = pd.Series({0: ["split", "once"], 1: ["split", expected_val]})
    tm.assert_series_equal(expected, result, check_index_type=False)


def test_split_to_dataframe_no_splits(any_string_dtype):
    s = pd.Series(["nosplit", "alsonosplit"], dtype=any_string_dtype)
    result = s.str.split("_", expand=True)
    exp = pd.DataFrame(
        {0: pd.Series(["nosplit", "alsonosplit"], dtype=any_string_dtype)}
    )
    tm.assert_frame_equal(result, exp)


def test_split_to_dataframe(any_string_dtype):
    s = pd.Series(["some_equal_splits", "with_no_nans"], dtype=any_string_dtype)
    result = s.str.split("_", expand=True)
    exp = pd.DataFrame(
        {0: ["some", "with"], 1: ["equal", "no"], 2: ["splits", "nans"]},
        dtype=any_string_dtype,
    )
    tm.assert_frame_equal(result, exp)


def test_split_to_dataframe_unequal_splits(any_string_dtype):
    s = pd.Series(
        ["some_unequal_splits", "one_of_these_things_is_not"], dtype=any_string_dtype
    )
    result = s.str.split("_", expand=True)
    exp = pd.DataFrame(
        {
            0: ["some", "one"],
            1: ["unequal", "of"],
            2: ["splits", "these"],
            3: [None, "things"],
            4: [None, "is"],
            5: [None, "not"],
        },
        dtype=any_string_dtype,
    )
    tm.assert_frame_equal(result, exp)


def test_split_to_dataframe_with_index(any_string_dtype):
    s = pd.Series(
        ["some_splits", "with_index"], index=["preserve", "me"], dtype=any_string_dtype
    )
    result = s.str.split("_", expand=True)
    exp = pd.DataFrame(
        {0: ["some", "with"], 1: ["splits", "index"]},
        index=["preserve", "me"],
        dtype=any_string_dtype,
    )
    tm.assert_frame_equal(result, exp)

    with pytest.raises(ValueError, match="expand must be"):
        s.str.split("_", expand="not_a_boolean")


def test_split_to_multiindex_expand_no_splits():
    # https://github.com/pandas-dev/pandas/issues/23677

    idx = pd.Index(["nosplit", "alsonosplit", np.nan])
    result = idx.str.split("_", expand=True)
    exp = idx
    tm.assert_index_equal(result, exp)
    assert result.nlevels == 1


def test_split_to_multiindex_expand():
    idx = pd.Index(["some_equal_splits", "with_no_nans", np.nan, None])
    result = idx.str.split("_", expand=True)
    exp = pd.MultiIndex.from_tuples(
        [
            ("some", "equal", "splits"),
            ("with", "no", "nans"),
            [np.nan, np.nan, np.nan],
            [None, None, None],
        ]
    )
    tm.assert_index_equal(result, exp)
    assert result.nlevels == 3


def test_split_to_multiindex_expand_unequal_splits():
    idx = pd.Index(["some_unequal_splits", "one_of_these_things_is_not", np.nan, None])
    result = idx.str.split("_", expand=True)
    exp = pd.MultiIndex.from_tuples(
        [
            ("some", "unequal", "splits", np.nan, np.nan, np.nan),
            ("one", "of", "these", "things", "is", "not"),
            (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan),
            (None, None, None, None, None, None),
        ]
    )
    tm.assert_index_equal(result, exp)
    assert result.nlevels == 6

    with pytest.raises(ValueError, match="expand must be"):
        idx.str.split("_", expand="not_a_boolean")


@pytest.mark.parametrize(
    "pat, expected_data",
    [
        (r"a(?=b)", [["aa"], ["", "b"], ["ba"], ["bb"]]),
        (r"(?<=a)b", [["aa"], ["a", ""], ["ba"], ["bb"]]),
        (r"a(?!b)", [["", "", ""], ["ab"], ["b", ""], ["bb"]]),
        (r"(?<!b)a", [["", "", ""], ["", "b"], ["ba"], ["bb"]]),
        ("ab", [["aa"], ["", ""], ["ba"], ["bb"]]),
    ],
)
def test_split_lookarounds(any_string_dtype, pat, expected_data):
    # https://github.com/pandas-dev/pandas/issues/60833
    ser = pd.Series(["aa", "ab", "ba", "bb", None], dtype=any_string_dtype)
    result = ser.str.split(pat, regex=True)
    if any_string_dtype == "object":
        null_result = None
    elif any_string_dtype == "str":
        null_result = np.nan
    elif any_string_dtype == "string":
        null_result = pd.NA
    else:
        raise ValueError(f"Unrecognized dtype: {any_string_dtype}")
    expected = pd.Series([*expected_data, null_result])
    tm.assert_series_equal(result, expected)


def test_split_regex_end_of_string(any_string_dtype):
    # https://github.com/pandas-dev/pandas/pull/63613
    ser = pd.Series(["baz", "bar", "bars", "bar\n"], dtype=any_string_dtype)

    # with dollar sign
    result = ser.str.split("r$", regex=True)
    expected = pd.Series([["baz"], ["ba", ""], ["bars"], ["ba", "\n"]], dtype=object)
    tm.assert_series_equal(result, expected)

    # with \Z (ensure this is translated to \z for pyarrow)
    result = ser.str.split(r"r\Z", regex=True)
    expected = pd.Series([["baz"], ["ba", ""], ["bars"], ["bar\n"]], dtype=object)
    tm.assert_series_equal(result, expected)

    # ensure finding a literal \Z still works
    ser = pd.Series([r"bar\Z", "bar", r"bar\Zs", "bar\n"], dtype=any_string_dtype)
    result = ser.str.split(r"r\\Z", regex=True)
    expected = pd.Series([["ba", ""], ["bar"], ["ba", "s"], ["bar\n"]], dtype=object)
    tm.assert_series_equal(result, expected)


def test_rsplit_to_dataframe_expand_no_splits(any_string_dtype):
    s = pd.Series(["nosplit", "alsonosplit"], dtype=any_string_dtype)
    result = s.str.rsplit("_", expand=True)
    exp = pd.DataFrame(
        {0: pd.Series(["nosplit", "alsonosplit"])}, dtype=any_string_dtype
    )
    tm.assert_frame_equal(result, exp)


def test_rsplit_to_dataframe_expand(any_string_dtype):
    s = pd.Series(["some_equal_splits", "with_no_nans"], dtype=any_string_dtype)
    result = s.str.rsplit("_", expand=True)
    exp = pd.DataFrame(
        {0: ["some", "with"], 1: ["equal", "no"], 2: ["splits", "nans"]},
        dtype=any_string_dtype,
    )
    tm.assert_frame_equal(result, exp)

    result = s.str.rsplit("_", expand=True, n=2)
    exp = pd.DataFrame(
        {0: ["some", "with"], 1: ["equal", "no"], 2: ["splits", "nans"]},
        dtype=any_string_dtype,
    )
    tm.assert_frame_equal(result, exp)

    result = s.str.rsplit("_", expand=True, n=1)
    exp = pd.DataFrame(
        {0: ["some_equal", "with_no"], 1: ["splits", "nans"]}, dtype=any_string_dtype
    )
    tm.assert_frame_equal(result, exp)


def test_rsplit_to_dataframe_expand_with_index(any_string_dtype):
    s = pd.Series(
        ["some_splits", "with_index"], index=["preserve", "me"], dtype=any_string_dtype
    )
    result = s.str.rsplit("_", expand=True)
    exp = pd.DataFrame(
        {0: ["some", "with"], 1: ["splits", "index"]},
        index=["preserve", "me"],
        dtype=any_string_dtype,
    )
    tm.assert_frame_equal(result, exp)


def test_rsplit_to_multiindex_expand_no_split():
    idx = pd.Index(["nosplit", "alsonosplit"])
    result = idx.str.rsplit("_", expand=True)
    exp = idx
    tm.assert_index_equal(result, exp)
    assert result.nlevels == 1


def test_rsplit_to_multiindex_expand():
    idx = pd.Index(["some_equal_splits", "with_no_nans"])
    result = idx.str.rsplit("_", expand=True)
    exp = pd.MultiIndex.from_tuples(
        [("some", "equal", "splits"), ("with", "no", "nans")]
    )
    tm.assert_index_equal(result, exp)
    assert result.nlevels == 3


def test_rsplit_to_multiindex_expand_n():
    idx = pd.Index(["some_equal_splits", "with_no_nans"])
    result = idx.str.rsplit("_", expand=True, n=1)
    exp = pd.MultiIndex.from_tuples([("some_equal", "splits"), ("with_no", "nans")])
    tm.assert_index_equal(result, exp)
    assert result.nlevels == 2


def test_split_nan_expand(any_string_dtype):
    # gh-18450
    s = pd.Series(["foo,bar,baz", np.nan], dtype=any_string_dtype)
    result = s.str.split(",", expand=True)
    exp = pd.DataFrame(
        [["foo", "bar", "baz"], [np.nan, np.nan, np.nan]], dtype=any_string_dtype
    )
    tm.assert_frame_equal(result, exp)

    # check that these are actually np.nan/pd.NA and not None
    # TODO see GH 18463
    # tm.assert_frame_equal does not differentiate
    if is_object_or_nan_string_dtype(any_string_dtype):
        assert all(np.isnan(x) for x in result.iloc[1])
    else:
        assert all(x is pd.NA for x in result.iloc[1])


def test_split_with_name_series(any_string_dtype):
    # GH 12617

    # should preserve name
    s = pd.Series(["a,b", "c,d"], name="xxx", dtype=any_string_dtype)
    res = s.str.split(",")
    exp = pd.Series([["a", "b"], ["c", "d"]], name="xxx")
    tm.assert_series_equal(res, exp)

    res = s.str.split(",", expand=True)
    exp = pd.DataFrame([["a", "b"], ["c", "d"]], dtype=any_string_dtype)
    tm.assert_frame_equal(res, exp)


def test_split_with_name_index():
    # GH 12617
    idx = pd.Index(["a,b", "c,d"], name="xxx")
    res = idx.str.split(",")
    exp = pd.Index([["a", "b"], ["c", "d"]], name="xxx")
    assert res.nlevels == 1
    tm.assert_index_equal(res, exp)

    res = idx.str.split(",", expand=True)
    exp = pd.MultiIndex.from_tuples([("a", "b"), ("c", "d")])
    assert res.nlevels == 2
    tm.assert_index_equal(res, exp)


@pytest.mark.parametrize(
    "method, exp",
    [
        [
            "partition",
            [
                ("a", "__", "b__c"),
                ("c", "__", "d__e"),
                np.nan,
                ("f", "__", "g__h"),
                None,
            ],
        ],
        [
            "rpartition",
            [
                ("a__b", "__", "c"),
                ("c__d", "__", "e"),
                np.nan,
                ("f__g", "__", "h"),
                None,
            ],
        ],
    ],
)
def test_partition_series_more_than_one_char(method, exp, any_string_dtype):
    # https://github.com/pandas-dev/pandas/issues/23558
    # more than one char
    s = pd.Series(
        ["a__b__c", "c__d__e", np.nan, "f__g__h", None], dtype=any_string_dtype
    )
    result = getattr(s.str, method)("__", expand=False)
    expected = pd.Series(exp)
    expected = _convert_na_value(s, expected)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "method, exp",
    [
        [
            "partition",
            [("a", " ", "b c"), ("c", " ", "d e"), np.nan, ("f", " ", "g h"), None],
        ],
        [
            "rpartition",
            [("a b", " ", "c"), ("c d", " ", "e"), np.nan, ("f g", " ", "h"), None],
        ],
    ],
)
def test_partition_series_none(any_string_dtype, method, exp):
    # https://github.com/pandas-dev/pandas/issues/23558
    # None
    s = pd.Series(["a b c", "c d e", np.nan, "f g h", None], dtype=any_string_dtype)
    result = getattr(s.str, method)(expand=False)
    expected = pd.Series(exp)
    expected = _convert_na_value(s, expected)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "method, exp",
    [
        [
            "partition",
            [("abc", "", ""), ("cde", "", ""), np.nan, ("fgh", "", ""), None],
        ],
        [
            "rpartition",
            [("", "", "abc"), ("", "", "cde"), np.nan, ("", "", "fgh"), None],
        ],
    ],
)
def test_partition_series_not_split(any_string_dtype, method, exp):
    # https://github.com/pandas-dev/pandas/issues/23558
    # Not split
    s = pd.Series(["abc", "cde", np.nan, "fgh", None], dtype=any_string_dtype)
    result = getattr(s.str, method)("_", expand=False)
    expected = pd.Series(exp)
    expected = _convert_na_value(s, expected)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "method, exp",
    [
        [
            "partition",
            [("a", "_", "b_c"), ("c", "_", "d_e"), np.nan, ("f", "_", "g_h")],
        ],
        [
            "rpartition",
            [("a_b", "_", "c"), ("c_d", "_", "e"), np.nan, ("f_g", "_", "h")],
        ],
    ],
)
def test_partition_series_unicode(any_string_dtype, method, exp):
    # https://github.com/pandas-dev/pandas/issues/23558
    # unicode
    s = pd.Series(["a_b_c", "c_d_e", np.nan, "f_g_h"], dtype=any_string_dtype)

    result = getattr(s.str, method)("_", expand=False)
    expected = pd.Series(exp)
    expected = _convert_na_value(s, expected)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("method", ["partition", "rpartition"])
def test_partition_series_stdlib(any_string_dtype, method):
    # https://github.com/pandas-dev/pandas/issues/23558
    # compare to standard lib
    s = pd.Series(["A_B_C", "B_C_D", "E_F_G", "EFGHEF"], dtype=any_string_dtype)
    result = getattr(s.str, method)("_", expand=False).tolist()
    assert result == [getattr(v, method)("_") for v in s]


@pytest.mark.parametrize(
    "method, exp",
    [
        [
            "partition",
            [("a", "_", "b_c"), ("c", "_", "d_e"), ("f", "_", "g_h"), np.nan, None],
        ],
        [
            "rpartition",
            [("a_b", "_", "c"), ("c_d", "_", "e"), ("f_g", "_", "h"), np.nan, None],
        ],
    ],
)
def test_partition_index(method, exp):
    # https://github.com/pandas-dev/pandas/issues/23558

    values = pd.Index(["a_b_c", "c_d_e", "f_g_h", np.nan, None])

    result = getattr(values.str, method)("_", expand=False)
    exp = pd.Index(np.array(exp, dtype=object), dtype=object)
    tm.assert_index_equal(result, exp)
    assert result.nlevels == 1


@pytest.mark.parametrize(
    "method, exp",
    [
        [
            "partition",
            {
                0: ["a", "c", np.nan, "f", None],
                1: ["_", "_", np.nan, "_", None],
                2: ["b_c", "d_e", np.nan, "g_h", None],
            },
        ],
        [
            "rpartition",
            {
                0: ["a_b", "c_d", np.nan, "f_g", None],
                1: ["_", "_", np.nan, "_", None],
                2: ["c", "e", np.nan, "h", None],
            },
        ],
    ],
)
def test_partition_to_dataframe(any_string_dtype, method, exp):
    # https://github.com/pandas-dev/pandas/issues/23558

    s = pd.Series(["a_b_c", "c_d_e", np.nan, "f_g_h", None], dtype=any_string_dtype)
    result = getattr(s.str, method)("_")
    expected = pd.DataFrame(
        exp,
        dtype=any_string_dtype,
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "method, exp",
    [
        [
            "partition",
            {
                0: ["a", "c", np.nan, "f", None],
                1: ["_", "_", np.nan, "_", None],
                2: ["b_c", "d_e", np.nan, "g_h", None],
            },
        ],
        [
            "rpartition",
            {
                0: ["a_b", "c_d", np.nan, "f_g", None],
                1: ["_", "_", np.nan, "_", None],
                2: ["c", "e", np.nan, "h", None],
            },
        ],
    ],
)
def test_partition_to_dataframe_from_series(any_string_dtype, method, exp):
    # https://github.com/pandas-dev/pandas/issues/23558
    s = pd.Series(["a_b_c", "c_d_e", np.nan, "f_g_h", None], dtype=any_string_dtype)
    result = getattr(s.str, method)("_", expand=True)
    expected = pd.DataFrame(
        exp,
        dtype=any_string_dtype,
    )
    tm.assert_frame_equal(result, expected)


def test_partition_with_name(any_string_dtype):
    # GH 12617

    s = pd.Series(["a,b", "c,d"], name="xxx", dtype=any_string_dtype)
    result = s.str.partition(",")
    expected = pd.DataFrame(
        {0: ["a", "c"], 1: [",", ","], 2: ["b", "d"]}, dtype=any_string_dtype
    )
    tm.assert_frame_equal(result, expected)


def test_partition_with_name_expand(any_string_dtype):
    # GH 12617
    # should preserve name
    s = pd.Series(["a,b", "c,d"], name="xxx", dtype=any_string_dtype)
    result = s.str.partition(",", expand=False)
    expected = pd.Series([("a", ",", "b"), ("c", ",", "d")], name="xxx")
    tm.assert_series_equal(result, expected)


def test_partition_index_with_name():
    idx = pd.Index(["a,b", "c,d"], name="xxx")
    result = idx.str.partition(",")
    expected = pd.MultiIndex.from_tuples([("a", ",", "b"), ("c", ",", "d")])
    assert result.nlevels == 3
    tm.assert_index_equal(result, expected)


def test_partition_index_with_name_expand_false():
    idx = pd.Index(["a,b", "c,d"], name="xxx")
    # should preserve name
    result = idx.str.partition(",", expand=False)
    expected = pd.Index(np.array([("a", ",", "b"), ("c", ",", "d")]), name="xxx")
    assert result.nlevels == 1
    tm.assert_index_equal(result, expected)


@pytest.mark.parametrize(
    "method, exp",
    [
        ["partition", {0: ["abc", "d"], 1: ["", ""], 2: ["", ""]}],
        ["rpartition", {0: ["", ""], 1: ["", ""], 2: ["abc", "d"]}],
    ],
)
def test_partition_to_dataframe_not_split(any_string_dtype, method, exp):
    # GH#63602 sep is absent, so two of the three columns are empty strings
    s = pd.Series(["abc", "d"], dtype=any_string_dtype)
    result = getattr(s.str, method)("_")
    expected = pd.DataFrame(exp, dtype=any_string_dtype)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "method, exp",
    [
        ["partition", {0: ["é", "ü"], 1: ["_", "_"], 2: ["ü_é", "é"]}],
        ["rpartition", {0: ["é_ü", "ü"], 1: ["_", "_"], 2: ["é", "é"]}],
    ],
)
def test_partition_to_dataframe_multibyte(any_string_dtype, method, exp):
    # GH#63602 splitting must happen on character, not byte, boundaries
    s = pd.Series(["é_ü_é", "ü_é"], dtype=any_string_dtype)
    result = getattr(s.str, method)("_")
    expected = pd.DataFrame(exp, dtype=any_string_dtype)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("method", ["partition", "rpartition"])
@pytest.mark.parametrize("expand", [True, False])
def test_partition_empty_sep(any_string_dtype, method, expand):
    # GH#63602 every dtype reports an empty separator the way str.partition does
    s = pd.Series(["a_b"], dtype=any_string_dtype)
    with pytest.raises(ValueError, match="empty separator"):
        getattr(s.str, method)("", expand=expand)


@pytest.mark.parametrize("method", ["partition", "rpartition"])
def test_partition_to_dataframe_empty(any_string_dtype, method):
    # GH#63602 an empty Series expands to no columns at all
    s = pd.Series([], dtype=any_string_dtype)
    result = getattr(s.str, method)("_")
    expected = pd.DataFrame(index=pd.RangeIndex(0), columns=pd.RangeIndex(0))
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("method", ["partition", "rpartition"])
def test_partition_to_dataframe_all_na(any_string_dtype, method):
    # GH#63602 with no non-null row to take a width from, expand to one
    #  all-NA column, as the object-dtype implementation does
    s = pd.Series([None, None], dtype=any_string_dtype)
    result = getattr(s.str, method)("_")
    expected = pd.DataFrame({0: [None, None]}, dtype=any_string_dtype)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("method", ["partition", "rpartition"])
def test_partition_sep_kwarg(any_string_dtype, method):
    # GH 22676; depr kwarg "pat" in favor of "sep"
    s = pd.Series(["a_b_c", "c_d_e", np.nan, "f_g_h"], dtype=any_string_dtype)

    expected = getattr(s.str, method)(sep="_")
    result = getattr(s.str, method)("_")
    tm.assert_frame_equal(result, expected)


def test_get():
    ser = pd.Series(["a_b_c", "c_d_e", np.nan, "f_g_h"])
    result = ser.str.split("_").str.get(1)
    expected = pd.Series(["b", "d", np.nan, "g"], dtype=object)
    tm.assert_series_equal(result, expected)


def test_get_mixed_object():
    ser = pd.Series(
        ["a_b_c", np.nan, "c_d_e", True, datetime(2011, 1, 1), None, 1, 2.0]
    )
    result = ser.str.split("_").str.get(1)
    expected = pd.Series(
        ["b", np.nan, "d", np.nan, np.nan, None, np.nan, np.nan], dtype=object
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("idx", [2, -3])
def test_get_bounds(idx):
    ser = pd.Series(["1_2_3_4_5", "6_7_8_9_10", "11_12"])
    result = ser.str.split("_").str.get(idx)
    expected = pd.Series(["3", "8", np.nan], dtype=object)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "idx, exp", [[2, [3, 3, np.nan, "b"]], [-1, [3, 3, np.nan, np.nan]]]
)
def test_get_complex(idx, exp):
    # GH 20671, getting value not in dict raising `KeyError`
    ser = pd.Series([(1, 2, 3), [1, 2, 3], {1, 2, 3}, {1: "a", 2: "b", 3: "c"}])

    result = ser.str.get(idx)
    expected = pd.Series(exp)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("to_type", [tuple, list, np.array])
def test_get_complex_nested(to_type):
    ser = pd.Series([to_type([to_type([1, 2])])])

    result = ser.str.get(0)
    expected = pd.Series([to_type([1, 2])])
    tm.assert_series_equal(result, expected)

    result = ser.str.get(1)
    expected = pd.Series([np.nan])
    tm.assert_series_equal(result, expected)


def test_get_strings(any_string_dtype):
    ser = pd.Series(["a", "ab", np.nan, "abc"], dtype=any_string_dtype)
    result = ser.str.get(2)
    expected = pd.Series([np.nan, np.nan, np.nan, "c"], dtype=any_string_dtype)
    tm.assert_series_equal(result, expected)


def test_split_arrow_dtype_regex_inference():
    # GH#58321 - ArrowDtype string should infer regex for multi-char patterns
    pa = pytest.importorskip("pyarrow")

    ser = pd.Series(["230/270/270", "240-290-290"], dtype=pd.ArrowDtype(pa.string()))
    result = ser.str.split(r"/|-", expand=True)
    expected = pd.DataFrame(
        {"0": ["230", "240"], "1": ["270", "290"], "2": ["270", "290"]},
        dtype=pd.ArrowDtype(pa.string()),
    )
    expected.columns = pd.RangeIndex(3)
    tm.assert_frame_equal(result, expected)


def test_split_arrow_dtype_regex_true():
    # GH#58321 - explicit regex=True should use regex splitting
    pa = pytest.importorskip("pyarrow")

    ser = pd.Series(["230/270/270", "240-290-290"], dtype=pd.ArrowDtype(pa.string()))
    result = ser.str.split(r"/|-", expand=True, regex=True)
    expected = pd.DataFrame(
        {"0": ["230", "240"], "1": ["270", "290"], "2": ["270", "290"]},
        dtype=pd.ArrowDtype(pa.string()),
    )
    expected.columns = pd.RangeIndex(3)
    tm.assert_frame_equal(result, expected)


def test_split_arrow_dtype_single_char_literal():
    # GH#58321 - single-char pattern with regex=None should be literal
    pa = pytest.importorskip("pyarrow")

    ser = pd.Series(["a.b.c"], dtype=pd.ArrowDtype(pa.string()))
    result = ser.str.split(".", expand=True)
    expected = pd.DataFrame(
        {"0": ["a"], "1": ["b"], "2": ["c"]},
        dtype=pd.ArrowDtype(pa.string()),
    )
    expected.columns = pd.RangeIndex(3)
    tm.assert_frame_equal(result, expected)


def test_split_arrow_dtype_regex_false_multichar():
    # GH#58321 - regex=False with multi-char pattern: "|" is literal, not alternation
    pa = pytest.importorskip("pyarrow")

    ser = pd.Series(["a|b|c", "d/e/f"], dtype=pd.ArrowDtype(pa.string()))
    result = ser.str.split("|", regex=False, expand=True)
    expected = pd.DataFrame(
        {"0": ["a", "d/e/f"], "1": ["b", None], "2": ["c", None]},
        dtype=pd.ArrowDtype(pa.string()),
    )
    expected.columns = pd.RangeIndex(3)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("method", ["split", "rsplit", "partition", "rpartition"])
def test_split_category_returns_list_not_string(any_string_dtype, method):
    # GH#66341
    ser = pd.Series(
        [
            "this is a regular sentence",
            "https://docs.python.org/3/tutorial/index.html",
            np.nan,
        ],
        dtype=any_string_dtype,
    ).astype("category")

    if method in ("split", "rsplit"):
        result = getattr(ser.str, method)()
        expected = pd.Series(
            [
                ["this", "is", "a", "regular", "sentence"],
                ["https://docs.python.org/3/tutorial/index.html"],
                np.nan,
            ]
        )
    elif method == "partition":
        result = ser.str.partition(expand=False)
        expected = pd.Series(
            [
                ("this", " ", "is a regular sentence"),
                ("https://docs.python.org/3/tutorial/index.html", "", ""),
                np.nan,
            ]
        )
    else:
        result = ser.str.rpartition(expand=False)
        expected = pd.Series(
            [
                ("this is a regular", " ", "sentence"),
                ("", "", "https://docs.python.org/3/tutorial/index.html"),
                np.nan,
            ]
        )

    assert result.dtype == object
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("method", ["split", "rsplit", "partition", "rpartition"])
def test_split_category_expand(any_string_dtype, method):
    # GH#66341 - expand=True on categorical should return string columns
    ser = pd.Series(
        ["a b", "c d"],
        dtype=any_string_dtype,
    ).astype("category")

    result = getattr(ser.str, method)(expand=True)

    expected_dtype = ser.dtype.categories.dtype
    assert isinstance(result, pd.DataFrame)
    for col_dtype in result.dtypes:
        assert col_dtype == expected_dtype
