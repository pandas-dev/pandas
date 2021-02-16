from datetime import datetime

import numpy as np
import pytest

import pandas as pd
from pandas import (
    DataFrame,
    Index,
    MultiIndex,
    Series,
    _testing as tm,
)


def test_split():
    values = Series(["a_b_c", "c_d_e", np.nan, "f_g_h"])

    result = values.str.split("_")
    exp = Series([["a", "b", "c"], ["c", "d", "e"], np.nan, ["f", "g", "h"]])
    tm.assert_series_equal(result, exp)

    # more than one char
    values = Series(["a__b__c", "c__d__e", np.nan, "f__g__h"])
    result = values.str.split("__")
    tm.assert_series_equal(result, exp)

    result = values.str.split("__", expand=False)
    tm.assert_series_equal(result, exp)

    # mixed
    mixed = Series(["a_b_c", np.nan, "d_e_f", True, datetime.today(), None, 1, 2.0])
    result = mixed.str.split("_")
    exp = Series(
        [
            ["a", "b", "c"],
            np.nan,
            ["d", "e", "f"],
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    assert isinstance(result, Series)
    tm.assert_almost_equal(result, exp)

    result = mixed.str.split("_", expand=False)
    assert isinstance(result, Series)
    tm.assert_almost_equal(result, exp)

    # regex split
    values = Series(["a,b_c", "c_d,e", np.nan, "f,g,h"])
    result = values.str.split("[,_]")
    exp = Series([["a", "b", "c"], ["c", "d", "e"], np.nan, ["f", "g", "h"]])
    tm.assert_series_equal(result, exp)


@pytest.mark.parametrize("dtype", [object, "string"])
@pytest.mark.parametrize("method", ["split", "rsplit"])
def test_split_n(dtype, method):
    s = Series(["a b", pd.NA, "b c"], dtype=dtype)
    expected = Series([["a", "b"], pd.NA, ["b", "c"]])

    result = getattr(s.str, method)(" ", n=None)
    tm.assert_series_equal(result, expected)

    result = getattr(s.str, method)(" ", n=0)
    tm.assert_series_equal(result, expected)


def test_rsplit():
    values = Series(["a_b_c", "c_d_e", np.nan, "f_g_h"])
    result = values.str.rsplit("_")
    exp = Series([["a", "b", "c"], ["c", "d", "e"], np.nan, ["f", "g", "h"]])
    tm.assert_series_equal(result, exp)

    # more than one char
    values = Series(["a__b__c", "c__d__e", np.nan, "f__g__h"])
    result = values.str.rsplit("__")
    tm.assert_series_equal(result, exp)

    result = values.str.rsplit("__", expand=False)
    tm.assert_series_equal(result, exp)

    # mixed
    mixed = Series(["a_b_c", np.nan, "d_e_f", True, datetime.today(), None, 1, 2.0])
    result = mixed.str.rsplit("_")
    exp = Series(
        [
            ["a", "b", "c"],
            np.nan,
            ["d", "e", "f"],
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    assert isinstance(result, Series)
    tm.assert_almost_equal(result, exp)

    result = mixed.str.rsplit("_", expand=False)
    assert isinstance(result, Series)
    tm.assert_almost_equal(result, exp)

    # regex split is not supported by rsplit
    values = Series(["a,b_c", "c_d,e", np.nan, "f,g,h"])
    result = values.str.rsplit("[,_]")
    exp = Series([["a,b_c"], ["c_d,e"], np.nan, ["f,g,h"]])
    tm.assert_series_equal(result, exp)

    # setting max number of splits, make sure it's from reverse
    values = Series(["a_b_c", "c_d_e", np.nan, "f_g_h"])
    result = values.str.rsplit("_", n=1)
    exp = Series([["a_b", "c"], ["c_d", "e"], np.nan, ["f_g", "h"]])
    tm.assert_series_equal(result, exp)


def test_split_blank_string():
    # expand blank split GH 20067
    values = Series([""], name="test")
    result = values.str.split(expand=True)
    exp = DataFrame([[]])  # NOTE: this is NOT an empty DataFrame
    tm.assert_frame_equal(result, exp)

    values = Series(["a b c", "a b", "", " "], name="test")
    result = values.str.split(expand=True)
    exp = DataFrame(
        [
            ["a", "b", "c"],
            ["a", "b", np.nan],
            [np.nan, np.nan, np.nan],
            [np.nan, np.nan, np.nan],
        ]
    )
    tm.assert_frame_equal(result, exp)


def test_split_noargs():
    # #1859
    s = Series(["Wes McKinney", "Travis  Oliphant"])
    result = s.str.split()
    expected = ["Travis", "Oliphant"]
    assert result[1] == expected
    result = s.str.rsplit()
    assert result[1] == expected


def test_split_maxsplit():
    # re.split 0, str.split -1
    s = Series(["bd asdf jfg", "kjasdflqw asdfnfk"])

    result = s.str.split(n=-1)
    xp = s.str.split()
    tm.assert_series_equal(result, xp)

    result = s.str.split(n=0)
    tm.assert_series_equal(result, xp)

    xp = s.str.split("asdf")
    result = s.str.split("asdf", n=0)
    tm.assert_series_equal(result, xp)

    result = s.str.split("asdf", n=-1)
    tm.assert_series_equal(result, xp)


def test_split_no_pat_with_nonzero_n():
    s = Series(["split once", "split once too!"])
    result = s.str.split(n=1)
    expected = Series({0: ["split", "once"], 1: ["split", "once too!"]})
    tm.assert_series_equal(expected, result, check_index_type=False)


def test_split_to_dataframe():
    s = Series(["nosplit", "alsonosplit"])
    result = s.str.split("_", expand=True)
    exp = DataFrame({0: Series(["nosplit", "alsonosplit"])})
    tm.assert_frame_equal(result, exp)

    s = Series(["some_equal_splits", "with_no_nans"])
    result = s.str.split("_", expand=True)
    exp = DataFrame({0: ["some", "with"], 1: ["equal", "no"], 2: ["splits", "nans"]})
    tm.assert_frame_equal(result, exp)

    s = Series(["some_unequal_splits", "one_of_these_things_is_not"])
    result = s.str.split("_", expand=True)
    exp = DataFrame(
        {
            0: ["some", "one"],
            1: ["unequal", "of"],
            2: ["splits", "these"],
            3: [np.nan, "things"],
            4: [np.nan, "is"],
            5: [np.nan, "not"],
        }
    )
    tm.assert_frame_equal(result, exp)

    s = Series(["some_splits", "with_index"], index=["preserve", "me"])
    result = s.str.split("_", expand=True)
    exp = DataFrame(
        {0: ["some", "with"], 1: ["splits", "index"]}, index=["preserve", "me"]
    )
    tm.assert_frame_equal(result, exp)

    with pytest.raises(ValueError, match="expand must be"):
        s.str.split("_", expand="not_a_boolean")


def test_split_to_multiindex_expand():
    # https://github.com/pandas-dev/pandas/issues/23677

    idx = Index(["nosplit", "alsonosplit", np.nan])
    result = idx.str.split("_", expand=True)
    exp = idx
    tm.assert_index_equal(result, exp)
    assert result.nlevels == 1

    idx = Index(["some_equal_splits", "with_no_nans", np.nan, None])
    result = idx.str.split("_", expand=True)
    exp = MultiIndex.from_tuples(
        [
            ("some", "equal", "splits"),
            ("with", "no", "nans"),
            [np.nan, np.nan, np.nan],
            [None, None, None],
        ]
    )
    tm.assert_index_equal(result, exp)
    assert result.nlevels == 3

    idx = Index(["some_unequal_splits", "one_of_these_things_is_not", np.nan, None])
    result = idx.str.split("_", expand=True)
    exp = MultiIndex.from_tuples(
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


def test_rsplit_to_dataframe_expand():
    s = Series(["nosplit", "alsonosplit"])
    result = s.str.rsplit("_", expand=True)
    exp = DataFrame({0: Series(["nosplit", "alsonosplit"])})
    tm.assert_frame_equal(result, exp)

    s = Series(["some_equal_splits", "with_no_nans"])
    result = s.str.rsplit("_", expand=True)
    exp = DataFrame({0: ["some", "with"], 1: ["equal", "no"], 2: ["splits", "nans"]})
    tm.assert_frame_equal(result, exp)

    result = s.str.rsplit("_", expand=True, n=2)
    exp = DataFrame({0: ["some", "with"], 1: ["equal", "no"], 2: ["splits", "nans"]})
    tm.assert_frame_equal(result, exp)

    result = s.str.rsplit("_", expand=True, n=1)
    exp = DataFrame({0: ["some_equal", "with_no"], 1: ["splits", "nans"]})
    tm.assert_frame_equal(result, exp)

    s = Series(["some_splits", "with_index"], index=["preserve", "me"])
    result = s.str.rsplit("_", expand=True)
    exp = DataFrame(
        {0: ["some", "with"], 1: ["splits", "index"]}, index=["preserve", "me"]
    )
    tm.assert_frame_equal(result, exp)


def test_rsplit_to_multiindex_expand():
    idx = Index(["nosplit", "alsonosplit"])
    result = idx.str.rsplit("_", expand=True)
    exp = idx
    tm.assert_index_equal(result, exp)
    assert result.nlevels == 1

    idx = Index(["some_equal_splits", "with_no_nans"])
    result = idx.str.rsplit("_", expand=True)
    exp = MultiIndex.from_tuples([("some", "equal", "splits"), ("with", "no", "nans")])
    tm.assert_index_equal(result, exp)
    assert result.nlevels == 3

    idx = Index(["some_equal_splits", "with_no_nans"])
    result = idx.str.rsplit("_", expand=True, n=1)
    exp = MultiIndex.from_tuples([("some_equal", "splits"), ("with_no", "nans")])
    tm.assert_index_equal(result, exp)
    assert result.nlevels == 2


def test_split_nan_expand():
    # gh-18450
    s = Series(["foo,bar,baz", np.nan])
    result = s.str.split(",", expand=True)
    exp = DataFrame([["foo", "bar", "baz"], [np.nan, np.nan, np.nan]])
    tm.assert_frame_equal(result, exp)

    # check that these are actually np.nan and not None
    # TODO see GH 18463
    # tm.assert_frame_equal does not differentiate
    assert all(np.isnan(x) for x in result.iloc[1])


def test_split_with_name():
    # GH 12617

    # should preserve name
    s = Series(["a,b", "c,d"], name="xxx")
    res = s.str.split(",")
    exp = Series([["a", "b"], ["c", "d"]], name="xxx")
    tm.assert_series_equal(res, exp)

    res = s.str.split(",", expand=True)
    exp = DataFrame([["a", "b"], ["c", "d"]])
    tm.assert_frame_equal(res, exp)

    idx = Index(["a,b", "c,d"], name="xxx")
    res = idx.str.split(",")
    exp = Index([["a", "b"], ["c", "d"]], name="xxx")
    assert res.nlevels == 1
    tm.assert_index_equal(res, exp)

    res = idx.str.split(",", expand=True)
    exp = MultiIndex.from_tuples([("a", "b"), ("c", "d")])
    assert res.nlevels == 2
    tm.assert_index_equal(res, exp)


def test_partition_series():
    # https://github.com/pandas-dev/pandas/issues/23558

    values = Series(["a_b_c", "c_d_e", np.nan, "f_g_h", None])

    result = values.str.partition("_", expand=False)
    exp = Series(
        [("a", "_", "b_c"), ("c", "_", "d_e"), np.nan, ("f", "_", "g_h"), None]
    )
    tm.assert_series_equal(result, exp)

    result = values.str.rpartition("_", expand=False)
    exp = Series(
        [("a_b", "_", "c"), ("c_d", "_", "e"), np.nan, ("f_g", "_", "h"), None]
    )
    tm.assert_series_equal(result, exp)

    # more than one char
    values = Series(["a__b__c", "c__d__e", np.nan, "f__g__h", None])
    result = values.str.partition("__", expand=False)
    exp = Series(
        [
            ("a", "__", "b__c"),
            ("c", "__", "d__e"),
            np.nan,
            ("f", "__", "g__h"),
            None,
        ]
    )
    tm.assert_series_equal(result, exp)

    result = values.str.rpartition("__", expand=False)
    exp = Series(
        [
            ("a__b", "__", "c"),
            ("c__d", "__", "e"),
            np.nan,
            ("f__g", "__", "h"),
            None,
        ]
    )
    tm.assert_series_equal(result, exp)

    # None
    values = Series(["a b c", "c d e", np.nan, "f g h", None])
    result = values.str.partition(expand=False)
    exp = Series(
        [("a", " ", "b c"), ("c", " ", "d e"), np.nan, ("f", " ", "g h"), None]
    )
    tm.assert_series_equal(result, exp)

    result = values.str.rpartition(expand=False)
    exp = Series(
        [("a b", " ", "c"), ("c d", " ", "e"), np.nan, ("f g", " ", "h"), None]
    )
    tm.assert_series_equal(result, exp)

    # Not split
    values = Series(["abc", "cde", np.nan, "fgh", None])
    result = values.str.partition("_", expand=False)
    exp = Series([("abc", "", ""), ("cde", "", ""), np.nan, ("fgh", "", ""), None])
    tm.assert_series_equal(result, exp)

    result = values.str.rpartition("_", expand=False)
    exp = Series([("", "", "abc"), ("", "", "cde"), np.nan, ("", "", "fgh"), None])
    tm.assert_series_equal(result, exp)

    # unicode
    values = Series(["a_b_c", "c_d_e", np.nan, "f_g_h"])

    result = values.str.partition("_", expand=False)
    exp = Series([("a", "_", "b_c"), ("c", "_", "d_e"), np.nan, ("f", "_", "g_h")])
    tm.assert_series_equal(result, exp)

    result = values.str.rpartition("_", expand=False)
    exp = Series([("a_b", "_", "c"), ("c_d", "_", "e"), np.nan, ("f_g", "_", "h")])
    tm.assert_series_equal(result, exp)

    # compare to standard lib
    values = Series(["A_B_C", "B_C_D", "E_F_G", "EFGHEF"])
    result = values.str.partition("_", expand=False).tolist()
    assert result == [v.partition("_") for v in values]
    result = values.str.rpartition("_", expand=False).tolist()
    assert result == [v.rpartition("_") for v in values]


def test_partition_index():
    # https://github.com/pandas-dev/pandas/issues/23558

    values = Index(["a_b_c", "c_d_e", "f_g_h", np.nan, None])

    result = values.str.partition("_", expand=False)
    exp = Index(
        np.array(
            [("a", "_", "b_c"), ("c", "_", "d_e"), ("f", "_", "g_h"), np.nan, None],
            dtype=object,
        )
    )
    tm.assert_index_equal(result, exp)
    assert result.nlevels == 1

    result = values.str.rpartition("_", expand=False)
    exp = Index(
        np.array(
            [("a_b", "_", "c"), ("c_d", "_", "e"), ("f_g", "_", "h"), np.nan, None],
            dtype=object,
        )
    )
    tm.assert_index_equal(result, exp)
    assert result.nlevels == 1

    result = values.str.partition("_")
    exp = Index(
        [
            ("a", "_", "b_c"),
            ("c", "_", "d_e"),
            ("f", "_", "g_h"),
            (np.nan, np.nan, np.nan),
            (None, None, None),
        ]
    )
    tm.assert_index_equal(result, exp)
    assert isinstance(result, MultiIndex)
    assert result.nlevels == 3

    result = values.str.rpartition("_")
    exp = Index(
        [
            ("a_b", "_", "c"),
            ("c_d", "_", "e"),
            ("f_g", "_", "h"),
            (np.nan, np.nan, np.nan),
            (None, None, None),
        ]
    )
    tm.assert_index_equal(result, exp)
    assert isinstance(result, MultiIndex)
    assert result.nlevels == 3


def test_partition_to_dataframe():
    # https://github.com/pandas-dev/pandas/issues/23558

    values = Series(["a_b_c", "c_d_e", np.nan, "f_g_h", None])
    result = values.str.partition("_")
    exp = DataFrame(
        {
            0: ["a", "c", np.nan, "f", None],
            1: ["_", "_", np.nan, "_", None],
            2: ["b_c", "d_e", np.nan, "g_h", None],
        }
    )
    tm.assert_frame_equal(result, exp)

    result = values.str.rpartition("_")
    exp = DataFrame(
        {
            0: ["a_b", "c_d", np.nan, "f_g", None],
            1: ["_", "_", np.nan, "_", None],
            2: ["c", "e", np.nan, "h", None],
        }
    )
    tm.assert_frame_equal(result, exp)

    values = Series(["a_b_c", "c_d_e", np.nan, "f_g_h", None])
    result = values.str.partition("_", expand=True)
    exp = DataFrame(
        {
            0: ["a", "c", np.nan, "f", None],
            1: ["_", "_", np.nan, "_", None],
            2: ["b_c", "d_e", np.nan, "g_h", None],
        }
    )
    tm.assert_frame_equal(result, exp)

    result = values.str.rpartition("_", expand=True)
    exp = DataFrame(
        {
            0: ["a_b", "c_d", np.nan, "f_g", None],
            1: ["_", "_", np.nan, "_", None],
            2: ["c", "e", np.nan, "h", None],
        }
    )
    tm.assert_frame_equal(result, exp)


def test_partition_with_name():
    # GH 12617

    s = Series(["a,b", "c,d"], name="xxx")
    res = s.str.partition(",")
    exp = DataFrame({0: ["a", "c"], 1: [",", ","], 2: ["b", "d"]})
    tm.assert_frame_equal(res, exp)

    # should preserve name
    res = s.str.partition(",", expand=False)
    exp = Series([("a", ",", "b"), ("c", ",", "d")], name="xxx")
    tm.assert_series_equal(res, exp)

    idx = Index(["a,b", "c,d"], name="xxx")
    res = idx.str.partition(",")
    exp = MultiIndex.from_tuples([("a", ",", "b"), ("c", ",", "d")])
    assert res.nlevels == 3
    tm.assert_index_equal(res, exp)

    # should preserve name
    res = idx.str.partition(",", expand=False)
    exp = Index(np.array([("a", ",", "b"), ("c", ",", "d")]), name="xxx")
    assert res.nlevels == 1
    tm.assert_index_equal(res, exp)


def test_partition_sep_kwarg():
    # GH 22676; depr kwarg "pat" in favor of "sep"
    values = Series(["a_b_c", "c_d_e", np.nan, "f_g_h"])

    expected = values.str.partition(sep="_")
    result = values.str.partition("_")
    tm.assert_frame_equal(result, expected)

    expected = values.str.rpartition(sep="_")
    result = values.str.rpartition("_")
    tm.assert_frame_equal(result, expected)


def test_get():
    values = Series(["a_b_c", "c_d_e", np.nan, "f_g_h"])

    result = values.str.split("_").str.get(1)
    expected = Series(["b", "d", np.nan, "g"])
    tm.assert_series_equal(result, expected)

    # mixed
    mixed = Series(["a_b_c", np.nan, "c_d_e", True, datetime.today(), None, 1, 2.0])

    rs = Series(mixed).str.split("_").str.get(1)
    xp = Series(["b", np.nan, "d", np.nan, np.nan, np.nan, np.nan, np.nan])

    assert isinstance(rs, Series)
    tm.assert_almost_equal(rs, xp)

    # bounds testing
    values = Series(["1_2_3_4_5", "6_7_8_9_10", "11_12"])

    # positive index
    result = values.str.split("_").str.get(2)
    expected = Series(["3", "8", np.nan])
    tm.assert_series_equal(result, expected)

    # negative index
    result = values.str.split("_").str.get(-3)
    expected = Series(["3", "8", np.nan])
    tm.assert_series_equal(result, expected)


def test_get_complex():
    # GH 20671, getting value not in dict raising `KeyError`
    values = Series([(1, 2, 3), [1, 2, 3], {1, 2, 3}, {1: "a", 2: "b", 3: "c"}])

    result = values.str.get(1)
    expected = Series([2, 2, np.nan, "a"])
    tm.assert_series_equal(result, expected)

    result = values.str.get(-1)
    expected = Series([3, 3, np.nan, np.nan])
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("to_type", [tuple, list, np.array])
def test_get_complex_nested(to_type):
    values = Series([to_type([to_type([1, 2])])])

    result = values.str.get(0)
    expected = Series([to_type([1, 2])])
    tm.assert_series_equal(result, expected)

    result = values.str.get(1)
    expected = Series([np.nan])
    tm.assert_series_equal(result, expected)
