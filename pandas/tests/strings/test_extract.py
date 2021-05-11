from datetime import datetime
import re

import numpy as np
import pytest

from pandas import (
    DataFrame,
    Index,
    MultiIndex,
    Series,
    _testing as tm,
)


def test_extract_expand_kwarg_wrong_type_raises(any_string_dtype):
    # TODO: should this raise TypeError
    values = Series(["fooBAD__barBAD", np.nan, "foo"], dtype=any_string_dtype)
    with pytest.raises(ValueError, match="expand must be True or False"):
        values.str.extract(".*(BAD[_]+).*(BAD)", expand=None)


def test_extract_expand_kwarg(any_string_dtype):
    s = Series(["fooBAD__barBAD", np.nan, "foo"], dtype=any_string_dtype)
    expected = DataFrame(["BAD__", np.nan, np.nan], dtype=any_string_dtype)

    result = s.str.extract(".*(BAD[_]+).*")
    tm.assert_frame_equal(result, expected)

    result = s.str.extract(".*(BAD[_]+).*", expand=True)
    tm.assert_frame_equal(result, expected)

    expected = DataFrame(
        [["BAD__", "BAD"], [np.nan, np.nan], [np.nan, np.nan]], dtype=any_string_dtype
    )
    result = s.str.extract(".*(BAD[_]+).*(BAD)", expand=False)
    tm.assert_frame_equal(result, expected)


def test_extract_expand_False_mixed_object():
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

    result = Series(mixed).str.extract(".*(BAD[_]+).*(BAD)", expand=False)
    er = [np.nan, np.nan]  # empty row
    expected = DataFrame([["BAD_", "BAD"], er, ["BAD_", "BAD"], er, er, er, er, er, er])
    tm.assert_frame_equal(result, expected)


def test_extract_expand_index_raises():
    # GH9980
    # Index only works with one regex group since
    # multi-group would expand to a frame
    idx = Index(["A1", "A2", "A3", "A4", "B5"])
    msg = "only one regex group is supported with Index"
    with pytest.raises(ValueError, match=msg):
        idx.str.extract("([AB])([123])", expand=False)


def test_extract_expand_no_capture_groups_raises(index_or_series, any_string_dtype):
    s_or_idx = index_or_series(["A1", "B2", "C3"], dtype=any_string_dtype)
    msg = "pattern contains no capture groups"

    # no groups
    with pytest.raises(ValueError, match=msg):
        s_or_idx.str.extract("[ABC][123]", expand=False)

    # only non-capturing groups
    with pytest.raises(ValueError, match=msg):
        s_or_idx.str.extract("(?:[AB]).*", expand=False)


def test_extract_expand_single_capture_group(index_or_series, any_string_dtype):
    # single group renames series/index properly
    s_or_idx = index_or_series(["A1", "A2"], dtype=any_string_dtype)
    result = s_or_idx.str.extract(r"(?P<uno>A)\d", expand=False)

    expected = index_or_series(["A", "A"], name="uno", dtype=any_string_dtype)
    if index_or_series == Series:
        tm.assert_series_equal(result, expected)
    else:
        tm.assert_index_equal(result, expected)


def test_extract_expand_capture_groups(any_string_dtype):
    s = Series(["A1", "B2", "C3"], dtype=any_string_dtype)
    # one group, no matches
    result = s.str.extract("(_)", expand=False)
    expected = Series([np.nan, np.nan, np.nan], dtype=any_string_dtype)
    tm.assert_series_equal(result, expected)

    # two groups, no matches
    result = s.str.extract("(_)(_)", expand=False)
    expected = DataFrame(
        [[np.nan, np.nan], [np.nan, np.nan], [np.nan, np.nan]], dtype=any_string_dtype
    )
    tm.assert_frame_equal(result, expected)

    # one group, some matches
    result = s.str.extract("([AB])[123]", expand=False)
    expected = Series(["A", "B", np.nan], dtype=any_string_dtype)
    tm.assert_series_equal(result, expected)

    # two groups, some matches
    result = s.str.extract("([AB])([123])", expand=False)
    expected = DataFrame(
        [["A", "1"], ["B", "2"], [np.nan, np.nan]], dtype=any_string_dtype
    )
    tm.assert_frame_equal(result, expected)

    # one named group
    result = s.str.extract("(?P<letter>[AB])", expand=False)
    expected = Series(["A", "B", np.nan], name="letter", dtype=any_string_dtype)
    tm.assert_series_equal(result, expected)

    # two named groups
    result = s.str.extract("(?P<letter>[AB])(?P<number>[123])", expand=False)
    expected = DataFrame(
        [["A", "1"], ["B", "2"], [np.nan, np.nan]],
        columns=["letter", "number"],
        dtype=any_string_dtype,
    )
    tm.assert_frame_equal(result, expected)

    # mix named and unnamed groups
    result = s.str.extract("([AB])(?P<number>[123])", expand=False)
    expected = DataFrame(
        [["A", "1"], ["B", "2"], [np.nan, np.nan]],
        columns=[0, "number"],
        dtype=any_string_dtype,
    )
    tm.assert_frame_equal(result, expected)

    # one normal group, one non-capturing group
    result = s.str.extract("([AB])(?:[123])", expand=False)
    expected = Series(["A", "B", np.nan], dtype=any_string_dtype)
    tm.assert_series_equal(result, expected)

    # two normal groups, one non-capturing group
    s = Series(["A11", "B22", "C33"], dtype=any_string_dtype)
    result = s.str.extract("([AB])([123])(?:[123])", expand=False)
    expected = DataFrame(
        [["A", "1"], ["B", "2"], [np.nan, np.nan]], dtype=any_string_dtype
    )
    tm.assert_frame_equal(result, expected)

    # one optional group followed by one normal group
    s = Series(["A1", "B2", "3"], dtype=any_string_dtype)
    result = s.str.extract("(?P<letter>[AB])?(?P<number>[123])", expand=False)
    expected = DataFrame(
        [["A", "1"], ["B", "2"], [np.nan, "3"]],
        columns=["letter", "number"],
        dtype=any_string_dtype,
    )
    tm.assert_frame_equal(result, expected)

    # one normal group followed by one optional group
    s = Series(["A1", "B2", "C"], dtype=any_string_dtype)
    result = s.str.extract("(?P<letter>[ABC])(?P<number>[123])?", expand=False)
    expected = DataFrame(
        [["A", "1"], ["B", "2"], ["C", np.nan]],
        columns=["letter", "number"],
        dtype=any_string_dtype,
    )
    tm.assert_frame_equal(result, expected)


def test_extract_expand_capture_groups_index(index, any_string_dtype):
    # https://github.com/pandas-dev/pandas/issues/6348
    # not passing index to the extractor
    data = ["A1", "B2", "C"]

    if len(index) < len(data):
        pytest.skip("Index too short")

    index = index[: len(data)]
    s = Series(data, index=index, dtype=any_string_dtype)

    result = s.str.extract(r"(\d)", expand=False)
    expected = Series(["1", "2", np.nan], index=index, dtype=any_string_dtype)
    tm.assert_series_equal(result, expected)

    result = s.str.extract(r"(?P<letter>\D)(?P<number>\d)?", expand=False)
    expected = DataFrame(
        [["A", "1"], ["B", "2"], ["C", np.nan]],
        columns=["letter", "number"],
        index=index,
        dtype=any_string_dtype,
    )
    tm.assert_frame_equal(result, expected)


def test_extract_single_series_name_is_preserved(any_string_dtype):
    s = Series(["a3", "b3", "c2"], name="bob", dtype=any_string_dtype)
    result = s.str.extract(r"(?P<sue>[a-z])", expand=False)
    expected = Series(["a", "b", "c"], name="sue", dtype=any_string_dtype)
    tm.assert_series_equal(result, expected)


def test_extract_expand_True(any_string_dtype):
    # Contains tests like those in test_match and some others.
    s = Series(["fooBAD__barBAD", np.nan, "foo"], dtype=any_string_dtype)

    result = s.str.extract(".*(BAD[_]+).*(BAD)", expand=True)
    expected = DataFrame(
        [["BAD__", "BAD"], [np.nan, np.nan], [np.nan, np.nan]], dtype=any_string_dtype
    )
    tm.assert_frame_equal(result, expected)


def test_extract_expand_True_mixed_object():
    er = [np.nan, np.nan]  # empty row
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

    result = mixed.str.extract(".*(BAD[_]+).*(BAD)", expand=True)
    expected = DataFrame([["BAD_", "BAD"], er, ["BAD_", "BAD"], er, er, er, er, er, er])
    tm.assert_frame_equal(result, expected)


def test_extract_expand_True_single_capture_group_raises(
    index_or_series, any_string_dtype
):
    # these should work for both Series and Index
    # no groups
    s_or_idx = index_or_series(["A1", "B2", "C3"], dtype=any_string_dtype)
    msg = "pattern contains no capture groups"
    with pytest.raises(ValueError, match=msg):
        s_or_idx.str.extract("[ABC][123]", expand=True)

    # only non-capturing groups
    with pytest.raises(ValueError, match=msg):
        s_or_idx.str.extract("(?:[AB]).*", expand=True)


def test_extract_expand_True_single_capture_group(index_or_series, any_string_dtype):
    # single group renames series/index properly
    s_or_idx = index_or_series(["A1", "A2"], dtype=any_string_dtype)
    result = s_or_idx.str.extract(r"(?P<uno>A)\d", expand=True)
    expected_dtype = "object" if index_or_series is Index else any_string_dtype
    expected = DataFrame({"uno": ["A", "A"]}, dtype=expected_dtype)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("name", [None, "series_name"])
def test_extract_series(name, any_string_dtype):
    # extract should give the same result whether or not the series has a name.
    s = Series(["A1", "B2", "C3"], name=name, dtype=any_string_dtype)

    # one group, no matches
    result = s.str.extract("(_)", expand=True)
    expected = DataFrame([np.nan, np.nan, np.nan], dtype=any_string_dtype)
    tm.assert_frame_equal(result, expected)

    # two groups, no matches
    result = s.str.extract("(_)(_)", expand=True)
    expected = DataFrame(
        [[np.nan, np.nan], [np.nan, np.nan], [np.nan, np.nan]], dtype=any_string_dtype
    )
    tm.assert_frame_equal(result, expected)

    # one group, some matches
    result = s.str.extract("([AB])[123]", expand=True)
    expected = DataFrame(["A", "B", np.nan], dtype=any_string_dtype)
    tm.assert_frame_equal(result, expected)

    # two groups, some matches
    result = s.str.extract("([AB])([123])", expand=True)
    expected = DataFrame(
        [["A", "1"], ["B", "2"], [np.nan, np.nan]], dtype=any_string_dtype
    )
    tm.assert_frame_equal(result, expected)

    # one named group
    result = s.str.extract("(?P<letter>[AB])", expand=True)
    expected = DataFrame({"letter": ["A", "B", np.nan]}, dtype=any_string_dtype)
    tm.assert_frame_equal(result, expected)

    # two named groups
    result = s.str.extract("(?P<letter>[AB])(?P<number>[123])", expand=True)
    expected = DataFrame(
        [["A", "1"], ["B", "2"], [np.nan, np.nan]],
        columns=["letter", "number"],
        dtype=any_string_dtype,
    )
    tm.assert_frame_equal(result, expected)

    # mix named and unnamed groups
    result = s.str.extract("([AB])(?P<number>[123])", expand=True)
    expected = DataFrame(
        [["A", "1"], ["B", "2"], [np.nan, np.nan]],
        columns=[0, "number"],
        dtype=any_string_dtype,
    )
    tm.assert_frame_equal(result, expected)

    # one normal group, one non-capturing group
    result = s.str.extract("([AB])(?:[123])", expand=True)
    expected = DataFrame(["A", "B", np.nan], dtype=any_string_dtype)
    tm.assert_frame_equal(result, expected)


def test_extract_optional_groups(any_string_dtype):

    # two normal groups, one non-capturing group
    s = Series(["A11", "B22", "C33"], dtype=any_string_dtype)
    result = s.str.extract("([AB])([123])(?:[123])", expand=True)
    expected = DataFrame(
        [["A", "1"], ["B", "2"], [np.nan, np.nan]], dtype=any_string_dtype
    )
    tm.assert_frame_equal(result, expected)

    # one optional group followed by one normal group
    s = Series(["A1", "B2", "3"], dtype=any_string_dtype)
    result = s.str.extract("(?P<letter>[AB])?(?P<number>[123])", expand=True)
    expected = DataFrame(
        [["A", "1"], ["B", "2"], [np.nan, "3"]],
        columns=["letter", "number"],
        dtype=any_string_dtype,
    )
    tm.assert_frame_equal(result, expected)

    # one normal group followed by one optional group
    s = Series(["A1", "B2", "C"], dtype=any_string_dtype)
    result = s.str.extract("(?P<letter>[ABC])(?P<number>[123])?", expand=True)
    expected = DataFrame(
        [["A", "1"], ["B", "2"], ["C", np.nan]],
        columns=["letter", "number"],
        dtype=any_string_dtype,
    )
    tm.assert_frame_equal(result, expected)


def test_extract_dataframe_capture_groups_index(index, any_string_dtype):
    # GH6348
    # not passing index to the extractor

    data = ["A1", "B2", "C"]

    if len(index) < len(data):
        pytest.skip("Index too short")

    index = index[: len(data)]
    s = Series(data, index=index, dtype=any_string_dtype)

    result = s.str.extract(r"(\d)", expand=True)
    expected = DataFrame(["1", "2", np.nan], index=index, dtype=any_string_dtype)
    tm.assert_frame_equal(result, expected)

    result = s.str.extract(r"(?P<letter>\D)(?P<number>\d)?", expand=True)
    expected = DataFrame(
        [["A", "1"], ["B", "2"], ["C", np.nan]],
        columns=["letter", "number"],
        index=index,
        dtype=any_string_dtype,
    )
    tm.assert_frame_equal(result, expected)


def test_extract_single_group_returns_frame(any_string_dtype):
    # GH11386 extract should always return DataFrame, even when
    # there is only one group. Prior to v0.18.0, extract returned
    # Series when there was only one group in the regex.
    s = Series(["a3", "b3", "c2"], name="series_name", dtype=any_string_dtype)
    result = s.str.extract(r"(?P<letter>[a-z])", expand=True)
    expected = DataFrame({"letter": ["a", "b", "c"]}, dtype=any_string_dtype)
    tm.assert_frame_equal(result, expected)


def test_extractall():
    subject_list = [
        "dave@google.com",
        "tdhock5@gmail.com",
        "maudelaperriere@gmail.com",
        "rob@gmail.com some text steve@gmail.com",
        "a@b.com some text c@d.com and e@f.com",
        np.nan,
        "",
    ]
    expected_tuples = [
        ("dave", "google", "com"),
        ("tdhock5", "gmail", "com"),
        ("maudelaperriere", "gmail", "com"),
        ("rob", "gmail", "com"),
        ("steve", "gmail", "com"),
        ("a", "b", "com"),
        ("c", "d", "com"),
        ("e", "f", "com"),
    ]
    named_pattern = r"""
    (?P<user>[a-z0-9]+)
    @
    (?P<domain>[a-z]+)
    \.
    (?P<tld>[a-z]{2,4})
    """
    expected_columns = ["user", "domain", "tld"]
    S = Series(subject_list)
    # extractall should return a DataFrame with one row for each
    # match, indexed by the subject from which the match came.
    expected_index = MultiIndex.from_tuples(
        [(0, 0), (1, 0), (2, 0), (3, 0), (3, 1), (4, 0), (4, 1), (4, 2)],
        names=(None, "match"),
    )
    expected_df = DataFrame(expected_tuples, expected_index, expected_columns)
    computed_df = S.str.extractall(named_pattern, re.VERBOSE)
    tm.assert_frame_equal(computed_df, expected_df)

    # The index of the input Series should be used to construct
    # the index of the output DataFrame:
    series_index = MultiIndex.from_tuples(
        [
            ("single", "Dave"),
            ("single", "Toby"),
            ("single", "Maude"),
            ("multiple", "robAndSteve"),
            ("multiple", "abcdef"),
            ("none", "missing"),
            ("none", "empty"),
        ]
    )
    Si = Series(subject_list, series_index)
    expected_index = MultiIndex.from_tuples(
        [
            ("single", "Dave", 0),
            ("single", "Toby", 0),
            ("single", "Maude", 0),
            ("multiple", "robAndSteve", 0),
            ("multiple", "robAndSteve", 1),
            ("multiple", "abcdef", 0),
            ("multiple", "abcdef", 1),
            ("multiple", "abcdef", 2),
        ],
        names=(None, None, "match"),
    )
    expected_df = DataFrame(expected_tuples, expected_index, expected_columns)
    computed_df = Si.str.extractall(named_pattern, re.VERBOSE)
    tm.assert_frame_equal(computed_df, expected_df)

    # MultiIndexed subject with names.
    Sn = Series(subject_list, series_index)
    Sn.index.names = ("matches", "description")
    expected_index.names = ("matches", "description", "match")
    expected_df = DataFrame(expected_tuples, expected_index, expected_columns)
    computed_df = Sn.str.extractall(named_pattern, re.VERBOSE)
    tm.assert_frame_equal(computed_df, expected_df)

    # optional groups.
    subject_list = ["", "A1", "32"]
    named_pattern = "(?P<letter>[AB])?(?P<number>[123])"
    computed_df = Series(subject_list).str.extractall(named_pattern)
    expected_index = MultiIndex.from_tuples(
        [(1, 0), (2, 0), (2, 1)], names=(None, "match")
    )
    expected_df = DataFrame(
        [("A", "1"), (np.nan, "3"), (np.nan, "2")],
        expected_index,
        columns=["letter", "number"],
    )
    tm.assert_frame_equal(computed_df, expected_df)

    # only one of two groups has a name.
    pattern = "([AB])?(?P<number>[123])"
    computed_df = Series(subject_list).str.extractall(pattern)
    expected_df = DataFrame(
        [("A", "1"), (np.nan, "3"), (np.nan, "2")],
        expected_index,
        columns=[0, "number"],
    )
    tm.assert_frame_equal(computed_df, expected_df)


def test_extractall_single_group():
    # extractall(one named group) returns DataFrame with one named
    # column.
    s = Series(["a3", "b3", "d4c2"], name="series_name")
    r = s.str.extractall(r"(?P<letter>[a-z])")
    i = MultiIndex.from_tuples([(0, 0), (1, 0), (2, 0), (2, 1)], names=(None, "match"))
    e = DataFrame({"letter": ["a", "b", "d", "c"]}, i)
    tm.assert_frame_equal(r, e)

    # extractall(one un-named group) returns DataFrame with one
    # un-named column.
    r = s.str.extractall(r"([a-z])")
    e = DataFrame(["a", "b", "d", "c"], i)
    tm.assert_frame_equal(r, e)


def test_extractall_single_group_with_quantifier():
    # extractall(one un-named group with quantifier) returns
    # DataFrame with one un-named column (GH13382).
    s = Series(["ab3", "abc3", "d4cd2"], name="series_name")
    r = s.str.extractall(r"([a-z]+)")
    i = MultiIndex.from_tuples([(0, 0), (1, 0), (2, 0), (2, 1)], names=(None, "match"))
    e = DataFrame(["ab", "abc", "d", "cd"], i)
    tm.assert_frame_equal(r, e)


@pytest.mark.parametrize(
    "data, names",
    [
        ([], (None,)),
        ([], ("i1",)),
        ([], (None, "i2")),
        ([], ("i1", "i2")),
        (["a3", "b3", "d4c2"], (None,)),
        (["a3", "b3", "d4c2"], ("i1", "i2")),
        (["a3", "b3", "d4c2"], (None, "i2")),
        (["a3", "b3", "d4c2"], ("i1", "i2")),
    ],
)
def test_extractall_no_matches(data, names):
    # GH19075 extractall with no matches should return a valid MultiIndex
    n = len(data)
    if len(names) == 1:
        i = Index(range(n), name=names[0])
    else:
        a = (tuple([i] * (n - 1)) for i in range(n))
        i = MultiIndex.from_tuples(a, names=names)
    s = Series(data, name="series_name", index=i, dtype="object")
    ei = MultiIndex.from_tuples([], names=(names + ("match",)))

    # one un-named group.
    r = s.str.extractall("(z)")
    e = DataFrame(columns=[0], index=ei)
    tm.assert_frame_equal(r, e)

    # two un-named groups.
    r = s.str.extractall("(z)(z)")
    e = DataFrame(columns=[0, 1], index=ei)
    tm.assert_frame_equal(r, e)

    # one named group.
    r = s.str.extractall("(?P<first>z)")
    e = DataFrame(columns=["first"], index=ei)
    tm.assert_frame_equal(r, e)

    # two named groups.
    r = s.str.extractall("(?P<first>z)(?P<second>z)")
    e = DataFrame(columns=["first", "second"], index=ei)
    tm.assert_frame_equal(r, e)

    # one named, one un-named.
    r = s.str.extractall("(z)(?P<second>z)")
    e = DataFrame(columns=[0, "second"], index=ei)
    tm.assert_frame_equal(r, e)


def test_extractall_stringindex():
    s = Series(["a1a2", "b1", "c1"], name="xxx")
    res = s.str.extractall(r"[ab](?P<digit>\d)")
    exp_idx = MultiIndex.from_tuples([(0, 0), (0, 1), (1, 0)], names=[None, "match"])
    exp = DataFrame({"digit": ["1", "2", "1"]}, index=exp_idx)
    tm.assert_frame_equal(res, exp)

    # index should return the same result as the default index without name
    # thus index.name doesn't affect to the result
    for idx in [
        Index(["a1a2", "b1", "c1"]),
        Index(["a1a2", "b1", "c1"], name="xxx"),
    ]:

        res = idx.str.extractall(r"[ab](?P<digit>\d)")
        tm.assert_frame_equal(res, exp)

    s = Series(
        ["a1a2", "b1", "c1"],
        name="s_name",
        index=Index(["XX", "yy", "zz"], name="idx_name"),
    )
    res = s.str.extractall(r"[ab](?P<digit>\d)")
    exp_idx = MultiIndex.from_tuples(
        [("XX", 0), ("XX", 1), ("yy", 0)], names=["idx_name", "match"]
    )
    exp = DataFrame({"digit": ["1", "2", "1"]}, index=exp_idx)
    tm.assert_frame_equal(res, exp)


def test_extractall_errors():
    # Does not make sense to use extractall with a regex that has
    # no capture groups. (it returns DataFrame with one column for
    # each capture group)
    s = Series(["a3", "b3", "d4c2"], name="series_name")
    with pytest.raises(ValueError, match="no capture groups"):
        s.str.extractall(r"[a-z]")


def test_extract_index_one_two_groups():
    s = Series(["a3", "b3", "d4c2"], index=["A3", "B3", "D4"], name="series_name")
    r = s.index.str.extract(r"([A-Z])", expand=True)
    e = DataFrame(["A", "B", "D"])
    tm.assert_frame_equal(r, e)

    # Prior to v0.18.0, index.str.extract(regex with one group)
    # returned Index. With more than one group, extract raised an
    # error (GH9980). Now extract always returns DataFrame.
    r = s.index.str.extract(r"(?P<letter>[A-Z])(?P<digit>[0-9])", expand=True)
    e_list = [("A", "3"), ("B", "3"), ("D", "4")]
    e = DataFrame(e_list, columns=["letter", "digit"])
    tm.assert_frame_equal(r, e)


def test_extractall_same_as_extract():
    s = Series(["a3", "b3", "c2"], name="series_name")

    pattern_two_noname = r"([a-z])([0-9])"
    extract_two_noname = s.str.extract(pattern_two_noname, expand=True)
    has_multi_index = s.str.extractall(pattern_two_noname)
    no_multi_index = has_multi_index.xs(0, level="match")
    tm.assert_frame_equal(extract_two_noname, no_multi_index)

    pattern_two_named = r"(?P<letter>[a-z])(?P<digit>[0-9])"
    extract_two_named = s.str.extract(pattern_two_named, expand=True)
    has_multi_index = s.str.extractall(pattern_two_named)
    no_multi_index = has_multi_index.xs(0, level="match")
    tm.assert_frame_equal(extract_two_named, no_multi_index)

    pattern_one_named = r"(?P<group_name>[a-z])"
    extract_one_named = s.str.extract(pattern_one_named, expand=True)
    has_multi_index = s.str.extractall(pattern_one_named)
    no_multi_index = has_multi_index.xs(0, level="match")
    tm.assert_frame_equal(extract_one_named, no_multi_index)

    pattern_one_noname = r"([a-z])"
    extract_one_noname = s.str.extract(pattern_one_noname, expand=True)
    has_multi_index = s.str.extractall(pattern_one_noname)
    no_multi_index = has_multi_index.xs(0, level="match")
    tm.assert_frame_equal(extract_one_noname, no_multi_index)


def test_extractall_same_as_extract_subject_index():
    # same as above tests, but s has an MultiIndex.
    i = MultiIndex.from_tuples(
        [("A", "first"), ("B", "second"), ("C", "third")],
        names=("capital", "ordinal"),
    )
    s = Series(["a3", "b3", "c2"], i, name="series_name")

    pattern_two_noname = r"([a-z])([0-9])"
    extract_two_noname = s.str.extract(pattern_two_noname, expand=True)
    has_match_index = s.str.extractall(pattern_two_noname)
    no_match_index = has_match_index.xs(0, level="match")
    tm.assert_frame_equal(extract_two_noname, no_match_index)

    pattern_two_named = r"(?P<letter>[a-z])(?P<digit>[0-9])"
    extract_two_named = s.str.extract(pattern_two_named, expand=True)
    has_match_index = s.str.extractall(pattern_two_named)
    no_match_index = has_match_index.xs(0, level="match")
    tm.assert_frame_equal(extract_two_named, no_match_index)

    pattern_one_named = r"(?P<group_name>[a-z])"
    extract_one_named = s.str.extract(pattern_one_named, expand=True)
    has_match_index = s.str.extractall(pattern_one_named)
    no_match_index = has_match_index.xs(0, level="match")
    tm.assert_frame_equal(extract_one_named, no_match_index)

    pattern_one_noname = r"([a-z])"
    extract_one_noname = s.str.extract(pattern_one_noname, expand=True)
    has_match_index = s.str.extractall(pattern_one_noname)
    no_match_index = has_match_index.xs(0, level="match")
    tm.assert_frame_equal(extract_one_noname, no_match_index)
