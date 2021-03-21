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


def test_extract_expand_None():
    values = Series(["fooBAD__barBAD", np.nan, "foo"])
    with pytest.raises(ValueError, match="expand must be True or False"):
        values.str.extract(".*(BAD[_]+).*(BAD)", expand=None)


def test_extract_expand_unspecified():
    values = Series(["fooBAD__barBAD", np.nan, "foo"])
    result_unspecified = values.str.extract(".*(BAD[_]+).*")
    assert isinstance(result_unspecified, DataFrame)
    result_true = values.str.extract(".*(BAD[_]+).*", expand=True)
    tm.assert_frame_equal(result_unspecified, result_true)


def test_extract_expand_False():
    # Contains tests like those in test_match and some others.
    values = Series(["fooBAD__barBAD", np.nan, "foo"])
    er = [np.nan, np.nan]  # empty row

    result = values.str.extract(".*(BAD[_]+).*(BAD)", expand=False)
    exp = DataFrame([["BAD__", "BAD"], er, er])
    tm.assert_frame_equal(result, exp)

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

    rs = Series(mixed).str.extract(".*(BAD[_]+).*(BAD)", expand=False)
    exp = DataFrame([["BAD_", "BAD"], er, ["BAD_", "BAD"], er, er, er, er, er, er])
    tm.assert_frame_equal(rs, exp)

    # unicode
    values = Series(["fooBAD__barBAD", np.nan, "foo"])

    result = values.str.extract(".*(BAD[_]+).*(BAD)", expand=False)
    exp = DataFrame([["BAD__", "BAD"], er, er])
    tm.assert_frame_equal(result, exp)

    # GH9980
    # Index only works with one regex group since
    # multi-group would expand to a frame
    idx = Index(["A1", "A2", "A3", "A4", "B5"])
    with pytest.raises(ValueError, match="supported"):
        idx.str.extract("([AB])([123])", expand=False)

    # these should work for both Series and Index
    for klass in [Series, Index]:
        # no groups
        s_or_idx = klass(["A1", "B2", "C3"])
        msg = "pattern contains no capture groups"
        with pytest.raises(ValueError, match=msg):
            s_or_idx.str.extract("[ABC][123]", expand=False)

        # only non-capturing groups
        with pytest.raises(ValueError, match=msg):
            s_or_idx.str.extract("(?:[AB]).*", expand=False)

        # single group renames series/index properly
        s_or_idx = klass(["A1", "A2"])
        result = s_or_idx.str.extract(r"(?P<uno>A)\d", expand=False)
        assert result.name == "uno"

        exp = klass(["A", "A"], name="uno")
        if klass == Series:
            tm.assert_series_equal(result, exp)
        else:
            tm.assert_index_equal(result, exp)

    s = Series(["A1", "B2", "C3"])
    # one group, no matches
    result = s.str.extract("(_)", expand=False)
    exp = Series([np.nan, np.nan, np.nan], dtype=object)
    tm.assert_series_equal(result, exp)

    # two groups, no matches
    result = s.str.extract("(_)(_)", expand=False)
    exp = DataFrame(
        [[np.nan, np.nan], [np.nan, np.nan], [np.nan, np.nan]], dtype=object
    )
    tm.assert_frame_equal(result, exp)

    # one group, some matches
    result = s.str.extract("([AB])[123]", expand=False)
    exp = Series(["A", "B", np.nan])
    tm.assert_series_equal(result, exp)

    # two groups, some matches
    result = s.str.extract("([AB])([123])", expand=False)
    exp = DataFrame([["A", "1"], ["B", "2"], [np.nan, np.nan]])
    tm.assert_frame_equal(result, exp)

    # one named group
    result = s.str.extract("(?P<letter>[AB])", expand=False)
    exp = Series(["A", "B", np.nan], name="letter")
    tm.assert_series_equal(result, exp)

    # two named groups
    result = s.str.extract("(?P<letter>[AB])(?P<number>[123])", expand=False)
    exp = DataFrame(
        [["A", "1"], ["B", "2"], [np.nan, np.nan]], columns=["letter", "number"]
    )
    tm.assert_frame_equal(result, exp)

    # mix named and unnamed groups
    result = s.str.extract("([AB])(?P<number>[123])", expand=False)
    exp = DataFrame([["A", "1"], ["B", "2"], [np.nan, np.nan]], columns=[0, "number"])
    tm.assert_frame_equal(result, exp)

    # one normal group, one non-capturing group
    result = s.str.extract("([AB])(?:[123])", expand=False)
    exp = Series(["A", "B", np.nan])
    tm.assert_series_equal(result, exp)

    # two normal groups, one non-capturing group
    result = Series(["A11", "B22", "C33"]).str.extract(
        "([AB])([123])(?:[123])", expand=False
    )
    exp = DataFrame([["A", "1"], ["B", "2"], [np.nan, np.nan]])
    tm.assert_frame_equal(result, exp)

    # one optional group followed by one normal group
    result = Series(["A1", "B2", "3"]).str.extract(
        "(?P<letter>[AB])?(?P<number>[123])", expand=False
    )
    exp = DataFrame(
        [["A", "1"], ["B", "2"], [np.nan, "3"]], columns=["letter", "number"]
    )
    tm.assert_frame_equal(result, exp)

    # one normal group followed by one optional group
    result = Series(["A1", "B2", "C"]).str.extract(
        "(?P<letter>[ABC])(?P<number>[123])?", expand=False
    )
    exp = DataFrame(
        [["A", "1"], ["B", "2"], ["C", np.nan]], columns=["letter", "number"]
    )
    tm.assert_frame_equal(result, exp)

    # GH6348
    # not passing index to the extractor
    def check_index(index):
        data = ["A1", "B2", "C"]
        index = index[: len(data)]
        s = Series(data, index=index)
        result = s.str.extract(r"(\d)", expand=False)
        exp = Series(["1", "2", np.nan], index=index)
        tm.assert_series_equal(result, exp)

        result = Series(data, index=index).str.extract(
            r"(?P<letter>\D)(?P<number>\d)?", expand=False
        )
        e_list = [["A", "1"], ["B", "2"], ["C", np.nan]]
        exp = DataFrame(e_list, columns=["letter", "number"], index=index)
        tm.assert_frame_equal(result, exp)

    i_funs = [
        tm.makeStringIndex,
        tm.makeUnicodeIndex,
        tm.makeIntIndex,
        tm.makeDateIndex,
        tm.makePeriodIndex,
        tm.makeRangeIndex,
    ]
    for index in i_funs:
        check_index(index())

    # single_series_name_is_preserved.
    s = Series(["a3", "b3", "c2"], name="bob")
    r = s.str.extract(r"(?P<sue>[a-z])", expand=False)
    e = Series(["a", "b", "c"], name="sue")
    tm.assert_series_equal(r, e)
    assert r.name == e.name


def test_extract_expand_True():
    # Contains tests like those in test_match and some others.
    values = Series(["fooBAD__barBAD", np.nan, "foo"])
    er = [np.nan, np.nan]  # empty row

    result = values.str.extract(".*(BAD[_]+).*(BAD)", expand=True)
    exp = DataFrame([["BAD__", "BAD"], er, er])
    tm.assert_frame_equal(result, exp)

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

    rs = Series(mixed).str.extract(".*(BAD[_]+).*(BAD)", expand=True)
    exp = DataFrame([["BAD_", "BAD"], er, ["BAD_", "BAD"], er, er, er, er, er, er])
    tm.assert_frame_equal(rs, exp)

    # these should work for both Series and Index
    for klass in [Series, Index]:
        # no groups
        s_or_idx = klass(["A1", "B2", "C3"])
        msg = "pattern contains no capture groups"
        with pytest.raises(ValueError, match=msg):
            s_or_idx.str.extract("[ABC][123]", expand=True)

        # only non-capturing groups
        with pytest.raises(ValueError, match=msg):
            s_or_idx.str.extract("(?:[AB]).*", expand=True)

        # single group renames series/index properly
        s_or_idx = klass(["A1", "A2"])
        result_df = s_or_idx.str.extract(r"(?P<uno>A)\d", expand=True)
        assert isinstance(result_df, DataFrame)
        result_series = result_df["uno"]
        tm.assert_series_equal(result_series, Series(["A", "A"], name="uno"))


def test_extract_series():
    # extract should give the same result whether or not the
    # series has a name.
    for series_name in None, "series_name":
        s = Series(["A1", "B2", "C3"], name=series_name)
        # one group, no matches
        result = s.str.extract("(_)", expand=True)
        exp = DataFrame([np.nan, np.nan, np.nan], dtype=object)
        tm.assert_frame_equal(result, exp)

        # two groups, no matches
        result = s.str.extract("(_)(_)", expand=True)
        exp = DataFrame(
            [[np.nan, np.nan], [np.nan, np.nan], [np.nan, np.nan]], dtype=object
        )
        tm.assert_frame_equal(result, exp)

        # one group, some matches
        result = s.str.extract("([AB])[123]", expand=True)
        exp = DataFrame(["A", "B", np.nan])
        tm.assert_frame_equal(result, exp)

        # two groups, some matches
        result = s.str.extract("([AB])([123])", expand=True)
        exp = DataFrame([["A", "1"], ["B", "2"], [np.nan, np.nan]])
        tm.assert_frame_equal(result, exp)

        # one named group
        result = s.str.extract("(?P<letter>[AB])", expand=True)
        exp = DataFrame({"letter": ["A", "B", np.nan]})
        tm.assert_frame_equal(result, exp)

        # two named groups
        result = s.str.extract("(?P<letter>[AB])(?P<number>[123])", expand=True)
        e_list = [["A", "1"], ["B", "2"], [np.nan, np.nan]]
        exp = DataFrame(e_list, columns=["letter", "number"])
        tm.assert_frame_equal(result, exp)

        # mix named and unnamed groups
        result = s.str.extract("([AB])(?P<number>[123])", expand=True)
        exp = DataFrame(e_list, columns=[0, "number"])
        tm.assert_frame_equal(result, exp)

        # one normal group, one non-capturing group
        result = s.str.extract("([AB])(?:[123])", expand=True)
        exp = DataFrame(["A", "B", np.nan])
        tm.assert_frame_equal(result, exp)


def test_extract_optional_groups():

    # two normal groups, one non-capturing group
    result = Series(["A11", "B22", "C33"]).str.extract(
        "([AB])([123])(?:[123])", expand=True
    )
    exp = DataFrame([["A", "1"], ["B", "2"], [np.nan, np.nan]])
    tm.assert_frame_equal(result, exp)

    # one optional group followed by one normal group
    result = Series(["A1", "B2", "3"]).str.extract(
        "(?P<letter>[AB])?(?P<number>[123])", expand=True
    )
    e_list = [["A", "1"], ["B", "2"], [np.nan, "3"]]
    exp = DataFrame(e_list, columns=["letter", "number"])
    tm.assert_frame_equal(result, exp)

    # one normal group followed by one optional group
    result = Series(["A1", "B2", "C"]).str.extract(
        "(?P<letter>[ABC])(?P<number>[123])?", expand=True
    )
    e_list = [["A", "1"], ["B", "2"], ["C", np.nan]]
    exp = DataFrame(e_list, columns=["letter", "number"])
    tm.assert_frame_equal(result, exp)

    # GH6348
    # not passing index to the extractor
    def check_index(index):
        data = ["A1", "B2", "C"]
        index = index[: len(data)]
        result = Series(data, index=index).str.extract(r"(\d)", expand=True)
        exp = DataFrame(["1", "2", np.nan], index=index)
        tm.assert_frame_equal(result, exp)

        result = Series(data, index=index).str.extract(
            r"(?P<letter>\D)(?P<number>\d)?", expand=True
        )
        e_list = [["A", "1"], ["B", "2"], ["C", np.nan]]
        exp = DataFrame(e_list, columns=["letter", "number"], index=index)
        tm.assert_frame_equal(result, exp)

    i_funs = [
        tm.makeStringIndex,
        tm.makeUnicodeIndex,
        tm.makeIntIndex,
        tm.makeDateIndex,
        tm.makePeriodIndex,
        tm.makeRangeIndex,
    ]
    for index in i_funs:
        check_index(index())


def test_extract_single_group_returns_frame():
    # GH11386 extract should always return DataFrame, even when
    # there is only one group. Prior to v0.18.0, extract returned
    # Series when there was only one group in the regex.
    s = Series(["a3", "b3", "c2"], name="series_name")
    r = s.str.extract(r"(?P<letter>[a-z])", expand=True)
    e = DataFrame({"letter": ["a", "b", "c"]})
    tm.assert_frame_equal(r, e)


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
