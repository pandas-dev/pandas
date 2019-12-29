from datetime import datetime
from io import StringIO
import re
from typing import Dict, List, Union

import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Index, Series, Timestamp, date_range
import pandas.util.testing as tm


@pytest.fixture
def mix_ab() -> Dict[str, List[Union[int, str]]]:
    return {"a": list(range(4)), "b": list("ab..")}


@pytest.fixture
def mix_abc() -> Dict[str, List[Union[float, str]]]:
    return {"a": list(range(4)), "b": list("ab.."), "c": ["a", "b", np.nan, "d"]}


class TestDataFrameReplace:
    def test_replace_inplace(self, datetime_frame, float_string_frame):
        datetime_frame["A"][:5] = np.nan
        datetime_frame["A"][-5:] = np.nan

        tsframe = datetime_frame.copy()
        tsframe.replace(np.nan, 0, inplace=True)
        tm.assert_frame_equal(tsframe, datetime_frame.fillna(0))

        # mixed type
        mf = float_string_frame
        mf.iloc[5:20, mf.columns.get_loc("foo")] = np.nan
        mf.iloc[-10:, mf.columns.get_loc("A")] = np.nan

        result = float_string_frame.replace(np.nan, 0)
        expected = float_string_frame.fillna(value=0)
        tm.assert_frame_equal(result, expected)

        tsframe = datetime_frame.copy()
        tsframe.replace([np.nan], [0], inplace=True)
        tm.assert_frame_equal(tsframe, datetime_frame.fillna(0))

    def test_regex_replace_scalar(self, mix_ab):
        obj = {"a": list("ab.."), "b": list("efgh")}
        dfobj = DataFrame(obj)
        dfmix = DataFrame(mix_ab)

        # simplest cases
        # regex -> value
        # obj frame
        res = dfobj.replace(r"\s*\.\s*", np.nan, regex=True)
        tm.assert_frame_equal(dfobj, res.fillna("."))

        # mixed
        res = dfmix.replace(r"\s*\.\s*", np.nan, regex=True)
        tm.assert_frame_equal(dfmix, res.fillna("."))

        # regex -> regex
        # obj frame
        res = dfobj.replace(r"\s*(\.)\s*", r"\1\1\1", regex=True)
        objc = obj.copy()
        objc["a"] = ["a", "b", "...", "..."]
        expec = DataFrame(objc)
        tm.assert_frame_equal(res, expec)

        # with mixed
        res = dfmix.replace(r"\s*(\.)\s*", r"\1\1\1", regex=True)
        mixc = mix_ab.copy()
        mixc["b"] = ["a", "b", "...", "..."]
        expec = DataFrame(mixc)
        tm.assert_frame_equal(res, expec)

        # everything with compiled regexs as well
        res = dfobj.replace(re.compile(r"\s*\.\s*"), np.nan, regex=True)
        tm.assert_frame_equal(dfobj, res.fillna("."))

        # mixed
        res = dfmix.replace(re.compile(r"\s*\.\s*"), np.nan, regex=True)
        tm.assert_frame_equal(dfmix, res.fillna("."))

        # regex -> regex
        # obj frame
        res = dfobj.replace(re.compile(r"\s*(\.)\s*"), r"\1\1\1")
        objc = obj.copy()
        objc["a"] = ["a", "b", "...", "..."]
        expec = DataFrame(objc)
        tm.assert_frame_equal(res, expec)

        # with mixed
        res = dfmix.replace(re.compile(r"\s*(\.)\s*"), r"\1\1\1")
        mixc = mix_ab.copy()
        mixc["b"] = ["a", "b", "...", "..."]
        expec = DataFrame(mixc)
        tm.assert_frame_equal(res, expec)

        res = dfmix.replace(regex=re.compile(r"\s*(\.)\s*"), value=r"\1\1\1")
        mixc = mix_ab.copy()
        mixc["b"] = ["a", "b", "...", "..."]
        expec = DataFrame(mixc)
        tm.assert_frame_equal(res, expec)

        res = dfmix.replace(regex=r"\s*(\.)\s*", value=r"\1\1\1")
        mixc = mix_ab.copy()
        mixc["b"] = ["a", "b", "...", "..."]
        expec = DataFrame(mixc)
        tm.assert_frame_equal(res, expec)

    def test_regex_replace_scalar_inplace(self, mix_ab):
        obj = {"a": list("ab.."), "b": list("efgh")}
        dfobj = DataFrame(obj)
        dfmix = DataFrame(mix_ab)

        # simplest cases
        # regex -> value
        # obj frame
        res = dfobj.copy()
        res.replace(r"\s*\.\s*", np.nan, regex=True, inplace=True)
        tm.assert_frame_equal(dfobj, res.fillna("."))

        # mixed
        res = dfmix.copy()
        res.replace(r"\s*\.\s*", np.nan, regex=True, inplace=True)
        tm.assert_frame_equal(dfmix, res.fillna("."))

        # regex -> regex
        # obj frame
        res = dfobj.copy()
        res.replace(r"\s*(\.)\s*", r"\1\1\1", regex=True, inplace=True)
        objc = obj.copy()
        objc["a"] = ["a", "b", "...", "..."]
        expec = DataFrame(objc)
        tm.assert_frame_equal(res, expec)

        # with mixed
        res = dfmix.copy()
        res.replace(r"\s*(\.)\s*", r"\1\1\1", regex=True, inplace=True)
        mixc = mix_ab.copy()
        mixc["b"] = ["a", "b", "...", "..."]
        expec = DataFrame(mixc)
        tm.assert_frame_equal(res, expec)

        # everything with compiled regexs as well
        res = dfobj.copy()
        res.replace(re.compile(r"\s*\.\s*"), np.nan, regex=True, inplace=True)
        tm.assert_frame_equal(dfobj, res.fillna("."))

        # mixed
        res = dfmix.copy()
        res.replace(re.compile(r"\s*\.\s*"), np.nan, regex=True, inplace=True)
        tm.assert_frame_equal(dfmix, res.fillna("."))

        # regex -> regex
        # obj frame
        res = dfobj.copy()
        res.replace(re.compile(r"\s*(\.)\s*"), r"\1\1\1", regex=True, inplace=True)
        objc = obj.copy()
        objc["a"] = ["a", "b", "...", "..."]
        expec = DataFrame(objc)
        tm.assert_frame_equal(res, expec)

        # with mixed
        res = dfmix.copy()
        res.replace(re.compile(r"\s*(\.)\s*"), r"\1\1\1", regex=True, inplace=True)
        mixc = mix_ab.copy()
        mixc["b"] = ["a", "b", "...", "..."]
        expec = DataFrame(mixc)
        tm.assert_frame_equal(res, expec)

        res = dfobj.copy()
        res.replace(regex=r"\s*\.\s*", value=np.nan, inplace=True)
        tm.assert_frame_equal(dfobj, res.fillna("."))

        # mixed
        res = dfmix.copy()
        res.replace(regex=r"\s*\.\s*", value=np.nan, inplace=True)
        tm.assert_frame_equal(dfmix, res.fillna("."))

        # regex -> regex
        # obj frame
        res = dfobj.copy()
        res.replace(regex=r"\s*(\.)\s*", value=r"\1\1\1", inplace=True)
        objc = obj.copy()
        objc["a"] = ["a", "b", "...", "..."]
        expec = DataFrame(objc)
        tm.assert_frame_equal(res, expec)

        # with mixed
        res = dfmix.copy()
        res.replace(regex=r"\s*(\.)\s*", value=r"\1\1\1", inplace=True)
        mixc = mix_ab.copy()
        mixc["b"] = ["a", "b", "...", "..."]
        expec = DataFrame(mixc)
        tm.assert_frame_equal(res, expec)

        # everything with compiled regexs as well
        res = dfobj.copy()
        res.replace(regex=re.compile(r"\s*\.\s*"), value=np.nan, inplace=True)
        tm.assert_frame_equal(dfobj, res.fillna("."))

        # mixed
        res = dfmix.copy()
        res.replace(regex=re.compile(r"\s*\.\s*"), value=np.nan, inplace=True)
        tm.assert_frame_equal(dfmix, res.fillna("."))

        # regex -> regex
        # obj frame
        res = dfobj.copy()
        res.replace(regex=re.compile(r"\s*(\.)\s*"), value=r"\1\1\1", inplace=True)
        objc = obj.copy()
        objc["a"] = ["a", "b", "...", "..."]
        expec = DataFrame(objc)
        tm.assert_frame_equal(res, expec)

        # with mixed
        res = dfmix.copy()
        res.replace(regex=re.compile(r"\s*(\.)\s*"), value=r"\1\1\1", inplace=True)
        mixc = mix_ab.copy()
        mixc["b"] = ["a", "b", "...", "..."]
        expec = DataFrame(mixc)
        tm.assert_frame_equal(res, expec)

    def test_regex_replace_list_obj(self):
        obj = {"a": list("ab.."), "b": list("efgh"), "c": list("helo")}
        dfobj = DataFrame(obj)

        # lists of regexes and values
        # list of [re1, re2, ..., reN] -> [v1, v2, ..., vN]
        to_replace_res = [r"\s*\.\s*", r"e|f|g"]
        values = [np.nan, "crap"]
        res = dfobj.replace(to_replace_res, values, regex=True)
        expec = DataFrame(
            {
                "a": ["a", "b", np.nan, np.nan],
                "b": ["crap"] * 3 + ["h"],
                "c": ["h", "crap", "l", "o"],
            }
        )
        tm.assert_frame_equal(res, expec)

        # list of [re1, re2, ..., reN] -> [re1, re2, .., reN]
        to_replace_res = [r"\s*(\.)\s*", r"(e|f|g)"]
        values = [r"\1\1", r"\1_crap"]
        res = dfobj.replace(to_replace_res, values, regex=True)
        expec = DataFrame(
            {
                "a": ["a", "b", "..", ".."],
                "b": ["e_crap", "f_crap", "g_crap", "h"],
                "c": ["h", "e_crap", "l", "o"],
            }
        )
        tm.assert_frame_equal(res, expec)

        # list of [re1, re2, ..., reN] -> [(re1 or v1), (re2 or v2), ..., (reN
        # or vN)]
        to_replace_res = [r"\s*(\.)\s*", r"e"]
        values = [r"\1\1", r"crap"]
        res = dfobj.replace(to_replace_res, values, regex=True)
        expec = DataFrame(
            {
                "a": ["a", "b", "..", ".."],
                "b": ["crap", "f", "g", "h"],
                "c": ["h", "crap", "l", "o"],
            }
        )
        tm.assert_frame_equal(res, expec)

        to_replace_res = [r"\s*(\.)\s*", r"e"]
        values = [r"\1\1", r"crap"]
        res = dfobj.replace(value=values, regex=to_replace_res)
        expec = DataFrame(
            {
                "a": ["a", "b", "..", ".."],
                "b": ["crap", "f", "g", "h"],
                "c": ["h", "crap", "l", "o"],
            }
        )
        tm.assert_frame_equal(res, expec)

    def test_regex_replace_list_obj_inplace(self):
        # same as above with inplace=True
        # lists of regexes and values
        obj = {"a": list("ab.."), "b": list("efgh"), "c": list("helo")}
        dfobj = DataFrame(obj)

        # lists of regexes and values
        # list of [re1, re2, ..., reN] -> [v1, v2, ..., vN]
        to_replace_res = [r"\s*\.\s*", r"e|f|g"]
        values = [np.nan, "crap"]
        res = dfobj.copy()
        res.replace(to_replace_res, values, inplace=True, regex=True)
        expec = DataFrame(
            {
                "a": ["a", "b", np.nan, np.nan],
                "b": ["crap"] * 3 + ["h"],
                "c": ["h", "crap", "l", "o"],
            }
        )
        tm.assert_frame_equal(res, expec)

        # list of [re1, re2, ..., reN] -> [re1, re2, .., reN]
        to_replace_res = [r"\s*(\.)\s*", r"(e|f|g)"]
        values = [r"\1\1", r"\1_crap"]
        res = dfobj.copy()
        res.replace(to_replace_res, values, inplace=True, regex=True)
        expec = DataFrame(
            {
                "a": ["a", "b", "..", ".."],
                "b": ["e_crap", "f_crap", "g_crap", "h"],
                "c": ["h", "e_crap", "l", "o"],
            }
        )
        tm.assert_frame_equal(res, expec)

        # list of [re1, re2, ..., reN] -> [(re1 or v1), (re2 or v2), ..., (reN
        # or vN)]
        to_replace_res = [r"\s*(\.)\s*", r"e"]
        values = [r"\1\1", r"crap"]
        res = dfobj.copy()
        res.replace(to_replace_res, values, inplace=True, regex=True)
        expec = DataFrame(
            {
                "a": ["a", "b", "..", ".."],
                "b": ["crap", "f", "g", "h"],
                "c": ["h", "crap", "l", "o"],
            }
        )
        tm.assert_frame_equal(res, expec)

        to_replace_res = [r"\s*(\.)\s*", r"e"]
        values = [r"\1\1", r"crap"]
        res = dfobj.copy()
        res.replace(value=values, regex=to_replace_res, inplace=True)
        expec = DataFrame(
            {
                "a": ["a", "b", "..", ".."],
                "b": ["crap", "f", "g", "h"],
                "c": ["h", "crap", "l", "o"],
            }
        )
        tm.assert_frame_equal(res, expec)

    def test_regex_replace_list_mixed(self, mix_ab):
        # mixed frame to make sure this doesn't break things
        dfmix = DataFrame(mix_ab)

        # lists of regexes and values
        # list of [re1, re2, ..., reN] -> [v1, v2, ..., vN]
        to_replace_res = [r"\s*\.\s*", r"a"]
        values = [np.nan, "crap"]
        mix2 = {"a": list(range(4)), "b": list("ab.."), "c": list("halo")}
        dfmix2 = DataFrame(mix2)
        res = dfmix2.replace(to_replace_res, values, regex=True)
        expec = DataFrame(
            {
                "a": mix2["a"],
                "b": ["crap", "b", np.nan, np.nan],
                "c": ["h", "crap", "l", "o"],
            }
        )
        tm.assert_frame_equal(res, expec)

        # list of [re1, re2, ..., reN] -> [re1, re2, .., reN]
        to_replace_res = [r"\s*(\.)\s*", r"(a|b)"]
        values = [r"\1\1", r"\1_crap"]
        res = dfmix.replace(to_replace_res, values, regex=True)
        expec = DataFrame({"a": mix_ab["a"], "b": ["a_crap", "b_crap", "..", ".."]})
        tm.assert_frame_equal(res, expec)

        # list of [re1, re2, ..., reN] -> [(re1 or v1), (re2 or v2), ..., (reN
        # or vN)]
        to_replace_res = [r"\s*(\.)\s*", r"a", r"(b)"]
        values = [r"\1\1", r"crap", r"\1_crap"]
        res = dfmix.replace(to_replace_res, values, regex=True)
        expec = DataFrame({"a": mix_ab["a"], "b": ["crap", "b_crap", "..", ".."]})
        tm.assert_frame_equal(res, expec)

        to_replace_res = [r"\s*(\.)\s*", r"a", r"(b)"]
        values = [r"\1\1", r"crap", r"\1_crap"]
        res = dfmix.replace(regex=to_replace_res, value=values)
        expec = DataFrame({"a": mix_ab["a"], "b": ["crap", "b_crap", "..", ".."]})
        tm.assert_frame_equal(res, expec)

    def test_regex_replace_list_mixed_inplace(self, mix_ab):
        dfmix = DataFrame(mix_ab)
        # the same inplace
        # lists of regexes and values
        # list of [re1, re2, ..., reN] -> [v1, v2, ..., vN]
        to_replace_res = [r"\s*\.\s*", r"a"]
        values = [np.nan, "crap"]
        res = dfmix.copy()
        res.replace(to_replace_res, values, inplace=True, regex=True)
        expec = DataFrame({"a": mix_ab["a"], "b": ["crap", "b", np.nan, np.nan]})
        tm.assert_frame_equal(res, expec)

        # list of [re1, re2, ..., reN] -> [re1, re2, .., reN]
        to_replace_res = [r"\s*(\.)\s*", r"(a|b)"]
        values = [r"\1\1", r"\1_crap"]
        res = dfmix.copy()
        res.replace(to_replace_res, values, inplace=True, regex=True)
        expec = DataFrame({"a": mix_ab["a"], "b": ["a_crap", "b_crap", "..", ".."]})
        tm.assert_frame_equal(res, expec)

        # list of [re1, re2, ..., reN] -> [(re1 or v1), (re2 or v2), ..., (reN
        # or vN)]
        to_replace_res = [r"\s*(\.)\s*", r"a", r"(b)"]
        values = [r"\1\1", r"crap", r"\1_crap"]
        res = dfmix.copy()
        res.replace(to_replace_res, values, inplace=True, regex=True)
        expec = DataFrame({"a": mix_ab["a"], "b": ["crap", "b_crap", "..", ".."]})
        tm.assert_frame_equal(res, expec)

        to_replace_res = [r"\s*(\.)\s*", r"a", r"(b)"]
        values = [r"\1\1", r"crap", r"\1_crap"]
        res = dfmix.copy()
        res.replace(regex=to_replace_res, value=values, inplace=True)
        expec = DataFrame({"a": mix_ab["a"], "b": ["crap", "b_crap", "..", ".."]})
        tm.assert_frame_equal(res, expec)

    def test_regex_replace_dict_mixed(self, mix_abc):
        dfmix = DataFrame(mix_abc)

        # dicts
        # single dict {re1: v1}, search the whole frame
        # need test for this...

        # list of dicts {re1: v1, re2: v2, ..., re3: v3}, search the whole
        # frame
        res = dfmix.replace({"b": r"\s*\.\s*"}, {"b": np.nan}, regex=True)
        res2 = dfmix.copy()
        res2.replace({"b": r"\s*\.\s*"}, {"b": np.nan}, inplace=True, regex=True)
        expec = DataFrame(
            {"a": mix_abc["a"], "b": ["a", "b", np.nan, np.nan], "c": mix_abc["c"]}
        )
        tm.assert_frame_equal(res, expec)
        tm.assert_frame_equal(res2, expec)

        # list of dicts {re1: re11, re2: re12, ..., reN: re1N}, search the
        # whole frame
        res = dfmix.replace({"b": r"\s*(\.)\s*"}, {"b": r"\1ty"}, regex=True)
        res2 = dfmix.copy()
        res2.replace({"b": r"\s*(\.)\s*"}, {"b": r"\1ty"}, inplace=True, regex=True)
        expec = DataFrame(
            {"a": mix_abc["a"], "b": ["a", "b", ".ty", ".ty"], "c": mix_abc["c"]}
        )
        tm.assert_frame_equal(res, expec)
        tm.assert_frame_equal(res2, expec)

        res = dfmix.replace(regex={"b": r"\s*(\.)\s*"}, value={"b": r"\1ty"})
        res2 = dfmix.copy()
        res2.replace(regex={"b": r"\s*(\.)\s*"}, value={"b": r"\1ty"}, inplace=True)
        expec = DataFrame(
            {"a": mix_abc["a"], "b": ["a", "b", ".ty", ".ty"], "c": mix_abc["c"]}
        )
        tm.assert_frame_equal(res, expec)
        tm.assert_frame_equal(res2, expec)

        # scalar -> dict
        # to_replace regex, {value: value}
        expec = DataFrame(
            {"a": mix_abc["a"], "b": [np.nan, "b", ".", "."], "c": mix_abc["c"]}
        )
        res = dfmix.replace("a", {"b": np.nan}, regex=True)
        res2 = dfmix.copy()
        res2.replace("a", {"b": np.nan}, regex=True, inplace=True)
        tm.assert_frame_equal(res, expec)
        tm.assert_frame_equal(res2, expec)

        res = dfmix.replace("a", {"b": np.nan}, regex=True)
        res2 = dfmix.copy()
        res2.replace(regex="a", value={"b": np.nan}, inplace=True)
        expec = DataFrame(
            {"a": mix_abc["a"], "b": [np.nan, "b", ".", "."], "c": mix_abc["c"]}
        )
        tm.assert_frame_equal(res, expec)
        tm.assert_frame_equal(res2, expec)

    def test_regex_replace_dict_nested(self, mix_abc):
        # nested dicts will not work until this is implemented for Series
        dfmix = DataFrame(mix_abc)
        res = dfmix.replace({"b": {r"\s*\.\s*": np.nan}}, regex=True)
        res2 = dfmix.copy()
        res4 = dfmix.copy()
        res2.replace({"b": {r"\s*\.\s*": np.nan}}, inplace=True, regex=True)
        res3 = dfmix.replace(regex={"b": {r"\s*\.\s*": np.nan}})
        res4.replace(regex={"b": {r"\s*\.\s*": np.nan}}, inplace=True)
        expec = DataFrame(
            {"a": mix_abc["a"], "b": ["a", "b", np.nan, np.nan], "c": mix_abc["c"]}
        )
        tm.assert_frame_equal(res, expec)
        tm.assert_frame_equal(res2, expec)
        tm.assert_frame_equal(res3, expec)
        tm.assert_frame_equal(res4, expec)

    def test_regex_replace_dict_nested_non_first_character(self):
        # GH 25259
        df = pd.DataFrame({"first": ["abc", "bca", "cab"]})
        expected = pd.DataFrame({"first": [".bc", "bc.", "c.b"]})
        result = df.replace({"a": "."}, regex=True)
        tm.assert_frame_equal(result, expected)

    def test_regex_replace_dict_nested_gh4115(self):
        df = pd.DataFrame({"Type": ["Q", "T", "Q", "Q", "T"], "tmp": 2})
        expected = DataFrame({"Type": [0, 1, 0, 0, 1], "tmp": 2})
        result = df.replace({"Type": {"Q": 0, "T": 1}})
        tm.assert_frame_equal(result, expected)

    def test_regex_replace_list_to_scalar(self, mix_abc):
        df = DataFrame(mix_abc)
        expec = DataFrame(
            {
                "a": mix_abc["a"],
                "b": np.array([np.nan] * 4),
                "c": [np.nan, np.nan, np.nan, "d"],
            }
        )
        res = df.replace([r"\s*\.\s*", "a|b"], np.nan, regex=True)
        res2 = df.copy()
        res3 = df.copy()
        res2.replace([r"\s*\.\s*", "a|b"], np.nan, regex=True, inplace=True)
        res3.replace(regex=[r"\s*\.\s*", "a|b"], value=np.nan, inplace=True)
        tm.assert_frame_equal(res, expec)
        tm.assert_frame_equal(res2, expec)
        tm.assert_frame_equal(res3, expec)

    def test_regex_replace_str_to_numeric(self, mix_abc):
        # what happens when you try to replace a numeric value with a regex?
        df = DataFrame(mix_abc)
        res = df.replace(r"\s*\.\s*", 0, regex=True)
        res2 = df.copy()
        res2.replace(r"\s*\.\s*", 0, inplace=True, regex=True)
        res3 = df.copy()
        res3.replace(regex=r"\s*\.\s*", value=0, inplace=True)
        expec = DataFrame({"a": mix_abc["a"], "b": ["a", "b", 0, 0], "c": mix_abc["c"]})
        tm.assert_frame_equal(res, expec)
        tm.assert_frame_equal(res2, expec)
        tm.assert_frame_equal(res3, expec)

    def test_regex_replace_regex_list_to_numeric(self, mix_abc):
        df = DataFrame(mix_abc)
        res = df.replace([r"\s*\.\s*", "b"], 0, regex=True)
        res2 = df.copy()
        res2.replace([r"\s*\.\s*", "b"], 0, regex=True, inplace=True)
        res3 = df.copy()
        res3.replace(regex=[r"\s*\.\s*", "b"], value=0, inplace=True)
        expec = DataFrame(
            {"a": mix_abc["a"], "b": ["a", 0, 0, 0], "c": ["a", 0, np.nan, "d"]}
        )
        tm.assert_frame_equal(res, expec)
        tm.assert_frame_equal(res2, expec)
        tm.assert_frame_equal(res3, expec)

    def test_regex_replace_series_of_regexes(self, mix_abc):
        df = DataFrame(mix_abc)
        s1 = Series({"b": r"\s*\.\s*"})
        s2 = Series({"b": np.nan})
        res = df.replace(s1, s2, regex=True)
        res2 = df.copy()
        res2.replace(s1, s2, inplace=True, regex=True)
        res3 = df.copy()
        res3.replace(regex=s1, value=s2, inplace=True)
        expec = DataFrame(
            {"a": mix_abc["a"], "b": ["a", "b", np.nan, np.nan], "c": mix_abc["c"]}
        )
        tm.assert_frame_equal(res, expec)
        tm.assert_frame_equal(res2, expec)
        tm.assert_frame_equal(res3, expec)

    def test_regex_replace_numeric_to_object_conversion(self, mix_abc):
        df = DataFrame(mix_abc)
        expec = DataFrame({"a": ["a", 1, 2, 3], "b": mix_abc["b"], "c": mix_abc["c"]})
        res = df.replace(0, "a")
        tm.assert_frame_equal(res, expec)
        assert res.a.dtype == np.object_

    @pytest.mark.parametrize("metachar", ["[]", "()", r"\d", r"\w", r"\s"])
    def test_replace_regex_metachar(self, metachar):
        df = DataFrame({"a": [metachar, "else"]})
        result = df.replace({"a": {metachar: "paren"}})
        expected = DataFrame({"a": ["paren", "else"]})
        tm.assert_frame_equal(result, expected)

    def test_replace(self, datetime_frame):
        datetime_frame["A"][:5] = np.nan
        datetime_frame["A"][-5:] = np.nan

        zero_filled = datetime_frame.replace(np.nan, -1e8)
        tm.assert_frame_equal(zero_filled, datetime_frame.fillna(-1e8))
        tm.assert_frame_equal(zero_filled.replace(-1e8, np.nan), datetime_frame)

        datetime_frame["A"][:5] = np.nan
        datetime_frame["A"][-5:] = np.nan
        datetime_frame["B"][:5] = -1e8

        # empty
        df = DataFrame(index=["a", "b"])
        tm.assert_frame_equal(df, df.replace(5, 7))

        # GH 11698
        # test for mixed data types.
        df = pd.DataFrame(
            [("-", pd.to_datetime("20150101")), ("a", pd.to_datetime("20150102"))]
        )
        df1 = df.replace("-", np.nan)
        expected_df = pd.DataFrame(
            [(np.nan, pd.to_datetime("20150101")), ("a", pd.to_datetime("20150102"))]
        )
        tm.assert_frame_equal(df1, expected_df)

    def test_replace_list(self):
        obj = {"a": list("ab.."), "b": list("efgh"), "c": list("helo")}
        dfobj = DataFrame(obj)

        # lists of regexes and values
        # list of [v1, v2, ..., vN] -> [v1, v2, ..., vN]
        to_replace_res = [r".", r"e"]
        values = [np.nan, "crap"]
        res = dfobj.replace(to_replace_res, values)
        expec = DataFrame(
            {
                "a": ["a", "b", np.nan, np.nan],
                "b": ["crap", "f", "g", "h"],
                "c": ["h", "crap", "l", "o"],
            }
        )
        tm.assert_frame_equal(res, expec)

        # list of [v1, v2, ..., vN] -> [v1, v2, .., vN]
        to_replace_res = [r".", r"f"]
        values = [r"..", r"crap"]
        res = dfobj.replace(to_replace_res, values)
        expec = DataFrame(
            {
                "a": ["a", "b", "..", ".."],
                "b": ["e", "crap", "g", "h"],
                "c": ["h", "e", "l", "o"],
            }
        )
        tm.assert_frame_equal(res, expec)

    def test_replace_with_empty_list(self):
        # GH 21977
        s = pd.Series([["a", "b"], [], np.nan, [1]])
        df = pd.DataFrame({"col": s})
        expected = df
        result = df.replace([], np.nan)
        tm.assert_frame_equal(result, expected)

        # GH 19266
        with pytest.raises(ValueError, match="cannot assign mismatch"):
            df.replace({np.nan: []})
        with pytest.raises(ValueError, match="cannot assign mismatch"):
            df.replace({np.nan: ["dummy", "alt"]})

    def test_replace_series_dict(self):
        # from GH 3064
        df = DataFrame({"zero": {"a": 0.0, "b": 1}, "one": {"a": 2.0, "b": 0}})
        result = df.replace(0, {"zero": 0.5, "one": 1.0})
        expected = DataFrame({"zero": {"a": 0.5, "b": 1}, "one": {"a": 2.0, "b": 1.0}})
        tm.assert_frame_equal(result, expected)

        result = df.replace(0, df.mean())
        tm.assert_frame_equal(result, expected)

        # series to series/dict
        df = DataFrame({"zero": {"a": 0.0, "b": 1}, "one": {"a": 2.0, "b": 0}})
        s = Series({"zero": 0.0, "one": 2.0})
        result = df.replace(s, {"zero": 0.5, "one": 1.0})
        expected = DataFrame({"zero": {"a": 0.5, "b": 1}, "one": {"a": 1.0, "b": 0.0}})
        tm.assert_frame_equal(result, expected)

        result = df.replace(s, df.mean())
        tm.assert_frame_equal(result, expected)

    def test_replace_convert(self):
        # gh 3907
        df = DataFrame([["foo", "bar", "bah"], ["bar", "foo", "bah"]])
        m = {"foo": 1, "bar": 2, "bah": 3}
        rep = df.replace(m)
        expec = Series([np.int64] * 3)
        res = rep.dtypes
        tm.assert_series_equal(expec, res)

    def test_replace_mixed(self, float_string_frame):
        mf = float_string_frame
        mf.iloc[5:20, mf.columns.get_loc("foo")] = np.nan
        mf.iloc[-10:, mf.columns.get_loc("A")] = np.nan

        result = float_string_frame.replace(np.nan, -18)
        expected = float_string_frame.fillna(value=-18)
        tm.assert_frame_equal(result, expected)
        tm.assert_frame_equal(result.replace(-18, np.nan), float_string_frame)

        result = float_string_frame.replace(np.nan, -1e8)
        expected = float_string_frame.fillna(value=-1e8)
        tm.assert_frame_equal(result, expected)
        tm.assert_frame_equal(result.replace(-1e8, np.nan), float_string_frame)

        # int block upcasting
        df = DataFrame(
            {
                "A": Series([1.0, 2.0], dtype="float64"),
                "B": Series([0, 1], dtype="int64"),
            }
        )
        expected = DataFrame(
            {
                "A": Series([1.0, 2.0], dtype="float64"),
                "B": Series([0.5, 1], dtype="float64"),
            }
        )
        result = df.replace(0, 0.5)
        tm.assert_frame_equal(result, expected)

        df.replace(0, 0.5, inplace=True)
        tm.assert_frame_equal(df, expected)

        # int block splitting
        df = DataFrame(
            {
                "A": Series([1.0, 2.0], dtype="float64"),
                "B": Series([0, 1], dtype="int64"),
                "C": Series([1, 2], dtype="int64"),
            }
        )
        expected = DataFrame(
            {
                "A": Series([1.0, 2.0], dtype="float64"),
                "B": Series([0.5, 1], dtype="float64"),
                "C": Series([1, 2], dtype="int64"),
            }
        )
        result = df.replace(0, 0.5)
        tm.assert_frame_equal(result, expected)

        # to object block upcasting
        df = DataFrame(
            {
                "A": Series([1.0, 2.0], dtype="float64"),
                "B": Series([0, 1], dtype="int64"),
            }
        )
        expected = DataFrame(
            {
                "A": Series([1, "foo"], dtype="object"),
                "B": Series([0, 1], dtype="int64"),
            }
        )
        result = df.replace(2, "foo")
        tm.assert_frame_equal(result, expected)

        expected = DataFrame(
            {
                "A": Series(["foo", "bar"], dtype="object"),
                "B": Series([0, "foo"], dtype="object"),
            }
        )
        result = df.replace([1, 2], ["foo", "bar"])
        tm.assert_frame_equal(result, expected)

        # test case from
        df = DataFrame(
            {"A": Series([3, 0], dtype="int64"), "B": Series([0, 3], dtype="int64")}
        )
        result = df.replace(3, df.mean().to_dict())
        expected = df.copy().astype("float64")
        m = df.mean()
        expected.iloc[0, 0] = m[0]
        expected.iloc[1, 1] = m[1]
        tm.assert_frame_equal(result, expected)

    def test_replace_simple_nested_dict(self):
        df = DataFrame({"col": range(1, 5)})
        expected = DataFrame({"col": ["a", 2, 3, "b"]})

        result = df.replace({"col": {1: "a", 4: "b"}})
        tm.assert_frame_equal(expected, result)

        # in this case, should be the same as the not nested version
        result = df.replace({1: "a", 4: "b"})
        tm.assert_frame_equal(expected, result)

    def test_replace_simple_nested_dict_with_nonexistent_value(self):
        df = DataFrame({"col": range(1, 5)})
        expected = DataFrame({"col": ["a", 2, 3, "b"]})

        result = df.replace({-1: "-", 1: "a", 4: "b"})
        tm.assert_frame_equal(expected, result)

        result = df.replace({"col": {-1: "-", 1: "a", 4: "b"}})
        tm.assert_frame_equal(expected, result)

    def test_replace_value_is_none(self, datetime_frame):
        orig_value = datetime_frame.iloc[0, 0]
        orig2 = datetime_frame.iloc[1, 0]

        datetime_frame.iloc[0, 0] = np.nan
        datetime_frame.iloc[1, 0] = 1

        result = datetime_frame.replace(to_replace={np.nan: 0})
        expected = datetime_frame.T.replace(to_replace={np.nan: 0}).T
        tm.assert_frame_equal(result, expected)

        result = datetime_frame.replace(to_replace={np.nan: 0, 1: -1e8})
        tsframe = datetime_frame.copy()
        tsframe.iloc[0, 0] = 0
        tsframe.iloc[1, 0] = -1e8
        expected = tsframe
        tm.assert_frame_equal(expected, result)
        datetime_frame.iloc[0, 0] = orig_value
        datetime_frame.iloc[1, 0] = orig2

    def test_replace_for_new_dtypes(self, datetime_frame):

        # dtypes
        tsframe = datetime_frame.copy().astype(np.float32)
        tsframe["A"][:5] = np.nan
        tsframe["A"][-5:] = np.nan

        zero_filled = tsframe.replace(np.nan, -1e8)
        tm.assert_frame_equal(zero_filled, tsframe.fillna(-1e8))
        tm.assert_frame_equal(zero_filled.replace(-1e8, np.nan), tsframe)

        tsframe["A"][:5] = np.nan
        tsframe["A"][-5:] = np.nan
        tsframe["B"][:5] = -1e8

        b = tsframe["B"]
        b[b == -1e8] = np.nan
        tsframe["B"] = b
        result = tsframe.fillna(method="bfill")
        tm.assert_frame_equal(result, tsframe.fillna(method="bfill"))

    @pytest.mark.parametrize(
        "frame, to_replace, value, expected",
        [
            (DataFrame({"ints": [1, 2, 3]}), 1, 0, DataFrame({"ints": [0, 2, 3]})),
            (
                DataFrame({"ints": [1, 2, 3]}, dtype=np.int32),
                1,
                0,
                DataFrame({"ints": [0, 2, 3]}, dtype=np.int32),
            ),
            (
                DataFrame({"ints": [1, 2, 3]}, dtype=np.int16),
                1,
                0,
                DataFrame({"ints": [0, 2, 3]}, dtype=np.int16),
            ),
            (
                DataFrame({"bools": [True, False, True]}),
                False,
                True,
                DataFrame({"bools": [True, True, True]}),
            ),
            (
                DataFrame({"complex": [1j, 2j, 3j]}),
                1j,
                0,
                DataFrame({"complex": [0j, 2j, 3j]}),
            ),
            (
                DataFrame(
                    {
                        "datetime64": Index(
                            [
                                datetime(2018, 5, 28),
                                datetime(2018, 7, 28),
                                datetime(2018, 5, 28),
                            ]
                        )
                    }
                ),
                datetime(2018, 5, 28),
                datetime(2018, 7, 28),
                DataFrame({"datetime64": Index([datetime(2018, 7, 28)] * 3)}),
            ),
            # GH 20380
            (
                DataFrame({"dt": [datetime(3017, 12, 20)], "str": ["foo"]}),
                "foo",
                "bar",
                DataFrame({"dt": [datetime(3017, 12, 20)], "str": ["bar"]}),
            ),
            (
                DataFrame(
                    {
                        "A": date_range("20130101", periods=3, tz="US/Eastern"),
                        "B": [0, np.nan, 2],
                    }
                ),
                Timestamp("20130102", tz="US/Eastern"),
                Timestamp("20130104", tz="US/Eastern"),
                DataFrame(
                    {
                        "A": [
                            Timestamp("20130101", tz="US/Eastern"),
                            Timestamp("20130104", tz="US/Eastern"),
                            Timestamp("20130103", tz="US/Eastern"),
                        ],
                        "B": [0, np.nan, 2],
                    }
                ),
            ),
        ],
    )
    def test_replace_dtypes(self, frame, to_replace, value, expected):
        result = getattr(frame, "replace")(to_replace, value)
        tm.assert_frame_equal(result, expected)

    def test_replace_input_formats_listlike(self):
        # both dicts
        to_rep = {"A": np.nan, "B": 0, "C": ""}
        values = {"A": 0, "B": -1, "C": "missing"}
        df = DataFrame(
            {"A": [np.nan, 0, np.inf], "B": [0, 2, 5], "C": ["", "asdf", "fd"]}
        )
        filled = df.replace(to_rep, values)
        expected = {k: v.replace(to_rep[k], values[k]) for k, v in df.items()}
        tm.assert_frame_equal(filled, DataFrame(expected))

        result = df.replace([0, 2, 5], [5, 2, 0])
        expected = DataFrame(
            {"A": [np.nan, 5, np.inf], "B": [5, 2, 0], "C": ["", "asdf", "fd"]}
        )
        tm.assert_frame_equal(result, expected)

        # scalar to dict
        values = {"A": 0, "B": -1, "C": "missing"}
        df = DataFrame(
            {"A": [np.nan, 0, np.nan], "B": [0, 2, 5], "C": ["", "asdf", "fd"]}
        )
        filled = df.replace(np.nan, values)
        expected = {k: v.replace(np.nan, values[k]) for k, v in df.items()}
        tm.assert_frame_equal(filled, DataFrame(expected))

        # list to list
        to_rep = [np.nan, 0, ""]
        values = [-2, -1, "missing"]
        result = df.replace(to_rep, values)
        expected = df.copy()
        for i in range(len(to_rep)):
            expected.replace(to_rep[i], values[i], inplace=True)
        tm.assert_frame_equal(result, expected)

        msg = r"Replacement lists must match in length\. Expecting 3 got 2"
        with pytest.raises(ValueError, match=msg):
            df.replace(to_rep, values[1:])

    def test_replace_input_formats_scalar(self):
        df = DataFrame(
            {"A": [np.nan, 0, np.inf], "B": [0, 2, 5], "C": ["", "asdf", "fd"]}
        )

        # dict to scalar
        to_rep = {"A": np.nan, "B": 0, "C": ""}
        filled = df.replace(to_rep, 0)
        expected = {k: v.replace(to_rep[k], 0) for k, v in df.items()}
        tm.assert_frame_equal(filled, DataFrame(expected))

        msg = "value argument must be scalar, dict, or Series"
        with pytest.raises(TypeError, match=msg):
            df.replace(to_rep, [np.nan, 0, ""])

        # list to scalar
        to_rep = [np.nan, 0, ""]
        result = df.replace(to_rep, -1)
        expected = df.copy()
        for i in range(len(to_rep)):
            expected.replace(to_rep[i], -1, inplace=True)
        tm.assert_frame_equal(result, expected)

    def test_replace_limit(self):
        pass

    def test_replace_dict_no_regex(self):
        answer = Series(
            {
                0: "Strongly Agree",
                1: "Agree",
                2: "Neutral",
                3: "Disagree",
                4: "Strongly Disagree",
            }
        )
        weights = {
            "Agree": 4,
            "Disagree": 2,
            "Neutral": 3,
            "Strongly Agree": 5,
            "Strongly Disagree": 1,
        }
        expected = Series({0: 5, 1: 4, 2: 3, 3: 2, 4: 1})
        result = answer.replace(weights)
        tm.assert_series_equal(result, expected)

    def test_replace_series_no_regex(self):
        answer = Series(
            {
                0: "Strongly Agree",
                1: "Agree",
                2: "Neutral",
                3: "Disagree",
                4: "Strongly Disagree",
            }
        )
        weights = Series(
            {
                "Agree": 4,
                "Disagree": 2,
                "Neutral": 3,
                "Strongly Agree": 5,
                "Strongly Disagree": 1,
            }
        )
        expected = Series({0: 5, 1: 4, 2: 3, 3: 2, 4: 1})
        result = answer.replace(weights)
        tm.assert_series_equal(result, expected)

    def test_replace_dict_tuple_list_ordering_remains_the_same(self):
        df = DataFrame(dict(A=[np.nan, 1]))
        res1 = df.replace(to_replace={np.nan: 0, 1: -1e8})
        res2 = df.replace(to_replace=(1, np.nan), value=[-1e8, 0])
        res3 = df.replace(to_replace=[1, np.nan], value=[-1e8, 0])

        expected = DataFrame({"A": [0, -1e8]})
        tm.assert_frame_equal(res1, res2)
        tm.assert_frame_equal(res2, res3)
        tm.assert_frame_equal(res3, expected)

    def test_replace_doesnt_replace_without_regex(self):
        raw = """fol T_opp T_Dir T_Enh
        0    1     0     0    vo
        1    2    vr     0     0
        2    2     0     0     0
        3    3     0    bt     0"""
        df = pd.read_csv(StringIO(raw), sep=r"\s+")
        res = df.replace({r"\D": 1})
        tm.assert_frame_equal(df, res)

    def test_replace_bool_with_string(self):
        df = DataFrame({"a": [True, False], "b": list("ab")})
        result = df.replace(True, "a")
        expected = DataFrame({"a": ["a", False], "b": df.b})
        tm.assert_frame_equal(result, expected)

    def test_replace_pure_bool_with_string_no_op(self):
        df = DataFrame(np.random.rand(2, 2) > 0.5)
        result = df.replace("asdf", "fdsa")
        tm.assert_frame_equal(df, result)

    def test_replace_bool_with_bool(self):
        df = DataFrame(np.random.rand(2, 2) > 0.5)
        result = df.replace(False, True)
        expected = DataFrame(np.ones((2, 2), dtype=bool))
        tm.assert_frame_equal(result, expected)

    def test_replace_with_dict_with_bool_keys(self):
        df = DataFrame({0: [True, False], 1: [False, True]})
        with pytest.raises(TypeError, match="Cannot compare types .+"):
            df.replace({"asdf": "asdb", True: "yes"})

    def test_replace_truthy(self):
        df = DataFrame({"a": [True, True]})
        r = df.replace([np.inf, -np.inf], np.nan)
        e = df
        tm.assert_frame_equal(r, e)

    def test_nested_dict_overlapping_keys_replace_int(self):
        # GH 27660 keep behaviour consistent for simple dictionary and
        # nested dictionary replacement
        df = DataFrame({"a": list(range(1, 5))})

        result = df.replace({"a": dict(zip(range(1, 5), range(2, 6)))})
        expected = df.replace(dict(zip(range(1, 5), range(2, 6))))
        tm.assert_frame_equal(result, expected)

    def test_nested_dict_overlapping_keys_replace_str(self):
        # GH 27660
        a = np.arange(1, 5)
        astr = a.astype(str)
        bstr = np.arange(2, 6).astype(str)
        df = DataFrame({"a": astr})
        result = df.replace(dict(zip(astr, bstr)))
        expected = df.replace({"a": dict(zip(astr, bstr))})
        tm.assert_frame_equal(result, expected)

    def test_replace_swapping_bug(self):
        df = pd.DataFrame({"a": [True, False, True]})
        res = df.replace({"a": {True: "Y", False: "N"}})
        expect = pd.DataFrame({"a": ["Y", "N", "Y"]})
        tm.assert_frame_equal(res, expect)

        df = pd.DataFrame({"a": [0, 1, 0]})
        res = df.replace({"a": {0: "Y", 1: "N"}})
        expect = pd.DataFrame({"a": ["Y", "N", "Y"]})
        tm.assert_frame_equal(res, expect)

    def test_replace_period(self):
        d = {
            "fname": {
                "out_augmented_AUG_2011.json": pd.Period(year=2011, month=8, freq="M"),
                "out_augmented_JAN_2011.json": pd.Period(year=2011, month=1, freq="M"),
                "out_augmented_MAY_2012.json": pd.Period(year=2012, month=5, freq="M"),
                "out_augmented_SUBSIDY_WEEK.json": pd.Period(
                    year=2011, month=4, freq="M"
                ),
                "out_augmented_AUG_2012.json": pd.Period(year=2012, month=8, freq="M"),
                "out_augmented_MAY_2011.json": pd.Period(year=2011, month=5, freq="M"),
                "out_augmented_SEP_2013.json": pd.Period(year=2013, month=9, freq="M"),
            }
        }

        df = pd.DataFrame(
            [
                "out_augmented_AUG_2012.json",
                "out_augmented_SEP_2013.json",
                "out_augmented_SUBSIDY_WEEK.json",
                "out_augmented_MAY_2012.json",
                "out_augmented_MAY_2011.json",
                "out_augmented_AUG_2011.json",
                "out_augmented_JAN_2011.json",
            ],
            columns=["fname"],
        )
        assert set(df.fname.values) == set(d["fname"].keys())
        # We don't support converting object -> specialized EA in
        # replace yet.
        expected = DataFrame(
            {"fname": [d["fname"][k] for k in df.fname.values]}, dtype=object
        )
        result = df.replace(d)
        tm.assert_frame_equal(result, expected)

    def test_replace_datetime(self):
        d = {
            "fname": {
                "out_augmented_AUG_2011.json": pd.Timestamp("2011-08"),
                "out_augmented_JAN_2011.json": pd.Timestamp("2011-01"),
                "out_augmented_MAY_2012.json": pd.Timestamp("2012-05"),
                "out_augmented_SUBSIDY_WEEK.json": pd.Timestamp("2011-04"),
                "out_augmented_AUG_2012.json": pd.Timestamp("2012-08"),
                "out_augmented_MAY_2011.json": pd.Timestamp("2011-05"),
                "out_augmented_SEP_2013.json": pd.Timestamp("2013-09"),
            }
        }

        df = pd.DataFrame(
            [
                "out_augmented_AUG_2012.json",
                "out_augmented_SEP_2013.json",
                "out_augmented_SUBSIDY_WEEK.json",
                "out_augmented_MAY_2012.json",
                "out_augmented_MAY_2011.json",
                "out_augmented_AUG_2011.json",
                "out_augmented_JAN_2011.json",
            ],
            columns=["fname"],
        )
        assert set(df.fname.values) == set(d["fname"].keys())
        expected = DataFrame({"fname": [d["fname"][k] for k in df.fname.values]})
        result = df.replace(d)
        tm.assert_frame_equal(result, expected)

    def test_replace_datetimetz(self):

        # GH 11326
        # behaving poorly when presented with a datetime64[ns, tz]
        df = DataFrame(
            {
                "A": date_range("20130101", periods=3, tz="US/Eastern"),
                "B": [0, np.nan, 2],
            }
        )
        result = df.replace(np.nan, 1)
        expected = DataFrame(
            {
                "A": date_range("20130101", periods=3, tz="US/Eastern"),
                "B": Series([0, 1, 2], dtype="float64"),
            }
        )
        tm.assert_frame_equal(result, expected)

        result = df.fillna(1)
        tm.assert_frame_equal(result, expected)

        result = df.replace(0, np.nan)
        expected = DataFrame(
            {
                "A": date_range("20130101", periods=3, tz="US/Eastern"),
                "B": [np.nan, np.nan, 2],
            }
        )
        tm.assert_frame_equal(result, expected)

        result = df.replace(
            Timestamp("20130102", tz="US/Eastern"),
            Timestamp("20130104", tz="US/Eastern"),
        )
        expected = DataFrame(
            {
                "A": [
                    Timestamp("20130101", tz="US/Eastern"),
                    Timestamp("20130104", tz="US/Eastern"),
                    Timestamp("20130103", tz="US/Eastern"),
                ],
                "B": [0, np.nan, 2],
            }
        )
        tm.assert_frame_equal(result, expected)

        result = df.copy()
        result.iloc[1, 0] = np.nan
        result = result.replace({"A": pd.NaT}, Timestamp("20130104", tz="US/Eastern"))
        tm.assert_frame_equal(result, expected)

        # coerce to object
        result = df.copy()
        result.iloc[1, 0] = np.nan
        result = result.replace({"A": pd.NaT}, Timestamp("20130104", tz="US/Pacific"))
        expected = DataFrame(
            {
                "A": [
                    Timestamp("20130101", tz="US/Eastern"),
                    Timestamp("20130104", tz="US/Pacific"),
                    Timestamp("20130103", tz="US/Eastern"),
                ],
                "B": [0, np.nan, 2],
            }
        )
        tm.assert_frame_equal(result, expected)

        result = df.copy()
        result.iloc[1, 0] = np.nan
        result = result.replace({"A": np.nan}, Timestamp("20130104"))
        expected = DataFrame(
            {
                "A": [
                    Timestamp("20130101", tz="US/Eastern"),
                    Timestamp("20130104"),
                    Timestamp("20130103", tz="US/Eastern"),
                ],
                "B": [0, np.nan, 2],
            }
        )
        tm.assert_frame_equal(result, expected)

    def test_replace_with_empty_dictlike(self, mix_abc):
        # GH 15289
        df = DataFrame(mix_abc)
        tm.assert_frame_equal(df, df.replace({}))
        tm.assert_frame_equal(df, df.replace(Series([], dtype=object)))

        tm.assert_frame_equal(df, df.replace({"b": {}}))
        tm.assert_frame_equal(df, df.replace(Series({"b": {}})))

    @pytest.mark.parametrize(
        "to_replace, method, expected",
        [
            (0, "bfill", {"A": [1, 1, 2], "B": [5, np.nan, 7], "C": ["a", "b", "c"]}),
            (
                np.nan,
                "bfill",
                {"A": [0, 1, 2], "B": [5.0, 7.0, 7.0], "C": ["a", "b", "c"]},
            ),
            ("d", "ffill", {"A": [0, 1, 2], "B": [5, np.nan, 7], "C": ["a", "b", "c"]}),
            (
                [0, 2],
                "bfill",
                {"A": [1, 1, 2], "B": [5, np.nan, 7], "C": ["a", "b", "c"]},
            ),
            (
                [1, 2],
                "pad",
                {"A": [0, 0, 0], "B": [5, np.nan, 7], "C": ["a", "b", "c"]},
            ),
            (
                (1, 2),
                "bfill",
                {"A": [0, 2, 2], "B": [5, np.nan, 7], "C": ["a", "b", "c"]},
            ),
            (
                ["b", "c"],
                "ffill",
                {"A": [0, 1, 2], "B": [5, np.nan, 7], "C": ["a", "a", "a"]},
            ),
        ],
    )
    def test_replace_method(self, to_replace, method, expected):
        # GH 19632
        df = DataFrame({"A": [0, 1, 2], "B": [5, np.nan, 7], "C": ["a", "b", "c"]})

        result = df.replace(to_replace=to_replace, value=None, method=method)
        expected = DataFrame(expected)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "replace_dict, final_data",
        [({"a": 1, "b": 1}, [[3, 3], [2, 2]]), ({"a": 1, "b": 2}, [[3, 1], [2, 3]])],
    )
    def test_categorical_replace_with_dict(self, replace_dict, final_data):
        # GH 26988
        df = DataFrame([[1, 1], [2, 2]], columns=["a", "b"], dtype="category")
        expected = DataFrame(final_data, columns=["a", "b"], dtype="category")
        expected["a"] = expected["a"].cat.set_categories([1, 2, 3])
        expected["b"] = expected["b"].cat.set_categories([1, 2, 3])
        result = df.replace(replace_dict, 3)
        tm.assert_frame_equal(result, expected)
        with pytest.raises(AssertionError):
            # ensure non-inplace call does not affect original
            tm.assert_frame_equal(df, expected)
        df.replace(replace_dict, 3, inplace=True)
        tm.assert_frame_equal(df, expected)

    @pytest.mark.parametrize(
        "df, to_replace, exp",
        [
            (
                {"col1": [1, 2, 3], "col2": [4, 5, 6]},
                {4: 5, 5: 6, 6: 7},
                {"col1": [1, 2, 3], "col2": [5, 6, 7]},
            ),
            (
                {"col1": [1, 2, 3], "col2": ["4", "5", "6"]},
                {"4": "5", "5": "6", "6": "7"},
                {"col1": [1, 2, 3], "col2": ["5", "6", "7"]},
            ),
        ],
    )
    def test_replace_commutative(self, df, to_replace, exp):
        # GH 16051
        # DataFrame.replace() overwrites when values are non-numeric
        # also added to data frame whilst issue was for series

        df = pd.DataFrame(df)

        expected = pd.DataFrame(exp)
        result = df.replace(to_replace)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "replacer",
        [
            pd.Timestamp("20170827"),
            np.int8(1),
            np.int16(1),
            np.float32(1),
            np.float64(1),
        ],
    )
    def test_replace_replacer_dtype(self, replacer):
        # GH26632
        df = pd.DataFrame(["a"])
        result = df.replace({"a": replacer, "b": replacer})
        expected = pd.DataFrame([replacer])
        tm.assert_frame_equal(result, expected)
