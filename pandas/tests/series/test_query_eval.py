import numpy as np
import pytest

import pandas.util._test_decorators as td

from pandas import Series, IndexMulti, Index
import pandas._testing as tm

PARSERS = "python", "pandas"
ENGINES = "python", pytest.param("numexpr", marks=td.skip_if_no_ne)


@pytest.fixture(params=PARSERS, ids=lambda x: x)
def parser(request):
    return request.param


@pytest.fixture(params=ENGINES, ids=lambda x: x)
def engine(request):
    return request.param


def skip_if_no_pandas_parser(parser):
    if parser != "pandas":
        pytest.skip(f"cannot evaluate with parser {repr(parser)}")


class TestSeriesEval:
    # smaller hits python, larger hits numexpr
    @pytest.mark.parametrize("n", [4, 4000])
    @pytest.mark.parametrize(
        "op_str,op,rop",
        [
            ("+", "__add__", "__radd__"),
            ("-", "__sub__", "__rsub__"),
            ("*", "__mul__", "__rmul__"),
            ("/", "__truediv__", "__rtruediv__"),
        ],
    )
    def test_ops(self, op_str, op, rop, n):
        # tst ops and reversed ops in evaluation
        # GH7198
        series = Series(1, index=range(n))
        series.iloc[0] = 2
        m = series.mean()

        base = Series(np.tile(m, n)) #noqa

        expected = eval(f"base {op_str} series")

        # ops as strings
        result = eval(f"m {op_str} series")
        tm.assert_series_equal(result, expected)

        # these are commutative
        if op in ["+", "*"]:
            result = getattr(series, op)(m)
            tm.assert_series_equal(result, expected)

        # these are not
        elif op in ["-", "/"]:
            result = getattr(series, rop)(m)
            tm.assert_series_equal(result, expected)

    def test_series_sub_numexpr_path(self):
        # GH7192: Note we need a large number of rows to ensure this
        # goes through the numexpr path
        series = Series(np.random.randn(25000))
        series.iloc[0:5] = np.nan
        expected = 1 - np.isnan(series.iloc[0:25])
        result = (1 - np.isnan(series)).iloc[0:25]
        tm.assert_series_equal(result, expected)

    def test_query_non_str(self):
        # GH 11485
        series = Series({"A": [1, 2, 3]})

        msg = "expr must be a string to be evaluated"
        with pytest.raises(ValueError, match=msg):
            series.query(lambda x: x.A == 1)

        with pytest.raises(ValueError, match=msg):
            series.query(111)

    def test_query_empty_string(self):
        # GH 13139
        series = Series({"A": [1, 2, 3]})

        msg = "expr cannot be an empty string"
        with pytest.raises(ValueError, match=msg):
            series.query("")

    def test_eval_resolvers_as_list(self):
        # GH 14095
        series = Series(np.random.randn(10))
        dict1 = {"a": 1}
        dict2 = {"b": 2}
        assert series.eval("a + b", resolvers=[dict1, dict2]) == dict1["a"] + dict2["b"]
        assert eval("a + b", resolvers=[dict1, dict2]) == dict1["a"] + dict2["b"]


class TestSeriesEvalWithSeries:
    def setup_method(self, method):
        self.index = Index(data=[2000, 2001, 2002], name="year")
        self.series = Series(np.random.randn(3), index=self.index)

    def teardown_method(self, method):
        del self.index
        del self.series

    @pytest.mark.parametrize("op", ["<", "<=", ">", ">=", "==", "!="])
    def test_bool_expr(self, op, parser, engine):
        res = self.series.eval(f"year {op} 2001", engine=engine, parser=parser)
        data1 = eval(f"2000 {op} 2001", engine=engine, parser=parser)
        data2 = eval(f"2001 {op} 2001", engine=engine, parser=parser)
        data3 = eval(f"2002 {op} 2001", engine=engine, parser=parser)
        expect = Series(data=[data1, data2, data3], index=self.index)
        # names are not checked due to different results based on engine
        # (python vs numexpr)
        tm.assert_series_equal(res, expect, check_names=False)

    def test_and_bitwise_operator(self, parser, engine):
        res = self.series.eval("(year < 2001) & (year != 2000)", engine=engine,
            parser=parser)
        expect = Series(data=[False, False, False], index=self.index)
        # names are not checked due to different results based on engine
        # (python vs numexpr)
        tm.assert_series_equal(res, expect, check_names=False)

    def test_or_bitwise_operator(self, parser, engine):
        res = self.series.eval("(year > 2001) | (year == 2000)", engine=engine,
            parser=parser)
        expect = Series(data=[True, False, True], index=self.index)
        # names are not checked due to different results based on engine (python
        # vs numexpr)
        tm.assert_series_equal(res, expect, check_names=False)


class TestSeriesQueryByIndexMethods:
    def setup_method(self, method):
        self.series = Series(np.random.randn(10), index=list(range(10)))
        self.frame = self.series.to_frame()

    def teardown_method(self, method):
        del self.series
        del self.frame

    # test the boolean operands
    def test_bool_operands(self):
        for op in ["<", "<=", ">", ">=", "==", "!="]:
            test = "index " + op + " 5"
            result = self.series.query(test)
            expected = self.frame.query(test)[0]
            # names are not checked since computation/eval.py adds name to
            # evaluated Series
            tm.assert_series_equal(result, expected, check_names=False)

    # test the operands that can join queries
    def test_and_bitwise_operator(self):
        test = "(index > 2) & (index < 8)"
        result = self.series.query(test)
        expected = self.frame.query(test)[0]
        # names are not checked since computation/eval.py adds name to
        # evaluated Series
        tm.assert_series_equal(result, expected, check_names=False)

    def test_or_bitwise_operator(self):
        test = "(index < 3) | (index > 7)"
        result = self.series.query(test)
        expected = self.frame.query(test)[0]
        # names are not checked since computation/eval.py adds name to
        # evaluated Series
        tm.assert_series_equal(result, expected, check_names=False)


class TestSeriesQueryByMultiIndex:
    def setup_method(self, method):
        self.series = Series(np.random.randn(10),
            index=[["a"] * 5 + ["b"] * 5, list(range(10))])
        self.frame = self.series.to_frame()

    def teardown_method(self, method):
        del self.series
        del self.frame

    # check against first level
    def test_query_first_level(self):
        test = "ilevel_0 == 'b'"
        result = self.series.query(test)
        expected = self.frame.query(test)[0]
        # names are not checked since computation/eval.py adds name to evaluated Series
        tm.assert_series_equal(result, expected, check_names=False)

    # check against not first level
    def test_query_not_first_level(self):
        test = "ilevel_1 > 4"
        result = self.series.query(test)
        expected = self.frame.query(test)[0]
        # names are not checked since computation/eval.py adds name to evaluated Series
        tm.assert_series_equal(result, expected, check_names=False)

    def test_both_levels(self):
        for op in ["&", "|"]:
            test = f"(ilevel_0 == 'b') " + op + " ((ilevel_1 % 2) == 0)"
            result = self.series.query(test)
            expected = self.frame.query(test)[0]
            # names are not checked since computation/eval.py adds name to evaluated
            # Series
            tm.assert_series_equal(result, expected, check_names=False)


def run_test(series, test):
    frame = series.to_frame()
    result = series.query(test)
    expected = frame.query(test)[0]
    tm.assert_series_equal(result, expected, check_names=False)


class TestSeriesQueryByIndex:
    def setup_method(self, method):
        self.series = Series(np.random.randn(10), index=list(range(10)))
        self.frame = self.series.to_frame()

    def teardown_method(self, method):
        del self.series
        del self.frame

    # test the boolean operands
    def test_bool_operands(self):
        for op in ["<", "<=", ">", ">=", "==", "!="]:
            run_test(self.series, "index " + op + " 5")
            run_test(self.series, "5 " + op + " index")

    # test list equality
    def test_list_equality(self):
        for op in ["==", "!="]:
            run_test(self.series, "index " + op + " [5]")
            run_test(self.series, "[5] " + op + " index")

    # test the operands that can join queries
    def test_and_bitwise_operator(self):
        run_test(self.series, "(index > 2) & (index < 8)")

    def test_or_bitwise_operator(self):
        run_test(self.series, "(index < 3) | (index > 7)")

    # test in and not in
    def test_in(self):
        run_test(self.series, "'a' in index")
        run_test(self.series, "['a'] in index")

    def test_not_in(self):
        run_test(self.series, "'a' not in index")
        run_test(self.series, "['a'] not in index")


class TestSeriesQueryBacktickQuoting:
    def test_single_backtick(self):
        series = Series(np.random.randn(10), index=Index(list(range(10)), name="B B"))
        run_test(series, "1 < `B B`")

    def test_double_backtick(self):
        series = Series(np.random.randn(10), index=MultiIndex.from_arrays(
            [list(range(10)), list(range(10))], names=["B B", "C C"]))
        run_test(series, "1 < `B B` and 4 < `C C`")

    def test_already_underscore(self):
        series = Series(np.random.randn(10), index=Index(list(range(10)), name = "B_B"))
        run_test(series, "1 < `B_B`")

    def test_same_name_but_underscore(self):
        series = Series(np.random.randn(10), index=MultiIndex.from_arrays(
            [list(range(10)), list(range(10))], names=["C_C", "C C"]))
        run_test(series, "1 < `C_C` and 4 < `C C`")

    def test_underscore_and_spaces(self):
        series = Series(np.random.randn(10), index=Index(list(range(10)), name="B_B B"))
        run_test(series, "1 < `B_B B`")

    def test_special_character_dot(self):
        series = Series(np.random.randn(10), index=Index(list(range(10)), name="B.B"))
        run_test(series, "1 < `B.B`")

    def test_special_character_hyphen(self):
        series = Series(np.random.randn(10), index=Index(list(range(10)), name="B-B"))
        run_test(series, "1 < `B-B`")

    def test_start_with_digit(self):
        series = Series(np.random.randn(10), index=Index(list(range(10)), name="1e1"))
        run_test(series, "1 < `1e1`")

    def test_keyword(self):
        series = Series(np.random.randn(10), index=Index(list(range(10)), name="def"))
        run_test(series, "1 < `def`")

    def test_empty_string(self):
        series = Series(np.random.randn(10), index=Index(list(range(10)), name=""))
        run_test(series, "1 < ``")

    def test_spaces(self):
        series = Series(np.random.randn(10), index=Index(list(range(10)), name="  "))
        run_test(series, "1 < `  `")

    def test_parenthesis(self):
        series = Series(np.random.randn(10), index=Index(list(range(10)), name="(xyz)"))
        run_test(series, "1 < `(xyz)`")

    def test_many_symbols(self):
        series = Series(np.random.randn(10), index=Index(list(range(10)), name="  &^ :!€$?(} >    <++*''  "))
        run_test(series, "1 < `  &^ :!€$?(} >    <++*''  `")

    def test_failing_character_outside_range(self):
        series = Series(np.random.randn(10), index=Index(list(range(10)), name="☺"))
        with pytest.raises(SyntaxError):
            series.query("`☺` > 4")

    def test_failing_hashtag(self):
        series = Series(np.random.randn(10), index=Index(list(range(10)), name="foo#bar"))
        with pytest.raises(SyntaxError):
            series.query("`foo#bar` > 4")


class TestSeriesQueryWithMultiIndex:
    def setup_method(self, method):
        multiIndex = MultiIndex.from_arrays([["a"] * 5 + ["b"] * 5,
            [str(x) for x in range(10)]], names=["alpha", "num"])
        self.series = Series(np.random.randn(10), index=multiIndex)

    def teardown_method(self, method):
        del self.series

    # test against level 0
    def test_equality(self):
        for op in ["==", "!="]:
            run_test(self.series, "alpha " + op + " 'b'")
            run_test(self.series, "'b' " + op + " alpha")

    def test_list_equality(self):
        for op in ["==", "!="]:
            run_test(self.series, "alpha " + op + " ['b']")
            run_test(self.series, "['b'] " + op + " alpha")

    def test_in_operators(self):
        for op in ["in", "not in"]:
            run_test(self.series, "['b'] " + op + " alpha")
            run_test(self.series, "'b' " + op + " alpha")

    # test against level 1
    def test_equality_level1(self):
        for op in ["==", "!="]:
            run_test(self.series, "num " + op + " '3'")
            run_test(self.series, "'3' " + op + " num")

    def test_list_equality_level1(self):
        for op in ["==", "!="]:
            run_test(self.series, "num " + op + " ['3']")
            run_test(self.series, "['3'] " + op + " num")

    def test_in_operator(self):
        for op in ["in", "not in"]:
            run_test(self.series, "['3'] " + op + " num")
            run_test(self.series, "'3' " + op + " num")


class TestSeriesQueryByUnamedMultiIndex:
    def setup_method(self, method):
        self.series = Series(np.random.randn(10), index=[["a"] * 5 + ["b"] * 5,
            list(range(10))])

    def teardown_method(self, method):
        del self.series

    # check against first level
    def test_query_first_level(self):
        run_test(self.series, "ilevel_0 == 'b'")
        run_test(self.series, "'b' == ilevel_0")

    # check against not first level
    def test_query_not_first_level(self):
        run_test(self.series, "ilevel_1 > 4")
        run_test(self.series, "4 > ilevel_1")

    def test_both_levels(self):
        for op in ["&", "|"]:
            run_test(self.series, "(ilevel_0 == 'b') " + op + " ((ilevel_1 % 2) == 0)")
            run_test(self.series, "((ilevel_1 % 2) == 0) " + op + " (ilevel_0 == 'b')")

    def test_levels_equality(self):
        index = [np.random.randint(5, size=100), np.random.randint(5, size=100)]
        series = Series(np.random.randn(100), index=index)

        # test equality
        run_test(series, "ilevel_0 == ilevel_1")
        run_test(series, "ilevel_1 == ilevel_0")

        # test inequality
        run_test(series, "ilevel_0 != ilevel_1")
        run_test(series, "ilevel_1 != ilevel_0")
