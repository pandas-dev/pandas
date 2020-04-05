import pytest
import numpy as np
from pandas import Series, eval, Index
import pandas._testing as tm
import pandas.util._test_decorators as td

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

        base = Series(np.tile(m, n))

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

    def test_bool_expr(self, parser, engine):
        res = self.series.eval("year == 2001", engine=engine, parser=parser)
        expect = Series(data=[False, True, False], index=self.index)
        # names are not checked due to different results based on engine (python vs numexpr)
        tm.assert_series_equal(res, expect, check_names=False)

    
class TestSeriesQueryByIndex:
    def setup_method(self, method):
        self.series = Series(np.random.randn(10), index=[x for x in range(10)])
        self.frame = self.series.to_frame()
    
    def teardown_method(self, method):
        del self.series
        del self.frame

    # test the boolean operands
    @pytest.mark.parametrize("op", ["<", "<=", ">", ">=", "==", "!="])
    def test_bool_operands(self, op):
        test = f"index {op} 5"
        result = self.series.query(test)
        expected = self.frame.query(test)[0]
        tm.assert_series_equal(result, expected)

    # test the operands that can join queries
    def test_and_bitwise_operator(self, op):
        test = f"(index > 2) & (index < 8)"
        result = self.series.query(test)
        expected = self.frame.query(test)[0]
        tm.assert_series_equal(result, expected)

    def test_or_bitwise_operator(self, op):
        test = f"(index < 3) | (index > 7)"
        result = self.series.query(test)
        expected = self.frame.query(test)[0]
        tm.assert_series_equal(result, expected)

class TestSeriesQueryByMultiIndex:
    def setup_method(self, method):
        self.series = Series(np.random.randn(10),
             index=[["a"]*5 + ["b"]*5, [x for x in range(10)]])
        self.frame = self.series.to_frame()
    
    def teardown_method(self, method):
        del self.series
        del self.frame

    # check against first level
    def test_query_first_level(self):
        test = "ilevel_0 == 'b'"
        result = self.series.query(test)
        expected = self.frame.query(test)[0]
        tm.assert_series_equal(result, expected)

    # check against not first level
    def test_query_not_first_level(self):
        test = "ilevel_1 > 4"
        result = self.series.query(test)
        expected = self.frame.query(test)[0]
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("op", ["&", "|"])
    def test_both_levels(self, op):
        test = f"(ilevel_0 == 'b') {op} ((ilevel_1 % 2) == 0)"
        result = self.series.query(test)
        expected = self.frame.query(test)[0]
        tm.assert_series_equal(result, expected)


