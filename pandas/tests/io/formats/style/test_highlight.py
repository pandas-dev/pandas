import numpy as np
import pytest

from pandas import (
    DataFrame,
    IndexSlice,
)

pytest.importorskip("jinja2")


class TestStylerHighlight:
    def setup_method(self, method):
        np.random.seed(24)
        self.s = DataFrame({"A": np.random.permutation(range(6))})
        self.df = DataFrame({"A": [0, 1], "B": np.random.randn(2)})

    def test_highlight_null(self):
        df = DataFrame({"A": [0, np.nan]})
        result = df.style.highlight_null()._compute().ctx
        expected = {(1, 0): [("background-color", "red")]}
        assert result == expected

    def test_highlight_null_subset(self):
        # GH 31345
        df = DataFrame({"A": [0, np.nan], "B": [0, np.nan]})
        result = (
            df.style.highlight_null(null_color="red", subset=["A"])
            .highlight_null(null_color="green", subset=["B"])
            ._compute()
            .ctx
        )
        expected = {
            (1, 0): [("background-color", "red")],
            (1, 1): [("background-color", "green")],
        }
        assert result == expected

    @pytest.mark.parametrize("f", ["highlight_min", "highlight_max"])
    def test_highlight_minmax_basic(self, f):
        expected = {
            (0, 0): [("background-color", "red")],
            (1, 0): [("background-color", "red")],
        }
        if f == "highlight_min":
            df = -self.df
        else:
            df = self.df
        result = getattr(df.style, f)(axis=1, color="red")._compute().ctx
        assert result == expected

    @pytest.mark.parametrize("f", ["highlight_min", "highlight_max"])
    @pytest.mark.parametrize(
        "kwargs",
        [
            {"axis": None, "color": "red"},  # test axis
            {"axis": 0, "subset": ["A"], "color": "red"},  # test subset
            {"axis": None, "props": "background-color: red"},  # test props
        ],
    )
    def test_highlight_minmax_ext(self, f, kwargs):
        expected = {(1, 0): [("background-color", "red")]}
        if f == "highlight_min":
            df = -self.df
        else:
            df = self.df
        result = getattr(df.style, f)(**kwargs)._compute().ctx
        assert result == expected

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"start": 0, "stop": 1},  # test basic range
            {"start": 0, "stop": 1, "props": "background-color: yellow"},  # test props
            {"start": -9, "stop": 9, "subset": ["A"]},  # test subset effective
            {"start": 0},  # test no stop
            {"stop": 1, "subset": ["A"]},  # test no start
            {"start": [0, 1], "axis": 0},  # test start as sequence
            {"start": DataFrame([[0, 1], [1, 1]]), "axis": None},  # test axis with seq
            {"start": 0, "stop": [0, 1], "axis": 0},  # test sequence stop
        ],
    )
    def test_highlight_range(self, kwargs):
        expected = {
            (0, 0): [("background-color", "yellow")],
            (1, 0): [("background-color", "yellow")],
        }
        result = self.df.style.highlight_range(**kwargs)._compute().ctx
        assert result == expected

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"q_low": 0.5, "q_high": 1, "axis": 1},  # test basic range
            {"q_low": 0.5, "q_high": 1, "axis": None},  # test axis
            {"q_low": 0, "q_high": 1, "subset": ["A"]},  # test subset
            {"q_low": 0.5, "axis": 1},  # test no high
            {"q_high": 1, "subset": ["A"], "axis": 0},  # test no low
            {"q_low": 0.5, "axis": 1, "props": "background-color: yellow"},  # tst props
        ],
    )
    def test_highlight_quantile(self, kwargs):
        expected = {
            (0, 0): [("background-color", "yellow")],
            (1, 0): [("background-color", "yellow")],
        }
        result = self.df.style.highlight_quantile(**kwargs)._compute().ctx
        assert result == expected

    @pytest.mark.skipif(
        np.__version__[:4] in ["1.16", "1.17"], reason="Numpy Issue #14831"
    )
    @pytest.mark.parametrize(
        "f,kwargs",
        [
            ("highlight_min", {"axis": 1, "subset": IndexSlice[1, :]}),
            ("highlight_max", {"axis": 0, "subset": [0]}),
            ("highlight_quantile", {"axis": None, "q_low": 0.6, "q_high": 0.8}),
            ("highlight_range", {"subset": [0]}),
        ],
    )
    @pytest.mark.parametrize(
        "df",
        [
            DataFrame([[0, 1], [2, 3]], dtype=int),
            DataFrame([[0, 1], [2, 3]], dtype=float),
            DataFrame([[0, 1], [2, 3]], dtype="datetime64[ns]"),
            DataFrame([[0, 1], [2, 3]], dtype=str),
            DataFrame([[0, 1], [2, 3]], dtype="timedelta64[ns]"),
        ],
    )
    def test_all_highlight_dtypes(self, f, kwargs, df):
        if f == "highlight_quantile" and isinstance(df.iloc[0, 0], str):
            return None  # quantile incompatible with str
        elif f == "highlight_range":
            kwargs["start"] = df.iloc[1, 0]  # set the range low for testing

        expected = {(1, 0): [("background-color", "yellow")]}
        result = getattr(df.style, f)(**kwargs)._compute().ctx
        assert result == expected
