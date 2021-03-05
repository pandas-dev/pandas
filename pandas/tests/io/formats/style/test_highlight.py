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

    def test_highlight_max(self):
        df = DataFrame([[1, 2], [3, 4]], columns=["A", "B"])
        css_seq = [("background-color", "yellow")]
        # max(df) = min(-df)
        for max_ in [True, False]:
            if max_:
                attr = "highlight_max"
            else:
                df = -df
                attr = "highlight_min"
            result = getattr(df.style, attr)()._compute().ctx
            assert result[(1, 1)] == css_seq

            result = getattr(df.style, attr)(color="green")._compute().ctx
            assert result[(1, 1)] == [("background-color", "green")]

            result = getattr(df.style, attr)(subset="A")._compute().ctx
            assert result[(1, 0)] == css_seq

            result = getattr(df.style, attr)(axis=0)._compute().ctx
            expected = {
                (1, 0): css_seq,
                (1, 1): css_seq,
            }
            assert result == expected

            result = getattr(df.style, attr)(axis=1)._compute().ctx
            expected = {
                (0, 1): css_seq,
                (1, 1): css_seq,
            }
            assert result == expected

        # separate since we can't negate the strs
        df["C"] = ["a", "b"]
        result = df.style.highlight_max()._compute().ctx
        expected = {(1, 1): css_seq}

        result = df.style.highlight_min()._compute().ctx
        expected = {(0, 0): css_seq}

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
            ("highlight_between", {"subset": [0]}),
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
        elif f == "highlight_between":
            kwargs["left"] = df.iloc[1, 0]  # set the range low for testing

        expected = {(1, 0): [("background-color", "yellow")]}
        result = getattr(df.style, f)(**kwargs)._compute().ctx
        assert result == expected
