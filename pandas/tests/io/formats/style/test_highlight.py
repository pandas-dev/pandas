import numpy as np
import pytest

from pandas import DataFrame

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
