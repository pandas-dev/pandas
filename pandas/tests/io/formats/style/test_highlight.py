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

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"left": 0, "right": 1},  # test basic range
            {"left": 0, "right": 1, "props": "background-color: yellow"},  # test props
            {"left": -9, "right": 9, "subset": ["A"]},  # test subset effective
            {"left": 0},  # test no right
            {"right": 1, "subset": ["A"]},  # test no left
            {"left": [0, 1], "axis": 0},  # test left as sequence
            {"left": DataFrame([[0, 1], [1, 1]]), "axis": None},  # test axis with seq
            {"left": 0, "right": [0, 1], "axis": 0},  # test sequence right
        ],
    )
    def test_highlight_between(self, kwargs):
        expected = {
            (0, 0): [("background-color", "yellow")],
            (1, 0): [("background-color", "yellow")],
        }
        result = self.df.style.highlight_between(**kwargs)._compute().ctx
        assert result == expected

    @pytest.mark.parametrize(
        "arg, map, axis",
        [
            ("left", [1, 2, 3], 0),
            ("left", [1, 2], 1),
            ("left", np.array([[1, 2], [1, 2]]), None),
            ("right", [1, 2, 3], 0),
            ("right", [1, 2], 1),
            ("right", np.array([[1, 2], [1, 2]]), None),
        ],
    )
    def test_highlight_between_raises(self, arg, map, axis):
        df = DataFrame([[1, 2, 3], [1, 2, 3]])
        msg = f"supplied '{arg}' is not correct shape"
        with pytest.raises(ValueError, match=msg):
            df.style.highlight_between(**{arg: map, "axis": axis})._compute()

    def test_highlight_between_raises2(self):
        with pytest.raises(ValueError, match="as string must be 'left' or 'right'"):
            self.df.style.highlight_between(inclusive="badstring")._compute()

        with pytest.raises(ValueError, match="'inclusive' must be boolean or string"):
            self.df.style.highlight_between(inclusive=1)._compute()

    def test_highlight_between_inclusive(self):
        kwargs = {"left": 0, "right": 1, "subset": ["A"]}
        result = self.df.style.highlight_between(**kwargs, inclusive=True)._compute()
        assert result.ctx == {
            (0, 0): [("background-color", "yellow")],
            (1, 0): [("background-color", "yellow")],
        }
        result = self.df.style.highlight_between(**kwargs, inclusive=False)._compute()
        assert result.ctx == {}
        result = self.df.style.highlight_between(**kwargs, inclusive="left")._compute()
        assert result.ctx == {(0, 0): [("background-color", "yellow")]}
        result = self.df.style.highlight_between(**kwargs, inclusive="right")._compute()
        assert result.ctx == {(1, 0): [("background-color", "yellow")]}
