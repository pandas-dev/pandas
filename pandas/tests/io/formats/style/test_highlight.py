import numpy as np
import pytest

from pandas import (
    DataFrame,
    IndexSlice,
    Timedelta,
    Timestamp,
)

pytest.importorskip("jinja2")


class TestStylerHighlight:
    def setup_method(self, method):
        self.s = DataFrame({"A": [4, 5, 1, 0, 3, 2]})
        self.df = DataFrame({"A": [0, 1], "B": [-0.609, -1.228]})

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
            {"left": DataFrame({"A": [0, 1], "B": [1, 1]}), "axis": None},  # test axis
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
            ("left", [1, 2, 3], 0),  # 0 axis has 2 elements not 3
            ("left", [1, 2], 1),  # 1 axis has 3 elements not 2
            ("left", np.array([[1, 2], [1, 2]]), None),  # df is (2,3) not (2,2)
            ("right", [1, 2, 3], 0),  # same tests as above for 'right' not 'left'
            ("right", [1, 2], 1),  # ..
            ("right", np.array([[1, 2], [1, 2]]), None),  # ..
        ],
    )
    def test_highlight_between_raises(self, arg, map, axis):
        df = DataFrame([[1, 2, 3], [1, 2, 3]])
        msg = f"supplied '{arg}' is not correct shape"
        with pytest.raises(ValueError, match=msg):
            df.style.highlight_between(**{arg: map, "axis": axis})._compute()

    def test_highlight_between_raises2(self):
        msg = "values can be 'both', 'left', 'right', or 'neither'"
        with pytest.raises(ValueError, match=msg):
            self.df.style.highlight_between(inclusive="badstring")._compute()

        with pytest.raises(ValueError, match=msg):
            self.df.style.highlight_between(inclusive=1)._compute()

    def test_highlight_between_inclusive(self):
        kwargs = {"left": 0, "right": 1, "subset": ["A"]}
        result = self.df.style.highlight_between(**kwargs, inclusive="both")._compute()
        assert result.ctx == {
            (0, 0): [("background-color", "yellow")],
            (1, 0): [("background-color", "yellow")],
        }
        result = self.df.style.highlight_between(
            **kwargs, inclusive="neither"
        )._compute()
        assert result.ctx == {}
        result = self.df.style.highlight_between(**kwargs, inclusive="left")._compute()
        assert result.ctx == {(0, 0): [("background-color", "yellow")]}
        result = self.df.style.highlight_between(**kwargs, inclusive="right")._compute()
        assert result.ctx == {(1, 0): [("background-color", "yellow")]}

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"q_left": 0.5, "q_right": 1, "axis": 1},  # base case
            {"q_left": 0.5, "q_right": 1, "axis": None},  # test axis
            {"q_left": 0, "q_right": 1, "subset": ["A"]},  # test subset
            {"q_left": 0.5, "axis": 1},  # test no high
            {"q_right": 1, "subset": ["A"], "axis": 0},  # test no low
            {"q_left": 0.5, "axis": 1, "props": "background-color: yellow"},  # tst prop
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
            ("highlight_quantile", {"axis": None, "q_left": 0.6, "q_right": 0.8}),
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
        if f == "highlight_quantile" and isinstance(
            df.iloc[0, 0], (str, Timestamp, Timedelta)
        ):
            return None  # quantile incompatible with str, datetime64, timedelta64
        elif f == "highlight_between":
            kwargs["left"] = df.iloc[1, 0]  # set the range low for testing

        expected = {(1, 0): [("background-color", "yellow")]}
        result = getattr(df.style, f)(**kwargs)._compute().ctx
        assert result == expected
