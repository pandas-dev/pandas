import numpy as np
import pytest

from pandas import DataFrame

pytest.importorskip("jinja2")

from pandas.io.formats.style import Styler


@pytest.fixture
def df():
    return DataFrame({"A": [0, np.nan, 10], "B": [1, None, 2]})


@pytest.fixture
def styler(df):
    return Styler(df, uuid_len=0)


def test_highlight_null(styler):
    result = styler.highlight_null()._compute().ctx
    expected = {
        (1, 0): [("background-color", "red")],
        (1, 1): [("background-color", "red")],
    }
    assert result == expected


def test_highlight_null_subset(styler):
    # GH 31345
    result = (
        styler.highlight_null(null_color="red", subset=["A"])
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
def test_highlight_minmax_basic(df, f):
    expected = {
        (0, 1): [("background-color", "red")],
        # ignores NaN row,
        (2, 0): [("background-color", "red")],
    }
    if f == "highlight_min":
        df = -df
    result = getattr(df.style, f)(axis=1, color="red")._compute().ctx
    assert result == expected


@pytest.mark.parametrize("f", ["highlight_min", "highlight_max"])
@pytest.mark.parametrize(
    "kwargs",
    [
        {"axis": None, "color": "red"},  # test axis
        {"axis": 0, "subset": ["A"], "color": "red"},  # test subset and ignores NaN
        {"axis": None, "props": "background-color: red"},  # test props
    ],
)
def test_highlight_minmax_ext(df, f, kwargs):
    expected = {(2, 0): [("background-color", "red")]}
    if f == "highlight_min":
        df = -df
    result = getattr(df.style, f)(**kwargs)._compute().ctx
    assert result == expected
