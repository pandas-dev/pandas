import numpy as np
import pytest

from pandas import (
    DataFrame,
    IndexSlice,
)

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


@pytest.mark.parametrize(
    "kwargs",
    [
        {"left": 0, "right": 1},  # test basic range
        {"left": 0, "right": 1, "props": "background-color: yellow"},  # test props
        {"left": -100, "right": 100, "subset": IndexSlice[[0, 1], :]},  # test subset
        {"left": 0, "subset": IndexSlice[[0, 1], :]},  # test no right
        {"right": 1},  # test no left
        {"left": [0, 0, 11], "axis": 0},  # test left as sequence
        {"left": DataFrame({"A": [0, 0, 11], "B": [1, 1, 11]}), "axis": None},  # axis
        {"left": 0, "right": [0, 1], "axis": 1},  # test sequence right
    ],
)
def test_highlight_between(styler, kwargs):
    expected = {
        (0, 0): [("background-color", "yellow")],
        (0, 1): [("background-color", "yellow")],
    }
    result = styler.highlight_between(**kwargs)._compute().ctx
    assert result == expected


@pytest.mark.parametrize(
    "arg, map, axis",
    [
        ("left", [1, 2], 0),  # 0 axis has 3 elements not 2
        ("left", [1, 2, 3], 1),  # 1 axis has 2 elements not 3
        ("left", np.array([[1, 2], [1, 2]]), None),  # df is (2,3) not (2,2)
        ("right", [1, 2], 0),  # same tests as above for 'right' not 'left'
        ("right", [1, 2, 3], 1),  # ..
        ("right", np.array([[1, 2], [1, 2]]), None),  # ..
    ],
)
def test_highlight_between_raises(arg, styler, map, axis):
    msg = f"supplied '{arg}' is not correct shape"
    with pytest.raises(ValueError, match=msg):
        styler.highlight_between(**{arg: map, "axis": axis})._compute()


def test_highlight_between_raises2(styler):
    msg = "values can be 'both', 'left', 'right', or 'neither'"
    with pytest.raises(ValueError, match=msg):
        styler.highlight_between(inclusive="badstring")._compute()

    with pytest.raises(ValueError, match=msg):
        styler.highlight_between(inclusive=1)._compute()


@pytest.mark.parametrize(
    "inclusive, expected",
    [
        (
            "both",
            {
                (0, 0): [("background-color", "yellow")],
                (0, 1): [("background-color", "yellow")],
            },
        ),
        ("neither", {}),
        ("left", {(0, 0): [("background-color", "yellow")]}),
        ("right", {(0, 1): [("background-color", "yellow")]}),
    ],
)
def test_highlight_between_inclusive(styler, inclusive, expected):
    kwargs = {"left": 0, "right": 1, "subset": IndexSlice[[0, 1], :]}
    result = styler.highlight_between(**kwargs, inclusive=inclusive)._compute()
    assert result.ctx == expected
