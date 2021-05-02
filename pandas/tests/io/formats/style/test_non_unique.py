import pytest

from pandas import (
    DataFrame,
    IndexSlice,
)

pytest.importorskip("jinja2")

from pandas.io.formats.style import Styler


@pytest.fixture
def df():
    return DataFrame(
        [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
        index=["i", "j", "j"],
        columns=["c", "d", "d"],
        dtype=float,
    )


@pytest.fixture
def styler(df):
    return Styler(df, uuid_len=0)


def test_format_non_unique(df):
    # GH 41269

    # test dict
    html = df.style.format({"d": "{:.1f}"}).render()
    for val in ["1.000000<", "4.000000<", "7.000000<"]:
        assert val in html
    for val in ["2.0<", "3.0<", "5.0<", "6.0<", "8.0<", "9.0<"]:
        assert val in html

    # test subset
    html = df.style.format(precision=1, subset=IndexSlice["j", "d"]).render()
    for val in ["1.000000<", "4.000000<", "7.000000<", "2.000000<", "3.000000<"]:
        assert val in html
    for val in ["5.0<", "6.0<", "8.0<", "9.0<"]:
        assert val in html


@pytest.mark.parametrize("func", ["apply", "applymap"])
def test_apply_applymap_non_unique_raises(df, func):
    # GH 41269
    if func == "apply":
        op = lambda s: ["color: red;"] * len(s)
    else:
        op = lambda v: "color: red;"

    with pytest.raises(KeyError, match="`Styler.apply` and `.applymap` are not"):
        # slice is non-unique on columns
        getattr(df.style, func)(op, subset=("i", "d"))._compute()

    with pytest.raises(KeyError, match="`Styler.apply` and `.applymap` are not"):
        # slice is non-unique on rows
        getattr(df.style, func)(op, subset=("j", "c"))._compute()

    # unique subset OK
    getattr(df.style, func)(op, subset=("i", "c"))._compute()


def test_table_styles_dict_non_unique_index(styler):
    styles = styler.set_table_styles(
        {"j": [{"selector": "td", "props": "a: v;"}]}, axis=1
    ).table_styles
    assert styles == [
        {"selector": "td.row1", "props": [("a", "v")]},
        {"selector": "td.row2", "props": [("a", "v")]},
    ]


def test_table_styles_dict_non_unique_columns(styler):
    styles = styler.set_table_styles(
        {"d": [{"selector": "td", "props": "a: v;"}]}, axis=0
    ).table_styles
    assert styles == [
        {"selector": "td.col1", "props": [("a", "v")]},
        {"selector": "td.col2", "props": [("a", "v")]},
    ]


def test_maybe_convert_css_raises(styler):
    with pytest.raises(ValueError, match="Styles supplied as string must follow CSS"):
        styler.applymap(lambda x: "bad-css;")._compute()
