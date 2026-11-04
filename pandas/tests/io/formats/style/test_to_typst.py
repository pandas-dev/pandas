from textwrap import dedent

import pytest

import pandas as pd

pytest.importorskip("jinja2")
from pandas.io.formats.style import Styler


@pytest.fixture
def df():
    return pd.DataFrame(
        {"A": [0, 1], "B": [-0.61, -1.22], "C": pd.Series(["ab", "cd"], dtype=object)}
    )


@pytest.fixture
def styler(df):
    return Styler(df, uuid_len=0, precision=2)


def test_basic_table(styler):
    result = styler.to_typst()
    expected = dedent(
        """\
    #table(
      columns: 4,
      [], [A], [B], [C],

      [0], [0], [-0.61], [ab],
      [1], [1], [-1.22], [cd],
    )"""
    )
    assert result == expected


def test_concat(styler):
    result = styler.concat(styler.data.agg(["sum"]).style).to_typst()
    expected = dedent(
        """\
    #table(
      columns: 4,
      [], [A], [B], [C],

      [0], [0], [-0.61], [ab],
      [1], [1], [-1.22], [cd],
      [sum], [1], [-1.830000], [abcd],
    )"""
    )
    assert result == expected


def test_concat_recursion(styler):
    df = styler.data
    styler1 = styler
    styler2 = Styler(df.agg(["sum"]), uuid_len=0, precision=3)
    styler3 = Styler(df.agg(["sum"]), uuid_len=0, precision=4)
    result = styler1.concat(styler2.concat(styler3)).to_typst()
    expected = dedent(
        """\
    #table(
      columns: 4,
      [], [A], [B], [C],

      [0], [0], [-0.61], [ab],
      [1], [1], [-1.22], [cd],
      [sum], [1], [-1.830], [abcd],
      [sum], [1], [-1.8300], [abcd],
    )"""
    )
    assert result == expected


def test_concat_chain(styler):
    df = styler.data
    styler1 = styler
    styler2 = Styler(df.agg(["sum"]), uuid_len=0, precision=3)
    styler3 = Styler(df.agg(["sum"]), uuid_len=0, precision=4)
    result = styler1.concat(styler2).concat(styler3).to_typst()
    expected = dedent(
        """\
    #table(
      columns: 4,
      [], [A], [B], [C],

      [0], [0], [-0.61], [ab],
      [1], [1], [-1.22], [cd],
      [sum], [1], [-1.830], [abcd],
      [sum], [1], [-1.8300], [abcd],
    )"""
    )
    assert result == expected


@pytest.mark.parametrize("hide_index", [True, False])
@pytest.mark.parametrize("hide_column", [True, False])
def test_hide(styler, hide_index, hide_column):
    # GH 64663
    if hide_index:
        styler.hide(axis="index")
    if hide_column:
        styler.hide(subset=["B"], axis="columns")
    result = styler.to_typst()
    expected = dedent(
        f"""\
    #table(
      columns: {4 - hide_index - hide_column},
      {"" if hide_index else "[], "}[A], {"" if hide_column else "[B], "}[C],

      {"" if hide_index else "[0], "}[0], {"" if hide_column else "[-0.61], "}[ab],
      {"" if hide_index else "[1], "}[1], {"" if hide_column else "[-1.22], "}[cd],
    )"""
    )
    assert result == expected


@pytest.mark.parametrize("hide_headers", [True, False])
def test_hide_all_columns(styler, hide_headers):
    # GH 64663
    styler.hide(axis="index").hide(styler.columns, axis="columns")
    if hide_headers:
        styler.hide(axis="columns")
    result = styler.to_typst()
    # Typst requires a positive column count even when there are no cells.
    expected = "#table( columns: 1, )"
    assert result.split() == expected.split()


@pytest.mark.parametrize("sparse_index", [True, False])
@pytest.mark.parametrize("level", [None, 0, 1, [0, 1]])
def test_hide_multiindex(level, sparse_index):
    # GH 64663
    df = pd.DataFrame(
        {"A": [1, 2]}, index=pd.MultiIndex.from_tuples([("i", 0), ("i", 1)])
    )
    styler = df.style
    if level is not None:
        styler.hide(axis="index", level=level)
    result = styler.to_typst(sparse_index=sparse_index)
    index_cells = {
        "None": (
            "[], [], ",
            "[i], [0], ",
            f"{'[]' if sparse_index else '[i]'}, [1], ",
            3,
        ),
        "0": ("[], ", "[0], ", "[1], ", 2),
        "1": ("[], ", "[i], ", f"{'[]' if sparse_index else '[i]'}, ", 2),
        "[0, 1]": ("", "", "", 1),
    }
    head, first, second, columns = index_cells[str(level)]
    expected = dedent(
        f"""\
    #table(
      columns: {columns},
      {head}[A],

      {first}[1],
      {second}[2],
    )"""
    )
    assert result == expected


def test_hide_named_multiindex():
    # GH 64663
    df = pd.DataFrame(
        {"A": [1, 2], "B": [3, 4]},
        index=pd.MultiIndex.from_tuples(
            [("i", 0, "x"), ("i", 1, "y")], names=["outer", "middle", "inner"]
        ),
    )
    result = df.style.hide(level=[0, 2]).hide(["A"], axis="columns").to_typst()
    expected = dedent(
        """\
    #table(
      columns: 2,
      [], [B],
      [middle], [],

      [0], [3],
      [1], [4],
    )"""
    )
    assert result == expected


@pytest.mark.parametrize("sparse_columns", [True, False])
@pytest.mark.parametrize("hide_column", [True, False])
def test_hide_multiindex_columns(sparse_columns, hide_column):
    # GH 64663
    columns = pd.MultiIndex.from_tuples([("X", "a"), ("X", "b"), ("Y", "c")])
    styler = pd.DataFrame([[1, 2, 3]], columns=columns).style.hide(axis="index")
    if hide_column:
        styler.hide([("Y", "c")], axis="columns")
    result = styler.to_typst(sparse_columns=sparse_columns)
    expected = dedent(
        f"""\
    #table(
      columns: {2 if hide_column else 3},
      [X], {"[]" if sparse_columns else "[X]"},{" " + "[Y]," if not hide_column else ""}
      [a], [b],{" " + "[c]," if not hide_column else ""}

      [1], [2],{" " + "[3]," if not hide_column else ""}
    )"""
    )
    assert result == expected


@pytest.mark.parametrize("hide_index", [True, False])
def test_hide_column_headers(styler, hide_index):
    # GH 64663
    styler.hide(axis="columns")
    if hide_index:
        styler.hide(axis="index")
    result = styler.to_typst()
    expected = dedent(
        f"""\
    #table(
      columns: {3 if hide_index else 4},

      {"" if hide_index else "[0], "}[0], [-0.61], [ab],
      {"" if hide_index else "[1], "}[1], [-1.22], [cd],
    )"""
    )
    assert result == expected


def test_hide_trimmed(styler):
    # GH 64663
    result = (
        styler.hide(axis="index")
        .hide(["B"], axis="columns")
        .to_typst(max_rows=1, max_columns=1)
    )
    expected = dedent(
        """\
    #table(
      columns: 2,
      [A], [...],

      [0], [...],
      [...], [...],
    )"""
    )
    assert result == expected


def test_hide_concat(styler):
    # GH 64663
    result = (
        styler.concat(styler.data.agg(["sum"]).style)
        .hide(axis="index")
        .hide(["A", "B"], axis="columns")
        .to_typst()
    )
    expected = dedent(
        """\
    #table(
      columns: 1,
      [C],

      [ab],
      [cd],
      [abcd],
    )"""
    )
    assert result == expected
