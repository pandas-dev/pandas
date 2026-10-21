import numpy as np
import pytest

import pandas as pd

pytest.importorskip("jinja2")
from pandas.io.formats.style import Styler


@pytest.fixture
def df():
    return pd.DataFrame(
        data=[[0, 1, 2], [3, 4, 5], [6, 7, 8]],
        columns=["A", "B", "C"],
        index=["x", "y", "z"],
    )


@pytest.fixture
def styler(df):
    return Styler(df, uuid_len=0)


@pytest.mark.parametrize(
    "data, columns, index",
    [
        # Test basic reindex and ignoring blank
        ([["Min", "Max"], [np.nan, ""]], ["A", "C"], ["x", "y"]),
        # Test non-referenced columns, reversed col names, short index
        ([["Max", "Min", "Bad-Col"]], ["C", "A", "D"], ["x"]),
    ],
)
def test_tooltip_render(data, columns, index, styler):
    ttips = pd.DataFrame(data=data, columns=columns, index=index)

    # GH 21266
    result = styler.set_tooltips(ttips).to_html()

    # test tooltip table level class
    assert "#T_ .pd-t {\n  visibility: hidden;\n" in result

    # test 'Min' tooltip added
    assert "#T_ #T__row0_col0:hover .pd-t {\n  visibility: visible;\n}" in result
    assert '#T_ #T__row0_col0 .pd-t::after {\n  content: "Min";\n}' in result
    assert 'class="data row0 col0" >0<span class="pd-t"></span></td>' in result

    # test 'Max' tooltip added
    assert "#T_ #T__row0_col2:hover .pd-t {\n  visibility: visible;\n}" in result
    assert '#T_ #T__row0_col2 .pd-t::after {\n  content: "Max";\n}' in result
    assert 'class="data row0 col2" >2<span class="pd-t"></span></td>' in result

    # test Nan, empty string and bad column ignored
    assert "#T_ #T__row1_col0:hover .pd-t {\n  visibility: visible;\n}" not in result
    assert "#T_ #T__row1_col1:hover .pd-t {\n  visibility: visible;\n}" not in result
    assert "#T_ #T__row0_col1:hover .pd-t {\n  visibility: visible;\n}" not in result
    assert "#T_ #T__row1_col2:hover .pd-t {\n  visibility: visible;\n}" not in result
    assert "Bad-Col" not in result


def test_tooltip_ignored(styler):
    # GH 21266
    result = styler.to_html()  # no set_tooltips() creates no <span>
    assert '<style type="text/css">\n</style>' in result
    assert '<span class="pd-t"></span>' not in result
    assert 'title="' not in result


def test_tooltip_css_class(styler):
    # GH 21266
    result = styler.set_tooltips(
        pd.DataFrame([["tooltip"]], index=["x"], columns=["A"]),
        css_class="other-class",
        props=[("color", "green")],
    ).to_html()
    assert "#T_ .other-class {\n  color: green;\n" in result
    assert '#T_ #T__row0_col0 .other-class::after {\n  content: "tooltip";\n' in result

    # GH 39563
    result = styler.set_tooltips(  # set_tooltips overwrites previous
        pd.DataFrame([["tooltip"]], index=["x"], columns=["A"]),
        css_class="another-class",
        props="color:green;color:red;",
    ).to_html()
    assert "#T_ .another-class {\n  color: green;\n  color: red;\n}" in result


@pytest.mark.parametrize(
    "data, columns, index",
    [
        # Test basic reindex and ignoring blank
        ([["Min", "Max"], [np.nan, ""]], ["A", "C"], ["x", "y"]),
        # Test non-referenced columns, reversed col names, short index
        ([["Max", "Min", "Bad-Col"]], ["C", "A", "D"], ["x"]),
    ],
)
def test_tooltip_render_as_title(data, columns, index, styler):
    ttips = pd.DataFrame(data=data, columns=columns, index=index)
    # GH 56605
    result = styler.set_tooltips(ttips, as_title_attribute=True).to_html()

    # test css not added
    assert "#T_ .pd-t {\n  visibility: hidden;\n" not in result

    # test 'Min' tooltip added as title attribute and css does not exist
    assert "#T_ #T__row0_col0:hover .pd-t {\n  visibility: visible;\n}" not in result
    assert '#T_ #T__row0_col0 .pd-t::after {\n  content: "Min";\n}' not in result
    assert 'class="data row0 col0"  title="Min">0</td>' in result

    # test 'Max' tooltip added as title attribute and css does not exist
    assert "#T_ #T__row0_col2:hover .pd-t {\n  visibility: visible;\n}" not in result
    assert '#T_ #T__row0_col2 .pd-t::after {\n  content: "Max";\n}' not in result
    assert 'class="data row0 col2"  title="Max">2</td>' in result

    # test Nan, empty string and bad column ignored
    assert "#T_ #T__row1_col0:hover .pd-t {\n  visibility: visible;\n}" not in result
    assert "#T_ #T__row1_col1:hover .pd-t {\n  visibility: visible;\n}" not in result
    assert "#T_ #T__row0_col1:hover .pd-t {\n  visibility: visible;\n}" not in result
    assert "#T_ #T__row1_col2:hover .pd-t {\n  visibility: visible;\n}" not in result
    assert "Bad-Col" not in result
    assert 'class="data row0 col1" >1</td>' in result
    assert 'class="data row1 col0" >3</td>' in result
    assert 'class="data row1 col1" >4</td>' in result
    assert 'class="data row1 col2" >5</td>' in result
    assert 'class="data row2 col0" >6</td>' in result
    assert 'class="data row2 col1" >7</td>' in result
    assert 'class="data row2 col2" >8</td>' in result


def test_tooltip_render_as_title_with_hidden_index_level():
    df = pd.DataFrame(
        data=[[0, 1, 2], [3, 4, 5], [6, 7, 8]],
        columns=["A", "B", "C"],
        index=pd.MultiIndex.from_arrays(
            [["x", "y", "z"], [1, 2, 3], ["aa", "bb", "cc"]],
            names=["alpha", "num", "char"],
        ),
    )
    ttips = pd.DataFrame(
        # Test basic reindex and ignoring blank, and hide level 2 (num) from index
        data=[["Min", "Max"], [np.nan, ""]],
        columns=["A", "C"],
        index=pd.MultiIndex.from_arrays(
            [["x", "y"], [1, 2], ["aa", "bb"]], names=["alpha", "num", "char"]
        ),
    )
    styler = Styler(df, uuid_len=0)
    styler = styler.hide(axis=0, level=-1, names=True)
    # GH 56605
    result = styler.set_tooltips(ttips, as_title_attribute=True).to_html()

    # test css not added
    assert "#T_ .pd-t {\n  visibility: hidden;\n" not in result

    # test 'Min' tooltip added as title attribute and css does not exist
    assert "#T_ #T__row0_col0:hover .pd-t {\n  visibility: visible;\n}" not in result
    assert '#T_ #T__row0_col0 .pd-t::after {\n  content: "Min";\n}' not in result
    assert 'class="data row0 col0"  title="Min">0</td>' in result

    # test 'Max' tooltip added as title attribute and css does not exist
    assert "#T_ #T__row0_col2:hover .pd-t {\n  visibility: visible;\n}" not in result
    assert '#T_ #T__row0_col2 .pd-t::after {\n  content: "Max";\n}' not in result
    assert 'class="data row0 col2"  title="Max">2</td>' in result

    # test Nan, empty string and bad column ignored
    assert "#T_ #T__row1_col0:hover .pd-t {\n  visibility: visible;\n}" not in result
    assert "#T_ #T__row1_col1:hover .pd-t {\n  visibility: visible;\n}" not in result
    assert "#T_ #T__row0_col1:hover .pd-t {\n  visibility: visible;\n}" not in result
    assert "#T_ #T__row1_col2:hover .pd-t {\n  visibility: visible;\n}" not in result
    assert "Bad-Col" not in result
    assert 'class="data row0 col1" >1</td>' in result
    assert 'class="data row1 col0" >3</td>' in result
    assert 'class="data row1 col1" >4</td>' in result
    assert 'class="data row1 col2" >5</td>' in result
    assert 'class="data row2 col0" >6</td>' in result
    assert 'class="data row2 col1" >7</td>' in result
    assert 'class="data row2 col2" >8</td>' in result


def _titles(html):
    return [line.strip() for line in html.splitlines() if "title=" in line]


@pytest.fixture
def hide_df():
    df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}, index=["a", "b", "c"])
    tt = pd.DataFrame({"A": ["ta", "tb", ""], "B": ["", "", ""]}, index=df.index)
    return df, tt


def test_tooltip_title_hidden_row(hide_df):
    # GH#71281 the title belongs to the cell of the visible row, not to the cell
    # at the same position of the body, and empty tooltips are not rendered
    df, tt = hide_df
    result = (
        Styler(df, uuid_len=0)
        .hide(["a"])
        .set_tooltips(tt, as_title_attribute=True)
        .to_html()
    )
    assert _titles(result) == [
        '<td id="T__row1_col0" class="data row1 col0"  title="tb">2</td>'
    ]
    assert 'title="ta"' not in result


def test_tooltip_title_hidden_last_row(hide_df):
    # GH#71281
    df, tt = hide_df
    result = (
        Styler(df, uuid_len=0)
        .hide(["c"])
        .set_tooltips(tt, as_title_attribute=True)
        .to_html()
    )
    assert _titles(result) == [
        '<td id="T__row0_col0" class="data row0 col0"  title="ta">1</td>',
        '<td id="T__row1_col0" class="data row1 col0"  title="tb">2</td>',
    ]


def test_tooltip_title_hidden_column(hide_df):
    # hidden cells never get a title
    df, tt = hide_df
    result = (
        Styler(df, uuid_len=0)
        .hide(["A"], axis=1)
        .set_tooltips(tt, as_title_attribute=True)
        .to_html()
    )
    assert _titles(result) == []


@pytest.mark.parametrize("kwargs", [{"max_rows": 2}, {"max_columns": 1}])
def test_tooltip_title_trimmed(kwargs):
    # GH#71281
    df = pd.DataFrame(np.arange(15).reshape(5, 3), columns=["A", "B", "C"])
    tt = df.astype(str)
    styler = Styler(df, uuid_len=0).set_tooltips(tt, as_title_attribute=True)
    result = styler.to_html(**kwargs)
    n_rows = kwargs.get("max_rows", 5)
    n_cols = kwargs.get("max_columns", 3)
    expected = [
        f'<td id="T__row{i}_col{j}" class="data row{i} col{j}"  '
        f'title="{tt.iloc[i, j]}">{df.iloc[i, j]}</td>'
        for i in range(n_rows)
        for j in range(n_cols)
    ]
    assert _titles(result) == expected


def test_tooltip_title_same_cells_as_span_mode(hide_df):
    # GH#71281
    df, tt = hide_df
    styler = Styler(df, uuid_len=0).hide(["a"])
    span = styler.set_tooltips(tt).to_html()
    title = styler.set_tooltips(tt, as_title_attribute=True).to_html()
    assert '#T__row1_col0 .pd-t::after {\n  content: "tb";\n}' in span
    assert 'content: "ta"' not in span
    assert 'class="data row1 col0"  title="tb">2</td>' in title
    assert 'title="ta"' not in title
