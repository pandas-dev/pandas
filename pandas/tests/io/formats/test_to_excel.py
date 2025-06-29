"""Tests formatting as writer-agnostic ExcelCells"""

import string

import numpy as np
import pytest

pytest.importorskip("jinja2")

from pandas.errors import CSSWarning

from pandas import (
    DataFrame,
    Index,
    MultiIndex,
    NaT,
    Timestamp,
)
import pandas._testing as tm

from pandas.io.formats.excel import (
    CssExcelCell,
    CSSToExcelConverter,
    ExcelFormatter,
)
from pandas.io.formats.style import Styler


@pytest.mark.parametrize(
    "css,expected",
    [
        # FONT
        # - name
        ("font-family: foo,bar", {"font": {"name": "foo"}}),
        ('font-family: "foo bar",baz', {"font": {"name": "foo bar"}}),
        ("font-family: foo,\nbar", {"font": {"name": "foo"}}),
        ("font-family: foo, bar,    baz", {"font": {"name": "foo"}}),
        ("font-family: bar, foo", {"font": {"name": "bar"}}),
        ("font-family: 'foo bar', baz", {"font": {"name": "foo bar"}}),
        ("font-family: 'foo \\'bar', baz", {"font": {"name": "foo 'bar"}}),
        ('font-family: "foo \\"bar", baz', {"font": {"name": 'foo "bar'}}),
        ('font-family: "foo ,bar", baz', {"font": {"name": "foo ,bar"}}),
        # - family
        ("font-family: serif", {"font": {"name": "serif", "family": 1}}),
        ("font-family: Serif", {"font": {"name": "serif", "family": 1}}),
        ("font-family: roman, serif", {"font": {"name": "roman", "family": 1}}),
        ("font-family: roman, sans-serif", {"font": {"name": "roman", "family": 2}}),
        ("font-family: roman, sans serif", {"font": {"name": "roman"}}),
        ("font-family: roman, sansserif", {"font": {"name": "roman"}}),
        ("font-family: roman, cursive", {"font": {"name": "roman", "family": 4}}),
        ("font-family: roman, fantasy", {"font": {"name": "roman", "family": 5}}),
        # - size
        ("font-size: 1em", {"font": {"size": 12}}),
        ("font-size: xx-small", {"font": {"size": 6}}),
        ("font-size: x-small", {"font": {"size": 7.5}}),
        ("font-size: small", {"font": {"size": 9.6}}),
        ("font-size: medium", {"font": {"size": 12}}),
        ("font-size: large", {"font": {"size": 13.5}}),
        ("font-size: x-large", {"font": {"size": 18}}),
        ("font-size: xx-large", {"font": {"size": 24}}),
        ("font-size: 50%", {"font": {"size": 6}}),
        # - bold
        ("font-weight: 100", {"font": {"bold": False}}),
        ("font-weight: 200", {"font": {"bold": False}}),
        ("font-weight: 300", {"font": {"bold": False}}),
        ("font-weight: 400", {"font": {"bold": False}}),
        ("font-weight: normal", {"font": {"bold": False}}),
        ("font-weight: lighter", {"font": {"bold": False}}),
        ("font-weight: bold", {"font": {"bold": True}}),
        ("font-weight: bolder", {"font": {"bold": True}}),
        ("font-weight: 700", {"font": {"bold": True}}),
        ("font-weight: 800", {"font": {"bold": True}}),
        ("font-weight: 900", {"font": {"bold": True}}),
        # - italic
        ("font-style: italic", {"font": {"italic": True}}),
        ("font-style: oblique", {"font": {"italic": True}}),
        # - underline
        ("text-decoration: underline", {"font": {"underline": "single"}}),
        ("text-decoration: overline", {}),
        ("text-decoration: none", {}),
        # - strike
        ("text-decoration: line-through", {"font": {"strike": True}}),
        (
            "text-decoration: underline line-through",
            {"font": {"strike": True, "underline": "single"}},
        ),
        (
            "text-decoration: underline; text-decoration: line-through",
            {"font": {"strike": True}},
        ),
        # - color
        ("color: red", {"font": {"color": "FF0000"}}),
        ("color: #ff0000", {"font": {"color": "FF0000"}}),
        ("color: #f0a", {"font": {"color": "FF00AA"}}),
        # - shadow
        ("text-shadow: none", {"font": {"shadow": False}}),
        ("text-shadow: 0px -0em 0px #CCC", {"font": {"shadow": False}}),
        ("text-shadow: 0px -0em 0px #999", {"font": {"shadow": False}}),
        ("text-shadow: 0px -0em 0px", {"font": {"shadow": False}}),
        ("text-shadow: 2px -0em 0px #CCC", {"font": {"shadow": True}}),
        ("text-shadow: 0px -2em 0px #CCC", {"font": {"shadow": True}}),
        ("text-shadow: 0px -0em 2px #CCC", {"font": {"shadow": True}}),
        ("text-shadow: 0px -0em 2px", {"font": {"shadow": True}}),
        ("text-shadow: 0px -2em", {"font": {"shadow": True}}),
        # FILL
        # - color, fillType
        (
            "background-color: red",
            {"fill": {"fgColor": "FF0000", "patternType": "solid"}},
        ),
        (
            "background-color: #ff0000",
            {"fill": {"fgColor": "FF0000", "patternType": "solid"}},
        ),
        (
            "background-color: #f0a",
            {"fill": {"fgColor": "FF00AA", "patternType": "solid"}},
        ),
        # BORDER
        # - style
        (
            "border-style: solid",
            {
                "border": {
                    "top": {"style": "medium"},
                    "bottom": {"style": "medium"},
                    "left": {"style": "medium"},
                    "right": {"style": "medium"},
                }
            },
        ),
        (
            "border-style: solid; border-width: thin",
            {
                "border": {
                    "top": {"style": "thin"},
                    "bottom": {"style": "thin"},
                    "left": {"style": "thin"},
                    "right": {"style": "thin"},
                }
            },
        ),
        (
            "border-top-style: solid; border-top-width: thin",
            {"border": {"top": {"style": "thin"}}},
        ),
        (
            "border-top-style: solid; border-top-width: 1pt",
            {"border": {"top": {"style": "thin"}}},
        ),
        ("border-top-style: solid", {"border": {"top": {"style": "medium"}}}),
        (
            "border-top-style: solid; border-top-width: medium",
            {"border": {"top": {"style": "medium"}}},
        ),
        (
            "border-top-style: solid; border-top-width: 2pt",
            {"border": {"top": {"style": "medium"}}},
        ),
        (
            "border-top-style: solid; border-top-width: thick",
            {"border": {"top": {"style": "thick"}}},
        ),
        (
            "border-top-style: solid; border-top-width: 4pt",
            {"border": {"top": {"style": "thick"}}},
        ),
        (
            "border-top-style: solid; border-top-width: none",
            {"border": {"top": {"style": "none"}}},
        ),
        (
            "border-top-style: solid; border-top-width: 0.000001pt",
            {"border": {"top": {"style": "none"}}},
        ),
        (
            "border-top-style: dotted",
            {"border": {"top": {"style": "mediumDashDotDot"}}},
        ),
        (
            "border-top-style: dotted; border-top-width: thin",
            {"border": {"top": {"style": "dotted"}}},
        ),
        ("border-top-style: dashed", {"border": {"top": {"style": "mediumDashed"}}}),
        (
            "border-top-style: dashed; border-top-width: thin",
            {"border": {"top": {"style": "dashed"}}},
        ),
        ("border-top-style: double", {"border": {"top": {"style": "double"}}}),
        # - color
        (
            "border-style: solid; border-color: #0000ff",
            {
                "border": {
                    "top": {"style": "medium", "color": "0000FF"},
                    "right": {"style": "medium", "color": "0000FF"},
                    "bottom": {"style": "medium", "color": "0000FF"},
                    "left": {"style": "medium", "color": "0000FF"},
                }
            },
        ),
        (
            "border-top-style: double; border-top-color: blue",
            {"border": {"top": {"style": "double", "color": "0000FF"}}},
        ),
        (
            "border-top-style: solid; border-top-color: #06c",
            {"border": {"top": {"style": "medium", "color": "0066CC"}}},
        ),
        (
            "border-top-color: blue",
            {"border": {"top": {"color": "0000FF", "style": "none"}}},
        ),
        (
            "border-top-style: slantDashDot; border-top-color: blue",
            {"border": {"top": {"style": "slantDashDot", "color": "0000FF"}}},
        ),
        # ALIGNMENT
        # - horizontal
        ("text-align: center", {"alignment": {"horizontal": "center"}}),
        ("text-align: left", {"alignment": {"horizontal": "left"}}),
        ("text-align: right", {"alignment": {"horizontal": "right"}}),
        ("text-align: justify", {"alignment": {"horizontal": "justify"}}),
        # - vertical
        ("vertical-align: top", {"alignment": {"vertical": "top"}}),
        ("vertical-align: text-top", {"alignment": {"vertical": "top"}}),
        ("vertical-align: middle", {"alignment": {"vertical": "center"}}),
        ("vertical-align: bottom", {"alignment": {"vertical": "bottom"}}),
        ("vertical-align: text-bottom", {"alignment": {"vertical": "bottom"}}),
        # - wrap_text
        ("white-space: nowrap", {"alignment": {"wrap_text": False}}),
        ("white-space: pre", {"alignment": {"wrap_text": False}}),
        ("white-space: pre-line", {"alignment": {"wrap_text": False}}),
        ("white-space: normal", {"alignment": {"wrap_text": True}}),
        # NUMBER FORMAT
        ("number-format: 0%", {"number_format": {"format_code": "0%"}}),
        (
            "number-format: 0§[Red](0)§-§@;",
            {"number_format": {"format_code": "0;[red](0);-;@"}},  # GH 46152
        ),
    ],
)
def test_css_to_excel(css, expected):
    convert = CSSToExcelConverter()
    assert expected == convert(css)


def test_css_to_excel_multiple():
    convert = CSSToExcelConverter()
    actual = convert(
        """
        font-weight: bold;
        text-decoration: underline;
        color: red;
        border-width: thin;
        text-align: center;
        vertical-align: top;
        unused: something;
    """
    )
    assert {
        "font": {"bold": True, "underline": "single", "color": "FF0000"},
        "border": {
            "top": {"style": "thin"},
            "right": {"style": "thin"},
            "bottom": {"style": "thin"},
            "left": {"style": "thin"},
        },
        "alignment": {"horizontal": "center", "vertical": "top"},
    } == actual


@pytest.mark.parametrize(
    "css",
    [
        "border-top-style: unhandled-border-style",
        "border-style: another-unhandled-style",
    ],
)
def test_css_to_excel_unhandled_border_style_warns(css):
    """Test that unhandled border styles raise a CSSWarning."""
    convert = CSSToExcelConverter()
    with tm.assert_produces_warning(CSSWarning, match="Unhandled border style format"):
        convert(css)


@pytest.mark.parametrize(
    "css,inherited,expected",
    [
        ("font-weight: bold", "", {"font": {"bold": True}}),
        ("", "font-weight: bold", {"font": {"bold": True}}),
        (
            "font-weight: bold",
            "font-style: italic",
            {"font": {"bold": True, "italic": True}},
        ),
        ("font-style: normal", "font-style: italic", {"font": {"italic": False}}),
        ("font-style: inherit", "", {}),
        (
            "font-style: normal; font-style: inherit",
            "font-style: italic",
            {"font": {"italic": True}},
        ),
    ],
)
def test_css_to_excel_inherited(css, inherited, expected):
    convert = CSSToExcelConverter(inherited)
    assert expected == convert(css)


@pytest.mark.parametrize(
    "input_color,output_color",
    (
        list(CSSToExcelConverter.NAMED_COLORS.items())
        + [("#" + rgb, rgb) for rgb in CSSToExcelConverter.NAMED_COLORS.values()]
        + [("#F0F", "FF00FF"), ("#ABC", "AABBCC")]
    ),
)
def test_css_to_excel_good_colors(input_color, output_color):
    # see gh-18392
    css = (
        f"border-top-color: {input_color}; "
        f"border-right-color: {input_color}; "
        f"border-bottom-color: {input_color}; "
        f"border-left-color: {input_color}; "
        f"background-color: {input_color}; "
        f"color: {input_color}"
    )

    expected = {}

    expected["fill"] = {"patternType": "solid", "fgColor": output_color}

    expected["font"] = {"color": output_color}

    expected["border"] = {
        k: {"color": output_color, "style": "none"}
        for k in ("top", "right", "bottom", "left")
    }

    with tm.assert_produces_warning(None):
        convert = CSSToExcelConverter()
        assert expected == convert(css)


@pytest.mark.parametrize("input_color", [None, "not-a-color"])
def test_css_to_excel_bad_colors(input_color):
    # see gh-18392
    css = (
        f"border-top-color: {input_color}; "
        f"border-right-color: {input_color}; "
        f"border-bottom-color: {input_color}; "
        f"border-left-color: {input_color}; "
        f"background-color: {input_color}; "
        f"color: {input_color}"
    )

    expected = {}

    if input_color is not None:
        expected["fill"] = {"patternType": "solid"}

    with tm.assert_produces_warning(CSSWarning, match="Unhandled color format"):
        convert = CSSToExcelConverter()
        assert expected == convert(css)


@pytest.mark.parametrize("input_color", ["#", "#1234567"])
def test_css_to_excel_invalid_color_raises(input_color):
    """Test that invalid colors raise a ValueError."""
    css = (
        f"border-top-color: {input_color}; "
        f"border-right-color: {input_color}; "
        f"border-bottom-color: {input_color}; "
        f"border-left-color: {input_color}; "
        f"background-color: {input_color}; "
        f"color: {input_color}"
    )

    convert = CSSToExcelConverter()
    with pytest.raises(ValueError, match=f"Unexpected color {input_color}"):
        convert(css)


def tests_css_named_colors_valid():
    upper_hexs = set(map(str.upper, string.hexdigits))
    for color in CSSToExcelConverter.NAMED_COLORS.values():
        assert len(color) == 6 and all(c in upper_hexs for c in color)


def test_css_named_colors_from_mpl_present():
    mpl_colors = pytest.importorskip("matplotlib.colors")

    pd_colors = CSSToExcelConverter.NAMED_COLORS
    for name, color in mpl_colors.CSS4_COLORS.items():
        assert name in pd_colors and pd_colors[name] == color[1:]


@pytest.mark.parametrize(
    "styles,expected",
    [
        ([("color", "green"), ("color", "red")], "color: red;"),
        ([("font-weight", "bold"), ("font-weight", "normal")], "font-weight: normal;"),
        ([("text-align", "center"), ("TEXT-ALIGN", "right")], "text-align: right;"),
    ],
)
def test_css_excel_cell_precedence(styles, expected):
    """It applies favors latter declarations over former declarations"""
    # See GH 47371
    converter = CSSToExcelConverter()
    converter._call_cached.cache_clear()
    css_styles = {(0, 0): styles}
    cell = CssExcelCell(
        row=0,
        col=0,
        val="",
        style=None,
        css_styles=css_styles,
        css_row=0,
        css_col=0,
        css_converter=converter,
    )
    converter._call_cached.cache_clear()

    assert cell.style == converter(expected)


@pytest.mark.parametrize(
    "styles,cache_hits,cache_misses",
    [
        ([[("color", "green"), ("color", "red"), ("color", "green")]], 0, 1),
        (
            [
                [("font-weight", "bold")],
                [("font-weight", "normal"), ("font-weight", "bold")],
            ],
            1,
            1,
        ),
        ([[("text-align", "center")], [("TEXT-ALIGN", "center")]], 1, 1),
        (
            [
                [("font-weight", "bold"), ("text-align", "center")],
                [("font-weight", "bold"), ("text-align", "left")],
            ],
            0,
            2,
        ),
        (
            [
                [("font-weight", "bold"), ("text-align", "center")],
                [("font-weight", "bold"), ("text-align", "left")],
                [("font-weight", "bold"), ("text-align", "center")],
            ],
            1,
            2,
        ),
    ],
)
def test_css_excel_cell_cache(styles, cache_hits, cache_misses):
    """It caches unique cell styles"""
    # See GH 47371
    converter = CSSToExcelConverter()
    converter._call_cached.cache_clear()

    css_styles = {(0, i): _style for i, _style in enumerate(styles)}
    for css_row, css_col in css_styles:
        CssExcelCell(
            row=0,
            col=0,
            val="",
            style=None,
            css_styles=css_styles,
            css_row=css_row,
            css_col=css_col,
            css_converter=converter,
        )
    cache_info = converter._call_cached.cache_info()
    converter._call_cached.cache_clear()

    assert cache_info.hits == cache_hits
    assert cache_info.misses == cache_misses


class TestExcelFormatter:
    def test_excel_sheet_size_error(self, tmp_path):
        """Testing that we throw an error when the sheet input is too large"""
        breaking_row_count = 2**20 + 1
        breaking_col_count = 2**14 + 1

        # Test for too many rows
        row_df = DataFrame(np.zeros(shape=(breaking_row_count, 1)))
        formatter_rows = ExcelFormatter(row_df)

        # Test for too many columns
        col_df = DataFrame(np.zeros(shape=(1, breaking_col_count)))
        formatter_cols = ExcelFormatter(col_df)

        msg = "This sheet is too large!"
        dummy_path = tmp_path / "dummy.xlsx"

        with pytest.raises(ValueError, match=msg):
            formatter_rows.write(writer=dummy_path)

        with pytest.raises(ValueError, match=msg):
            formatter_cols.write(writer=dummy_path)

    @pytest.mark.parametrize(
        "formatter_kwargs",
        [
            {},
            {"style_converter": CSSToExcelConverter()},
        ],
    )
    def test_styler_input_recognized(self, formatter_kwargs):
        # GH 48567
        """
        Testing that when the df input is a Styler, it is recognized and
        written to self.style_converter.
        """
        df = DataFrame({"A": [1, 2]})
        styler = Styler(df)
        formatter = ExcelFormatter(styler, **formatter_kwargs)

        assert formatter.styler is styler
        assert formatter.df is df
        assert isinstance(formatter.style_converter, CSSToExcelConverter)

    @pytest.mark.parametrize("columns", [
        ["A","C"],
        ["A","B","C"]
    ])
    def test_column_missmatch(self, columns):
        """Testing that we throw an error when the specified columns are not found"""
        df = DataFrame({"A": [1, 2], "B": [3, 4]})
        msg="Not all names specified in 'columns' are found"
        with pytest.raises(KeyError, match=msg):
            ExcelFormatter(df, cols=columns)

    @pytest.mark.parametrize("columns",[
    ["C"],
    ["D"],
    ["C","D"]
    ])
    def test_column_full_miss(self,columns):
        """Testing that we throw an error when all of the columns are not in"""
        df = DataFrame({"A": [1, 2], "B": [3, 4]})
        msg="Passed columns are not ALL present dataframe"
        with pytest.raises(KeyError, match=msg):
            ExcelFormatter(df, cols=columns)

    @pytest.mark.parametrize("merge_cells", [
        "1",
        "invalid",
        1,
        0
    ])
    def test_merge_cells_not_bool_throws(self, merge_cells):
        """Testing that we throw an error when merge_cells is not a boolean"""
        df = DataFrame({"A": [1, 2], "B": [3, 4]})
        msg = f"Unexpected value for {merge_cells=}."
        with pytest.raises(ValueError, match=msg):
            ExcelFormatter(df, merge_cells=merge_cells)


    @pytest.mark.parametrize(
        "val, na_rep, expected",
        [
            (np.nan, "missing", "missing"),
            (None, "missing", "missing"),
            (np.nan, "", ""),
        ],
    )
    def test_format_value_handles_na(self, val, na_rep, expected):
        """
        Test that _format_value correctly handles scalar missing values.
        """
        df = DataFrame()
        formatter = ExcelFormatter(df, na_rep=na_rep)
        result = formatter._format_value(val)
        assert result == expected

    @pytest.mark.parametrize(
        "val, float_format, inf_rep, expected",
        [
            (1.12345, "%.2f", "inf", 1.12),
            (np.inf, None, "inf_string", "inf_string"),
            (-np.inf, None, "inf_string", "-inf_string"),
            (1.123, None, "inf", 1.123),
        ],
    )
    def test_format_value_handles_float(self, val, float_format, inf_rep, expected):
        """
        Test that _format_value correctly handles float values.
        """
        df = DataFrame()
        formatter = ExcelFormatter(
            df, float_format=float_format, inf_rep=inf_rep
        )
        result = formatter._format_value(val)
        assert result == expected

    def test_format_value_throws_for_tz_aware_dt(self):
        """
        Test that _format_value raises ValueError for tz-aware datetimes.
        """
        df = DataFrame()
        formatter = ExcelFormatter(df)
        val = Timestamp("2025-06-27 10:00", tz="UTC")
        msg = (
            "Excel does not support datetimes with "
            "timezones. Please ensure that datetimes "
            "are timezone unaware before writing to Excel."
        )
        with pytest.raises(ValueError, match=msg):
            formatter._format_value(val)

    def test_formated_header_mi_multi_index_throws(self):
        header = [
            ("Cool", "A"), ("Cool", "B"),
            ("Amazing", "C"), ("Amazing", "D")
        ]
        columns = MultiIndex.from_tuples(header)
        rng = np.random.default_rng()
        df = DataFrame(
            rng.standard_normal((4, 4)),
            columns=columns
        )
        formatter = ExcelFormatter(df, index=False)
        assert(formatter.columns.nlevels > 1 and not formatter.index)
        with pytest.raises(NotImplementedError):
            list(formatter._format_header_mi())

    def test_returns_none_no_header(self):
        df = DataFrame({"A": [1,2], "B": [3,4]})
        formatter = ExcelFormatter(df,header=False)
        assert(not formatter._has_aliases)
        assert(list(formatter._format_header_mi()) == [])

    @pytest.mark.parametrize(
        "df, merge_cells, expected_cells",
        [
            # Case 1: MultiIndex columns, merge_cells=True
            (
                DataFrame(
                    np.zeros((1, 4)),
                    columns=MultiIndex.from_product([["A", "B"], ["C", "D"]]),
                ),
                True,
                [
                    (0, 0, None, None, None),
                    (1, 0, None, None, None),
                    (0, 1, "A", 0, 2),  # row, col, val, mergestart, mergeend
                    (0, 3, "B", 0, 4),
                    (1, 1, "C", None, None),
                    (1, 2, "D", None, None),
                    (1, 3, "C", None, None),
                    (1, 4, "D", None, None),
                ],
            ),
            # Case 2: MultiIndex columns, merge_cells=False
            (
                DataFrame(
                    np.zeros((1, 4)),
                    columns=MultiIndex.from_product([["A", "B"], ["C", "D"]]),
                ),
                False,
                [
                    (0, 0, None, None, None),
                    (1, 0, None, None, None),
                    (0, 1, "A", None, None),
                    (0, 2, "A", None, None),
                    (0, 3, "B", None, None),
                    (0, 4, "B", None, None),
                    (1, 1, "C", None, None),
                    (1, 2, "D", None, None),
                    (1, 3, "C", None, None),
                    (1, 4, "D", None, None),
                ],
            ),
            # Case 3: MultiIndex on both index and columns
            (
                DataFrame(
                    np.zeros((1, 1)),
                    columns=MultiIndex.from_product([["A"], ["B"]]),
                    index=MultiIndex.from_product([["X"], ["Y"]]),
                ),
                True,
                [
                    (0, 1, None, None, None),
                    (1, 1, None, None, None),
                    (0, 2, "A", None, None),
                    (1, 2, "B", None, None),
                ],
            ),
        ],
    )
    def test_format_header_mi_general(self, df, merge_cells, expected_cells):
        """Test general behavior of _format_header_mi."""
        formatter = ExcelFormatter(df, merge_cells=merge_cells)
        result = list(formatter._format_header_mi())

        # Extract relevant fields for comparison
        result_tuples = [
            (cell.row, cell.col, cell.val, cell.mergestart, cell.mergeend)
            for cell in result
        ]
        assert result_tuples == expected_cells


    @pytest.mark.parametrize(
        "df, formatter_options, expected_cells",
        [
            # Case 1: Default index=True
            (
                DataFrame({"A": [1], "B": [2]}),
                {},
                [
                    (0, 1, "A", None, None),
                    (0, 2, "B", None, None),
                ],
            ),
            # Case 2: index=False
            (
                DataFrame({"A": [1], "B": [2]}),
                {"index": False},
                [
                    (0, 0, "A", None, None),
                    (0, 1, "B", None, None),
                ],
            ),
            # Case 3: With MultiIndex on rows
            (
                DataFrame(
                    {"A": [1], "B": [2]},
                    index=MultiIndex.from_product([["X"], ["Y"]]),
                ),
                {},
                [
                    (0, 2, "A", None, None),
                    (0, 3, "B", None, None),
                ],
            ),
            # Case 4: With custom header aliases
            (
                DataFrame({"A": [1], "B": [2]}),
                {"header": ["C", "D"]},
                [
                    (0, 1, "C", None, None),
                    (0, 2, "D", None, None),
                ],
            ),
        ],
    )
    def test_format_header_regular_general(self, df, formatter_options, expected_cells):
        """Test general behavior of _format_header_regular."""
        formatter = ExcelFormatter(df, **formatter_options)
        result = list(formatter._format_header_regular())

        # Extract relevant fields for comparison
        result_tuples = [
            (cell.row, cell.col, cell.val, cell.mergestart, cell.mergeend)
            for cell in result
        ]
        assert result_tuples == expected_cells

    @pytest.mark.parametrize(
        "df, formatter_options, expected_cells",
        [
            # Case 1: Default index=True, with index name
            (
                DataFrame({"A": [10]}, index=Index(["r1"], name="idx")),
                {},
                [
                    (1, 0, "idx", None, None),  # Index label
                    (2, 0, "r1", None, None),   # Index value
                    (2, 1, 10, None, None),   # Body value
                ],
            ),
            # Case 2: index=False
            (
                DataFrame({"A": [10]}),
                {"index": False},
                [
                    (2, 0, 10, None, None),  # Body value, coloffset=0
                ],
            ),
            # Case 3: With string index_label
            (
                DataFrame({"A": [10]}, index=Index(["r1"], name="idx")),
                {"index_label": "custom_idx"},
                [
                    (1, 0, "custom_idx", None, None),
                    (2, 0, "r1", None, None),
                    (2, 1, 10, None, None),
                ],
            ),
            # Case 4: With list index_label
            (
                DataFrame({"A": [10]}, index=Index(["r1"], name="idx")),
                {"index_label": ["custom_idx", "ignored"]},
                [
                    (1, 0, "custom_idx", None, None),
                    (2, 0, "r1", None, None),
                    (2, 1, 10, None, None),
                ],
            ),
            # Case 5: With MultiIndex columns
            (
                DataFrame(
                    [[10]],
                    index=Index(["r1"], name="idx"),
                    columns=MultiIndex.from_product([["A"], ["B"]]),
                ),
                {},
                [
                    (3, 0, "idx", None, None),  # Index label, row is pushed down
                    (4, 0, "r1", None, None),
                    (4, 1, 10, None, None),
                ],
            ),
        ],
    )
    def test_format_regular_rows_general(self, df, formatter_options, expected_cells):
        """Test general behavior of _format_regular_rows."""
        formatter = ExcelFormatter(df, **formatter_options)
        # Simulate header writing to set rowcounter correctly
        if formatter.header:
            formatter.rowcounter = len(formatter.df.columns.names)

        result = list(formatter._format_regular_rows())

        result_tuples = [
            (cell.row, cell.col, cell.val, cell.mergestart, cell.mergeend)
            for cell in result
        ]
        assert result_tuples == expected_cells

    @pytest.mark.parametrize(
        "df, formatter_options, expected_cells",
        [
            # Case 1: merge_cells=True (default)
            (
                DataFrame(
                    {"A": [10, 20]},
                    index=MultiIndex.from_tuples(
                        [("L1", "a"), ("L1", "b")], names=["idx1", "idx2"]
                    ),
                ),
                {"merge_cells": True},
                [
                    (0, 0, "idx1", None, None),
                    (0, 1, "idx2", None, None),
                    (1, 0, "L1", 2, 0),  # Merged cell
                    (1, 1, "a", None, None),
                    (2, 1, "b", None, None),
                    (1, 2, 10, None, None),  # Body
                    (2, 2, 20, None, None),
                ],
            ),
            # Case 2: merge_cells=False
            (
                DataFrame(
                    {"A": [10, 20]},
                    index=MultiIndex.from_tuples(
                        [("L1", "a"), ("L2", "b")], names=["idx1", "idx2"]
                    ),
                ),
                {"merge_cells": False},
                [
                    (0, 0, "idx1", None, None),
                    (0, 1, "idx2", None, None),
                    (1, 0, "L1", None, None),  # Non-merged
                    (2, 0, "L2", None, None),
                    (1, 1, "a", None, None),
                    (2, 1, "b", None, None),
                    (1, 2, 10, None, None),
                    (2, 2, 20, None, None),
                ],
            ),
            # Case 3: With custom index_label
            (
                DataFrame(
                    {"A": [10]},
                    index=MultiIndex.from_tuples([("L1", "a")], names=["idx1", "idx2"]),
                ),
                {"index_label": ["Custom1", "Custom2"]},
                [
                    (0, 0, "Custom1", None, None),
                    (0, 1, "Custom2", None, None),
                    (1, 0, "L1", None, None),
                    (1, 1, "a", None, None),
                    (1, 2, 10, None, None),
                ],
            ),
            # Case 4: With MultiIndex on columns
            (
                DataFrame(
                    [[10]],
                    index=MultiIndex.from_tuples([("L1", "a")], names=["idx1", "idx2"]),
                    columns=MultiIndex.from_product([["A"], ["B"]]),
                ),
                {},
                [
                    (2, 0, "idx1", None, None),  # Row pushed down
                    (2, 1, "idx2", None, None),
                    (3, 0, "L1", None, None),
                    (3, 1, "a", None, None),
                    (3, 2, 10, None, None),
                ],
            ),
        ],
    )
    def test_format_hierarchical_rows_general(
        self, df, formatter_options, expected_cells
    ):
        """Test general behavior of _format_hierarchical_rows."""
        formatter = ExcelFormatter(df, **formatter_options)
        # Simulate header writing to set rowcounter correctly
        if formatter.header:
            formatter.rowcounter = df.columns.nlevels - 1

        result = list(formatter._format_hierarchical_rows())

        result_tuples = [
            (cell.row, cell.col, cell.val, cell.mergestart, cell.mergeend)
            for cell in result
        ]
        assert result_tuples == expected_cells

    @pytest.mark.parametrize(
        "df, coloffset, rowcounter, expected_cells",
        [
            # Case 1: No offsets
            (
                DataFrame({"A": [10, 20], "B": [30, 40]}),
                0,
                0,
                [
                    (0, 0, 10), (1, 0, 20),  # Column A
                    (0, 1, 30), (1, 1, 40),  # Column B
                ],
            ),
            # Case 2: With offsets
            (
                DataFrame({"A": [10, 20], "B": [30, 40]}),
                1,  # Simulates index=True
                1,  # Simulates header=True
                [
                    (1, 1, 10), (2, 1, 20),  # Column A
                    (1, 2, 30), (2, 2, 40),  # Column B
                ],
            ),
        ],
    )
    def test_generate_body_general(self, df, coloffset, rowcounter, expected_cells):
        """Test general behavior of _generate_body."""
        formatter = ExcelFormatter(df)
        formatter.rowcounter = rowcounter  # Simulate state after header formatting
        result = list(formatter._generate_body(coloffset))

        # Only checking row, col, and val for this internal method
        result_tuples = [(cell.row, cell.col, cell.val) for cell in result]
        assert result_tuples == expected_cells

    @pytest.mark.parametrize(
        "df, formatter_options, expected_values",
        [
            # Case 1: Regular DF with various formats
            (
                DataFrame({"A": [np.nan, 1.2345]}, index=Index([0, 1], name="idx")),
                {
                    "na_rep": "NULL",
                    "float_format": "%.2f",
                    "inf_rep": "Infinity",
                },
                [
                    "A",        # Column header
                    "idx",      # Index header
                    0,          # Index value
                    1,          # Index value
                    "NULL",     # Formatted NaN
                    1.23,       # Formatted float
                ],
            ),
            # Case 2: No index or header
            (
                DataFrame([np.inf]),
                {"index": False, "header": False, "inf_rep": "INF"},
                ["INF"],  # Only the formatted body value
            ),
            # Case 3: MultiIndex on columns and rows
            (
                DataFrame(
                    [[np.nan]],
                    index=MultiIndex.from_tuples([(1, NaT)], names=["i1", "i2"]),
                    columns=MultiIndex.from_tuples([("A", "a")], names=["c1", "c2"]),
                ),
                {"na_rep": "Missing"},
                [
                    "c1", "c2",             # Column headers
                    "A", "a",
                    "i1", "i2",             # Index headers
                    1, "Missing",           # Index values (NaT formatted)
                    "Missing",              # Body value (NaN formatted)
                ],
            ),
        ],
    )
    def test_get_formatted_cells_integration(self, df, formatter_options,
                                             expected_values):
        """
        Integration test for get_formatted_cells to ensure it chains all
        formatting methods and applies final value formatting correctly.
        """
        formatter = ExcelFormatter(df, **formatter_options)
        result = list(formatter.get_formatted_cells())
        result_values = [cell.val for cell in result]

        assert result_values == expected_values
        # Flatten the structure and check against expected formatted values
