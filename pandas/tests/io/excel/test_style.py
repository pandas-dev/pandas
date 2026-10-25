import contextlib
import uuid

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm

from pandas.io.excel import ExcelWriter
from pandas.io.formats.excel import ExcelFormatter

pytest.importorskip("jinja2")
# jinja2 is currently required for Styler.__init__(). Technically Styler.to_excel
# could compute styles and render to excel without jinja2, since there is no
# 'template' file, but this needs the import error to delayed until render time.


@pytest.fixture
def tmp_excel(tmp_path):
    tmp = tmp_path / f"{uuid.uuid4()}.xlsx"
    tmp.touch()
    return str(tmp)


def assert_equal_cell_styles(cell1, cell2):
    # TODO: should find a better way to check equality
    assert cell1.alignment.__dict__ == cell2.alignment.__dict__
    assert cell1.border.__dict__ == cell2.border.__dict__
    assert cell1.fill.__dict__ == cell2.fill.__dict__
    assert cell1.font.__dict__ == cell2.font.__dict__
    assert cell1.number_format == cell2.number_format
    assert cell1.protection.__dict__ == cell2.protection.__dict__


def test_styler_default_values(tmp_excel):
    # GH 54154
    openpyxl = pytest.importorskip("openpyxl")
    df = pd.DataFrame([{"A": 1, "B": 2, "C": 3}, {"A": 1, "B": 2, "C": 3}])

    with ExcelWriter(tmp_excel, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="custom")

    with contextlib.closing(openpyxl.load_workbook(tmp_excel)) as wb:
        # Check font, spacing, indentation
        assert wb["custom"].cell(1, 1).font.bold is False
        assert wb["custom"].cell(1, 1).alignment.horizontal is None
        assert wb["custom"].cell(1, 1).alignment.vertical is None

        # Check border
        assert wb["custom"].cell(1, 1).border.bottom.color is None
        assert wb["custom"].cell(1, 1).border.top.color is None
        assert wb["custom"].cell(1, 1).border.left.color is None
        assert wb["custom"].cell(1, 1).border.right.color is None


@pytest.mark.parametrize("engine", ["xlsxwriter", "openpyxl"])
def test_styler_to_excel_unstyled(engine, tmp_excel):
    # compare DataFrame.to_excel and Styler.to_excel when no styles applied
    pytest.importorskip(engine)
    df = pd.DataFrame(np.random.default_rng(2).standard_normal((2, 2)))
    with ExcelWriter(tmp_excel, engine=engine) as writer:
        df.to_excel(writer, sheet_name="dataframe")
        df.style.to_excel(writer, sheet_name="unstyled")

    openpyxl = pytest.importorskip("openpyxl")  # test loading only with openpyxl
    with contextlib.closing(openpyxl.load_workbook(tmp_excel)) as wb:
        for col1, col2 in zip(
            wb["dataframe"].columns,
            wb["unstyled"].columns,
            strict=True,
        ):
            assert len(col1) == len(col2)
            for cell1, cell2 in zip(col1, col2, strict=True):
                assert cell1.value == cell2.value
                assert_equal_cell_styles(cell1, cell2)


shared_style_params = [
    (
        "background-color: #111222",
        ["fill", "fgColor", "rgb"],
        {"xlsxwriter": "FF111222", "openpyxl": "00111222"},
    ),
    (
        "color: #111222",
        ["font", "color", "value"],
        {"xlsxwriter": "FF111222", "openpyxl": "00111222"},
    ),
    ("font-family: Arial;", ["font", "name"], "arial"),
    ("font-weight: bold;", ["font", "b"], True),
    ("font-style: italic;", ["font", "i"], True),
    ("text-decoration: underline;", ["font", "u"], "single"),
    ("number-format: $??,???.00;", ["number_format"], "$??,???.00"),
    ("text-align: left;", ["alignment", "horizontal"], "left"),
    (
        "vertical-align: bottom;",
        ["alignment", "vertical"],
        {"xlsxwriter": None, "openpyxl": "bottom"},  # xlsxwriter Fails
    ),
    ("vertical-align: middle;", ["alignment", "vertical"], "center"),
    # Border widths
    ("border-left: 2pt solid red", ["border", "left", "style"], "medium"),
    ("border-left: 1pt dotted red", ["border", "left", "style"], "dotted"),
    ("border-left: 2pt dotted red", ["border", "left", "style"], "mediumDashDotDot"),
    ("border-left: 1pt dashed red", ["border", "left", "style"], "dashed"),
    ("border-left: 2pt dashed red", ["border", "left", "style"], "mediumDashed"),
    ("border-left: 1pt solid red", ["border", "left", "style"], "thin"),
    ("border-left: 3pt solid red", ["border", "left", "style"], "thick"),
    # Border expansion
    (
        "border-left: 2pt solid #111222",
        ["border", "left", "color", "rgb"],
        {"xlsxwriter": "FF111222", "openpyxl": "00111222"},
    ),
    ("border: 1pt solid red", ["border", "top", "style"], "thin"),
    (
        "border: 1pt solid #111222",
        ["border", "top", "color", "rgb"],
        {"xlsxwriter": "FF111222", "openpyxl": "00111222"},
    ),
    ("border: 1pt solid red", ["border", "right", "style"], "thin"),
    (
        "border: 1pt solid #111222",
        ["border", "right", "color", "rgb"],
        {"xlsxwriter": "FF111222", "openpyxl": "00111222"},
    ),
    ("border: 1pt solid red", ["border", "bottom", "style"], "thin"),
    (
        "border: 1pt solid #111222",
        ["border", "bottom", "color", "rgb"],
        {"xlsxwriter": "FF111222", "openpyxl": "00111222"},
    ),
    ("border: 1pt solid red", ["border", "left", "style"], "thin"),
    (
        "border: 1pt solid #111222",
        ["border", "left", "color", "rgb"],
        {"xlsxwriter": "FF111222", "openpyxl": "00111222"},
    ),
    # Border styles
    (
        "border-left-style: hair; border-left-color: black",
        ["border", "left", "style"],
        "hair",
    ),
]


def test_styler_custom_style(tmp_excel):
    # GH 54154
    css_style = "background-color: #111222"
    openpyxl = pytest.importorskip("openpyxl")
    df = pd.DataFrame([{"A": 1, "B": 2}, {"A": 1, "B": 2}])

    with ExcelWriter(tmp_excel, engine="openpyxl") as writer:
        styler = df.style.map(lambda x: css_style)
        styler.to_excel(writer, sheet_name="custom", index=False)

    with contextlib.closing(openpyxl.load_workbook(tmp_excel)) as wb:
        # Check font, spacing, indentation
        assert wb["custom"].cell(1, 1).font.bold is False
        assert wb["custom"].cell(1, 1).alignment.horizontal is None
        assert wb["custom"].cell(1, 1).alignment.vertical is None

        # Check border
        assert wb["custom"].cell(1, 1).border.bottom.color is None
        assert wb["custom"].cell(1, 1).border.top.color is None
        assert wb["custom"].cell(1, 1).border.left.color is None
        assert wb["custom"].cell(1, 1).border.right.color is None

        # Check background color
        assert wb["custom"].cell(2, 1).fill.fgColor.index == "00111222"
        assert wb["custom"].cell(3, 1).fill.fgColor.index == "00111222"
        assert wb["custom"].cell(2, 2).fill.fgColor.index == "00111222"
        assert wb["custom"].cell(3, 2).fill.fgColor.index == "00111222"


@pytest.mark.parametrize("engine", ["xlsxwriter", "openpyxl"])
@pytest.mark.parametrize("css, attrs, expected", shared_style_params)
def test_styler_to_excel_basic(engine, css, attrs, expected, tmp_excel):
    pytest.importorskip(engine)
    df = pd.DataFrame(np.random.default_rng(2).standard_normal((1, 1)))
    styler = df.style.map(lambda x: css)

    with ExcelWriter(tmp_excel, engine=engine) as writer:
        df.to_excel(writer, sheet_name="dataframe")
        styler.to_excel(writer, sheet_name="styled")

    openpyxl = pytest.importorskip("openpyxl")  # test loading only with openpyxl
    with contextlib.closing(openpyxl.load_workbook(tmp_excel)) as wb:
        # test unstyled data cell does not have expected styles
        # test styled cell has expected styles
        u_cell, s_cell = wb["dataframe"].cell(2, 2), wb["styled"].cell(2, 2)
    for attr in attrs:
        u_cell, s_cell = getattr(u_cell, attr, None), getattr(s_cell, attr)

    if isinstance(expected, dict):
        assert u_cell is None or u_cell != expected[engine]
        assert s_cell == expected[engine]
    else:
        assert u_cell is None or u_cell != expected
        assert s_cell == expected


@pytest.mark.parametrize("engine", ["xlsxwriter", "openpyxl"])
@pytest.mark.parametrize("css, attrs, expected", shared_style_params)
def test_styler_to_excel_basic_indexes(engine, css, attrs, expected, tmp_excel):
    pytest.importorskip(engine)
    df = pd.DataFrame(np.random.default_rng(2).standard_normal((1, 1)))

    styler = df.style
    styler.map_index(lambda x: css, axis=0)
    styler.map_index(lambda x: css, axis=1)

    null_styler = df.style
    null_styler.map(lambda x: "null: css;")
    null_styler.map_index(lambda x: "null: css;", axis=0)
    null_styler.map_index(lambda x: "null: css;", axis=1)

    with ExcelWriter(tmp_excel, engine=engine) as writer:
        null_styler.to_excel(writer, sheet_name="null_styled")
        styler.to_excel(writer, sheet_name="styled")

    openpyxl = pytest.importorskip("openpyxl")  # test loading only with openpyxl
    with contextlib.closing(openpyxl.load_workbook(tmp_excel)) as wb:
        # test null styled index cells does not have expected styles
        # test styled cell has expected styles
        ui_cell, si_cell = wb["null_styled"].cell(2, 1), wb["styled"].cell(2, 1)
        uc_cell, sc_cell = wb["null_styled"].cell(1, 2), wb["styled"].cell(1, 2)
    for attr in attrs:
        ui_cell, si_cell = getattr(ui_cell, attr, None), getattr(si_cell, attr)
        uc_cell, sc_cell = getattr(uc_cell, attr, None), getattr(sc_cell, attr)

    if isinstance(expected, dict):
        assert ui_cell is None or ui_cell != expected[engine]
        assert si_cell == expected[engine]
        assert uc_cell is None or uc_cell != expected[engine]
        assert sc_cell == expected[engine]
    else:
        assert ui_cell is None or ui_cell != expected
        assert si_cell == expected
        assert uc_cell is None or uc_cell != expected
        assert sc_cell == expected


# From https://openpyxl.readthedocs.io/en/stable/api/openpyxl.styles.borders.html
# Note: Leaving behavior of "width"-type styles undefined; user should use border-width
# instead
excel_border_styles = [
    # "thin",
    "dashed",
    "mediumDashDot",
    "dashDotDot",
    "hair",
    "dotted",
    "mediumDashDotDot",
    # "medium",
    "double",
    "dashDot",
    "slantDashDot",
    # "thick",
    "mediumDashed",
]


@pytest.mark.parametrize("engine", ["xlsxwriter", "openpyxl"])
@pytest.mark.parametrize("border_style", excel_border_styles)
def test_styler_to_excel_border_style(engine, border_style, tmp_excel):
    css = f"border-left: {border_style} black thin"
    attrs = ["border", "left", "style"]
    expected = border_style

    pytest.importorskip(engine)
    df = pd.DataFrame(np.random.default_rng(2).standard_normal((1, 1)))
    styler = df.style.map(lambda x: css)

    with ExcelWriter(tmp_excel, engine=engine) as writer:
        df.to_excel(writer, sheet_name="dataframe")
        styler.to_excel(writer, sheet_name="styled")

    openpyxl = pytest.importorskip("openpyxl")  # test loading only with openpyxl
    with contextlib.closing(openpyxl.load_workbook(tmp_excel)) as wb:
        # test unstyled data cell does not have expected styles
        # test styled cell has expected styles
        u_cell, s_cell = wb["dataframe"].cell(2, 2), wb["styled"].cell(2, 2)
    for attr in attrs:
        u_cell, s_cell = getattr(u_cell, attr, None), getattr(s_cell, attr)

    if isinstance(expected, dict):
        assert u_cell is None or u_cell != expected[engine]
        assert s_cell == expected[engine]
    else:
        assert u_cell is None or u_cell != expected
        assert s_cell == expected


def test_styler_custom_converter(tmp_excel):
    openpyxl = pytest.importorskip("openpyxl")

    def custom_converter(css):
        return {"font": {"color": {"rgb": "111222"}}}

    df = pd.DataFrame(np.random.default_rng(2).standard_normal((1, 1)))
    styler = df.style.map(lambda x: "color: #888999")
    with ExcelWriter(tmp_excel, engine="openpyxl") as writer:
        ExcelFormatter(styler, style_converter=custom_converter).write(
            writer, sheet_name="custom"
        )

    with contextlib.closing(openpyxl.load_workbook(tmp_excel)) as wb:
        assert wb["custom"].cell(2, 2).font.color.value == "00111222"


@pytest.mark.single_cpu
@pytest.mark.filterwarnings(
    "ignore:The default engine for reading:pandas.errors.Pandas4Warning"
)
def test_styler_to_s3(s3_bucket_public, s3so):
    # GH#46381
    mock_bucket_name = s3_bucket_public.name
    target_file = f"{uuid.uuid4()}.xlsx"
    df = pd.DataFrame({"x": [1, 2, 3], "y": [2, 4, 6]})
    styler = df.style.set_sticky(axis="index")
    uri = f"s3://{mock_bucket_name}/{target_file}"
    styler.to_excel(uri, storage_options=s3so)
    result = pd.read_excel(uri, index_col=0, storage_options=s3so)
    tm.assert_frame_equal(result, df)


@pytest.mark.parametrize("merge_cells", [True, False, "columns"])
def test_format_hierarchical_rows_periodindex(merge_cells):
    # GH#60099
    df = pd.DataFrame(
        {"A": [1, 2]},
        index=pd.MultiIndex.from_arrays(
            [
                pd.period_range(start="2006-10-06", end="2006-10-07", freq="D"),
                ["X", "Y"],
            ],
            names=["date", "category"],
        ),
    )
    formatter = ExcelFormatter(df, merge_cells=merge_cells)
    formatted_cells = formatter._format_hierarchical_rows()

    for cell in formatted_cells:
        if cell.row != 0 and cell.col == 0:
            assert isinstance(cell.val, pd.Timestamp), (
                "Period should be converted to Timestamp"
            )


@pytest.mark.parametrize("cols", [["B"], ["C", "A"], ["B", "B"], ["A", "B", "C"]])
def test_styler_to_excel_columns_subset(cols):
    # GH#42276 styles must follow their column when `columns` subsets or reorders
    df = pd.DataFrame({"A": [1, 2], "B": [3, 4], "C": [5, 6]})
    styler = df.style.map(lambda x: f"color: {x}").map_index(
        lambda x: f"color: {x}", axis=1
    )
    formatter = ExcelFormatter(styler, cols=cols, style_converter=dict)
    styles = {
        (cell.row, cell.col): cell.style for cell in formatter.get_formatted_cells()
    }

    for col_num, col in enumerate(cols, start=1):
        assert styles[0, col_num] == {"color": col}
        for row_num, val in enumerate(df[col], start=1):
            assert styles[row_num, col_num] == {"color": str(val)}


def test_styler_to_excel_columns_all_duplicate_labels():
    # selecting every column of a frame with duplicate labels is a no-op reindex
    df = pd.DataFrame([[1, 2, 3]], columns=["A", "A", "B"])
    styler = df.style.map_index(lambda x: f"color: {x}", axis=1)
    formatter = ExcelFormatter(styler, cols=["A", "A", "B"], style_converter=dict)
    cells = formatter.get_formatted_cells()
    styles = {cell.col: cell.style for cell in cells if cell.row == 0}
    assert styles == {1: {"color": "A"}, 2: {"color": "A"}, 3: {"color": "B"}}


@pytest.mark.parametrize("merge_cells", [True, False])
def test_styler_to_excel_columns_subset_multiindex(merge_cells):
    # GH#42276
    columns = pd.MultiIndex.from_tuples([("a", "x"), ("b", "y")])
    df = pd.DataFrame([[1, 2]], columns=columns)
    styler = df.style.map_index(lambda x: f"color: {x}", axis=1)
    formatter = ExcelFormatter(
        styler, cols=[("b", "y")], merge_cells=merge_cells, style_converter=dict
    )
    styles = {
        (cell.row, cell.col): cell.style
        for cell in formatter.get_formatted_cells()
        if cell.row < 2 and cell.col == 1
    }
    assert styles == {(0, 1): {"color": "b"}, (1, 1): {"color": "y"}}


def _styled_cells(styler, **kwargs):
    # style every cell "color: <its value>" so misplaced styles are detectable
    styler = (
        styler.map(lambda x: f"color: {x}")
        .map_index(lambda x: f"color: {x}", axis=0)
        .map_index(lambda x: f"color: {x}", axis=1)
    )
    formatter = ExcelFormatter(styler, style_converter=dict, **kwargs)
    return {
        (cell.row, cell.col): (cell.val, cell.style)
        for cell in formatter.get_formatted_cells()
    }


def _expected_cells(grid):
    return {
        (row, col): (val, {"color": str(val)})
        for row, vals in enumerate(grid)
        for col, val in enumerate(vals)
        if val is not None
    }


@pytest.mark.parametrize(
    "cols, header", [(None, ["A", "C"]), (["C", "B", "A"], ["C", "A"])]
)
def test_styler_to_excel_hidden_rows_and_columns(cols, header):
    # GH#42276
    df = pd.DataFrame(
        {"A": [1, 2, 3], "B": [4, 5, 6], "C": [7, 8, 9]}, index=["x", "y", "z"]
    )
    styler = df.style.hide(["y"]).hide(["B"], axis=1)
    result = _styled_cells(styler, cols=cols)
    body = df.loc[["x", "z"], header]
    grid = [[None, *header]] + [[label, *row] for label, row in body.iterrows()]
    assert result == _expected_cells(grid)


@pytest.mark.parametrize("merge_cells", [True, False])
def test_styler_to_excel_hidden_levels(merge_cells):
    # GH#42276
    index = pd.MultiIndex.from_tuples([("a", "i", "p"), ("b", "j", "q")])
    columns = pd.MultiIndex.from_tuples([("X", "u"), ("Y", "v")])
    df = pd.DataFrame([[1, 2], [3, 4]], index=index, columns=columns)
    styler = df.style.hide(axis=0, level=1).hide(axis=1, level=0)
    result = _styled_cells(styler, merge_cells=merge_cells)
    grid = [[None, None, "u", "v"], ["a", "p", 1, 2], ["b", "q", 3, 4]]
    assert result == _expected_cells(grid)


def test_styler_to_excel_hidden_index_and_header():
    # GH#42276
    df = pd.DataFrame({"A": [1, 2]}, index=["x", "y"])
    result = _styled_cells(df.style.hide(axis=0).hide(axis=1))
    assert result == _expected_cells([[1], [2]])


def test_styler_to_excel_hidden_index_names():
    # GH#42276
    df = pd.DataFrame({"A": [1]}, index=pd.Index(["x"], name="name"))
    assert _styled_cells(df.style)[0, 0] == ("name", None)
    result = _styled_cells(df.style.hide(axis=0, names=True))
    assert result == _expected_cells([[None, "A"], ["x", 1]])


@pytest.mark.parametrize(
    "level, expected",
    [
        (1, [["x0", "x2", "A"], ["a", "p", 1]]),
        ([0, 1], [["x2", "A"], ["p", 1]]),
    ],
)
def test_styler_to_excel_hidden_levels_index_label(level, expected):
    # GH#42276 index_label names every level, including hidden ones
    index = pd.MultiIndex.from_tuples([("a", "i", "p")], names=["L0", "L1", "L2"])
    df = pd.DataFrame({"A": [1]}, index=index)
    styler = df.style.hide(axis=0, level=level)
    formatter = ExcelFormatter(styler, index_label=["x0", "x1", "x2"])
    result = {
        (cell.row, cell.col): cell.val for cell in formatter.get_formatted_cells()
    }
    assert result == {
        (row, col): val
        for row, vals in enumerate(expected)
        for col, val in enumerate(vals)
    }


def test_styler_to_excel_hidden_columns_header_aliases():
    # GH#42276 header aliases name every column, including hidden ones
    df = pd.DataFrame({"A": [1], "B": [2], "C": [3]})
    styler = df.style.hide(["B"], axis=1)
    formatter = ExcelFormatter(styler, header=["a", "b", "c"])
    result = {
        (cell.row, cell.col): cell.val for cell in formatter.get_formatted_cells()
    }
    assert result == {(0, 1): "a", (0, 2): "c", (1, 0): 0, (1, 1): 1, (1, 2): 3}

    msg = "Writing 3 cols but got 2 aliases"
    with pytest.raises(ValueError, match=msg):
        ExcelFormatter(styler, header=["a", "c"])


def test_styler_to_excel_hidden_column_names():
    # GH#42276
    columns = pd.MultiIndex.from_tuples([("X", "u")], names=["n0", "n1"])
    df = pd.DataFrame([[1]], index=["x"], columns=columns)
    assert _styled_cells(df.style)[0, 0] == ("n0", None)
    result = _styled_cells(df.style.hide(axis=1, names=True))
    # row 2 is left blank for the index names
    expected = _expected_cells([[None, "X"], [None, "u"], [], ["x", 1]])
    assert result == {**expected, (0, 0): ("", None), (1, 0): ("", None)}
