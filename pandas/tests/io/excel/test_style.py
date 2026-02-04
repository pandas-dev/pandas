import contextlib
import uuid

import numpy as np
import pytest

import pandas.util._test_decorators as td

from pandas import (
    DataFrame,
    MultiIndex,
    Timestamp,
    period_range,
    read_excel,
)
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
    df = DataFrame([{"A": 1, "B": 2, "C": 3}, {"A": 1, "B": 2, "C": 3}])

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
    df = DataFrame(np.random.default_rng(2).standard_normal((2, 2)))
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
    df = DataFrame([{"A": 1, "B": 2}, {"A": 1, "B": 2}])

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
    df = DataFrame(np.random.default_rng(2).standard_normal((1, 1)))
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
    df = DataFrame(np.random.default_rng(2).standard_normal((1, 1)))

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
    df = DataFrame(np.random.default_rng(2).standard_normal((1, 1)))
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

    df = DataFrame(np.random.default_rng(2).standard_normal((1, 1)))
    styler = df.style.map(lambda x: "color: #888999")
    with ExcelWriter(tmp_excel, engine="openpyxl") as writer:
        ExcelFormatter(styler, style_converter=custom_converter).write(
            writer, sheet_name="custom"
        )

    with contextlib.closing(openpyxl.load_workbook(tmp_excel)) as wb:
        assert wb["custom"].cell(2, 2).font.color.value == "00111222"


@pytest.mark.single_cpu
@td.skip_if_not_us_locale
def test_styler_to_s3(s3_bucket_public, s3so):
    # GH#46381
    mock_bucket_name = s3_bucket_public.name
    target_file = f"{uuid.uuid4()}.xlsx"
    df = DataFrame({"x": [1, 2, 3], "y": [2, 4, 6]})
    styler = df.style.set_sticky(axis="index")
    uri = f"s3://{mock_bucket_name}/{target_file}"
    styler.to_excel(uri, storage_options=s3so)
    result = read_excel(uri, index_col=0, storage_options=s3so)
    tm.assert_frame_equal(result, df)


@pytest.mark.parametrize("merge_cells", [True, False, "columns"])
def test_format_hierarchical_rows_periodindex(merge_cells):
    # GH#60099
    df = DataFrame(
        {"A": [1, 2]},
        index=MultiIndex.from_arrays(
            [
                period_range(start="2006-10-06", end="2006-10-07", freq="D"),
                ["X", "Y"],
            ],
            names=["date", "category"],
        ),
    )
    formatter = ExcelFormatter(df, merge_cells=merge_cells)
    formatted_cells = formatter._format_hierarchical_rows()

    for cell in formatted_cells:
        if cell.row != 0 and cell.col == 0:
            assert isinstance(cell.val, Timestamp), (
                "Period should be converted to Timestamp"
            )
