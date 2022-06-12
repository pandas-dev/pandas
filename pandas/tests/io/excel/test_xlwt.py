from datetime import (
    date,
    datetime,
)
import re

import numpy as np
import pytest

from pandas import (
    DataFrame,
    MultiIndex,
    options,
)
import pandas._testing as tm

from pandas.io.excel import (
    ExcelWriter,
    _XlwtWriter,
)

xlwt = pytest.importorskip("xlwt")

pytestmark = pytest.mark.parametrize("ext,", [".xls"])


def test_excel_raise_error_on_multiindex_columns_and_no_index(ext):
    # MultiIndex as columns is not yet implemented 9794
    cols = MultiIndex.from_tuples(
        [("site", ""), ("2014", "height"), ("2014", "weight")]
    )
    df = DataFrame(np.random.randn(10, 3), columns=cols)

    msg = (
        "Writing to Excel with MultiIndex columns and no index "
        "\\('index'=False\\) is not yet implemented."
    )
    with pytest.raises(NotImplementedError, match=msg):
        with tm.ensure_clean(ext) as path:
            df.to_excel(path, index=False)


def test_excel_multiindex_columns_and_index_true(ext):
    cols = MultiIndex.from_tuples(
        [("site", ""), ("2014", "height"), ("2014", "weight")]
    )
    df = DataFrame(np.random.randn(10, 3), columns=cols)
    with tm.ensure_clean(ext) as path:
        df.to_excel(path, index=True)


def test_excel_multiindex_index(ext):
    # MultiIndex as index works so assert no error #9794
    cols = MultiIndex.from_tuples(
        [("site", ""), ("2014", "height"), ("2014", "weight")]
    )
    df = DataFrame(np.random.randn(3, 10), index=cols)
    with tm.ensure_clean(ext) as path:
        df.to_excel(path, index=False)


def test_to_excel_styleconverter(ext):
    hstyle = {
        "font": {"bold": True},
        "borders": {"top": "thin", "right": "thin", "bottom": "thin", "left": "thin"},
        "alignment": {"horizontal": "center", "vertical": "top"},
    }

    xls_style = _XlwtWriter._convert_to_style(hstyle)
    assert xls_style.font.bold
    assert xlwt.Borders.THIN == xls_style.borders.top
    assert xlwt.Borders.THIN == xls_style.borders.right
    assert xlwt.Borders.THIN == xls_style.borders.bottom
    assert xlwt.Borders.THIN == xls_style.borders.left
    assert xlwt.Alignment.HORZ_CENTER == xls_style.alignment.horz
    assert xlwt.Alignment.VERT_TOP == xls_style.alignment.vert


def test_write_append_mode_raises(ext):
    msg = "Append mode is not supported with xlwt!"

    with tm.ensure_clean(ext) as f:
        with pytest.raises(ValueError, match=msg):
            ExcelWriter(f, engine="xlwt", mode="a")


def test_to_excel_xlwt_warning(ext):
    # GH 26552
    df = DataFrame(np.random.randn(3, 10))
    with tm.ensure_clean(ext) as path:
        with tm.assert_produces_warning(
            FutureWarning,
            match="As the xlwt package is no longer maintained",
        ):
            df.to_excel(path)


def test_option_xls_writer_deprecated(ext):
    # GH 26552
    with tm.assert_produces_warning(
        FutureWarning,
        match="As the xlwt package is no longer maintained",
        check_stacklevel=False,
    ):
        options.io.excel.xls.writer = "xlwt"


@pytest.mark.parametrize("style_compression", [0, 2])
def test_kwargs(ext, style_compression):
    # GH 42286
    kwargs = {"style_compression": style_compression}
    with tm.ensure_clean(ext) as f:
        msg = re.escape("Use of **kwargs is deprecated")
        with tm.assert_produces_warning(FutureWarning, match=msg):
            with ExcelWriter(f, engine="xlwt", **kwargs) as writer:
                assert (
                    writer.book._Workbook__styles.style_compression == style_compression
                )
                # xlwt won't allow us to close without writing something
                DataFrame().to_excel(writer)


@pytest.mark.parametrize("style_compression", [0, 2])
def test_engine_kwargs(ext, style_compression):
    # GH 42286
    engine_kwargs = {"style_compression": style_compression}
    with tm.ensure_clean(ext) as f:
        with ExcelWriter(f, engine="xlwt", engine_kwargs=engine_kwargs) as writer:
            assert writer.book._Workbook__styles.style_compression == style_compression
            # xlwt won't allow us to close without writing something
            DataFrame().to_excel(writer)


def test_book_and_sheets_consistent(ext):
    # GH#45687 - Ensure sheets is updated if user modifies book
    with tm.ensure_clean(ext) as f:
        with ExcelWriter(f) as writer:
            assert writer.sheets == {}
            sheet = writer.book.add_sheet("test_name")
            assert writer.sheets == {"test_name": sheet}


@pytest.mark.parametrize("attr", ["fm_date", "fm_datetime"])
def test_deprecated_attr(ext, attr):
    # GH#45572
    with tm.ensure_clean(ext) as path:
        with ExcelWriter(path, engine="xlwt") as writer:
            msg = f"{attr} is not part of the public API"
            with tm.assert_produces_warning(FutureWarning, match=msg):
                getattr(writer, attr)


def test_write_date_datetime_format(ext):
    # see gh-44284
    #
    # Test that custom date/datetime formats are respected
    # by inspecting formatting info in written file.
    xlrd = pytest.importorskip("xlrd")

    df = DataFrame(
        [
            [date(2014, 1, 31), datetime(1998, 5, 26, 23, 33, 4)],
            [date(1999, 9, 24), datetime(2014, 2, 28, 13, 5, 13)],
        ],
        index=["X", "Y"],
        columns=["DATE", "DATETIME"],
    )

    with tm.ensure_clean(ext) as f:
        with ExcelWriter(
            f,
            engine="xlwt",
            date_format="DD.MM.YYYY",
            datetime_format="DD.MM.YYYY HH-MM-SS",
        ) as writer:
            df.to_excel(writer, "test1")

        # formatting_info defaults to False
        # so have to use xlrd.open_workbook() directly
        with xlrd.open_workbook(f, formatting_info=True) as book:
            sh = book["test1"]
            xf_list = book.xf_list
            format_map = book.format_map

            date_cells = (sh[1, 1], sh[2, 1])
            assert all(
                format_map[xf_list[cell.xf_index].format_key].format_str == "DD.MM.YYYY"
                for cell in date_cells
            )

            datetime_cells = (sh[1, 2], sh[2, 2])
            assert all(
                format_map[xf_list[cell.xf_index].format_key].format_str
                == "DD.MM.YYYY HH-MM-SS"
                for cell in datetime_cells
            )
