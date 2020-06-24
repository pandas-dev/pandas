import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, MultiIndex
import pandas._testing as tm

from pandas.io.excel import ExcelWriter, _XlwtWriter

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
    df = pd.DataFrame(np.random.randn(10, 3), columns=cols)
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
