import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd
from pandas import DataFrame, MultiIndex
from pandas.util.testing import ensure_clean

from pandas.io.excel import ExcelWriter, _XlwtWriter
from pandas.tests.io.excel.writer.test_base import _WriterBase


@td.skip_if_no('xlwt')
@pytest.mark.parametrize("merge_cells,ext,engine", [
    (None, '.xls', 'xlwt')])
class TestXlwtTests(_WriterBase):

    def test_excel_raise_error_on_multiindex_columns_and_no_index(
            self, merge_cells, ext, engine):
        # MultiIndex as columns is not yet implemented 9794
        cols = MultiIndex.from_tuples([('site', ''),
                                       ('2014', 'height'),
                                       ('2014', 'weight')])
        df = DataFrame(np.random.randn(10, 3), columns=cols)
        with pytest.raises(NotImplementedError):
            with ensure_clean(ext) as path:
                df.to_excel(path, index=False)

    def test_excel_multiindex_columns_and_index_true(self, merge_cells, ext,
                                                     engine):
        cols = MultiIndex.from_tuples([('site', ''),
                                       ('2014', 'height'),
                                       ('2014', 'weight')])
        df = pd.DataFrame(np.random.randn(10, 3), columns=cols)
        with ensure_clean(ext) as path:
            df.to_excel(path, index=True)

    def test_excel_multiindex_index(self, merge_cells, ext, engine):
        # MultiIndex as index works so assert no error #9794
        cols = MultiIndex.from_tuples([('site', ''),
                                       ('2014', 'height'),
                                       ('2014', 'weight')])
        df = DataFrame(np.random.randn(3, 10), index=cols)
        with ensure_clean(ext) as path:
            df.to_excel(path, index=False)

    def test_to_excel_styleconverter(self, merge_cells, ext, engine):
        import xlwt

        hstyle = {"font": {"bold": True},
                  "borders": {"top": "thin",
                              "right": "thin",
                              "bottom": "thin",
                              "left": "thin"},
                  "alignment": {"horizontal": "center", "vertical": "top"}}

        xls_style = _XlwtWriter._convert_to_style(hstyle)
        assert xls_style.font.bold
        assert xlwt.Borders.THIN == xls_style.borders.top
        assert xlwt.Borders.THIN == xls_style.borders.right
        assert xlwt.Borders.THIN == xls_style.borders.bottom
        assert xlwt.Borders.THIN == xls_style.borders.left
        assert xlwt.Alignment.HORZ_CENTER == xls_style.alignment.horz
        assert xlwt.Alignment.VERT_TOP == xls_style.alignment.vert

    def test_write_append_mode_raises(self, merge_cells, ext, engine):
        msg = "Append mode is not supported with xlwt!"

        with ensure_clean(ext) as f:
            with pytest.raises(ValueError, match=msg):
                ExcelWriter(f, engine=engine, mode='a')
