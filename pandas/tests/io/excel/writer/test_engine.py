from distutils.version import LooseVersion
import os

import numpy as np
import pytest

from pandas.compat import PY36
import pandas.util._test_decorators as td

import pandas as pd
from pandas import DataFrame
import pandas.util.testing as tm
from pandas.util.testing import ensure_clean

from pandas.io.excel import (
    ExcelFile, ExcelWriter, _OpenpyxlWriter, _XlsxWriter, _XlwtWriter,
    register_writer)
from pandas.io.formats.excel import ExcelFormatter


class TestExcelWriterEngineTests(object):

    @pytest.mark.parametrize('klass,ext', [
        pytest.param(_XlsxWriter, '.xlsx', marks=pytest.mark.skipif(
            not td.safe_import('xlsxwriter'), reason='No xlsxwriter')),
        pytest.param(_OpenpyxlWriter, '.xlsx', marks=pytest.mark.skipif(
            not td.safe_import('openpyxl'), reason='No openpyxl')),
        pytest.param(_XlwtWriter, '.xls', marks=pytest.mark.skipif(
            not td.safe_import('xlwt'), reason='No xlwt'))
    ])
    def test_ExcelWriter_dispatch(self, klass, ext):
        with ensure_clean(ext) as path:
            writer = ExcelWriter(path)
            if ext == '.xlsx' and td.safe_import('xlsxwriter'):
                # xlsxwriter has preference over openpyxl if both installed
                assert isinstance(writer, _XlsxWriter)
            else:
                assert isinstance(writer, klass)

    def test_ExcelWriter_dispatch_raises(self):
        with pytest.raises(ValueError, match='No engine'):
            ExcelWriter('nothing')

    @pytest.mark.filterwarnings("ignore:\\nPanel:FutureWarning")
    def test_register_writer(self):
        # some awkward mocking to test out dispatch and such actually works
        called_save = []
        called_write_cells = []

        class DummyClass(ExcelWriter):
            called_save = False
            called_write_cells = False
            supported_extensions = ['test', 'xlsx', 'xls']
            engine = 'dummy'

            def save(self):
                called_save.append(True)

            def write_cells(self, *args, **kwargs):
                called_write_cells.append(True)

        def check_called(func):
            func()
            assert len(called_save) >= 1
            assert len(called_write_cells) >= 1
            del called_save[:]
            del called_write_cells[:]

        with pd.option_context('io.excel.xlsx.writer', 'dummy'):
            register_writer(DummyClass)
            writer = ExcelWriter('something.test')
            assert isinstance(writer, DummyClass)
            df = tm.makeCustomDataframe(1, 1)

            func = lambda: df.to_excel('something.test')
            check_called(func)
            check_called(lambda: df.to_excel('something.xlsx'))
            check_called(
                lambda: df.to_excel(
                    'something.xls', engine='dummy'))


@pytest.mark.parametrize('engine', [
    pytest.param('xlwt',
                 marks=pytest.mark.xfail(reason='xlwt does not support '
                                                'openpyxl-compatible '
                                                'style dicts')),
    'xlsxwriter',
    'openpyxl',
])
def test_styler_to_excel(engine):
    def style(df):
        # XXX: RGB colors not supported in xlwt
        return DataFrame([['font-weight: bold', '', ''],
                          ['', 'color: blue', ''],
                          ['', '', 'text-decoration: underline'],
                          ['border-style: solid', '', ''],
                          ['', 'font-style: italic', ''],
                          ['', '', 'text-align: right'],
                          ['background-color: red', '', ''],
                          ['number-format: 0%', '', ''],
                          ['', '', ''],
                          ['', '', ''],
                          ['', '', '']],
                         index=df.index, columns=df.columns)

    def assert_equal_style(cell1, cell2):
        # XXX: should find a better way to check equality
        assert cell1.alignment.__dict__ == cell2.alignment.__dict__
        assert cell1.border.__dict__ == cell2.border.__dict__
        assert cell1.fill.__dict__ == cell2.fill.__dict__
        assert cell1.font.__dict__ == cell2.font.__dict__
        assert cell1.number_format == cell2.number_format
        assert cell1.protection.__dict__ == cell2.protection.__dict__

    def custom_converter(css):
        # use bold iff there is custom style attached to the cell
        if css.strip(' \n;'):
            return {'font': {'bold': True}}
        return {}

    pytest.importorskip('jinja2')
    pytest.importorskip(engine)

    # Prepare spreadsheets

    df = DataFrame(np.random.randn(11, 3))
    with ensure_clean('.xlsx' if engine != 'xlwt' else '.xls') as path:
        writer = ExcelWriter(path, engine=engine)
        df.to_excel(writer, sheet_name='frame')
        df.style.to_excel(writer, sheet_name='unstyled')
        styled = df.style.apply(style, axis=None)
        styled.to_excel(writer, sheet_name='styled')
        ExcelFormatter(styled, style_converter=custom_converter).write(
            writer, sheet_name='custom')
        writer.save()

        if engine not in ('openpyxl', 'xlsxwriter'):
            # For other engines, we only smoke test
            return
        openpyxl = pytest.importorskip('openpyxl')
        wb = openpyxl.load_workbook(path)

        # (1) compare DataFrame.to_excel and Styler.to_excel when unstyled
        n_cells = 0
        for col1, col2 in zip(wb['frame'].columns,
                              wb['unstyled'].columns):
            assert len(col1) == len(col2)
            for cell1, cell2 in zip(col1, col2):
                assert cell1.value == cell2.value
                assert_equal_style(cell1, cell2)
                n_cells += 1

        # ensure iteration actually happened:
        assert n_cells == (11 + 1) * (3 + 1)

        # (2) check styling with default converter

        # XXX: openpyxl (as at 2.4) prefixes colors with 00, xlsxwriter with FF
        alpha = '00' if engine == 'openpyxl' else 'FF'

        n_cells = 0
        for col1, col2 in zip(wb['frame'].columns,
                              wb['styled'].columns):
            assert len(col1) == len(col2)
            for cell1, cell2 in zip(col1, col2):
                ref = '%s%d' % (cell2.column, cell2.row)
                # XXX: this isn't as strong a test as ideal; we should
                #      confirm that differences are exclusive
                if ref == 'B2':
                    assert not cell1.font.bold
                    assert cell2.font.bold
                elif ref == 'C3':
                    assert cell1.font.color.rgb != cell2.font.color.rgb
                    assert cell2.font.color.rgb == alpha + '0000FF'
                elif ref == 'D4':
                    # This fails with engine=xlsxwriter due to
                    # https://bitbucket.org/openpyxl/openpyxl/issues/800
                    if engine == 'xlsxwriter' \
                       and (LooseVersion(openpyxl.__version__) <
                            LooseVersion('2.4.6')):
                        pass
                    else:
                        assert cell1.font.underline != cell2.font.underline
                        assert cell2.font.underline == 'single'
                elif ref == 'B5':
                    assert not cell1.border.left.style
                    assert (cell2.border.top.style ==
                            cell2.border.right.style ==
                            cell2.border.bottom.style ==
                            cell2.border.left.style ==
                            'medium')
                elif ref == 'C6':
                    assert not cell1.font.italic
                    assert cell2.font.italic
                elif ref == 'D7':
                    assert (cell1.alignment.horizontal !=
                            cell2.alignment.horizontal)
                    assert cell2.alignment.horizontal == 'right'
                elif ref == 'B8':
                    assert cell1.fill.fgColor.rgb != cell2.fill.fgColor.rgb
                    assert cell1.fill.patternType != cell2.fill.patternType
                    assert cell2.fill.fgColor.rgb == alpha + 'FF0000'
                    assert cell2.fill.patternType == 'solid'
                elif ref == 'B9':
                    assert cell1.number_format == 'General'
                    assert cell2.number_format == '0%'
                else:
                    assert_equal_style(cell1, cell2)

                assert cell1.value == cell2.value
                n_cells += 1

        assert n_cells == (11 + 1) * (3 + 1)

        # (3) check styling with custom converter
        n_cells = 0
        for col1, col2 in zip(wb['frame'].columns,
                              wb['custom'].columns):
            assert len(col1) == len(col2)
            for cell1, cell2 in zip(col1, col2):
                ref = '%s%d' % (cell2.column, cell2.row)
                if ref in ('B2', 'C3', 'D4', 'B5', 'C6', 'D7', 'B8', 'B9'):
                    assert not cell1.font.bold
                    assert cell2.font.bold
                else:
                    assert_equal_style(cell1, cell2)

                assert cell1.value == cell2.value
                n_cells += 1

        assert n_cells == (11 + 1) * (3 + 1)


@td.skip_if_no('openpyxl')
@pytest.mark.skipif(not PY36, reason='requires fspath')
class TestFSPath(object):

    def test_excelfile_fspath(self):
        with tm.ensure_clean('foo.xlsx') as path:
            df = DataFrame({"A": [1, 2]})
            df.to_excel(path)
            xl = ExcelFile(path)
            result = os.fspath(xl)
            assert result == path

    def test_excelwriter_fspath(self):
        with tm.ensure_clean('foo.xlsx') as path:
            writer = ExcelWriter(path)
            assert os.fspath(writer) == str(path)
