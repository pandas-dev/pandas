"""Tests formatting as writer-agnostic ExcelCells

ExcelFormatter is tested implicitly in pandas/tests/io/test_excel.py
"""

import pytest

import numpy as np

import pandas as pd
from pandas.compat import openpyxl_compat
from pandas.io.excel import ExcelWriter
from pandas import DataFrame
from pandas.io.formats.excel import CSSToExcelConverter
from pandas.io.formats.excel import ExcelFormatter
from pandas.util.testing import ensure_clean


@pytest.mark.parametrize('css,expected', [
    # FONT
    # - name
    ('font-family: foo,bar', {'font': {'name': 'foo'}}),
    ('font-family: "foo bar",baz', {'font': {'name': 'foo bar'}}),
    ('font-family: foo,\nbar', {'font': {'name': 'foo'}}),
    ('font-family: foo, bar,    baz', {'font': {'name': 'foo'}}),
    ('font-family: bar, foo', {'font': {'name': 'bar'}}),
    ('font-family: \'foo bar\', baz', {'font': {'name': 'foo bar'}}),
    ('font-family: \'foo \\\'bar\', baz', {'font': {'name': 'foo \'bar'}}),
    ('font-family: "foo \\"bar", baz', {'font': {'name': 'foo "bar'}}),
    ('font-family: "foo ,bar", baz', {'font': {'name': 'foo ,bar'}}),
    # - family
    ('font-family: serif', {'font': {'name': 'serif', 'family': 1}}),
    ('font-family: Serif', {'font': {'name': 'serif', 'family': 1}}),
    ('font-family: roman, serif', {'font': {'name': 'roman', 'family': 1}}),
    ('font-family: roman, sans-serif', {'font': {'name': 'roman',
                                                 'family': 2}}),
    ('font-family: roman, sans serif', {'font': {'name': 'roman'}}),
    ('font-family: roman, sansserif', {'font': {'name': 'roman'}}),
    ('font-family: roman, cursive', {'font': {'name': 'roman', 'family': 4}}),
    ('font-family: roman, fantasy', {'font': {'name': 'roman', 'family': 5}}),
    # - size
    ('font-size: 1em', {'font': {'size': 12}}),
    # - bold
    ('font-weight: 100', {'font': {'bold': False}}),
    ('font-weight: 200', {'font': {'bold': False}}),
    ('font-weight: 300', {'font': {'bold': False}}),
    ('font-weight: 400', {'font': {'bold': False}}),
    ('font-weight: normal', {'font': {'bold': False}}),
    ('font-weight: lighter', {'font': {'bold': False}}),
    ('font-weight: bold', {'font': {'bold': True}}),
    ('font-weight: bolder', {'font': {'bold': True}}),
    ('font-weight: 700', {'font': {'bold': True}}),
    ('font-weight: 800', {'font': {'bold': True}}),
    ('font-weight: 900', {'font': {'bold': True}}),
    # - italic
    # - underline
    ('text-decoration: underline',
     {'font': {'underline': 'single'}}),
    ('text-decoration: overline',
     {}),
    ('text-decoration: none',
     {}),
    # - strike
    ('text-decoration: line-through',
     {'font': {'strike': True}}),
    ('text-decoration: underline line-through',
     {'font': {'strike': True, 'underline': 'single'}}),
    ('text-decoration: underline; text-decoration: line-through',
     {'font': {'strike': True}}),
    # - color
    ('color: red', {'font': {'color': 'FF0000'}}),
    ('color: #ff0000', {'font': {'color': 'FF0000'}}),
    ('color: #f0a', {'font': {'color': 'FF00AA'}}),
    # - shadow
    ('text-shadow: none', {'font': {'shadow': False}}),
    ('text-shadow: 0px -0em 0px #CCC', {'font': {'shadow': False}}),
    ('text-shadow: 0px -0em 0px #999', {'font': {'shadow': False}}),
    ('text-shadow: 0px -0em 0px', {'font': {'shadow': False}}),
    ('text-shadow: 2px -0em 0px #CCC', {'font': {'shadow': True}}),
    ('text-shadow: 0px -2em 0px #CCC', {'font': {'shadow': True}}),
    ('text-shadow: 0px -0em 2px #CCC', {'font': {'shadow': True}}),
    ('text-shadow: 0px -0em 2px', {'font': {'shadow': True}}),
    ('text-shadow: 0px -2em', {'font': {'shadow': True}}),

    # FILL
    # - color, fillType
    ('background-color: red', {'fill': {'fgColor': 'FF0000',
                                        'patternType': 'solid'}}),
    ('background-color: #ff0000', {'fill': {'fgColor': 'FF0000',
                                            'patternType': 'solid'}}),
    ('background-color: #f0a', {'fill': {'fgColor': 'FF00AA',
                                         'patternType': 'solid'}}),
    # BORDER
    # - style
    ('border-style: solid',
     {'border': {'top': {'style': 'medium'},
                 'bottom': {'style': 'medium'},
                 'left': {'style': 'medium'},
                 'right': {'style': 'medium'}}}),
    ('border-style: solid; border-width: thin',
     {'border': {'top': {'style': 'thin'},
                 'bottom': {'style': 'thin'},
                 'left': {'style': 'thin'},
                 'right': {'style': 'thin'}}}),

    ('border-top-style: solid; border-top-width: thin',
     {'border': {'top': {'style': 'thin'}}}),
    ('border-top-style: solid; border-top-width: 1pt',
     {'border': {'top': {'style': 'thin'}}}),
    ('border-top-style: solid',
     {'border': {'top': {'style': 'medium'}}}),
    ('border-top-style: solid; border-top-width: medium',
     {'border': {'top': {'style': 'medium'}}}),
    ('border-top-style: solid; border-top-width: 2pt',
     {'border': {'top': {'style': 'medium'}}}),
    ('border-top-style: solid; border-top-width: thick',
     {'border': {'top': {'style': 'thick'}}}),
    ('border-top-style: solid; border-top-width: 4pt',
     {'border': {'top': {'style': 'thick'}}}),

    ('border-top-style: dotted',
     {'border': {'top': {'style': 'mediumDashDotDot'}}}),
    ('border-top-style: dotted; border-top-width: thin',
     {'border': {'top': {'style': 'dotted'}}}),
    ('border-top-style: dashed',
     {'border': {'top': {'style': 'mediumDashed'}}}),
    ('border-top-style: dashed; border-top-width: thin',
     {'border': {'top': {'style': 'dashed'}}}),
    ('border-top-style: double',
     {'border': {'top': {'style': 'double'}}}),
    # - color
    ('border-style: solid; border-color: #0000ff',
     {'border': {'top': {'style': 'medium', 'color': '0000FF'},
                 'right': {'style': 'medium', 'color': '0000FF'},
                 'bottom': {'style': 'medium', 'color': '0000FF'},
                 'left': {'style': 'medium', 'color': '0000FF'}}}),
    ('border-top-style: double; border-top-color: blue',
     {'border': {'top': {'style': 'double', 'color': '0000FF'}}}),
    ('border-top-style: solid; border-top-color: #06c',
     {'border': {'top': {'style': 'medium', 'color': '0066CC'}}}),
    # ALIGNMENT
    # - horizontal
    ('text-align: center',
     {'alignment': {'horizontal': 'center'}}),
    ('text-align: left',
     {'alignment': {'horizontal': 'left'}}),
    ('text-align: right',
     {'alignment': {'horizontal': 'right'}}),
    ('text-align: justify',
     {'alignment': {'horizontal': 'justify'}}),
    # - vertical
    ('vertical-align: top',
     {'alignment': {'vertical': 'top'}}),
    ('vertical-align: text-top',
     {'alignment': {'vertical': 'top'}}),
    ('vertical-align: middle',
     {'alignment': {'vertical': 'center'}}),
    ('vertical-align: bottom',
     {'alignment': {'vertical': 'bottom'}}),
    ('vertical-align: text-bottom',
     {'alignment': {'vertical': 'bottom'}}),
    # - wrap_text
    ('white-space: nowrap',
     {'alignment': {'wrap_text': False}}),
    ('white-space: pre',
     {'alignment': {'wrap_text': False}}),
    ('white-space: pre-line',
     {'alignment': {'wrap_text': False}}),
    ('white-space: normal',
     {'alignment': {'wrap_text': True}}),
])
def test_css_to_excel(css, expected):
    pytest.importorskip('cssdecl')
    convert = CSSToExcelConverter()
    assert expected == convert(css)


def test_css_to_excel_multiple():
    pytest.importorskip('cssdecl')
    convert = CSSToExcelConverter()
    actual = convert('''
        font-weight: bold;
        text-decoration: underline;
        color: red;
        border-width: thin;
        text-align: center;
        vertical-align: top;
        unused: something;
    ''')
    assert {"font": {"bold": True, "underline": "single", "color": "FF0000"},
            "border": {"top": {"style": "thin"},
                       "right": {"style": "thin"},
                       "bottom": {"style": "thin"},
                       "left": {"style": "thin"}},
            "alignment": {"horizontal": "center",
                          "vertical": "top"}} == actual


@pytest.mark.parametrize('css,inherited,expected', [
    ('font-weight: bold', '',
     {'font': {'bold': True}}),
    ('', 'font-weight: bold',
     {'font': {'bold': True}}),
    ('font-weight: bold', 'font-style: italic',
     {'font': {'bold': True, 'italic': True}}),
    ('font-style: normal', 'font-style: italic',
     {'font': {'italic': False}}),
    ('font-style: inherit', '', {}),
    ('font-style: normal; font-style: inherit', 'font-style: italic',
     {'font': {'italic': True}}),
])
def test_css_to_excel_inherited(css, inherited, expected):
    pytest.importorskip('cssdecl')
    convert = CSSToExcelConverter(inherited)
    assert expected == convert(css)


@pytest.fixture
def styled_dataframe():
    def style(df):
        # XXX: RGB colors not supported in xlwt
        return DataFrame([['font-weight: bold', '', ''],
                          ['', 'color: blue', ''],
                          ['', '', 'text-decoration: underline'],
                          ['border-style: solid', '', ''],
                          ['', 'font-style: italic', ''],
                          ['', '', 'text-align: right'],
                          ['background-color: red', '', ''],
                          ['', '', ''],
                          ['', '', ''],
                          ['', '', '']],
                         index=df.index, columns=df.columns)

    pytest.importorskip('jinja2')
    df = DataFrame(np.random.randn(10, 3))
    return df.style.apply(style, axis=None)


def assert_equal_style(cell1, cell2):
    # XXX: should find a better way to check equality
    # Neither OpenPyXl's Cell, nor its style objects have __eq__ defined
    assert cell1.alignment.__dict__ == cell2.alignment.__dict__
    assert cell1.border.__dict__ == cell2.border.__dict__
    assert cell1.fill.__dict__ == cell2.fill.__dict__
    assert cell1.font.__dict__ == cell2.font.__dict__
    assert cell1.number_format == cell2.number_format
    assert cell1.protection.__dict__ == cell2.protection.__dict__


@pytest.mark.parametrize('engine', [
    pytest.mark.xfail('xlwt', reason='xlwt does not support '
                                     'openpyxl-compatible style dicts'),
    'xlsxwriter',
    'openpyxl',
])
def test_styler_to_excel(engine, styled_dataframe):
    def custom_converter(css):
        # use bold iff there is custom style attached to the cell
        if css.strip(' \n;'):
            return {'font': {'bold': True}}
        return {}

    pytest.importorskip('cssdecl')
    pytest.importorskip(engine)

    if engine == 'openpyxl' and openpyxl_compat.is_compat(major_ver=1):
        pytest.xfail('openpyxl1 does not support some openpyxl2-compatible '
                     'style dicts')

    styled = styled_dataframe
    df = styled.data

    # Prepare spreadsheets

    with ensure_clean('.xlsx' if engine != 'xlwt' else '.xls') as path:
        writer = ExcelWriter(path, engine=engine)
        df.to_excel(writer, sheet_name='frame')
        df.style.to_excel(writer, sheet_name='unstyled')
        styled.to_excel(writer, sheet_name='styled')
        ExcelFormatter(styled, style_converter=custom_converter).write(
            writer, sheet_name='custom')

    # For engines other than openpyxl 2, we only smoke test
    if engine != 'openpyxl':
        return
    if not openpyxl_compat.is_compat(major_ver=2):
        pytest.skip('incompatible openpyxl version')

    # (1) compare DataFrame.to_excel and Styler.to_excel when unstyled
    n_cells = 0
    for col1, col2 in zip(writer.sheets['frame'].columns,
                          writer.sheets['unstyled'].columns):
        assert len(col1) == len(col2)
        for cell1, cell2 in zip(col1, col2):
            assert cell1.value == cell2.value
            assert_equal_style(cell1, cell2)
            n_cells += 1

    # ensure iteration actually happened:
    assert n_cells == (10 + 1) * (3 + 1)

    # (2) check styling with default converter
    n_cells = 0
    for col1, col2 in zip(writer.sheets['frame'].columns,
                          writer.sheets['styled'].columns):
        assert len(col1) == len(col2)
        for cell1, cell2 in zip(col1, col2):
            ref = '%s%d' % (cell2.column, cell2.row)
            # XXX: this isn't as strong a test as ideal; we should
            #      differences are exclusive
            if ref == 'B2':
                assert not cell1.font.bold
                assert cell2.font.bold
            elif ref == 'C3':
                assert cell1.font.color.rgb != cell2.font.color.rgb
                assert cell2.font.color.rgb == '000000FF'
            elif ref == 'D4':
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
                assert cell2.fill.fgColor.rgb == '00FF0000'
                assert cell2.fill.patternType == 'solid'
            else:
                assert_equal_style(cell1, cell2)

            assert cell1.value == cell2.value
            n_cells += 1

    assert n_cells == (10 + 1) * (3 + 1)

    # (3) check styling with custom converter
    n_cells = 0
    for col1, col2 in zip(writer.sheets['frame'].columns,
                          writer.sheets['custom'].columns):
        assert len(col1) == len(col2)
        for cell1, cell2 in zip(col1, col2):
            ref = '%s%d' % (cell2.column, cell2.row)
            if ref in ('B2', 'C3', 'D4', 'B5', 'C6', 'D7', 'B8'):
                assert not cell1.font.bold
                assert cell2.font.bold
            else:
                assert_equal_style(cell1, cell2)

            assert cell1.value == cell2.value
            n_cells += 1

    assert n_cells == (10 + 1) * (3 + 1)


@pytest.mark.parametrize('engine', [
    pytest.mark.xfail('xlwt', reason='xlwt does not support '
                                     'openpyxl-compatible style dicts'),
    'xlsxwriter',
    'openpyxl',
])
def test_styler_to_excel_no_cssdecl(engine, styled_dataframe):
    def custom_converter(css):
        # use bold iff there is custom style attached to the cell
        if css.strip(' \n;'):
            return {'font': {'bold': True}}
        return {}

    try:
        import cssdecl  # noqa
    except ImportError:
        with pytest.raises(ImportError) as rec:
            pd.DataFrame({"A": [1, 2]}).style.to_excel(engine)
        assert rec.match("not installed")
    else:
        pytest.skip('Test only run if cssdecl not installed')

    pytest.importorskip(engine)

    styled = styled_dataframe
    df = styled.data

    # Prepare spreadsheets

    with ensure_clean('.xlsx' if engine != 'xlwt' else '.xls') as path:
        writer = ExcelWriter(path, engine=engine)
        df.to_excel(writer, sheet_name='frame')
        with pytest.raises(ImportError):
            # default style_converter requires cssdecl
            df.style.to_excel(writer, sheet_name='unstyled')
            styled.to_excel(writer, sheet_name='styled')

        ExcelFormatter(styled, style_converter=custom_converter).write(
            writer, sheet_name='custom')

        writer.save()

        openpyxl = pytest.importorskip('openpyxl')
        if openpyxl_compat.is_compat(major_ver=1):
            # smoke test only
            return

        wb = openpyxl.load_workbook(path)

        # check styling with custom converter
        n_cells = 0
        for col1, col2 in zip(wb['frame'].columns,
                              wb['custom'].columns):
            assert len(col1) == len(col2)
            for cell1, cell2 in zip(col1, col2):
                ref = '%s%d' % (cell2.column, cell2.row)
                if ref in ('B2', 'C3', 'D4', 'B5', 'C6', 'D7', 'B8'):
                    assert not cell1.font.bold
                    assert cell2.font.bold
                else:
                    assert_equal_style(cell1, cell2)

                assert cell1.value == cell2.value
                n_cells += 1

        assert n_cells == (10 + 1) * (3 + 1)
