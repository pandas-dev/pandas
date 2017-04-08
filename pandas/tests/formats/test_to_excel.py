"""Tests formatting as writer-agnostic ExcelCells

ExcelFormatter is tested implicitly in pandas/tests/io/test_excel.py
"""

import pytest

from pandas.formats.excel import CSSToExcelConverter


@pytest.mark.parametrize('css,expected', [
    # FONT
    # - name
    ('font-family: foo,bar', {'font': {'name': 'foo'}}),
    pytest.mark.xfail(('font-family: "foo bar",baz',
                       {'font': {'name': 'foo bar'}})),
    ('font-family: foo,\nbar', {'font': {'name': 'foo'}}),
    ('font-family: foo, bar,    baz', {'font': {'name': 'foo'}}),
    ('font-family: bar, foo', {'font': {'name': 'bar'}}),
    # - family
    ('font-family: serif', {'font': {'name': 'serif', 'family': 1}}),
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
     {'font': {'underline': 'single', 'strike': False}}),
    ('text-decoration: overline',
     {'font': {'strike': False}}),
    ('text-decoration: none',
     {'font': {'strike': False}}),
    # - strike
    ('text-decoration: line-through',
     {'font': {'strike': True, 'underline': False}}),
    ('text-decoration: underline line-through',
     {'font': {'strike': True, 'underline': True}}),
    ('text-decoration: underline; text-decoration: line-through',
     {'font': {'strike': True, 'underline': False}}),
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
    # text-shadow with color first not yet implemented
    pytest.mark.xfail(('text-shadow: #CCC 3px 3px 3px',
                       {'font': {'shadow': True}})),
    pytest.mark.xfail(('text-shadow: #999 0px 0px 0px',
                       {'font': {'shadow': False}})),
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
    # TODO: need to check this produces valid OpenXML without color
    ('border-style: solid',
     {'border': {'top': {'style': 'medium'},
                 'bottom': {'style': 'medium'},
                 'left': {'style': 'medium'},
                 'right': {'style': 'medium'}}}),
    ('border-top-style: solid',
     {'border': {'top': {'style': 'medium'}}}),
    ('border-top-style: dotted',
     {'border': {'top': {'style': 'mediumDashDotDot'}}}),
    ('border-top-style: dashed',
     {'border': {'top': {'style': 'mediumDashed'}}}),
    # - color
    # ALIGNMENT
    # - horizontal
    # - vertical
    # - wrap_text
])
def test_css_to_excel(css, expected):
    convert = CSSToExcelConverter()
    assert expected == convert(css)


def test_css_to_excel_multiple():
    convert = CSSToExcelConverter()
    actual = convert('''
        font-weight: bold;
        border-width: thin;
        text-align: center;
        vertical-align: top;
        unused: something;
    ''')
    assert {"font": {"bold": True},
            "border": {"top": {"style": "hair"},
                       "right": {"style": "hair"},
                       "bottom": {"style": "hair"},
                       "left": {"style": "hair"}},
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
    convert = CSSToExcelConverter(inherited)
    assert expected == convert(css)


@pytest.mark.xfail
def test_css_to_excel_warns_when_not_supported():
    pass
