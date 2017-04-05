"""Tests formatting as writer-agnostic ExcelCells

Most of the conversion to Excel is tested in pandas/tests/io/test_excel.py
"""

import pytest

from pandas.formats.format import CSSResolver, CSSWarning, CSSToExcelConverter


# Test parsing and normalising of CSS


def assert_resolves(css, props, inherited=None):
    resolve = CSSResolver()
    actual = resolve(css, inherited=inherited)
    assert props == actual


def test_css_parse_whitespace():
    pass  # TODO


def test_css_parse_case():
    pass  # TODO


def test_css_parse_empty():
    pass  # TODO


def test_css_parse_invalid():
    pass  # TODO


@pytest.mark.xfail
def test_css_parse_comments():
    pass  # TODO


@pytest.mark.xfail
def test_css_parse_strings():
    pass  # TODO


@pytest.mark.parametrize(
    'shorthand,expansions',
    [('margin', ['margin-top', 'margin-right',
                 'margin-bottom', 'margin-left']),
     ('padding', ['padding-top', 'padding-right',
                  'padding-bottom', 'padding-left']),
     ('border-width', ['border-top-width', 'border-right-width',
                       'border-bottom-width', 'border-left-width']),
     ('border-color', ['border-top-color', 'border-right-color',
                       'border-bottom-color', 'border-left-color']),
     ('border-style', ['border-top-style', 'border-right-style',
                       'border-bottom-style', 'border-left-style']),
     ])
def test_css_direction_shorthands(shorthand, expansions):
    top, right, bottom, left = expansions

    assert_resolves('%s: thin' % shorthand,
                    {top: 'thin', right: 'thin',
                     bottom: 'thin', left: 'thin'})

    assert_resolves('%s: thin thick' % shorthand,
                    {top: 'thin', right: 'thick',
                     bottom: 'thin', left: 'thick'})

    assert_resolves('%s: thin thick medium' % shorthand,
                    {top: 'thin', right: 'thick',
                     bottom: 'medium', left: 'thick'})

    assert_resolves('%s: thin thick medium none' % shorthand,
                    {top: 'thin', right: 'thick',
                     bottom: 'medium', left: 'none'})

    with pytest.warns(CSSWarning):
        assert_resolves('%s: thin thick medium none medium' % shorthand,
                        {})


@pytest.mark.xfail
@pytest.mark.parametrize('css,props', [
    ('font: italic bold 12pt helvetica,sans-serif',
     {'font-family': 'helvetica,sans-serif',
      'font-style': 'italic',
      'font-weight': 'bold',
      'font-size': '12pt'}),
    ('font: bold italic 12pt helvetica,sans-serif',
     {'font-family': 'helvetica,sans-serif',
      'font-style': 'italic',
      'font-weight': 'bold',
      'font-size': '12pt'}),
])
def test_css_font_shorthand(css, props):
    assert_resolves(css, props)


@pytest.mark.xfail
def test_css_background_shorthand():
    pass  # TODO


def test_css_override():
    pass  # TODO


def test_css_override_inherited():
    pass  # TODO


def test_css_default_inherited():
    pass  # TODO


def test_css_none_absent():
    pass  # TODO


def test_css_font_size():
    pass  # TODO


def test_css_font_size_invalid():
    pass  # TODO


# Test translation of CSS to ExcelCell.style values


@pytest.mark.parametrize('css,expected', [
    # FONT
    # - name
    # - family
    # - size
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
     {'font': {'underline': True, 'strike': False}}),
    ('text-decoration: overline',
     {'font': {'underline': False, 'strike': False}}),
    ('text-decoration: none',
     {'font': {'underline': False, 'strike': False}}),
    # - strike
    ('text-decoration: line-through',
     {'font': {'strike': True, 'underline': False}}),
    ('text-decoration: underline line-through',
     {'font': {'strike': True, 'underline': True}}),
    ('text-decoration: underline; text-decoration: line-through',
     {'font': {'strike': True, 'underline': False}}),
    # - color
    # - shadow
    # FILL
    # - color, fillType
    # BORDER
    # - style
    # - color
    # ALIGNMENT
    # - horizontal
    # - vertical
    # - wrap_text
])
def test_css_to_excel(css, expected):
    convert = CSSToExcelConverter()
    assert expected == convert(css)


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
