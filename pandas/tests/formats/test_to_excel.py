"""Tests formatting as writer-agnostic ExcelCells

Most of the conversion to Excel is tested in pandas/tests/io/test_excel.py
"""

import pytest

from pandas.formats.format import CSSResolver, CSSWarning, CSSToExcelConverter


# Tests for CSSResolver


def assert_resolves(css, props, inherited=None):
    resolve = CSSResolver()
    actual = resolve(css, inherited=inherited)
    assert props == actual


def assert_same_resolution(css1, css2, inherited=None):
    resolve = CSSResolver()
    resolved1 = resolve(css1, inherited=inherited)
    resolved2 = resolve(css2, inherited=inherited)
    assert resolved1 == resolved2


@pytest.mark.parametrize('name,norm,abnorm', [
    ('whitespace', 'hello: world; foo: bar',
     ' \t hello \t :\n  world \n  ;  \n foo: \tbar\n\n'),
    ('case', 'hello: world; foo: bar', 'Hello: WORLD; foO: bar'),
    ('empty-decl', 'hello: world; foo: bar',
     '; hello: world;; foo: bar;\n; ;'),
    ('empty-list', '', ';'),
])
def test_css_parse_normalisation(name, norm, abnorm):
    assert_same_resolution(norm, abnorm)


@pytest.mark.xfail
def test_css_parse_comments():
    assert_same_resolution('hello: world',
                           'hello/* foo */:/* bar \n */ world /*;not:here*/')


@pytest.mark.xfail
def test_css_parse_strings():
    # semicolons in strings
    assert_resolves('background-image: url(\'http://blah.com/foo?a;b=c\')',
                    {'background-image': 'url(\'http://blah.com/foo?a;b=c\')'})
    assert_resolves('background-image: url("http://blah.com/foo?a;b=c")',
                    {'background-image': 'url("http://blah.com/foo?a;b=c")'})


def test_css_parse_invalid():
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
def test_css_side_shorthands(shorthand, expansions):
    top, right, bottom, left = expansions

    assert_resolves('%s: 1pt' % shorthand,
                    {top: '1pt', right: '1pt',
                     bottom: '1pt', left: '1pt'})

    assert_resolves('%s: 1pt 4pt' % shorthand,
                    {top: '1pt', right: '4pt',
                     bottom: '1pt', left: '4pt'})

    assert_resolves('%s: 1pt 4pt 2pt' % shorthand,
                    {top: '1pt', right: '4pt',
                     bottom: '2pt', left: '4pt'})

    assert_resolves('%s: 1pt 4pt 2pt 0pt' % shorthand,
                    {top: '1pt', right: '4pt',
                     bottom: '2pt', left: '0pt'})

    with pytest.warns(CSSWarning):
        assert_resolves('%s: 1pt 1pt 1pt 1pt 1pt' % shorthand,
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
@pytest.mark.parametrize('css,props', [
    ('background: blue', {'background-color': 'blue'}),
    ('background: fixed blue',
     {'background-color': 'blue', 'background-attachment': 'fixed'}),
])
def test_css_background_shorthand(css, props):
    assert_resolves(css, props)


@pytest.mark.xfail
@pytest.mark.parametrize('style,equiv', [
    ('border: 1px solid red',
     'border-width: 1px; border-style: solid; border-color: red'),
    ('border: solid red 1px',
     'border-width: 1px; border-style: solid; border-color: red'),
    ('border: red solid',
     'border-style: solid; border-color: red'),
])
def test_css_border_shorthand(style, equiv):
    assert_same_resolution(style, equiv)


@pytest.mark.parametrize('style,inherited,equiv', [
    ('margin: 1px; margin: 2px', '',
     'margin: 2px'),
    ('margin: 1px', 'margin: 2px',
     'margin: 1px'),
    ('margin: 1px; margin: inherit', 'margin: 2px',
     'margin: 2px'),
    ('margin: 1px; margin-top: 2px', '',
     'margin-left: 1px; margin-right: 1px; ' +
     'margin-bottom: 1px; margin-top: 2px'),
    ('margin-top: 2px', 'margin: 1px',
     'margin: 1px; margin-top: 2px'),
    ('margin: 1px', 'margin-top: 2px',
     'margin: 1px'),
    ('margin: 1px; margin-top: inherit', 'margin: 2px',
     'margin: 1px; margin-top: 2px'),
])
def test_css_precedence(style, inherited, equiv):
    resolve = CSSResolver()
    inherited_props = resolve(inherited)
    style_props = resolve(style, inherited=inherited_props)
    equiv_props = resolve(equiv)
    assert style_props == equiv_props


@pytest.mark.parametrize('style,equiv', [
    ('margin: 1px; margin-top: inherit',
     'margin-bottom: 1px; margin-right: 1px; margin-left: 1px'),
    ('margin-top: inherit', ''),
    ('margin-top: initial', ''),
])
def test_css_none_absent(style, equiv):
    assert_same_resolution(style, equiv)


@pytest.mark.parametrize('size,resolved', [
    ('xx-small', '6pt'),
    ('x-small', '%fpt' % 7.5),
    ('small', '%fpt' % 9.6),
    ('medium', '12pt'),
    ('large', '%fpt' % 13.5),
    ('x-large', '18pt'),
    ('xx-large', '24pt'),

    ('8px', '6pt'),
    ('1.25pc', '15pt'),
    ('.25in', '18pt'),
    ('02.54cm', '72pt'),
    ('25.4mm', '72pt'),
    ('101.6q', '72pt'),
    ('101.6q', '72pt'),
])
@pytest.mark.parametrize('relative_to',  # invariant to inherited size
                         [None, '16pt'])
def test_css_absolute_font_size(size, relative_to, resolved):
    if relative_to is None:
        inherited = None
    else:
        inherited = {'font-size': relative_to}
    assert_resolves('font-size: %s' % size, {'font-size': resolved},
                    inherited=inherited)


@pytest.mark.parametrize('size,relative_to,resolved', [
    ('1em', None, '12pt'),
    ('1.0em', None, '12pt'),
    ('1.25em', None, '15pt'),
    ('1em', '16pt', '16pt'),
    ('1.0em', '16pt', '16pt'),
    ('1.25em', '16pt', '20pt'),
    ('1rem', '16pt', '12pt'),
    ('1.0rem', '16pt', '12pt'),
    ('1.25rem', '16pt', '15pt'),
    ('100%', None, '12pt'),
    ('125%', None, '15pt'),
    ('100%', '16pt', '16pt'),
    ('125%', '16pt', '20pt'),
    ('2ex', None, '12pt'),
    ('2.0ex', None, '12pt'),
    ('2.50ex', None, '15pt'),
    ('inherit', '16pt', '16pt'),
    # TODO: smaller, larger

    ('smaller', None, '10pt'),
    ('smaller', '18pt', '15pt'),
    ('larger', None, '%fpt' % 14.4),
    ('larger', '15pt', '18pt'),
])
def test_css_relative_font_size(size, relative_to, resolved):
    if relative_to is None:
        inherited = None
    else:
        inherited = {'font-size': relative_to}
    assert_resolves('font-size: %s' % size, {'font-size': resolved},
                    inherited=inherited)


def test_css_font_size_invalid():
    pass  # TODO


# Tests for CSSToExcelConverter


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
    ('font-color: red', {'font': {'color': 'FF0000'}}),
    ('font-color: #ff0000', {'font': {'color': 'FF0000'}}),
    ('font-color: #f0a', {'font': {'color': 'FF00AA'}}),
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
