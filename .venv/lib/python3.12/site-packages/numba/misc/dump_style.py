try:
    from pygments.styles.default import DefaultStyle
except ImportError:
    msg = "Please install pygments to see highlighted dumps"
    raise ImportError(msg)

import numba.core.config
from pygments.styles.manni import ManniStyle
from pygments.styles.monokai import MonokaiStyle
from pygments.styles.native import NativeStyle

from pygments.lexer import RegexLexer, include, bygroups, words
from pygments.token import Text, Name, String,  Punctuation, Keyword, \
    Operator, Number

from pygments.style import Style


class NumbaIRLexer(RegexLexer):
    """
    Pygments style lexer for Numba IR (for use with highlighting etc).
    """
    name = 'Numba_IR'
    aliases = ['numba_ir']
    filenames = ['*.numba_ir']

    identifier = r'\$[a-zA-Z0-9._]+'
    fun_or_var = r'([a-zA-Z_]+[a-zA-Z0-9]*)'

    tokens = {
        'root' : [
            (r'(label)(\ [0-9]+)(:)$',
                bygroups(Keyword, Name.Label, Punctuation)),

            (r' = ', Operator),
            include('whitespace'),
            include('keyword'),

            (identifier, Name.Variable),
            (fun_or_var + r'(\()',
                bygroups(Name.Function, Punctuation)),
            (fun_or_var + r'(\=)',
                bygroups(Name.Attribute, Punctuation)),
            (fun_or_var, Name.Constant),
            (r'[0-9]+', Number),

            # <built-in function some>
            (r'<[^>\n]*>', String),

            (r'[=<>{}\[\]()*.,!\':]|x\b', Punctuation)
        ],

        'keyword':[
            (words((
                'del', 'jump', 'call', 'branch',
            ), suffix=' '), Keyword),
        ],

        'whitespace': [
            (r'(\n|\s)', Text),
        ],
    }


def by_colorscheme():
    """
    Get appropriate style for highlighting according to
    NUMBA_COLOR_SCHEME setting
    """
    styles = DefaultStyle.styles.copy()
    styles.update({
        Name.Variable:        "#888888",
    })
    custom_default = type('CustomDefaultStyle', (Style, ), {'styles': styles})

    style_map = {
        'no_color' : custom_default,
        'dark_bg' : MonokaiStyle,
        'light_bg' : ManniStyle,
        'blue_bg' : NativeStyle,
        'jupyter_nb' : DefaultStyle,
    }

    return style_map[numba.core.config.COLOR_SCHEME]
