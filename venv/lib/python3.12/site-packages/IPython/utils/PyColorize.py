import keyword
import os
import sys
import token
import tokenize
import warnings
from io import StringIO
from typing import Any, Type, TypeAlias

import pygments
from pygments.formatters.terminal256 import Terminal256Formatter
from pygments.style import Style
from pygments.styles import get_style_by_name
from pygments.token import Token, _TokenType
from functools import cache

from typing import TypedDict


TokenStream: TypeAlias = list[tuple[_TokenType, str]]


__all__ = ["Parser", "Theme"]


class Symbols(TypedDict):
    top_line: str
    arrow_body: str
    arrow_head: str


_default_symbols: Symbols = Symbols(
    top_line="-",
    arrow_body="-",
    arrow_head=">",
)


class Theme:
    name: str
    base: str | None
    extra_style: dict[_TokenType, str]
    symbols: Symbols

    def __init__(
        self,
        name: str,
        base: str | None,
        extra_style: dict[_TokenType, str],
        *,
        symbols: Symbols | None = None,
    ) -> None:
        self.name = name
        self.base = base
        self.extra_style = extra_style
        s: Symbols = symbols if symbols is not None else _default_symbols
        self.symbols = {**_default_symbols, **s}
        self._formatter = Terminal256Formatter(style=self.as_pygments_style())

    @cache
    def as_pygments_style(self) -> Type[Style]:
        if self.base is not None:
            base_styles = get_style_by_name(self.base).styles
        else:
            base_styles = {}

        class MyStyle(Style):
            styles = {**base_styles, **self.extra_style}

        return MyStyle

    def format(self, stream: TokenStream) -> str:
        return pygments.format(stream, self._formatter)

    def make_arrow(self, width: int) -> str:
        """generate the leading arrow in front of traceback or debugger"""
        if width >= 2:
            return (
                self.symbols["arrow_body"] * (width - 2)
                + self.symbols["arrow_head"]
                + " "
            )
        elif width == 1:
            return self.symbols["arrow_head"]
        return ""


generate_tokens = tokenize.generate_tokens


#############################################################################
### Python Source Parser (does Highlighting)
#############################################################################

_KEYWORD = token.NT_OFFSET + 1
_TEXT = token.NT_OFFSET + 2

# ****************************************************************************

_pygment_token_mapping: dict[int, _TokenType] = {
    token.NUMBER: Token.Literal.Number,
    token.OP: Token.Operator,
    token.STRING: Token.Literal.String,
    token.COMMENT: Token.Comment,
    token.NAME: Token.Name,
    token.ERRORTOKEN: Token.Error,
    _KEYWORD: Token.Keyword,
    _TEXT: Token.Text,
}

# technically BW is not nocolor, we should have a no-style, style
nocolors_theme = Theme("nocolor", None, {})


linux_theme = Theme(
    "linux",
    "monokai",
    {
        Token.Header: "ansibrightred",
        Token.LinenoEm: "ansibrightgreen",
        Token.Lineno: "ansigreen",
        Token.ValEm: "ansibrightblue",
        Token.VName: "ansicyan",
        Token.Caret: "",
        Token.Filename: "ansibrightgreen",
        Token.ExcName: "ansibrightred",
        Token.Topline: "ansibrightred",
        Token.FilenameEm: "ansigreen",
        Token.Normal: "",
        Token.NormalEm: "ansibrightcyan",
        Token.Line: "ansiyellow",
        Token.TB.Name: "ansimagenta",
        Token.TB.NameEm: "ansibrightmagenta",
        Token.Breakpoint: "",
        Token.Breakpoint.Enabled: "ansibrightred",
        Token.Breakpoint.Disabled: "ansired",
        Token.Prompt: "ansibrightgreen",
        Token.PromptNum: "ansigreen bold",
        Token.OutPrompt: "ansibrightred",
        Token.OutPromptNum: "ansired bold",
    },
)

neutral_pygments_equiv = {
    Token.Header: "ansired",
    Token.LinenoEm: "ansigreen",
    Token.Lineno: "ansibrightgreen",
    Token.ValEm: "ansiblue",
    Token.VName: "ansicyan",
    Token.Caret: "",
    Token.Filename: "ansibrightgreen",
    Token.FilenameEm: "ansigreen",
    Token.ExcName: "ansired",
    Token.Topline: "ansired",
    Token.Normal: "",
    Token.NormalEm: "ansicyan",
    Token.Line: "ansired",
    Token.TB.Name: "ansibrightmagenta",
    Token.TB.NameEm: "ansimagenta",
    Token.Breakpoint: "",
    Token.Breakpoint.Enabled: "ansibrightred",
    Token.Breakpoint.Disabled: "ansired",
    ## specific override of pygments defaults for visibility
    Token.Number: "ansigreen",
    Token.Operator: "noinherit",
    Token.String: "ansiyellow",
    Token.Name.Function: "ansiblue",
    Token.Name.Class: "bold ansiblue",
    Token.Name.Namespace: "bold ansiblue",
    Token.Name.Variable.Magic: "ansiblue",
    Token.Prompt: "ansigreen",
    Token.OutPrompt: "ansired",
}


neutral_pygments_nt = {
    **neutral_pygments_equiv,
    Token.PromptNum: "ansigreen bold",
    Token.OutPromptNum: "ansired bold",
}
neutral_pygments_posix = {
    **neutral_pygments_equiv,
    Token.PromptNum: "ansibrightgreen bold",
    Token.OutPromptNum: "ansibrightred bold",
}


neutral_nt = Theme("neutral:nt", "default", neutral_pygments_nt)
neutral_posix = Theme("neutral:posix", "default", neutral_pygments_posix)


# Hack: the 'neutral' colours are not very visible on a dark background on
# Windows. Since Windows command prompts have a dark background by default, and
# relatively few users are likely to alter that, we will use the 'Linux' colours,
# designed for a dark background, as the default on Windows. Changing it here
# avoids affecting the prompt colours rendered by prompt_toolkit, where the
# neutral defaults do work OK.
if os.name == "nt":
    neutral_theme = neutral_nt
else:
    neutral_theme = neutral_posix


lightbg_theme = Theme(
    "lightbg",
    "pastie",
    {
        Token.Header: "ansired",
        Token.LinenoEm: "ansigreen",
        Token.Lineno: "ansibrightgreen",
        Token.ValEm: "ansiblue",
        Token.VName: "ansicyan",
        Token.Caret: "",
        Token.Filename: "ansigreen",
        Token.FilenameEm: "ansibrightgreen",
        Token.ExcName: "ansired",
        Token.Topline: "ansired",
        Token.Normal: "",
        Token.NormalEm: "ansicyan",
        Token.Line: "ansired",
        Token.TB.Name: "ansibrightmagenta",
        Token.TB.NameEm: "ansimagenta",
        Token.Breakpoint: "",
        Token.Breakpoint.Enabled: "ansibrightred",
        Token.Breakpoint.Disabled: "ansired",
        Token.Prompt: "ansibrightblue",
        Token.PromptNum: "ansiblue bold",
        Token.OutPrompt: "ansibrightred",
        Token.OutPromptNum: "ansired bold",
    },
)

PRIDE_RED = "#E40303"
PRIDE_ORANGE = "#FF8C00"
PRIDE_YELLOW = "#FFED00"
PRIDE_GREEN = "#008026"
PRIDE_INDIGO = "#004CFF"
PRIDE_VIOLET = "#732982"
pride_theme = Theme(
    "pride",
    "pastie",
    {
        Token.Header: PRIDE_INDIGO,
        Token.LinenoEm: f"{PRIDE_GREEN} italic",
        Token.Lineno: f"{PRIDE_GREEN} bold",
        Token.ValEm: f"{PRIDE_INDIGO} italic",
        Token.VName: "ansicyan",
        Token.Caret: "",
        Token.Filename: f"{PRIDE_YELLOW}",
        Token.FilenameEm: f"bg:{PRIDE_VIOLET}",
        Token.ExcName: f"{PRIDE_ORANGE}",
        Token.Topline: f"{PRIDE_RED}",
        Token.Normal: "",
        Token.NormalEm: "bold",
        Token.Line: "ansired",
        Token.TB.Name: "ansibrightmagenta",
        Token.TB.NameEm: "ansimagenta",
        Token.Breakpoint: "",
        Token.Breakpoint.Enabled: "ansibrightred",
        Token.Breakpoint.Disabled: "ansired",
        Token.Prompt: "ansibrightblue",
        Token.Prompt.Continuation.L1: f"ansiwhite bg:{PRIDE_RED}",
        Token.Prompt.Continuation.L2: f"ansiwhite bg:{PRIDE_ORANGE}",
        Token.Prompt.Continuation.L3: f"ansiblack bg:{PRIDE_YELLOW}",
        Token.Prompt.Continuation.L4: f"ansiwhite bg:{PRIDE_GREEN}",
        Token.Prompt.Continuation.L5: f"ansiwhite bg:{PRIDE_INDIGO}",
        Token.Prompt.Continuation.L6: f"ansiwhite bg:{PRIDE_VIOLET}",
        Token.PromptNum: "ansiblue bold",
        Token.OutPrompt: "ansibrightred",
        Token.OutPromptNum: "ansired bold",
    },
    symbols={"arrow_body": "\u2500", "arrow_head": "\u25b6", "top_line": "\u2500"},
)


C1 = "#D52D00"
C2 = "#EF7627"
C3 = "#FF9A56"
White = "#FFFFFF"
C5 = "#D162A4"
C6 = "#B55690"
C7 = "#A30262"

pl = {
    # Token.Whitespace: "#bbbbbb",
    Token.Comment: "#888888",
    Token.String: C5,
    Token.String.Escape: C1,
    Token.Keyword: f"italic {C2}",
    Token.Name.Class: C2,
    Token.Name.Exception: C1,
    Token.Name.Builtin: C3,
    Token.Name.Variable: C6,
    Token.Name.Constant: C7,
    Token.Name.Decorator: C2,
    Token.Number: C7,
    Token.Generic.Deleted: f"bg:{C1} #000000",
    Token.Generic.Emph: "italic",
    Token.Generic.Strong: "bold",
    Token.Generic.EmphStrong: "bold italic",
}

pridel_theme = Theme(
    "pride:l",
    None,
    {
        Token.Header: C3,
        Token.LinenoEm: C3,
        Token.Lineno: C2,
        Token.ValEm: C2,
        Token.VName: C2,
        Token.Caret: "",
        Token.Filename: C2,
        Token.FilenameEm: C3,
        Token.ExcName: C1,
        Token.Topline: C1,
        Token.Normal: "",
        Token.NormalEm: "bold",
        Token.Line: C2,
        Token.TB.Name: C6,
        Token.TB.NameEm: C7,
        Token.Breakpoint: "",
        Token.Breakpoint.Enabled: C1,
        Token.Breakpoint.Disabled: C7,
        Token.Prompt: C1,
        Token.PromptNum: C2,
        Token.Prompt.Continuation: C7,
        Token.Prompt.Continuation.L1: C2,
        Token.Prompt.Continuation.L2: C3,
        Token.Prompt.Continuation.L3: White,
        Token.Prompt.Continuation.L4: C5,
        Token.Prompt.Continuation.L5: C6,
        Token.Prompt.Continuation.L6: C7,
        Token.OutPrompt: C6,
        Token.OutPromptNum: C5,
        **pl,
    },
    symbols={"arrow_body": "\u2500", "arrow_head": "\u25b6", "top_line": "\u2500"},
)

GRUVBOX_VAL_EM = "#D79921"
GRUVBOX_V_NAME = "#83A598"
GRUVBOX_FILENAME = "#FBF1C7"
GRUVBOX_EXCEPTION_NAME = "#FB4934"
GRUVBOX_TOPLINE = "#CC241D"
GRUVBOX_BREAKPOINT_ENABLED = "#FB4934"
GRUVBOX_BREAKPOINT_DISABLED = "#CC241D"
GRUVBOX_PROMPT = "#689D6A"
GRUVBOX_PROMPT_NUM = "#8EC07C"
GRUVBOX_OUT_PROMPT = "#B16286"
GRUVBOX_OUT_PROMPT_NUM = "#D3869B"
gruvbox_dark_theme = Theme(
    "gruvbox-dark",
    "gruvbox-dark",
    {
        Token.Lineno: GRUVBOX_PROMPT_NUM,
        Token.LinenoEm: f"{GRUVBOX_PROMPT_NUM} bold",
        Token.ValEm: f"{GRUVBOX_VAL_EM} bold",
        Token.VName: GRUVBOX_V_NAME,
        Token.Caret: "",
        Token.Filename: GRUVBOX_FILENAME,
        Token.FilenameEm: f"{GRUVBOX_FILENAME} bold",
        Token.ExcName: f"{GRUVBOX_EXCEPTION_NAME} bold",
        Token.Topline: GRUVBOX_TOPLINE,
        Token.Breakpoint.Enabled: GRUVBOX_BREAKPOINT_ENABLED,
        Token.Breakpoint.Disabled: GRUVBOX_BREAKPOINT_DISABLED,
        Token.Prompt: GRUVBOX_PROMPT,
        Token.PromptNum: f"{GRUVBOX_PROMPT_NUM} bold",
        Token.OutPrompt: GRUVBOX_OUT_PROMPT,
        Token.OutPromptNum: f"{GRUVBOX_OUT_PROMPT_NUM} bold",
    },
    symbols={"arrow_body": "\u2500", "arrow_head": "\u25b6", "top_line": "\u2500"},
)

theme_table: dict[str, Theme] = {
    "nocolor": nocolors_theme,
    "linux": linux_theme,
    "neutral": neutral_theme,
    "neutral:nt": neutral_nt,
    "neutral:posix": neutral_posix,
    "lightbg": lightbg_theme,
    "pride": pride_theme,
    "pride:l": pridel_theme,
    "gruvbox-dark": gruvbox_dark_theme,
}


class Parser:
    """Format colored Python source."""

    _theme_name: str
    out: Any
    pos: int
    lines: list[int]
    raw: str

    def __init__(self, out: Any = sys.stdout, *, theme_name: str | None = None) -> None:
        """Create a parser with a specified color table and output channel.

        Call format() to process code.
        """

        assert theme_name is not None

        self.out = out
        self.pos = 0
        self.lines = []
        self.raw = ""
        if theme_name is not None:
            if theme_name in ["Linux", "LightBG", "Neutral", "NoColor"]:
                warnings.warn(
                    f"Theme names and color schemes are lowercase in IPython 9.0 use {theme_name.lower()} instead",
                    DeprecationWarning,
                    stacklevel=2,
                )
                theme_name = theme_name.lower()
        if not theme_name:
            self.theme_name = "nocolor"
        else:
            self.theme_name = theme_name

    @property
    def theme_name(self) -> str:
        return self._theme_name

    @theme_name.setter
    def theme_name(self, value: str) -> None:
        assert value == value.lower()
        self._theme_name = value

    @property
    def style(self) -> str:
        assert False
        return self._theme_name

    @style.setter
    def style(self, val: str) -> None:
        assert False
        assert val == val.lower()
        self._theme_name = val

    def format(self, raw: str, out: Any = None) -> str | None:
        return self.format2(raw, out)[0]

    def format2(self, raw: str, out: Any = None) -> tuple[str | None, bool]:
        """Parse and send the colored source.

        If out is not specified, the defaults (given to constructor) are used.

        out should be a file-type object. Optionally, out can be given as the
        string 'str' and the parser will automatically return the output in a
        string."""

        string_output = 0
        if out == "str" or self.out == "str" or isinstance(self.out, StringIO):
            # XXX - I don't really like this state handling logic, but at this
            # point I don't want to make major changes, so adding the
            # isinstance() check is the simplest I can do to ensure correct
            # behavior.
            out_old = self.out
            self.out = StringIO()
            string_output = 1
        elif out is not None:
            self.out = out
        else:
            raise ValueError(
                '`out` or `self.out` should be file-like or the value `"str"`'
            )

        # Fast return of the unmodified input for nocolor scheme
        # TODO:
        if self.theme_name == "nocolor":
            error = False
            self.out.write(raw)
            if string_output:
                return raw, error
            return None, error

        # local shorthands

        # Remove trailing whitespace and normalize tabs
        self.raw = raw.expandtabs().rstrip()

        # store line offsets in self.lines
        self.lines = [0, 0]
        pos = 0
        raw_find = self.raw.find
        lines_append = self.lines.append
        while True:
            pos = raw_find("\n", pos) + 1
            if not pos:
                break
            lines_append(pos)
        lines_append(len(self.raw))

        # parse the source and write it
        self.pos = 0
        text = StringIO(self.raw)

        error = False
        try:
            for atoken in generate_tokens(text.readline):
                self(*atoken)
        except tokenize.TokenError as ex:
            msg = ex.args[0]
            line = ex.args[1][0]
            self.out.write(
                theme_table[self.theme_name].format(
                    [
                        (Token, "\n\n"),
                        (
                            Token.Error,
                            f"*** ERROR: {msg}{self.raw[self.lines[line] :]}",
                        ),
                        (Token, "\n"),
                    ]
                )
            )
            error = True
        self.out.write(
            theme_table[self.theme_name].format(
                [
                    (Token, "\n"),
                ]
            )
        )

        if string_output:
            output = self.out.getvalue()
            self.out = out_old
            return (output, error)
        return (None, error)

    def _inner_call_(
        self, toktype: int, toktext: str, start_pos: tuple[int, int]
    ) -> str:
        """like call but write to a temporary buffer"""
        srow, scol = start_pos

        # calculate new positions
        oldpos = self.pos
        newpos = self.lines[srow] + scol
        self.pos = newpos + len(toktext)

        # send the original whitespace, if needed
        if newpos > oldpos:
            acc = self.raw[oldpos:newpos]
        else:
            acc = ""

        # skip indenting tokens
        if toktype in [token.INDENT, token.DEDENT]:
            self.pos = newpos
            return acc

        # map token type to a color group
        if token.LPAR <= toktype <= token.OP:
            toktype = token.OP
        elif toktype == token.NAME and keyword.iskeyword(toktext):
            toktype = _KEYWORD
        pyg_tok_type = _pygment_token_mapping.get(toktype, Token.Text)

        # send text, pygments should take care of splitting on newline and resending
        # the correct self.colors after the new line, which is necessary for pagers
        acc += theme_table[self.theme_name].format([(pyg_tok_type, toktext)])
        return acc

    def __call__(
        self,
        toktype: int,
        toktext: str,
        start_pos: tuple[int, int],
        end_pos: tuple[int, int],
        line: str,
    ) -> None:
        """Token handler, with syntax highlighting."""
        self.out.write(self._inner_call_(toktype, toktext, start_pos))
