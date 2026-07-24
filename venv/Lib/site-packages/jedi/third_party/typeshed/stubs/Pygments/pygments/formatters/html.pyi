from _typeshed import Incomplete, SupportsWrite
from collections.abc import Iterable
from typing import TypeVar

from pygments.formatter import Formatter
from pygments.token import _TokenType

_T = TypeVar("_T", str, bytes)

__all__ = ["HtmlFormatter"]

class HtmlFormatter(Formatter[_T]):
    title: Incomplete
    nowrap: Incomplete
    noclasses: Incomplete
    classprefix: Incomplete
    cssclass: Incomplete
    cssstyles: Incomplete
    prestyles: Incomplete
    cssfile: Incomplete
    noclobber_cssfile: Incomplete
    tagsfile: Incomplete
    tagurlformat: Incomplete
    filename: Incomplete
    wrapcode: Incomplete
    span_element_openers: Incomplete
    linenos: int
    linenostart: Incomplete
    linenostep: Incomplete
    linenospecial: Incomplete
    nobackground: Incomplete
    lineseparator: Incomplete
    lineanchors: Incomplete
    linespans: Incomplete
    anchorlinenos: Incomplete
    hl_lines: Incomplete
    def get_style_defs(self, arg=None): ...
    def get_token_style_defs(self, arg=None): ...
    def get_background_style_defs(self, arg=None): ...
    def get_linenos_style_defs(self): ...
    def get_css_prefix(self, arg): ...
    def wrap(self, source): ...
    def format_unencoded(self, tokensource: Iterable[tuple[_TokenType, str]], outfile: SupportsWrite[str]) -> None: ...
