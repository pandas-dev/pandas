from collections.abc import Callable, Container, Iterable, Mapping, Sequence
from typing import Any, Final, Literal, NamedTuple
from typing_extensions import Self, TypeAlias

__all__ = ["tabulate", "tabulate_formats", "simple_separated_format"]

__version__: Final[str]
# These constants are meant to be configurable
# https://github.com/astanin/python-tabulate#text-formatting
PRESERVE_WHITESPACE: bool
MIN_PADDING: int
# https://github.com/astanin/python-tabulate#wide-fullwidth-cjk-symbols
WIDE_CHARS_MODE: bool
SEPARATING_LINE: str

class Line(NamedTuple):
    begin: str
    hline: str
    sep: str
    end: str

class DataRow(NamedTuple):
    begin: str
    sep: str
    end: str

_TableFormatLine: TypeAlias = None | Line | Callable[[list[int], list[str]], str]
_TableFormatRow: TypeAlias = None | DataRow | Callable[[list[Any], list[int], list[str]], str]

class TableFormat(NamedTuple):
    lineabove: _TableFormatLine
    linebelowheader: _TableFormatLine
    linebetweenrows: _TableFormatLine
    linebelow: _TableFormatLine
    headerrow: _TableFormatRow
    datarow: _TableFormatRow
    padding: int
    with_header_hide: Container[str] | None

LATEX_ESCAPE_RULES: Final[dict[str, str]]
tabulate_formats: list[str]
multiline_formats: dict[str, str]

def simple_separated_format(separator: str) -> TableFormat: ...
def tabulate(
    # The key is converted using str().
    tabular_data: Mapping[Any, Iterable[Any]] | Iterable[Iterable[Any]],
    headers: str | dict[str, str] | Sequence[str] = (),
    tablefmt: str | TableFormat = "simple",
    floatfmt: str | Iterable[str] = "g",
    intfmt: str | Iterable[str] = "",
    numalign: str | None = "default",
    stralign: str | None = "default",
    missingval: str | Iterable[str] = "",
    showindex: str | bool | Iterable[Any] = "default",
    disable_numparse: bool | Iterable[int] = False,
    colglobalalign: Literal["right", "center", "decimal", "left"] | None = None,
    colalign: Iterable[str | None] | None = None,
    preserve_whitespace: bool = False,
    maxcolwidths: int | Iterable[int | None] | None = None,
    headersglobalalign: Literal["right", "center", "left"] | None = None,
    headersalign: Iterable[str | None] | None = None,
    rowalign: str | Iterable[str] | None = None,
    maxheadercolwidths: int | Iterable[int] | None = None,
    break_long_words: bool = True,
    break_on_hyphens: bool = True,
) -> str: ...

class JupyterHTMLStr(str):
    @property
    def str(self) -> Self: ...
