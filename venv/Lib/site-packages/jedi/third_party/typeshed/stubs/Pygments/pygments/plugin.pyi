import sys
from _typeshed import Incomplete
from collections.abc import Generator
from typing import Final

from pygments.filter import Filter
from pygments.formatter import Formatter
from pygments.lexer import Lexer
from pygments.style import Style

LEXER_ENTRY_POINT: Final = "pygments.lexers"
FORMATTER_ENTRY_POINT: Final = "pygments.formatters"
STYLE_ENTRY_POINT: Final = "pygments.styles"
FILTER_ENTRY_POINT: Final = "pygments.filters"

if sys.version_info >= (3, 10):
    from importlib.metadata import EntryPoints
    def iter_entry_points(group_name: str) -> EntryPoints: ...

else:
    from importlib.metadata import EntryPoint

    def iter_entry_points(group_name: str) -> tuple[EntryPoint, ...] | list[EntryPoint]: ...

def find_plugin_lexers() -> Generator[type[Lexer]]: ...
def find_plugin_formatters() -> Generator[tuple[str, type[Formatter[Incomplete]]]]: ...
def find_plugin_styles() -> Generator[tuple[str, type[Style]]]: ...
def find_plugin_filters() -> Generator[tuple[str, type[Filter]]]: ...
