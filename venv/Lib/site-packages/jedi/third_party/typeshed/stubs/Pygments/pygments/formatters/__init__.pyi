from _typeshed import Incomplete
from collections.abc import Generator

from ..formatter import Formatter
from .bbcode import BBCodeFormatter as BBCodeFormatter
from .groff import GroffFormatter as GroffFormatter
from .html import HtmlFormatter as HtmlFormatter
from .img import (
    BmpImageFormatter as BmpImageFormatter,
    GifImageFormatter as GifImageFormatter,
    ImageFormatter as ImageFormatter,
    JpgImageFormatter as JpgImageFormatter,
)
from .irc import IRCFormatter as IRCFormatter
from .latex import LatexFormatter as LatexFormatter
from .other import NullFormatter as NullFormatter, RawTokenFormatter as RawTokenFormatter, TestcaseFormatter as TestcaseFormatter
from .pangomarkup import PangoMarkupFormatter as PangoMarkupFormatter
from .rtf import RtfFormatter as RtfFormatter
from .svg import SvgFormatter as SvgFormatter
from .terminal import TerminalFormatter as TerminalFormatter
from .terminal256 import Terminal256Formatter as Terminal256Formatter, TerminalTrueColorFormatter as TerminalTrueColorFormatter

__all__ = [
    "get_formatter_by_name",
    "get_formatter_for_filename",
    "get_all_formatters",
    "load_formatter_from_file",
    "BBCodeFormatter",
    "BmpImageFormatter",
    "GifImageFormatter",
    "GroffFormatter",
    "HtmlFormatter",
    "IRCFormatter",
    "ImageFormatter",
    "JpgImageFormatter",
    "LatexFormatter",
    "NullFormatter",
    "PangoMarkupFormatter",
    "RawTokenFormatter",
    "RtfFormatter",
    "SvgFormatter",
    "Terminal256Formatter",
    "TerminalFormatter",
    "TerminalTrueColorFormatter",
    "TestcaseFormatter",
]

def get_all_formatters() -> Generator[type[Formatter[Incomplete]]]: ...
def get_formatter_by_name(_alias, **options): ...
def load_formatter_from_file(filename, formattername: str = "CustomFormatter", **options): ...
def get_formatter_for_filename(fn, **options): ...
