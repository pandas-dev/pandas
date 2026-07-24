from typing import Literal
from typing_extensions import TypeAlias

from docutils import parsers

_ParserName: TypeAlias = Literal["pycmark", "myst", "recommonmark"]

commonmark_parser_names: tuple[_ParserName, ...]
Parser: type[parsers.Parser]  # if Parser is None or parser_name is empty string, user cannot import current module
parser_name: _ParserName
