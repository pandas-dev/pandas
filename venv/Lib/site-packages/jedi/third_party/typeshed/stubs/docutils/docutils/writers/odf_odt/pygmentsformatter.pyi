from _typeshed import Incomplete
from typing import Any

# Formatter[str] from types-pygments
class _Formatter:
    name: Any
    aliases: Any
    filenames: Any
    unicodeoutput: bool
    style: Any
    full: Any
    title: Any
    encoding: Any
    options: Any
    def __init__(self, *, encoding: None = None, outencoding: None = None, **options) -> None: ...
    def get_style_defs(self, arg: str = ""): ...
    def format(self, tokensource, outfile): ...

class OdtPygmentsFormatter(_Formatter):
    rststyle_function: Incomplete
    escape_function: Incomplete
    def __init__(self, rststyle_function, escape_function) -> None: ...
    def rststyle(self, name, parameters=()): ...
    def get_style_defs(self, arg: str = ""): ...
    def format(self, tokensource, outfile): ...

class OdtPygmentsProgFormatter(OdtPygmentsFormatter):
    def format(self, tokensource, outfile) -> None: ...

class OdtPygmentsLaTeXFormatter(OdtPygmentsFormatter):
    def format(self, tokensource, outfile) -> None: ...
