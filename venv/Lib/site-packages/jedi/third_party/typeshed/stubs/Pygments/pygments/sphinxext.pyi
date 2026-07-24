from typing import Any, Final

from docutils.parsers.rst import Directive

MODULEDOC: Final[str]
LEXERDOC: Final[str]
FMTERDOC: Final[str]
FILTERDOC: Final[str]

class PygmentsDoc(Directive):
    filenames: set[str]
    def document_lexers_overview(self) -> str: ...
    def document_lexers(self) -> str: ...
    def document_formatters(self) -> str: ...
    def document_filters(self) -> str: ...

def setup(app: Any) -> None: ...  # Actual type of 'app' is sphinx.application.Sphinx
