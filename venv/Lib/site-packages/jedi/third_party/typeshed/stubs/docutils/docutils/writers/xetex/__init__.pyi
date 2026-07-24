from typing import ClassVar, Final

from docutils import nodes
from docutils.utils import Reporter
from docutils.writers import latex2e

__docformat__: Final = "reStructuredText"

class Writer(latex2e.Writer):
    default_template: ClassVar[str]
    default_preamble: ClassVar[str]
    config_section: ClassVar[str]
    config_section_dependencies: ClassVar[tuple[str, ...]]
    translator_class: type[XeLaTeXTranslator]

class Babel(latex2e.Babel):
    language_code: str
    reporter: Reporter
    language: str
    warn_msg: str  # type: ignore[misc]
    quote_index: int
    quotes: tuple[str, ...]
    literal_double_quote: str
    key: str
    def __init__(self, language_code: str, reporter: Reporter) -> None: ...

class XeLaTeXTranslator(latex2e.LaTeXTranslator):
    is_xetex: bool  # type: ignore[misc]
    def __init__(self, document: nodes.document) -> None: ...
    def to_latex_length(self, length_str: str, node: nodes.Node | None = None) -> str: ...
