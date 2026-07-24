from typing import ClassVar, Final

from docutils import parsers

__docformat__: Final = "reStructuredText"

class Parser(parsers.Parser):
    supported: ClassVar[tuple[str, ...]]
    config_section_dependencies: ClassVar[tuple[str, ...]]
