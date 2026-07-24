from typing import ClassVar, Final

from docutils import writers

__docformat__: Final = "reStructuredText"

class Writer(writers.UnfilteredWriter[str]):
    supported: ClassVar[tuple[str, ...]]
    config_section: ClassVar[str]
    config_section_dependencies: ClassVar[tuple[str]]
    def translate(self) -> None: ...
