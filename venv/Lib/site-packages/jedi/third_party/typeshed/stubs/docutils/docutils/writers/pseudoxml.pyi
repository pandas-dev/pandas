from typing import ClassVar, Final

from docutils import writers

__docformat__: Final = "reStructuredText"

class Writer(writers.Writer[str]):
    config_section: ClassVar[str]
    config_section_dependencies: ClassVar[tuple[str, ...]]
