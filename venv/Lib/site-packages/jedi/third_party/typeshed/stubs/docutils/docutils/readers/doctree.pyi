from typing import ClassVar, Final, TypeVar

from docutils import readers

_S = TypeVar("_S", bound=str | bytes)

__docformat__: Final = "reStructuredText"

class Reader(readers.ReReader[_S]):
    config_section_dependencies: ClassVar[tuple[str, ...]]
