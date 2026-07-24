from typing import ClassVar, Final, TypeVar

from docutils import nodes, readers

__docformat__: Final = "reStructuredText"

_S = TypeVar("_S", bound=str | bytes)

class Reader(readers.Reader[_S]):
    document: nodes.document | None  # type: ignore[assignment]
    config_section_dependencies: ClassVar[tuple[str, ...]]
