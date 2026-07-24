from _typeshed import Incomplete
from collections.abc import Generator, Iterable
from typing import ClassVar, Final, Literal

from docutils import nodes
from docutils.transforms import Transform

__docformat__: Final = "reStructuredText"

class Decorations(Transform):
    default_priority: ClassVar[int]
    def apply(self) -> None: ...
    def generate_header(self) -> None: ...
    def generate_footer(self) -> list[nodes.paragraph] | None: ...

class ExposeInternals(Transform):
    default_priority: ClassVar[int]
    def not_Text(self, node: object) -> bool: ...  # node passing to isinstance() method
    def apply(self) -> None: ...

class Messages(Transform):
    default_priority: ClassVar[int]
    def apply(self) -> None: ...

class FilterMessages(Transform):
    default_priority: ClassVar[int]
    def apply(self) -> None: ...

class TestMessages(Transform):
    __test__: bool
    default_priority: ClassVar[int]
    def apply(self) -> None: ...

class StripComments(Transform):
    default_priority: ClassVar[int]
    def apply(self) -> None: ...

class StripClassesAndElements(Transform):
    default_priority: ClassVar[int]
    strip_elements: set[Incomplete]
    def apply(self) -> None: ...
    def check_classes(self, node: object) -> bool: ...

class SmartQuotes(Transform):
    default_priority: ClassVar[int]
    nodes_to_skip: ClassVar[tuple[type[nodes.Node | nodes.Special], ...]]
    literal_nodes: ClassVar[tuple[type[nodes.Node | nodes.Body], ...]]
    smartquotes_action: ClassVar[str]
    unsupported_languages: set[str]
    def __init__(self, document: nodes.document, startnode: nodes.Node | None) -> None: ...
    def get_tokens(self, txtnodes: Iterable[nodes.Node]) -> Generator[tuple[Literal["literal", "plain"], str]]: ...
    def apply(self) -> None: ...

class Validate(Transform):
    default_priority: ClassVar[int]
    def apply(self) -> None: ...
