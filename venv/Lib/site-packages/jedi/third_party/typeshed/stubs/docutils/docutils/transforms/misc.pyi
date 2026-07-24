from typing import ClassVar, Final

from docutils import nodes
from docutils.transforms import Transform

__docformat__: Final = "reStructuredText"

class CallBack(Transform):
    default_priority: ClassVar[int]
    def apply(self) -> None: ...

class ClassAttribute(Transform):
    default_priority: ClassVar[int]
    def apply(self) -> None: ...

class Transitions(Transform):
    default_priority: ClassVar[int]
    def apply(self) -> None: ...
    def visit_transition(self, node: nodes.transition) -> None: ...
