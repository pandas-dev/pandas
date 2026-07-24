from _typeshed import Incomplete, Unused
from collections.abc import Iterable, Sequence
from typing import ClassVar, Final, NoReturn

from docutils import nodes
from docutils.transforms import Transform

__docformat__: Final = "reStructuredText"

class SectNum(Transform):
    default_priority: ClassVar[int]
    maxdepth: int
    startvalue: int
    prefix: str
    suffix: str
    def apply(self) -> None: ...
    def update_section_numbers(self, node: nodes.Element, prefix: Iterable[str] = (), depth: int = 0) -> None: ...

class Contents(Transform):
    default_priority: ClassVar[int]
    toc_id: Incomplete
    backlinks: Incomplete
    def apply(self) -> None: ...
    def build_contents(
        self, node: nodes.Element, level: int = 0
    ) -> nodes.bullet_list | list[None]: ...  # return empty list if entries is empty
    def copy_and_filter(self, node: nodes.Node) -> Sequence[nodes.Node]: ...

class ContentsFilter(nodes.TreeCopyVisitor):
    def get_entry_text(self) -> Sequence[nodes.Node]: ...
    def ignore_node_but_process_children(self, node: Unused) -> NoReturn: ...
    visit_problematic = ignore_node_but_process_children
    visit_reference = ignore_node_but_process_children
    visit_target = ignore_node_but_process_children
