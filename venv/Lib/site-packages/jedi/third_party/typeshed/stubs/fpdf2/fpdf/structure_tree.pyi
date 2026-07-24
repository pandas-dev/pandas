from _typeshed import Incomplete, Unused
from collections import defaultdict
from collections.abc import Iterable, Iterator

from .encryption import StandardSecurityHandler
from .syntax import PDFArray, PDFObject, PDFString

class NumberTree(PDFObject):
    __slots__ = ("_id", "nums")
    nums: defaultdict[Incomplete, list[Incomplete]]
    def __init__(self) -> None: ...
    def serialize(self, obj_dict: Unused = None, _security_handler: StandardSecurityHandler | None = None) -> str: ...

class StructTreeRoot(PDFObject):
    __slots__ = ("_id", "type", "parent_tree", "k")
    type: str
    parent_tree: NumberTree
    k: PDFArray[Incomplete]
    def __init__(self) -> None: ...

class StructElem(PDFObject):
    __slots__ = ("_id", "type", "s", "p", "k", "t", "alt", "pg", "_page_number")
    type: str
    s: str
    p: PDFObject
    k: PDFArray[Incomplete]
    t: PDFString | None
    alt: PDFString | None
    pg: Incomplete | None
    def __init__(
        self,
        struct_type: str,
        parent: PDFObject,
        kids: Iterable[int] | Iterable[StructElem],
        page_number: int | None = None,
        title: str | None = None,
        alt: str | None = None,
    ) -> None: ...
    def page_number(self) -> int | None: ...

class StructureTreeBuilder:
    struct_tree_root: Incomplete
    doc_struct_elem: Incomplete
    struct_elem_per_mc: Incomplete
    def __init__(self) -> None: ...
    def add_marked_content(
        self, page_number: int, struct_type: str, mcid: int | None = None, title: str | None = None, alt_text: str | None = None
    ) -> tuple[Incomplete, Incomplete]: ...
    def next_mcid_for_page(self, page_number: int) -> int: ...
    def empty(self) -> bool: ...
    def __iter__(self) -> Iterator[Incomplete]: ...
