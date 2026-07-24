from _typeshed import Incomplete
from collections.abc import Generator, Iterable
from dataclasses import dataclass

from .fonts import TextStyle
from .fpdf import FPDF
from .structure_tree import StructElem
from .syntax import Destination, PDFObject, PDFString

@dataclass
class OutlineSection:
    __slots__ = ("name", "level", "page_number", "dest", "struct_elem")
    name: str
    level: int
    page_number: int
    dest: Destination
    struct_elem: StructElem | None = None

class OutlineItemDictionary(PDFObject):
    __slots__ = ("_id", "title", "parent", "prev", "next", "first", "last", "count", "dest", "struct_elem")
    title: PDFString
    parent: Incomplete | None
    prev: Incomplete | None
    next: Incomplete | None
    first: Incomplete | None
    last: Incomplete | None
    count: int
    dest: Destination | None
    struct_elem: StructElem | None
    def __init__(self, title: str, dest: Destination | None = None, struct_elem: StructElem | None = None) -> None: ...

class OutlineDictionary(PDFObject):
    __slots__ = ("_id", "type", "first", "last", "count")
    type: str
    first: Incomplete | None
    last: Incomplete | None
    count: int
    def __init__(self) -> None: ...

def build_outline_objs(
    sections: Iterable[Incomplete],
) -> Generator[Incomplete, None, list[OutlineDictionary | OutlineItemDictionary]]: ...

class TableOfContents:
    text_style: TextStyle
    use_section_title_styles: bool
    level_indent: float
    line_spacing: float
    ignore_pages_before_toc: bool

    def __init__(
        self,
        text_style: TextStyle | None = None,
        use_section_title_styles: bool = False,
        level_indent: float = 7.5,
        line_spacing: float = 1.5,
        ignore_pages_before_toc: bool = True,
    ) -> None: ...
    def get_text_style(self, pdf: FPDF, item: OutlineSection) -> TextStyle: ...
    def render_toc_item(self, pdf: FPDF, item: OutlineSection) -> None: ...
    def render_toc(self, pdf: FPDF, outline: Iterable[OutlineSection]) -> None: ...
