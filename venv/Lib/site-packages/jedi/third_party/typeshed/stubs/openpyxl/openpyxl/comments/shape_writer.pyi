from _typeshed import Incomplete
from typing import Any
from typing_extensions import TypeAlias

from ..xml._functions_overloads import _lxml_Element, _ParentElement

_RootElement: TypeAlias = _ParentElement[Any] | _lxml_Element

vmlns: str
officens: str
excelns: str

class ShapeWriter:
    vml: Incomplete
    vml_path: Incomplete
    comments: Incomplete
    def __init__(self, comments) -> None: ...
    def add_comment_shapetype(self, root: _RootElement) -> None: ...
    def add_comment_shape(self, root: _RootElement, idx, coord, height, width) -> None: ...
    # Any object missing "findall" is replaced by an Element
    def write(self, root: _RootElement | None) -> str: ...
