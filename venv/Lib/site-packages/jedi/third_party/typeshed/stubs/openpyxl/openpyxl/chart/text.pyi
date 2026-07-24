from _typeshed import Incomplete, Unused
from typing import ClassVar, Literal

from openpyxl.chart.data_source import StrRef
from openpyxl.descriptors.base import Alias, Typed
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.drawing.text import ListStyle, RichTextProperties
from openpyxl.xml.functions import Element

class RichText(Serialisable):
    tagname: ClassVar[str]
    bodyPr: Typed[RichTextProperties, Literal[False]]
    properties: Alias
    lstStyle: Typed[ListStyle, Literal[True]]
    p: Incomplete
    paragraphs: Alias
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, bodyPr: RichTextProperties | None = None, lstStyle: ListStyle | None = None, p=None) -> None: ...

class Text(Serialisable):
    tagname: ClassVar[str]
    strRef: Typed[StrRef, Literal[True]]
    rich: Typed[RichText, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, strRef: StrRef | None = None, rich: RichText | None = None) -> None: ...
    def to_tree(self, tagname: str | None = None, idx: Unused = None, namespace: str | None = None) -> Element: ...
