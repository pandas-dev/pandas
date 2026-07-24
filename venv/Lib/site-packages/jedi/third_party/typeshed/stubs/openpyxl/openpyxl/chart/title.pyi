from _typeshed import Unused
from typing import ClassVar, Literal

from openpyxl.chart.layout import Layout
from openpyxl.chart.shapes import GraphicalProperties
from openpyxl.chart.text import RichText, Text
from openpyxl.descriptors import Strict, Typed
from openpyxl.descriptors.base import Alias, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.nested import NestedBool
from openpyxl.descriptors.serialisable import Serialisable

from ..xml._functions_overloads import _HasTagAndGet

class Title(Serialisable):
    tagname: ClassVar[str]
    tx: Typed[Text, Literal[True]]
    text: Alias
    layout: Typed[Layout, Literal[True]]
    overlay: NestedBool[Literal[True]]
    spPr: Typed[GraphicalProperties, Literal[True]]
    graphicalProperties: Alias
    txPr: Typed[RichText, Literal[True]]
    body: Alias
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        tx: Text | None = None,
        layout: Layout | None = None,
        overlay: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        spPr: GraphicalProperties | None = None,
        txPr: RichText | None = None,
        extLst: Unused = None,
    ) -> None: ...

def title_maker(text) -> Title: ...

class TitleDescriptor(Typed[Title, Literal[True]]):
    expected_type: type[Title]
    allow_none: Literal[True]
    def __set__(self, instance: Serialisable | Strict, value: str | Title | None) -> None: ...
