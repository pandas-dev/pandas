from _typeshed import Unused
from collections.abc import Iterable

from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.workbook.child import _WorkbookChild
from openpyxl.worksheet._read_only import ReadOnlyWorksheet

from .named_styles import NamedStyle
from .proxy import StyleProxy

class StyleDescriptor:
    collection: str
    key: str
    def __init__(self, collection: str, key: str) -> None: ...
    def __set__(self, instance: StyleableObject, value: Serialisable) -> None: ...
    def __get__(self, instance: StyleableObject, cls: Unused) -> StyleProxy: ...

class NumberFormatDescriptor:
    key: str
    collection: str
    def __set__(self, instance: StyleableObject, value: str) -> None: ...
    def __get__(self, instance: StyleableObject, cls: Unused) -> str: ...

class NamedStyleDescriptor:
    key: str
    collection: str
    def __set__(self, instance: StyleableObject, value: NamedStyle | str) -> None: ...
    def __get__(self, instance: StyleableObject, cls: Unused) -> str: ...

class StyleArrayDescriptor:
    key: str
    def __init__(self, key: str) -> None: ...
    def __set__(self, instance: StyleableObject, value: int) -> None: ...
    def __get__(self, instance: StyleableObject, cls: Unused) -> bool: ...

class StyleableObject:
    __slots__ = ("parent", "_style")
    font: StyleDescriptor
    fill: StyleDescriptor
    border: StyleDescriptor
    number_format: NumberFormatDescriptor
    protection: StyleDescriptor
    alignment: StyleDescriptor
    style: NamedStyleDescriptor
    quotePrefix: StyleArrayDescriptor
    pivotButton: StyleArrayDescriptor
    parent: _WorkbookChild | ReadOnlyWorksheet
    def __init__(
        self, sheet: _WorkbookChild | ReadOnlyWorksheet, style_array: bytes | bytearray | Iterable[int] | None = None
    ) -> None: ...
    @property
    def style_id(self) -> int: ...
    @property
    def has_style(self) -> bool: ...
