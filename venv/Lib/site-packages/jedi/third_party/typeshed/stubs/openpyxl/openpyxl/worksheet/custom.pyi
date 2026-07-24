from _typeshed import Incomplete
from typing import ClassVar, Literal

from openpyxl.descriptors.base import String
from openpyxl.descriptors.serialisable import Serialisable

class CustomProperty(Serialisable):
    tagname: ClassVar[str]
    name: String[Literal[False]]
    def __init__(self, name: str) -> None: ...

class CustomProperties(Serialisable):
    tagname: ClassVar[str]
    customPr: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, customPr=()) -> None: ...
