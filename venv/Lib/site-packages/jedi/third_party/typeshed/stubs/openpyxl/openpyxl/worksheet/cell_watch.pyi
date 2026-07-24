from _typeshed import Incomplete
from typing import ClassVar, Literal

from openpyxl.descriptors.base import String
from openpyxl.descriptors.serialisable import Serialisable

class CellWatch(Serialisable):
    tagname: ClassVar[str]
    r: String[Literal[True]]
    def __init__(self, r: str) -> None: ...

class CellWatches(Serialisable):
    tagname: ClassVar[str]
    cellWatch: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, cellWatch=()) -> None: ...
