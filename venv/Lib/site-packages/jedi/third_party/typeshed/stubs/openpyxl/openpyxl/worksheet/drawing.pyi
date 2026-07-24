from _typeshed import Incomplete
from typing import ClassVar

from openpyxl.descriptors.serialisable import Serialisable

class Drawing(Serialisable):
    tagname: ClassVar[str]
    id: Incomplete
    def __init__(self, id=None) -> None: ...
