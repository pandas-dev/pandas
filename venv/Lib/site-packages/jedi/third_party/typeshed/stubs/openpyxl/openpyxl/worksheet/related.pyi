from _typeshed import Incomplete, Unused

from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.xml.functions import Element

class Related(Serialisable):
    id: Incomplete
    def __init__(self, id=None) -> None: ...
    def to_tree(self, tagname: str | None, idx: Unused = None) -> Element: ...  # type: ignore[override]
