from typing import ClassVar, Literal

from openpyxl.descriptors.base import Bool, String, Typed, _ConvertibleToBool
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.styles.colors import Color

class ChartsheetProperties(Serialisable):
    tagname: ClassVar[str]
    published: Bool[Literal[True]]
    codeName: String[Literal[True]]
    tabColor: Typed[Color, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self, published: _ConvertibleToBool | None = None, codeName: str | None = None, tabColor: Color | None = None
    ) -> None: ...
