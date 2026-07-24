from _typeshed import Incomplete
from collections.abc import Iterator
from typing import ClassVar, Final, Literal
from typing_extensions import TypeAlias

from openpyxl.descriptors.base import Alias, Bool, NoneSet, Typed, _ConvertibleToBool
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.styles.colors import Color, ColorDescriptor

_SideStyle: TypeAlias = Literal[
    "dashDot",
    "dashDotDot",
    "dashed",
    "dotted",
    "double",
    "hair",
    "medium",
    "mediumDashDot",
    "mediumDashDotDot",
    "mediumDashed",
    "slantDashDot",
    "thick",
    "thin",
]

BORDER_NONE: Final = None
BORDER_DASHDOT: Final = "dashDot"
BORDER_DASHDOTDOT: Final = "dashDotDot"
BORDER_DASHED: Final = "dashed"
BORDER_DOTTED: Final = "dotted"
BORDER_DOUBLE: Final = "double"
BORDER_HAIR: Final = "hair"
BORDER_MEDIUM: Final = "medium"
BORDER_MEDIUMDASHDOT: Final = "mediumDashDot"
BORDER_MEDIUMDASHDOTDOT: Final = "mediumDashDotDot"
BORDER_MEDIUMDASHED: Final = "mediumDashed"
BORDER_SLANTDASHDOT: Final = "slantDashDot"
BORDER_THICK: Final = "thick"
BORDER_THIN: Final = "thin"

class Side(Serialisable):
    color: ColorDescriptor[Literal[True]]
    style: NoneSet[_SideStyle]
    border_style: Alias
    def __init__(
        self, style: _SideStyle | Literal["none"] | None = None, color: str | Color | None = None, border_style=None
    ) -> None: ...

class Border(Serialisable):
    tagname: ClassVar[str]
    __elements__: ClassVar[tuple[str, ...]]
    start: Typed[Side, Literal[True]]
    end: Typed[Side, Literal[True]]
    left: Typed[Side, Literal[True]]
    right: Typed[Side, Literal[True]]
    top: Typed[Side, Literal[True]]
    bottom: Typed[Side, Literal[True]]
    diagonal: Typed[Side, Literal[True]]
    vertical: Typed[Side, Literal[True]]
    horizontal: Typed[Side, Literal[True]]
    outline: Bool[Literal[False]]
    diagonalUp: Bool[Literal[False]]
    diagonalDown: Bool[Literal[False]]
    diagonal_direction: Incomplete
    def __init__(
        self,
        left: Side | None = None,
        right: Side | None = None,
        top: Side | None = None,
        bottom: Side | None = None,
        diagonal: Side | None = None,
        diagonal_direction=None,
        vertical: Side | None = None,
        horizontal: Side | None = None,
        diagonalUp: _ConvertibleToBool = False,
        diagonalDown: _ConvertibleToBool = False,
        outline: _ConvertibleToBool = True,
        start: Side | None = None,
        end: Side | None = None,
    ) -> None: ...
    def __iter__(self) -> Iterator[tuple[str, str]]: ...

DEFAULT_BORDER: Final[Border]
