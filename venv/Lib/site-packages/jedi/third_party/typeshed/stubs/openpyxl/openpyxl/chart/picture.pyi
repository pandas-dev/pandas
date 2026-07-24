from _typeshed import ConvertibleToFloat
from typing import ClassVar, Literal
from typing_extensions import TypeAlias

from openpyxl.descriptors.base import _ConvertibleToBool
from openpyxl.descriptors.nested import NestedBool, NestedFloat, NestedNoneSet, _NestedNoneSetParam
from openpyxl.descriptors.serialisable import Serialisable

from ..xml._functions_overloads import _HasTagAndGet

_PictureOptionsPictureFormat: TypeAlias = Literal["stretch", "stack", "stackScale"]

class PictureOptions(Serialisable):
    tagname: ClassVar[str]
    applyToFront: NestedBool[Literal[True]]
    applyToSides: NestedBool[Literal[True]]
    applyToEnd: NestedBool[Literal[True]]
    pictureFormat: NestedNoneSet[_PictureOptionsPictureFormat]
    pictureStackUnit: NestedFloat[Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        applyToFront: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        applyToSides: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        applyToEnd: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        pictureFormat: _NestedNoneSetParam[_PictureOptionsPictureFormat] = None,
        pictureStackUnit: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
    ) -> None: ...
