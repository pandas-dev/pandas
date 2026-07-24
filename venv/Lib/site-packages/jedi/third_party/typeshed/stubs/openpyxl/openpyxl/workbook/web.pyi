from _typeshed import ConvertibleToInt, Incomplete, Unused
from typing import ClassVar, Literal, overload
from typing_extensions import TypeAlias

from openpyxl.descriptors.base import Bool, Integer, NoneSet, String, _ConvertibleToBool
from openpyxl.descriptors.serialisable import Serialisable

_WebPublishingTargetScreenSize: TypeAlias = Literal[
    "544x376",
    "640x480",
    "720x512",
    "800x600",
    "1024x768",
    "1152x882",
    "1152x900",
    "1280x1024",
    "1600x1200",
    "1800x1440",
    "1920x1200",
]

class WebPublishObject(Serialisable):
    tagname: ClassVar[str]
    id: Integer[Literal[False]]
    divId: String[Literal[False]]
    sourceObject: String[Literal[True]]
    destinationFile: String[Literal[False]]
    title: String[Literal[True]]
    autoRepublish: Bool[Literal[True]]
    @overload
    def __init__(
        self,
        id: ConvertibleToInt,
        divId: str,
        sourceObject: str | None = None,
        *,
        destinationFile: str,
        title: str | None = None,
        autoRepublish: _ConvertibleToBool | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        id: ConvertibleToInt,
        divId: str,
        sourceObject: str | None,
        destinationFile: str,
        title: str | None = None,
        autoRepublish: _ConvertibleToBool | None = None,
    ) -> None: ...

class WebPublishObjectList(Serialisable):
    tagname: ClassVar[str]
    # Overwritten by property below
    # count: Integer
    webPublishObject: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, count: Unused = None, webPublishObject=()) -> None: ...
    @property
    def count(self) -> int: ...

class WebPublishing(Serialisable):
    tagname: ClassVar[str]
    css: Bool[Literal[True]]
    thicket: Bool[Literal[True]]
    longFileNames: Bool[Literal[True]]
    vml: Bool[Literal[True]]
    allowPng: Bool[Literal[True]]
    targetScreenSize: NoneSet[_WebPublishingTargetScreenSize]
    dpi: Integer[Literal[True]]
    codePage: Integer[Literal[True]]
    characterSet: String[Literal[True]]
    def __init__(
        self,
        css: _ConvertibleToBool | None = None,
        thicket: _ConvertibleToBool | None = None,
        longFileNames: _ConvertibleToBool | None = None,
        vml: _ConvertibleToBool | None = None,
        allowPng: _ConvertibleToBool | None = None,
        targetScreenSize: _WebPublishingTargetScreenSize | Literal["none"] | None = "800x600",
        dpi: ConvertibleToInt | None = None,
        codePage: ConvertibleToInt | None = None,
        characterSet: str | None = None,
    ) -> None: ...
