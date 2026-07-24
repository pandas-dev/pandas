from _typeshed import Incomplete
from typing import ClassVar, Literal
from typing_extensions import TypeAlias

from openpyxl.descriptors.base import Bool, NoneSet, String, _ConvertibleToBool
from openpyxl.descriptors.serialisable import Serialisable

_SmartTagPropertiesShow: TypeAlias = Literal["all", "noIndicator"]

class SmartTag(Serialisable):
    tagname: ClassVar[str]
    namespaceUri: String[Literal[True]]
    name: String[Literal[True]]
    url: String[Literal[True]]
    def __init__(self, namespaceUri: str | None = None, name: str | None = None, url: str | None = None) -> None: ...

class SmartTagList(Serialisable):
    tagname: ClassVar[str]
    smartTagType: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, smartTagType=()) -> None: ...

class SmartTagProperties(Serialisable):
    tagname: ClassVar[str]
    embed: Bool[Literal[True]]
    show: NoneSet[_SmartTagPropertiesShow]
    def __init__(
        self, embed: _ConvertibleToBool | None = None, show: _SmartTagPropertiesShow | Literal["none"] | None = None
    ) -> None: ...
