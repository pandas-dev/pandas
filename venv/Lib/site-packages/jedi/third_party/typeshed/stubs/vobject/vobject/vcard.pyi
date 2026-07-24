from typing import AnyStr

from .base import ContentLine
from .behavior import Behavior

class Name:
    family: str | list[str]
    given: str | list[str]
    additional: str | list[str]
    prefix: str | list[str]
    suffix: str | list[str]
    def __init__(
        self,
        family: str | list[str] = "",
        given: str | list[str] = "",
        additional: str | list[str] = "",
        prefix: str | list[str] = "",
        suffix: str | list[str] = "",
    ) -> None: ...
    @staticmethod
    def toString(val: str | list[str] | tuple[str, ...]) -> str: ...
    def __eq__(self, other: object) -> bool: ...

class Address:
    box: str | list[str]
    extended: str | list[str]
    street: str | list[str]
    city: str | list[str]
    region: str | list[str]
    code: str | list[str]
    country: str | list[str]
    def __init__(
        self,
        street: str | list[str] = "",
        city: str | list[str] = "",
        region: str | list[str] = "",
        code: str | list[str] = "",
        country: str | list[str] = "",
        box: str | list[str] = "",
        extended: str | list[str] = "",
    ) -> None: ...
    @staticmethod
    def toString(val: str | list[str] | tuple[str, ...], join_char: str = "\n") -> str: ...
    lines: tuple[str, ...]
    one_line: tuple[str, ...]
    def __eq__(self, other: object) -> bool: ...

class VCardTextBehavior(Behavior):
    allowGroup: bool
    base64string: str
    @classmethod
    def decode(cls, line: ContentLine) -> None: ...
    @classmethod
    def encode(cls, line: ContentLine) -> None: ...

class VCardBehavior(Behavior):
    allowGroup: bool
    defaultBehavior: type[VCardTextBehavior]

class VCard3_0(VCardBehavior):
    name: str
    description: str
    versionString: str
    isComponent: bool
    sortFirst: tuple[str, ...]
    @classmethod
    def generateImplicitParameters(cls, obj) -> None: ...

class FN(VCardTextBehavior):
    name: str
    description: str

class Label(VCardTextBehavior):
    name: str
    description: str

class GEO(VCardBehavior): ...

wacky_apple_photo_serialize: bool
REALLY_LARGE: float

class Photo(VCardTextBehavior):
    name: str
    description: str
    @classmethod
    def valueRepr(cls, line: ContentLine) -> str: ...
    @classmethod
    def serialize(cls, obj, buf, lineLength, validate, *args, **kwargs) -> None: ...  # type: ignore[override]

def toListOrString(string: str) -> str | list[str]: ...
def splitFields(string: str) -> list[str | list[str]]: ...
def toList(stringOrList: AnyStr | list[AnyStr]) -> list[AnyStr]: ...
def serializeFields(obj, order=None): ...

NAME_ORDER: tuple[str, ...]
ADDRESS_ORDER: tuple[str, ...]

class NameBehavior(VCardBehavior):
    hasNative: bool
    @staticmethod
    def transformToNative(obj): ...
    @staticmethod
    def transformFromNative(obj): ...

class AddressBehavior(VCardBehavior):
    hasNative: bool
    @staticmethod
    def transformToNative(obj): ...
    @staticmethod
    def transformFromNative(obj): ...

class OrgBehavior(VCardBehavior):
    hasNative: bool
    @staticmethod
    def transformToNative(obj): ...
    @staticmethod
    def transformFromNative(obj): ...
