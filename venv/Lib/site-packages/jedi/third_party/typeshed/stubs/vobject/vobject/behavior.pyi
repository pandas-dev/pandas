from _typeshed import Incomplete
from typing import ClassVar

class Behavior:
    name: str
    description: str
    versionString: str
    knownChildren: ClassVar[dict[str, tuple[int, int | None, int | None]]]
    quotedPrintable: bool
    defaultBehavior: Incomplete
    hasNative: bool
    isComponent: bool
    allowGroup: bool
    forceUTC: bool
    sortFirst: Incomplete
    @classmethod
    def validate(cls, obj, raiseException: bool = False, complainUnrecognized: bool = False) -> bool: ...
    @classmethod
    def lineValidate(cls, line, raiseException, complainUnrecognized): ...
    @classmethod
    def decode(cls, line) -> None: ...
    @classmethod
    def encode(cls, line) -> None: ...
    @classmethod
    def transformToNative(cls, obj): ...
    @classmethod
    def transformFromNative(cls, obj) -> None: ...
    @classmethod
    def generateImplicitParameters(cls, obj) -> None: ...
    @classmethod
    def serialize(cls, obj, buf, lineLength, validate: bool = True, *args, **kwargs): ...
    @classmethod
    def valueRepr(cls, line): ...
