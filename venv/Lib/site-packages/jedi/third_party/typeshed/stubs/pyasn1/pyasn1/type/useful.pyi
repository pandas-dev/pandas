import datetime

from pyasn1.type import char
from pyasn1.type.tag import TagSet

__all__ = ["ObjectDescriptor", "GeneralizedTime", "UTCTime"]

class ObjectDescriptor(char.GraphicString):
    tagSet: TagSet
    typeId: int

class TimeMixIn:
    class FixedOffset(datetime.tzinfo):
        def __init__(self, offset: int = 0, name: str = "UTC") -> None: ...
        def utcoffset(self, dt): ...
        def tzname(self, dt): ...
        def dst(self, dt): ...

    UTC: FixedOffset
    @property
    def asDateTime(self): ...
    @classmethod
    def fromDateTime(cls, dt): ...

class GeneralizedTime(char.VisibleString, TimeMixIn):
    tagSet: TagSet
    typeId: int

class UTCTime(char.VisibleString, TimeMixIn):
    tagSet: TagSet
    typeId: int
