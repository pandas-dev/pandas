from pyasn1.type import base
from pyasn1.type.tag import TagSet

__all__ = ["endOfOctets"]

class EndOfOctets(base.SimpleAsn1Type):
    defaultValue: int
    tagSet: TagSet
    def __new__(cls, *args, **kwargs): ...

endOfOctets: EndOfOctets
