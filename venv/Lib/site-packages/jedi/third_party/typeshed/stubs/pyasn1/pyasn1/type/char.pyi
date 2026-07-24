from pyasn1.type import univ
from pyasn1.type.tag import TagSet

__all__ = [
    "NumericString",
    "PrintableString",
    "TeletexString",
    "T61String",
    "VideotexString",
    "IA5String",
    "GraphicString",
    "VisibleString",
    "ISO646String",
    "GeneralString",
    "UniversalString",
    "BMPString",
    "UTF8String",
]

class AbstractCharacterString(univ.OctetString):
    def __bytes__(self) -> bytes: ...
    def prettyIn(self, value): ...
    def asOctets(self, padding: bool = True): ...
    def asNumbers(self, padding: bool = True): ...
    def prettyOut(self, value): ...
    def prettyPrint(self, scope: int = 0): ...
    def __reversed__(self): ...

class NumericString(AbstractCharacterString):
    tagSet: TagSet
    encoding: str
    typeId: int

class PrintableString(AbstractCharacterString):
    tagSet: TagSet
    encoding: str
    typeId: int

class TeletexString(AbstractCharacterString):
    tagSet: TagSet
    encoding: str
    typeId: int

class T61String(TeletexString):
    typeId: int

class VideotexString(AbstractCharacterString):
    tagSet: TagSet
    encoding: str
    typeId: int

class IA5String(AbstractCharacterString):
    tagSet: TagSet
    encoding: str
    typeId: int

class GraphicString(AbstractCharacterString):
    tagSet: TagSet
    encoding: str
    typeId: int

class VisibleString(AbstractCharacterString):
    tagSet: TagSet
    encoding: str
    typeId: int

class ISO646String(VisibleString):
    typeId: int

class GeneralString(AbstractCharacterString):
    tagSet: TagSet
    encoding: str
    typeId: int

class UniversalString(AbstractCharacterString):
    tagSet: TagSet
    encoding: str
    typeId: int

class BMPString(AbstractCharacterString):
    tagSet: TagSet
    encoding: str
    typeId: int

class UTF8String(AbstractCharacterString):
    tagSet: TagSet
    encoding: str
    typeId: int
