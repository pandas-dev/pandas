from _typeshed import Unused
from abc import abstractmethod

from pyasn1.type.base import Asn1Type
from pyasn1.type.tag import TagSet

__all__ = ["Encoder", "encode"]

class AbstractItemEncoder:
    supportIndefLenMode: bool
    eooIntegerSubstrate: tuple[int, int]
    eooOctetsSubstrate: bytes
    def encodeTag(self, singleTag, isConstructed): ...
    def encodeLength(self, length, defMode): ...
    @abstractmethod
    def encodeValue(self, value, asn1Spec, encodeFun, **options) -> None: ...
    def encode(self, value, asn1Spec: Asn1Type | None = None, encodeFun=None, **options): ...

class EndOfOctetsEncoder(AbstractItemEncoder):
    def encodeValue(self, value, asn1Spec, encodeFun, **options): ...

class BooleanEncoder(AbstractItemEncoder):
    supportIndefLenMode: bool
    def encodeValue(self, value, asn1Spec, encodeFun, **options): ...

class IntegerEncoder(AbstractItemEncoder):
    supportIndefLenMode: bool
    supportCompactZero: bool
    def encodeValue(self, value, asn1Spec, encodeFun, **options): ...

class BitStringEncoder(AbstractItemEncoder):
    def encodeValue(self, value, asn1Spec, encodeFun, **options): ...

class OctetStringEncoder(AbstractItemEncoder):
    def encodeValue(self, value, asn1Spec, encodeFun, **options): ...

class NullEncoder(AbstractItemEncoder):
    supportIndefLenMode: bool
    def encodeValue(self, value, asn1Spec, encodeFun, **options): ...

class ObjectIdentifierEncoder(AbstractItemEncoder):
    supportIndefLenMode: bool
    def encodeValue(self, value, asn1Spec, encodeFun, **options): ...

class RealEncoder(AbstractItemEncoder):
    supportIndefLenMode: bool
    binEncBase: int
    def encodeValue(self, value, asn1Spec, encodeFun, **options): ...

class SequenceEncoder(AbstractItemEncoder):
    omitEmptyOptionals: bool
    def encodeValue(self, value, asn1Spec, encodeFun, **options): ...

class SequenceOfEncoder(AbstractItemEncoder):
    def encodeValue(self, value, asn1Spec, encodeFun, **options): ...

class ChoiceEncoder(AbstractItemEncoder):
    def encodeValue(self, value, asn1Spec, encodeFun, **options): ...

class AnyEncoder(OctetStringEncoder):
    def encodeValue(self, value, asn1Spec, encodeFun, **options): ...

TAG_MAP: dict[TagSet, AbstractItemEncoder]
TYPE_MAP: dict[int, AbstractItemEncoder]
# deprecated aliases
tagMap = TAG_MAP
typeMap = TYPE_MAP

class SingleItemEncoder:
    fixedDefLengthMode: bool | None
    fixedChunkSize: int | None
    TAG_MAP: dict[TagSet, AbstractItemEncoder]
    TYPE_MAP: dict[int, AbstractItemEncoder]

    def __init__(self, tagMap=..., typeMap=..., **ignored: Unused) -> None: ...
    def __call__(self, value, asn1Spec: Asn1Type | None = None, **options): ...

class Encoder:
    SINGLE_ITEM_ENCODER: type[SingleItemEncoder]

    def __init__(self, tagMap=..., typeMap=..., **options: Unused) -> None: ...
    def __call__(self, pyObject, asn1Spec: Asn1Type | None = None, **options): ...

encode: Encoder
