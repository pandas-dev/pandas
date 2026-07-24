from _typeshed import Unused
from abc import abstractmethod
from collections import OrderedDict

from pyasn1.type.tag import TagSet

__all__ = ["encode"]

class AbstractItemEncoder:
    @abstractmethod
    def encode(self, value, encodeFun, **options) -> None: ...

class BooleanEncoder(AbstractItemEncoder):
    def encode(self, value, encodeFun, **options): ...

class IntegerEncoder(AbstractItemEncoder):
    def encode(self, value, encodeFun, **options): ...

class BitStringEncoder(AbstractItemEncoder):
    def encode(self, value, encodeFun, **options): ...

class OctetStringEncoder(AbstractItemEncoder):
    def encode(self, value, encodeFun, **options): ...

class TextStringEncoder(AbstractItemEncoder):
    def encode(self, value, encodeFun, **options): ...

class NullEncoder(AbstractItemEncoder):
    def encode(self, value, encodeFun, **options) -> None: ...

class ObjectIdentifierEncoder(AbstractItemEncoder):
    def encode(self, value, encodeFun, **options): ...

class RealEncoder(AbstractItemEncoder):
    def encode(self, value, encodeFun, **options): ...

class SetEncoder(AbstractItemEncoder):
    protoDict = dict
    def encode(self, value, encodeFun, **options): ...

class SequenceEncoder(SetEncoder):
    protoDict = OrderedDict

class SequenceOfEncoder(AbstractItemEncoder):
    def encode(self, value, encodeFun, **options): ...

class ChoiceEncoder(SequenceEncoder): ...

class AnyEncoder(AbstractItemEncoder):
    def encode(self, value, encodeFun, **options): ...

TAG_MAP: dict[TagSet, AbstractItemEncoder]
TYPE_MAP: dict[int, AbstractItemEncoder]
# deprecated aliases
tagMap = TAG_MAP
typeMap = TYPE_MAP

class SingleItemEncoder:
    TAG_MAP: dict[TagSet, AbstractItemEncoder]
    TYPE_MAP: dict[int, AbstractItemEncoder]

    def __init__(self, tagMap=..., typeMap=..., **ignored: Unused) -> None: ...
    def __call__(self, value, **options): ...

class Encoder:
    SINGLE_ITEM_ENCODER: type[SingleItemEncoder]

    def __init__(self, *, tagMap=..., typeMap=..., **options: Unused): ...
    def __call__(self, pyObject, asn1Spec=None, **options): ...

encode: SingleItemEncoder
