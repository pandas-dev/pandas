from typing import ClassVar

from pyasn1.codec.ber import encoder
from pyasn1.type.tag import TagSet

__all__ = ["Encoder", "encode"]

class BooleanEncoder(encoder.IntegerEncoder):
    def encodeValue(self, value, asn1Spec, encodeFun, **options): ...

class RealEncoder(encoder.RealEncoder): ...

class TimeEncoderMixIn:
    Z_CHAR: ClassVar[int]
    PLUS_CHAR: ClassVar[int]
    MINUS_CHAR: ClassVar[int]
    COMMA_CHAR: ClassVar[int]
    DOT_CHAR: ClassVar[int]
    ZERO_CHAR: ClassVar[int]
    MIN_LENGTH: ClassVar[int]
    MAX_LENGTH: ClassVar[int]
    def encodeValue(self, value, asn1Spec, encodeFun, **options): ...

class GeneralizedTimeEncoder(TimeEncoderMixIn, encoder.OctetStringEncoder): ...
class UTCTimeEncoder(TimeEncoderMixIn, encoder.OctetStringEncoder): ...

class SetOfEncoder(encoder.SequenceOfEncoder):
    def encodeValue(self, value, asn1Spec, encodeFun, **options): ...

class SequenceOfEncoder(encoder.SequenceOfEncoder):
    def encodeValue(self, value, asn1Spec, encodeFun, **options): ...

class SetEncoder(encoder.SequenceEncoder):
    def encodeValue(self, value, asn1Spec, encodeFun, **options): ...

class SequenceEncoder(encoder.SequenceEncoder):
    omitEmptyOptionals: bool

TAG_MAP: dict[TagSet, encoder.AbstractItemEncoder]
TYPE_MAP: dict[int, encoder.AbstractItemEncoder]
# deprecated aliases
tagMap = TAG_MAP
typeMap = TYPE_MAP

class SingleItemEncoder(encoder.SingleItemEncoder):
    fixedDefLengthMode: bool
    fixedChunkSize: int

    TAG_MAP: dict[TagSet, encoder.AbstractItemEncoder]
    TYPE_MAP: dict[int, encoder.AbstractItemEncoder]

class Encoder(encoder.Encoder):
    SINGLE_ITEM_ENCODER: type[SingleItemEncoder]

encode: Encoder
