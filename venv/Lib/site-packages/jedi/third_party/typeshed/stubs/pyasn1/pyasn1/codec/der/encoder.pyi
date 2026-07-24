from pyasn1.codec.ber.encoder import AbstractItemEncoder
from pyasn1.codec.cer import encoder
from pyasn1.type.tag import TagSet

__all__ = ["Encoder", "encode"]

class SetEncoder(encoder.SetEncoder): ...

TAG_MAP: dict[TagSet, AbstractItemEncoder]
TYPE_MAP: dict[int, AbstractItemEncoder]
# deprecated aliases
tagMap = TAG_MAP
typeMap = TYPE_MAP

class SingleItemEncoder(encoder.SingleItemEncoder):
    fixedDefLengthMode: bool
    fixedChunkSize: int

    TAG_MAP: dict[TagSet, AbstractItemEncoder]
    TYPE_MAP: dict[int, AbstractItemEncoder]

class Encoder(encoder.Encoder):
    SINGLE_ITEM_ENCODER: type[SingleItemEncoder]

encode: Encoder
