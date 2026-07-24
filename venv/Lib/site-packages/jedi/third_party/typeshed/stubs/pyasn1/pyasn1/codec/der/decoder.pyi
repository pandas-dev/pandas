from pyasn1.codec.ber.decoder import AbstractPayloadDecoder
from pyasn1.codec.cer import decoder
from pyasn1.type.tag import TagSet

__all__ = ["decode", "StreamingDecoder"]

class BitStringPayloadDecoder(decoder.BitStringPayloadDecoder):
    supportConstructedForm: bool

class OctetStringPayloadDecoder(decoder.OctetStringPayloadDecoder):
    supportConstructedForm: bool

RealPayloadDecoder = decoder.RealPayloadDecoder

TAG_MAP: dict[TagSet, AbstractPayloadDecoder]
TYPE_MAP: dict[int, AbstractPayloadDecoder]
# deprecated aliases
tagMap = TAG_MAP
typeMap = TYPE_MAP

class SingleItemDecoder(decoder.SingleItemDecoder):
    TAG_MAP: dict[TagSet, AbstractPayloadDecoder]
    TYPE_MAP: dict[int, AbstractPayloadDecoder]

    supportIndefLength: bool

class StreamingDecoder(decoder.StreamingDecoder):
    SINGLE_ITEM_DECODER: type[SingleItemDecoder]

class Decoder(decoder.Decoder):
    STREAMING_DECODER: type[StreamingDecoder]

decode: Decoder
