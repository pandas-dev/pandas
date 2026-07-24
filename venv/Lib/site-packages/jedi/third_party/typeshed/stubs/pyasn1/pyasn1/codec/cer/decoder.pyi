from _typeshed import Unused

from pyasn1.codec.ber import decoder
from pyasn1.type import univ
from pyasn1.type.tag import TagSet

__all__ = ["decode", "StreamingDecoder"]

class BooleanPayloadDecoder(decoder.AbstractSimplePayloadDecoder):
    protoComponent: univ.Boolean
    def valueDecoder(
        self,
        substrate,
        asn1Spec,
        tagSet: TagSet | None = None,
        length: int | None = None,
        state: Unused = None,
        decodeFun: Unused = None,
        substrateFun: Unused = None,
        **options,
    ): ...

BitStringPayloadDecoder = decoder.BitStringPayloadDecoder
OctetStringPayloadDecoder = decoder.OctetStringPayloadDecoder
RealPayloadDecoder = decoder.RealPayloadDecoder

TAG_MAP: dict[TagSet, decoder.AbstractPayloadDecoder]
TYPE_MAP: dict[int, decoder.AbstractPayloadDecoder]
# deprecated aliases
tagMap = TAG_MAP
typeMap = TYPE_MAP

class SingleItemDecoder(decoder.SingleItemDecoder):
    TAG_MAP: dict[TagSet, decoder.AbstractPayloadDecoder]
    TYPE_MAP: dict[int, decoder.AbstractPayloadDecoder]

class StreamingDecoder(decoder.StreamingDecoder):
    SINGLE_ITEM_DECODER: type[SingleItemDecoder]

class Decoder(decoder.Decoder):
    STREAMING_DECODER: type[StreamingDecoder]

decode: Decoder
