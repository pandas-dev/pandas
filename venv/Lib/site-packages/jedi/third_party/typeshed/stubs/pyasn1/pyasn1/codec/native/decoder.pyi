from _typeshed import Incomplete, Unused
from collections.abc import Callable

from pyasn1.type.tag import TagSet

__all__ = ["decode"]

class AbstractScalarPayloadDecoder:
    def __call__(self, pyObject, asn1Spec, decodeFun: Unused = None, **options): ...

class BitStringPayloadDecoder(AbstractScalarPayloadDecoder):
    def __call__(self, pyObject, asn1Spec, decodeFun: Unused = None, **options): ...

class SequenceOrSetPayloadDecoder:
    def __call__(self, pyObject, asn1Spec, decodeFun: Callable[..., Incomplete] | None = None, **options): ...

class SequenceOfOrSetOfPayloadDecoder:
    def __call__(self, pyObject, asn1Spec, decodeFun: Callable[..., Incomplete] | None = None, **options): ...

class ChoicePayloadDecoder:
    def __call__(self, pyObject, asn1Spec, decodeFun: Callable[..., Incomplete] | None = None, **options): ...

TAG_MAP: dict[TagSet, AbstractScalarPayloadDecoder | SequenceOrSetPayloadDecoder | ChoicePayloadDecoder]
TYPE_MAP: dict[int, AbstractScalarPayloadDecoder | SequenceOrSetPayloadDecoder | ChoicePayloadDecoder]
# deprecated aliases
tagMap = TAG_MAP
typeMap = TYPE_MAP

class SingleItemDecoder:
    TAG_MAP: dict[TagSet, AbstractScalarPayloadDecoder | SequenceOrSetPayloadDecoder | ChoicePayloadDecoder]
    TYPE_MAP: dict[int, AbstractScalarPayloadDecoder | SequenceOrSetPayloadDecoder | ChoicePayloadDecoder]

    def __init__(self, tagMap=..., typeMap=..., **ignored: Unused) -> None: ...
    def __call__(self, pyObject, asn1Spec, **options): ...

class Decoder:
    SINGLE_ITEM_DECODER: type[SingleItemDecoder]

    def __init__(self, *, tagMap=..., typeMap=..., **options: Unused) -> None: ...
    def __call__(self, pyObject, asn1Spec=None, **kwargs): ...

decode: Decoder
