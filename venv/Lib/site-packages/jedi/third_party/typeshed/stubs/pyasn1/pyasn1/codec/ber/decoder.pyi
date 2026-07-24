from _typeshed import Incomplete, Unused
from abc import ABCMeta, abstractmethod
from collections.abc import Callable

from pyasn1.type import base, char, univ, useful
from pyasn1.type.base import Asn1Type
from pyasn1.type.tag import TagSet

__all__ = ["StreamingDecoder", "Decoder", "decode"]

class AbstractPayloadDecoder:
    protoComponent: Asn1Type | None
    @abstractmethod
    def valueDecoder(
        self,
        substrate,
        asn1Spec,
        tagSet: TagSet | None = None,
        length: int | None = None,
        state=None,
        decodeFun: Callable[..., Incomplete] | None = None,
        substrateFun: Callable[..., Incomplete] | None = None,
        **options,
    ) -> None: ...
    # Abstract, but implementation is optional
    def indefLenValueDecoder(
        self,
        substrate,
        asn1Spec,
        tagSet: TagSet | None = None,
        length: int | None = None,
        state=None,
        decodeFun: Callable[..., Incomplete] | None = None,
        substrateFun: Callable[..., Incomplete] | None = None,
        **options,
    ) -> None: ...

class AbstractSimplePayloadDecoder(AbstractPayloadDecoder, metaclass=ABCMeta):
    @staticmethod
    def substrateCollector(asn1Object, substrate, length, options): ...

class RawPayloadDecoder(AbstractSimplePayloadDecoder):
    protoComponent: univ.Any
    def valueDecoder(
        self,
        substrate,
        asn1Spec,
        tagSet: TagSet | None = None,
        length: int | None = None,
        state: Unused = None,
        decodeFun: Callable[..., Incomplete] | None = None,
        substrateFun: Callable[..., Incomplete] | None = None,
        **options,
    ): ...
    def indefLenValueDecoder(
        self,
        substrate,
        asn1Spec,
        tagSet: TagSet | None = None,
        length: int | None = None,
        state: Unused = None,
        decodeFun: Callable[..., Incomplete] | None = None,
        substrateFun: Callable[..., Incomplete] | None = None,
        **options,
    ): ...

class IntegerPayloadDecoder(AbstractSimplePayloadDecoder):
    protoComponent: univ.Integer
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

class BooleanPayloadDecoder(IntegerPayloadDecoder):
    protoComponent: univ.Boolean

class BitStringPayloadDecoder(AbstractSimplePayloadDecoder):
    protoComponent: univ.BitString
    supportConstructedForm: bool
    def valueDecoder(
        self,
        substrate,
        asn1Spec,
        tagSet: TagSet | None = None,
        length: int | None = None,
        state: Unused = None,
        decodeFun: Callable[..., Incomplete] | None = None,
        substrateFun: Callable[..., Incomplete] | None = None,
        **options,
    ): ...
    def indefLenValueDecoder(
        self,
        substrate,
        asn1Spec,
        tagSet: TagSet | None = None,
        length: int | None = None,
        state: Unused = None,
        decodeFun: Callable[..., Incomplete] | None = None,
        substrateFun: Callable[..., Incomplete] | None = None,
        **options,
    ): ...

class OctetStringPayloadDecoder(AbstractSimplePayloadDecoder):
    protoComponent: univ.OctetString
    supportConstructedForm: bool
    def valueDecoder(
        self,
        substrate,
        asn1Spec,
        tagSet: TagSet | None = None,
        length: int | None = None,
        state: Unused = None,
        decodeFun: Callable[..., Incomplete] | None = None,
        substrateFun: Callable[..., Incomplete] | None = None,
        **options,
    ): ...
    def indefLenValueDecoder(
        self,
        substrate,
        asn1Spec,
        tagSet: TagSet | None = None,
        length: int | None = None,
        state: Unused = None,
        decodeFun: Callable[..., Incomplete] | None = None,
        substrateFun: Callable[..., Incomplete] | None = None,
        **options,
    ): ...

class NullPayloadDecoder(AbstractSimplePayloadDecoder):
    protoComponent: univ.Null
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

class ObjectIdentifierPayloadDecoder(AbstractSimplePayloadDecoder):
    protoComponent: univ.ObjectIdentifier
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

class RealPayloadDecoder(AbstractSimplePayloadDecoder):
    protoComponent: univ.Real
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

class AbstractConstructedPayloadDecoder(AbstractPayloadDecoder, metaclass=ABCMeta):
    protoComponent: base.ConstructedAsn1Type | None

class ConstructedPayloadDecoderBase(AbstractConstructedPayloadDecoder):
    protoRecordComponent: univ.SequenceAndSetBase | None
    protoSequenceComponent: univ.SequenceOfAndSetOfBase | None
    def valueDecoder(
        self,
        substrate,
        asn1Spec,
        tagSet: TagSet | None = None,
        length: int | None = None,
        state: Unused = None,
        decodeFun: Callable[..., Incomplete] | None = None,
        substrateFun: Callable[..., Incomplete] | None = None,
        **options,
    ): ...
    def indefLenValueDecoder(
        self,
        substrate,
        asn1Spec,
        tagSet: TagSet | None = None,
        length: int | None = None,
        state: Unused = None,
        decodeFun: Callable[..., Incomplete] | None = None,
        substrateFun: Callable[..., Incomplete] | None = None,
        **options,
    ): ...

class SequenceOrSequenceOfPayloadDecoder(ConstructedPayloadDecoderBase):
    protoRecordComponent: univ.Sequence
    protoSequenceComponent: univ.SequenceOf

class SequencePayloadDecoder(SequenceOrSequenceOfPayloadDecoder):
    protoComponent: univ.Sequence

class SequenceOfPayloadDecoder(SequenceOrSequenceOfPayloadDecoder):
    protoComponent: univ.SequenceOf

class SetOrSetOfPayloadDecoder(ConstructedPayloadDecoderBase):
    protoRecordComponent: univ.Set
    protoSequenceComponent: univ.SetOf

class SetPayloadDecoder(SetOrSetOfPayloadDecoder):
    protoComponent: univ.Set

class SetOfPayloadDecoder(SetOrSetOfPayloadDecoder):
    protoComponent: univ.SetOf

class ChoicePayloadDecoder(AbstractConstructedPayloadDecoder):
    protoComponent: univ.Choice
    def valueDecoder(
        self,
        substrate,
        asn1Spec,
        tagSet: TagSet | None = None,
        length: int | None = None,
        state=None,
        decodeFun: Callable[..., Incomplete] | None = None,
        substrateFun: Callable[..., Incomplete] | None = None,
        **options,
    ): ...
    def indefLenValueDecoder(
        self,
        substrate,
        asn1Spec,
        tagSet: TagSet | None = None,
        length: int | None = None,
        state=None,
        decodeFun: Callable[..., Incomplete] | None = None,
        substrateFun: Callable[..., Incomplete] | None = None,
        **options,
    ): ...

class AnyPayloadDecoder(AbstractSimplePayloadDecoder):
    protoComponent: univ.Any
    def valueDecoder(
        self,
        substrate,
        asn1Spec,
        tagSet: TagSet | None = None,
        length: int | None = None,
        state: Unused = None,
        decodeFun: Unused = None,
        substrateFun: Callable[..., Incomplete] | None = None,
        **options,
    ): ...
    def indefLenValueDecoder(
        self,
        substrate,
        asn1Spec,
        tagSet: TagSet | None = None,
        length: int | None = None,
        state: Unused = None,
        decodeFun: Callable[..., Incomplete] | None = None,
        substrateFun: Callable[..., Incomplete] | None = None,
        **options,
    ): ...

class UTF8StringPayloadDecoder(OctetStringPayloadDecoder):
    protoComponent: char.UTF8String

class NumericStringPayloadDecoder(OctetStringPayloadDecoder):
    protoComponent: char.NumericString

class PrintableStringPayloadDecoder(OctetStringPayloadDecoder):
    protoComponent: char.PrintableString

class TeletexStringPayloadDecoder(OctetStringPayloadDecoder):
    protoComponent: char.TeletexString

class VideotexStringPayloadDecoder(OctetStringPayloadDecoder):
    protoComponent: char.VideotexString

class IA5StringPayloadDecoder(OctetStringPayloadDecoder):
    protoComponent: char.IA5String

class GraphicStringPayloadDecoder(OctetStringPayloadDecoder):
    protoComponent: char.GraphicString

class VisibleStringPayloadDecoder(OctetStringPayloadDecoder):
    protoComponent: char.VisibleString

class GeneralStringPayloadDecoder(OctetStringPayloadDecoder):
    protoComponent: char.GeneralString

class UniversalStringPayloadDecoder(OctetStringPayloadDecoder):
    protoComponent: char.UniversalString

class BMPStringPayloadDecoder(OctetStringPayloadDecoder):
    protoComponent: char.BMPString

class ObjectDescriptorPayloadDecoder(OctetStringPayloadDecoder):
    protoComponent: useful.ObjectDescriptor

class GeneralizedTimePayloadDecoder(OctetStringPayloadDecoder):
    protoComponent: useful.GeneralizedTime

class UTCTimePayloadDecoder(OctetStringPayloadDecoder):
    protoComponent: useful.UTCTime

TAG_MAP: dict[TagSet, AbstractPayloadDecoder]
TYPE_MAP: dict[int, AbstractPayloadDecoder]
# deprecated aliases
tagMap = TAG_MAP
typeMap = TYPE_MAP

class SingleItemDecoder:
    defaultErrorState: int
    defaultRawDecoder: AnyPayloadDecoder
    supportIndefLength: bool
    TAG_MAP: dict[TagSet, AbstractPayloadDecoder]
    TYPE_MAP: dict[int, AbstractPayloadDecoder]
    def __init__(self, tagMap=..., typeMap=..., **ignored: Unused) -> None: ...
    def __call__(
        self,
        substrate,
        asn1Spec: Asn1Type | None = None,
        tagSet: TagSet | None = None,
        length: int | None = None,
        state=0,
        decodeFun: Unused = None,
        substrateFun: Callable[..., Incomplete] | None = None,
        **options,
    ): ...

decode: Decoder

class StreamingDecoder:
    SINGLE_ITEM_DECODER: type[SingleItemDecoder]

    def __init__(self, substrate, asn1Spec=None, *, tagMap=..., typeMap=..., **ignored: Unused) -> None: ...
    def __iter__(self): ...

class Decoder:
    STREAMING_DECODER: type[StreamingDecoder]

    @classmethod
    def __call__(cls, substrate, asn1Spec=None, *, tagMap=..., typeMap=..., **ignored: Unused): ...
