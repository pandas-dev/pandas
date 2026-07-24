from _typeshed import Incomplete, IndexableBuffer, SliceableBuffer, Unused
from collections.abc import Callable, Mapping
from typing import Any, Final, TypeVar, overload
from typing_extensions import TypeAlias

from pyasn1.codec.ber.encoder import AbstractItemEncoder
from pyasn1.type.tag import TagSet

# Use _typeshed._SupportsGetItemBuffer after PEP 688
_SupportsGetItemBuffer: TypeAlias = SliceableBuffer | IndexableBuffer
_R = TypeVar("_R")
_B = TypeVar("_B", bound=_SupportsGetItemBuffer)
# The possible return type is a union of all other decode methods, ie: AnyOf[Incomplete | bool]
_AllDecodersReturnType: TypeAlias = Any

CLASSES: Final[dict[tuple[bool, bool], int]]

class LDAPBooleanEncoder(AbstractItemEncoder):
    supportIndefLenMode: bool
    # Requires pyasn1 > 0.3.7
    def encodeValue(self, value: bool | int, asn1Spec: Unused, encodeFun: Unused, **options: Unused): ...

customTagMap: dict[TagSet, AbstractItemEncoder]
customTypeMap: dict[int, AbstractItemEncoder]

def compute_ber_size(data): ...
def decode_message_fast(message): ...
@overload
def decode_sequence(message: _B, start: int, stop: int, context_decoders: Mapping[int, Callable[[_B, int, int], _R]]) -> _R: ...
@overload
def decode_sequence(
    message: _SupportsGetItemBuffer, start: int, stop: int, context_decoders: None = None
) -> _AllDecodersReturnType: ...
def decode_integer(message, start: int, stop: int, context_decoders: Unused = None): ...
def decode_octet_string(message, start: int, stop: int, context_decoders: Unused = None): ...
def decode_boolean(message, start: int, stop: int, context_decoders: Unused = None): ...
def decode_bind_response(message, start: int, stop: int, context_decoders: Unused = None): ...
def decode_extended_response(message, start: int, stop: int, context_decoders: Unused = None): ...
def decode_intermediate_response(message, start: int, stop: int, context_decoders: Unused = None): ...
def decode_controls(message, start: int, stop: int, context_decoders: Unused = None): ...
def ldap_result_to_dict_fast(response): ...
def get_byte(x): ...
def get_bytes(x): ...

# The possible return type is a union of all other decode methods, ie: AnyOf[Incomplete | bool]
DECODERS: dict[tuple[int, int], Callable[..., _AllDecodersReturnType]]
BIND_RESPONSE_CONTEXT: dict[int, Callable[..., Incomplete]]
EXTENDED_RESPONSE_CONTEXT: dict[int, Callable[..., Incomplete]]
INTERMEDIATE_RESPONSE_CONTEXT: dict[int, Callable[..., Incomplete]]
LDAP_MESSAGE_CONTEXT: dict[int, Callable[..., Incomplete]]
CONTROLS_CONTEXT: dict[int, Callable[..., Incomplete]]
