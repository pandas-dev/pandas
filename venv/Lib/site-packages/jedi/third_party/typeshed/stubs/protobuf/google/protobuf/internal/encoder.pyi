from collections.abc import Callable
from typing_extensions import TypeAlias

from google.protobuf.descriptor import FieldDescriptor

_Sizer: TypeAlias = Callable[[int, bool, bool], int]

Int32Sizer: _Sizer
Int64Sizer: _Sizer
EnumSizer: _Sizer
UInt32Sizer: _Sizer
UInt64Sizer: _Sizer
SInt32Sizer: _Sizer
SInt64Sizer: _Sizer
Fixed32Sizer: _Sizer
SFixed32Sizer: _Sizer
FloatSizer: _Sizer
Fixed64Sizer: _Sizer
SFixed64Sizer: _Sizer
DoubleSizer: _Sizer
BoolSizer: _Sizer

def StringSizer(field_number: int, is_repeated: bool, is_packed: bool) -> _Sizer: ...
def BytesSizer(field_number: int, is_repeated: bool, is_packed: bool) -> _Sizer: ...
def GroupSizer(field_number: int, is_repeated: bool, is_packed: bool) -> _Sizer: ...
def MessageSizer(field_number: int, is_repeated: bool, is_packed: bool) -> _Sizer: ...
def MessageSetItemSizer(field_number: int) -> _Sizer: ...
def MapSizer(field_descriptor: FieldDescriptor, is_message_map: bool) -> _Sizer: ...
def TagBytes(field_number: int, wire_type: int) -> bytes: ...

_Encoder: TypeAlias = Callable[[Callable[[bytes], int], bytes, bool], int]

Int32Encoder: _Encoder
Int64Encoder: _Encoder
EnumEncoder: _Encoder
UInt32Encoder: _Encoder
UInt64Encoder: _Encoder
SInt32Encoder: _Encoder
SInt64Encoder: _Encoder
Fixed32Encoder: _Encoder
Fixed64Encoder: _Encoder
SFixed32Encoder: _Encoder
SFixed64Encoder: _Encoder
FloatEncoder: _Encoder
DoubleEncoder: _Encoder

def BoolEncoder(field_number: int, is_repeated: bool, is_packed: bool) -> _Encoder: ...
def StringEncoder(field_number: int, is_repeated: bool, is_packed: bool) -> _Encoder: ...
def BytesEncoder(field_number: int, is_repeated: bool, is_packed: bool) -> _Encoder: ...
def GroupEncoder(field_number: int, is_repeated: bool, is_packed: bool) -> _Encoder: ...
def MessageEncoder(field_number: int, is_repeated: bool, is_packed: bool) -> _Encoder: ...
def MessageSetItemEncoder(field_number: int) -> _Encoder: ...
def MapEncoder(field_descriptor: FieldDescriptor) -> _Encoder: ...
