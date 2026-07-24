from collections.abc import Mapping, Sequence
from typing import Any

from .descriptor_pb2 import (
    DescriptorProto,
    EnumDescriptorProto,
    EnumOptions,
    EnumValueOptions,
    FieldOptions,
    FileDescriptorProto,
    FileOptions,
    MessageOptions,
    MethodDescriptorProto,
    MethodOptions,
    OneofOptions,
    ServiceDescriptorProto,
    ServiceOptions,
)
from .descriptor_pool import DescriptorPool
from .message import Message

class Error(Exception): ...
class TypeTransformationError(Error): ...

class DescriptorMetaclass(type):
    def __instancecheck__(cls, obj: object) -> bool: ...

_internal_create_key: object
_USE_C_DESCRIPTORS: bool

class DescriptorBase(metaclass=DescriptorMetaclass):
    has_options: bool
    def __init__(self, file, options, serialized_options, options_class_name) -> None: ...
    def GetOptions(self) -> Any: ...  # Any: overridden with specific *Options in subclasses

class _NestedDescriptorBase(DescriptorBase):
    name: str
    full_name: str
    file: FileDescriptor
    containing_type: Descriptor | None
    def __init__(
        self,
        options,
        options_class_name,
        name,
        full_name,
        file,
        containing_type,
        serialized_start=None,
        serialized_end=None,
        serialized_options=None,
    ) -> None: ...
    def CopyToProto(self, proto: Any) -> None: ...  # Any: overridden with specific *Proto in subclasses

class Descriptor(_NestedDescriptorBase):
    fields: Sequence[FieldDescriptor]
    fields_by_number: Mapping[int, FieldDescriptor]
    fields_by_name: Mapping[str, FieldDescriptor]
    @property
    def fields_by_camelcase_name(self) -> Mapping[str, FieldDescriptor]: ...
    nested_types: Sequence[Descriptor]
    nested_types_by_name: Mapping[str, Descriptor]
    enum_types: Sequence[EnumDescriptor]
    enum_types_by_name: Mapping[str, EnumDescriptor]
    enum_values_by_name: dict[str, EnumValueDescriptor]
    extensions: Sequence[FieldDescriptor]
    extensions_by_name: Mapping[str, FieldDescriptor]
    is_extendable: bool
    extension_ranges: list[tuple[int, int]]
    oneofs: Sequence[OneofDescriptor]
    oneofs_by_name: Mapping[str, OneofDescriptor]
    def __init__(
        self,
        name: str,
        full_name: str,
        filename: str | None,
        containing_type: Descriptor | None,
        fields: list[FieldDescriptor],
        nested_types: list[FieldDescriptor],
        enum_types: list[EnumDescriptor],
        extensions: list[FieldDescriptor],
        options: MessageOptions | None = None,
        serialized_options: bytes | None = None,
        is_extendable: bool | None = True,
        extension_ranges: list[tuple[int, int]] | None = None,
        oneofs: list[OneofDescriptor] | None = None,
        file: FileDescriptor | None = None,
        serialized_start: int | None = None,
        serialized_end: int | None = None,
        syntax: str | None = None,
        is_map_entry: bool = False,
        create_key: object | None = None,
    ): ...
    def EnumValueName(self, enum: str, value: int) -> str: ...
    def CopyToProto(self, proto: DescriptorProto) -> None: ...
    def GetOptions(self) -> MessageOptions: ...

class FieldDescriptor(DescriptorBase):
    TYPE_DOUBLE: int
    TYPE_FLOAT: int
    TYPE_INT64: int
    TYPE_UINT64: int
    TYPE_INT32: int
    TYPE_FIXED64: int
    TYPE_FIXED32: int
    TYPE_BOOL: int
    TYPE_STRING: int
    TYPE_GROUP: int
    TYPE_MESSAGE: int
    TYPE_BYTES: int
    TYPE_UINT32: int
    TYPE_ENUM: int
    TYPE_SFIXED32: int
    TYPE_SFIXED64: int
    TYPE_SINT32: int
    TYPE_SINT64: int
    MAX_TYPE: int
    CPPTYPE_INT32: int
    CPPTYPE_INT64: int
    CPPTYPE_UINT32: int
    CPPTYPE_UINT64: int
    CPPTYPE_DOUBLE: int
    CPPTYPE_FLOAT: int
    CPPTYPE_BOOL: int
    CPPTYPE_ENUM: int
    CPPTYPE_STRING: int
    CPPTYPE_MESSAGE: int
    MAX_CPPTYPE: int
    LABEL_OPTIONAL: int
    LABEL_REQUIRED: int
    LABEL_REPEATED: int
    MAX_LABEL: int
    MAX_FIELD_NUMBER: int
    FIRST_RESERVED_FIELD_NUMBER: int
    LAST_RESERVED_FIELD_NUMBER: int
    def __new__(
        cls,
        name,
        full_name,
        index,
        number,
        type,
        cpp_type,
        label,
        default_value,
        message_type,
        enum_type,
        containing_type,
        is_extension,
        extension_scope,
        options=None,
        serialized_options=None,
        has_default_value=True,
        containing_oneof=None,
        json_name=None,
        file=None,
        create_key=None,
    ): ...
    name: str
    full_name: str
    index: int
    number: int
    type: int
    cpp_type: int
    @property
    def is_required(self) -> bool: ...
    @property
    def is_repeated(self) -> bool: ...
    @property
    def camelcase_name(self) -> str: ...
    @property
    def has_presence(self) -> bool: ...
    @property
    def is_packed(self) -> bool: ...
    has_default_value: bool
    default_value: Any  # Any: str, int, float, bytes, or bool
    containing_type: Descriptor | None
    message_type: Descriptor | None
    enum_type: EnumDescriptor | None
    is_extension: bool
    extension_scope: Descriptor | None
    containing_oneof: OneofDescriptor | None
    json_name: str
    def __init__(
        self,
        name,
        full_name,
        index,
        number,
        type,
        cpp_type,
        label,
        default_value,
        message_type,
        enum_type,
        containing_type,
        is_extension,
        extension_scope,
        options=None,
        serialized_options=None,
        has_default_value=True,
        containing_oneof=None,
        json_name=None,
        file=None,
        create_key=None,
    ) -> None: ...
    @staticmethod
    def ProtoTypeToCppProtoType(proto_type: int) -> int: ...
    def GetOptions(self) -> FieldOptions: ...

class EnumDescriptor(_NestedDescriptorBase):
    def __new__(
        cls,
        name,
        full_name,
        filename,
        values,
        containing_type=None,
        options=None,
        serialized_options=None,
        file=None,
        serialized_start=None,
        serialized_end=None,
        create_key=None,
    ): ...
    values: Sequence[EnumValueDescriptor]
    values_by_name: Mapping[str, EnumValueDescriptor]
    values_by_number: Mapping[int, EnumValueDescriptor]
    def __init__(
        self,
        name,
        full_name,
        filename,
        values,
        containing_type=None,
        options=None,
        serialized_options=None,
        file=None,
        serialized_start=None,
        serialized_end=None,
        create_key=None,
    ) -> None: ...
    @property
    def is_closed(self) -> bool: ...
    def CopyToProto(self, proto: EnumDescriptorProto) -> None: ...
    def GetOptions(self) -> EnumOptions: ...

class EnumValueDescriptor(DescriptorBase):
    def __new__(cls, name, index, number, type=None, options=None, serialized_options=None, create_key=None): ...
    name: str
    index: int
    number: int
    type: EnumDescriptor
    def __init__(self, name, index, number, type=None, options=None, serialized_options=None, create_key=None) -> None: ...
    def GetOptions(self) -> EnumValueOptions: ...

class OneofDescriptor(DescriptorBase):
    def __new__(cls, name, full_name, index, containing_type, fields, options=None, serialized_options=None, create_key=None): ...
    name: str
    full_name: str
    index: int
    containing_type: Descriptor
    fields: Sequence[FieldDescriptor]
    def __init__(
        self, name, full_name, index, containing_type, fields, options=None, serialized_options=None, create_key=None
    ) -> None: ...
    def GetOptions(self) -> OneofOptions: ...

class ServiceDescriptor(_NestedDescriptorBase):
    index: int
    methods: Sequence[MethodDescriptor]
    methods_by_name: Mapping[str, MethodDescriptor]
    def __init__(
        self,
        name: str,
        full_name: str,
        index: int,
        methods: list[MethodDescriptor],
        options: ServiceOptions | None = None,
        serialized_options: bytes | None = None,
        file: FileDescriptor | None = None,
        serialized_start: int | None = None,
        serialized_end: int | None = None,
        create_key: object | None = None,
    ): ...
    def FindMethodByName(self, name: str) -> MethodDescriptor: ...
    def CopyToProto(self, proto: ServiceDescriptorProto) -> None: ...
    def GetOptions(self) -> ServiceOptions: ...

class MethodDescriptor(DescriptorBase):
    def __new__(
        cls,
        name,
        full_name,
        index,
        containing_service,
        input_type,
        output_type,
        client_streaming=False,
        server_streaming=False,
        options=None,
        serialized_options=None,
        create_key=None,
    ): ...
    name: str
    full_name: str
    index: int
    containing_service: ServiceDescriptor
    input_type: Descriptor
    output_type: Descriptor
    client_streaming: bool
    server_streaming: bool
    def __init__(
        self,
        name,
        full_name,
        index,
        containing_service,
        input_type,
        output_type,
        client_streaming=False,
        server_streaming=False,
        options=None,
        serialized_options=None,
        create_key=None,
    ) -> None: ...
    def CopyToProto(self, proto: MethodDescriptorProto) -> None: ...
    def GetOptions(self) -> MethodOptions: ...

class FileDescriptor(DescriptorBase):
    def __new__(
        cls,
        name,
        package,
        options=None,
        serialized_options=None,
        serialized_pb=None,
        dependencies=None,
        public_dependencies=None,
        syntax=None,
        edition=None,
        pool=None,
        create_key=None,
    ): ...
    _options: Any
    pool: DescriptorPool
    message_types_by_name: Mapping[str, Descriptor]
    name: str
    package: str
    serialized_pb: bytes
    enum_types_by_name: Mapping[str, EnumDescriptor]
    extensions_by_name: Mapping[str, FieldDescriptor]
    services_by_name: Mapping[str, ServiceDescriptor]
    dependencies: Sequence[FileDescriptor]
    public_dependencies: Sequence[FileDescriptor]
    def __init__(
        self,
        name,
        package,
        options=None,
        serialized_options=None,
        serialized_pb=None,
        dependencies=None,
        public_dependencies=None,
        syntax=None,
        edition=None,
        pool=None,
        create_key=None,
    ) -> None: ...
    def CopyToProto(self, proto: FileDescriptorProto) -> None: ...
    def GetOptions(self) -> FileOptions: ...

def _ParseOptions(message: Message, string: bytes) -> Message: ...
def MakeDescriptor(
    desc_proto: DescriptorProto,
    package: str = "",
    build_file_if_cpp: bool = True,
    syntax: str | None = None,
    edition: str | None = None,
    file_desc: FileDescriptor | None = None,
) -> Descriptor: ...
