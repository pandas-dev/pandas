from _typeshed import Incomplete
from collections.abc import Iterable, Mapping
from typing import ClassVar, final

from google._upb._message import Descriptor, FileDescriptor, MessageMeta
from google.protobuf import message
from google.protobuf.internal import containers

DESCRIPTOR: FileDescriptor

@final
class ServerReflectionRequest(message.Message, metaclass=MessageMeta):
    HOST_FIELD_NUMBER: ClassVar[int]
    FILE_BY_FILENAME_FIELD_NUMBER: ClassVar[int]
    FILE_CONTAINING_SYMBOL_FIELD_NUMBER: ClassVar[int]
    FILE_CONTAINING_EXTENSION_FIELD_NUMBER: ClassVar[int]
    ALL_EXTENSION_NUMBERS_OF_TYPE_FIELD_NUMBER: ClassVar[int]
    LIST_SERVICES_FIELD_NUMBER: ClassVar[int]
    host: str
    file_by_filename: str
    file_containing_symbol: str
    file_containing_extension: ExtensionRequest
    all_extension_numbers_of_type: str
    list_services: str
    def __init__(
        self,
        host: str | None = ...,
        file_by_filename: str | None = ...,
        file_containing_symbol: str | None = ...,
        file_containing_extension: ExtensionRequest | Mapping[Incomplete, Incomplete] | None = ...,
        all_extension_numbers_of_type: str | None = ...,
        list_services: str | None = ...,
    ) -> None: ...
    DESCRIPTOR: Descriptor

@final
class ExtensionRequest(message.Message, metaclass=MessageMeta):
    CONTAINING_TYPE_FIELD_NUMBER: ClassVar[int]
    EXTENSION_NUMBER_FIELD_NUMBER: ClassVar[int]
    containing_type: str
    extension_number: int
    def __init__(self, containing_type: str | None = ..., extension_number: int | None = ...) -> None: ...
    DESCRIPTOR: Descriptor

@final
class ServerReflectionResponse(message.Message, metaclass=MessageMeta):
    VALID_HOST_FIELD_NUMBER: ClassVar[int]
    ORIGINAL_REQUEST_FIELD_NUMBER: ClassVar[int]
    FILE_DESCRIPTOR_RESPONSE_FIELD_NUMBER: ClassVar[int]
    ALL_EXTENSION_NUMBERS_RESPONSE_FIELD_NUMBER: ClassVar[int]
    LIST_SERVICES_RESPONSE_FIELD_NUMBER: ClassVar[int]
    ERROR_RESPONSE_FIELD_NUMBER: ClassVar[int]
    valid_host: str
    original_request: ServerReflectionRequest
    file_descriptor_response: FileDescriptorResponse
    all_extension_numbers_response: ExtensionNumberResponse
    list_services_response: ListServiceResponse
    error_response: ErrorResponse
    def __init__(
        self,
        valid_host: str | None = ...,
        original_request: ServerReflectionRequest | Mapping[Incomplete, Incomplete] | None = ...,
        file_descriptor_response: FileDescriptorResponse | Mapping[Incomplete, Incomplete] | None = ...,
        all_extension_numbers_response: ExtensionNumberResponse | Mapping[Incomplete, Incomplete] | None = ...,
        list_services_response: ListServiceResponse | Mapping[Incomplete, Incomplete] | None = ...,
        error_response: ErrorResponse | Mapping[Incomplete, Incomplete] | None = ...,
    ) -> None: ...
    DESCRIPTOR: Descriptor

@final
class FileDescriptorResponse(message.Message, metaclass=MessageMeta):
    FILE_DESCRIPTOR_PROTO_FIELD_NUMBER: ClassVar[int]
    file_descriptor_proto: containers.RepeatedScalarFieldContainer[bytes]
    def __init__(self, file_descriptor_proto: Iterable[bytes] | None = ...) -> None: ...
    DESCRIPTOR: Descriptor

@final
class ExtensionNumberResponse(message.Message, metaclass=MessageMeta):
    BASE_TYPE_NAME_FIELD_NUMBER: ClassVar[int]
    EXTENSION_NUMBER_FIELD_NUMBER: ClassVar[int]
    base_type_name: str
    extension_number: containers.RepeatedScalarFieldContainer[int]
    def __init__(self, base_type_name: str | None = ..., extension_number: Iterable[int] | None = ...) -> None: ...
    DESCRIPTOR: Descriptor

@final
class ListServiceResponse(message.Message, metaclass=MessageMeta):
    SERVICE_FIELD_NUMBER: ClassVar[int]
    service: containers.RepeatedCompositeFieldContainer[ServiceResponse]
    def __init__(self, service: Iterable[ServiceResponse | Mapping[Incomplete, Incomplete]] | None = ...) -> None: ...
    DESCRIPTOR: Descriptor

@final
class ServiceResponse(message.Message, metaclass=MessageMeta):
    NAME_FIELD_NUMBER: ClassVar[int]
    name: str
    def __init__(self, name: str | None = ...) -> None: ...
    DESCRIPTOR: Descriptor

@final
class ErrorResponse(message.Message, metaclass=MessageMeta):
    ERROR_CODE_FIELD_NUMBER: ClassVar[int]
    ERROR_MESSAGE_FIELD_NUMBER: ClassVar[int]
    error_code: int
    error_message: str
    def __init__(self, error_code: int | None = ..., error_message: str | None = ...) -> None: ...
    DESCRIPTOR: Descriptor
