from typing import Any
from typing_extensions import TypeAlias

_Data: TypeAlias = str | bytes | bytearray

class BaseSerializer:
    content_types: list[str] | None
    key: str | None
    def get_content_type(self) -> str: ...
    def loads(self, data: _Data) -> Any: ...
    def dumps(self, data: _Data) -> Any: ...

class JsonSerializer(BaseSerializer):
    content_types: list[str]
    key: str

class YamlSerializer(BaseSerializer):
    content_types: list[str]
    key: str

class Serializer:
    serializers: list[BaseSerializer]
    default: str
    def __init__(self, default: str | None = None, serializers: list[BaseSerializer] | None = None) -> None: ...
    def get_serializer(self, name: str | None = None, content_type: str | None = None) -> BaseSerializer: ...
    def loads(self, data: _Data, format: str | None = None) -> Any: ...
    def dumps(self, data: _Data, format: str | None = None) -> Any: ...
    def get_content_type(self, format: str | None = None) -> str: ...
