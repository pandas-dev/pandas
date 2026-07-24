from _typeshed import StrPath
from typing import Any

from .formats.base_formats import Format
from .resources import ModelResource

def get_resource_class(model_or_resource_class: str) -> ModelResource[Any]: ...

MIME_TYPE_FORMAT_MAPPING: dict[str, type[Format]]

def get_format_class(format_name: str, file_name: StrPath, encoding: str | None = None) -> Format: ...
def get_default_format_names() -> str: ...
