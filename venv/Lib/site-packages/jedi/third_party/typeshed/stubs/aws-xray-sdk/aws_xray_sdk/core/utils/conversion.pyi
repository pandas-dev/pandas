from logging import Logger
from typing import Any, TypeVar, overload

_K = TypeVar("_K")

log: Logger

@overload
def metadata_to_dict(obj: dict[_K, Any]) -> dict[_K, Any]: ...
@overload
def metadata_to_dict(obj: type) -> str: ...
@overload
def metadata_to_dict(obj: Any) -> Any: ...
