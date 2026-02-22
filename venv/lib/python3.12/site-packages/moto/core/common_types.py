from typing import Any, TypeVar, Union

TYPE_RESPONSE = tuple[int, dict[str, Any], Union[str, bytes]]
TYPE_IF_NONE = TypeVar("TYPE_IF_NONE")
