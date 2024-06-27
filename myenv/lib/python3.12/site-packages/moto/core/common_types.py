from typing import Any, Dict, Tuple, TypeVar, Union

TYPE_RESPONSE = Tuple[int, Dict[str, Any], Union[str, bytes]]
TYPE_IF_NONE = TypeVar("TYPE_IF_NONE")
