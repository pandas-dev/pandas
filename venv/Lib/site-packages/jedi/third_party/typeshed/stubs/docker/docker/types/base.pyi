from collections.abc import Mapping
from typing import TypeVar

_VT = TypeVar("_VT")

class DictType(dict[str, _VT]):
    def __init__(self, init: Mapping[str, _VT]) -> None: ...
