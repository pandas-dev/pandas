from collections import defaultdict
from typing import TypeVar

_KT = TypeVar("_KT")
_VT = TypeVar("_VT")

class BoundDictionary(defaultdict[_KT, _VT]):
    reference: str | None
    def __init__(self, reference: str | None = None, *args, **kw) -> None: ...
