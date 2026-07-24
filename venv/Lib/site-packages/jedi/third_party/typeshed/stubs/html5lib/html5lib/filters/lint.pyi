from _typeshed import Incomplete
from collections.abc import Iterable

from . import base

spaceCharacters: str

class Filter(base.Filter[dict[str, Incomplete]]):
    require_matching_tags: bool
    def __init__(self, source: Iterable[dict[str, Incomplete]], require_matching_tags: bool = True) -> None: ...
