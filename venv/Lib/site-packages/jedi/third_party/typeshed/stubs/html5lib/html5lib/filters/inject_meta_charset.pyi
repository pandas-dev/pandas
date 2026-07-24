from _typeshed import Incomplete
from collections.abc import Iterable

from . import base

class Filter(base.Filter[dict[str, Incomplete]]):
    encoding: str | None
    def __init__(self, source: Iterable[dict[str, Incomplete]], encoding: str | None) -> None: ...
