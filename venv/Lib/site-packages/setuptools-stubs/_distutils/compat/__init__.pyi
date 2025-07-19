from collections.abc import Iterable
from typing import TypeVar

_IterableT = TypeVar("_IterableT", bound=Iterable[str])

def consolidate_linker_args(args: _IterableT) -> _IterableT | str: ...
