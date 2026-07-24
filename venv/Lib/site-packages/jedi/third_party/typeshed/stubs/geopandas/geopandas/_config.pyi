from _typeshed import SupportsItems, Unused
from collections.abc import Callable
from typing import Any, NamedTuple

class Option(NamedTuple):
    key: str
    default_value: Any  # Can be "any" type
    doc: str
    validator: Callable[[object], Unused]
    callback: Callable[[str, object], Unused] | None

class Options:
    def __init__(self, options: SupportsItems[str, Option]) -> None: ...
    # Accept and return arbitrary values
    def __setattr__(self, key: str, value: Any) -> None: ...
    def __getattr__(self, key: str) -> Any: ...

display_precision: Option
use_pygeos: Option
io_engine: Option
options: Options
