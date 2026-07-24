from _typeshed import SupportsWrite
from collections.abc import Callable
from types import FrameType
from typing import Any

from .watch_element import WatchElement

class WatchPrint:
    # User-defined callbacks passed to `__init__` set as instance variables have arguments of type `Any` to be
    # compatible with more precisely-annotated signatures. These callbacks are passed from `watchpoints.watch.Watch`.

    custom_printer: Callable[[Any], None] | None
    file: str | SupportsWrite[str] | None
    stack_limit: int | None

    def __init__(
        self,
        file: str | SupportsWrite[str] | None = ...,
        stack_limit: int | None = None,
        custom_printer: Callable[[Any], None] | None = None,  # User-defined printing callback
    ) -> None: ...
    def __call__(self, frame: FrameType, elem: WatchElement, exec_info: tuple[str, str, int | None]) -> None: ...
    def getsourceline(self, exec_info: tuple[str, str, int | None]) -> str: ...
    def printer(self, obj: object) -> None: ...
