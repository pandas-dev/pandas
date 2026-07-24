from collections.abc import Callable
from queue import Queue
from threading import Lock, Thread
from typing import ClassVar, Literal
from typing_extensions import TypeAlias

from ._keyboard_event import KeyboardEvent
from ._mouse_event import _MouseEvent

_Event: TypeAlias = KeyboardEvent | _MouseEvent

class GenericListener:
    lock: ClassVar[Lock]
    handlers: list[Callable[[_Event], bool | None]]
    listening: bool
    queue: Queue[_Event]
    listening_thread: Thread | None
    processing_thread: Thread | None
    def invoke_handlers(self, event: _Event) -> Literal[1] | None: ...
    def start_if_necessary(self) -> None: ...
    def pre_process_event(self, event: _Event) -> None: ...
    def process(self) -> None: ...
    def add_handler(self, handler: Callable[[_Event], bool | None]) -> None: ...
    def remove_handler(self, handler: Callable[[_Event], bool | None]) -> None: ...
