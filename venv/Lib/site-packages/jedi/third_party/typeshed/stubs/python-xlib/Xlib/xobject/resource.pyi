from Xlib._typing import ErrorHandler
from Xlib.display import _BaseDisplay

class Resource:
    display: _BaseDisplay
    id: int
    owner: int
    def __init__(self, display: _BaseDisplay, rid: int, owner: int = 0) -> None: ...
    def __resource__(self) -> int: ...
    def kill_client(self, onerror: ErrorHandler[object] | None = None) -> None: ...
