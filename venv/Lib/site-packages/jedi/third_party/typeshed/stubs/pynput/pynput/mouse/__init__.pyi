from typing import Any

from pynput import _util

from ._base import Button as Button, Controller as Controller, Listener as Listener

class Events(_util.Events[Any, Listener]):
    class Move(_util.Events.Event):
        x: int
        y: int
        injected: bool
        def __init__(self, x: int, y: int, injected: bool) -> None: ...

    class Click(_util.Events.Event):
        x: int
        y: int
        button: Button
        pressed: bool
        injected: bool
        def __init__(self, x: int, y: int, button: Button, pressed: bool, injected: bool) -> None: ...

    class Scroll(_util.Events.Event):
        x: int
        y: int
        dx: int
        dy: int
        injected: bool
        def __init__(self, x: int, y: int, dx: int, dy: int, injected: bool) -> None: ...

    def __init__(self) -> None: ...
    def __next__(self) -> Move | Click | Scroll: ...
    def get(self, timeout: float | None = None) -> Move | Click | Scroll | None: ...
