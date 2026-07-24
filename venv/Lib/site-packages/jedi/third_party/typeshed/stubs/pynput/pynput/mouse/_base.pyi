import enum
import sys
from collections.abc import Callable
from types import TracebackType
from typing import Any
from typing_extensions import Self

from pynput._util import AbstractListener

class Button(enum.Enum):
    unknown = 0
    left = 1
    middle = 2
    right = 3
    if sys.platform == "linux":
        button8 = 8
        button9 = 9
        button10 = 10
        button11 = 11
        button12 = 12
        button13 = 13
        button14 = 14
        button15 = 15
        button16 = 16
        button17 = 17
        button18 = 18
        button19 = 19
        button20 = 20
        button21 = 21
        button22 = 22
        button23 = 23
        button24 = 24
        button25 = 25
        button26 = 26
        button27 = 27
        button28 = 28
        button29 = 29
        button30 = 30
        scroll_down = 5
        scroll_left = 6
        scroll_right = 7
        scroll_up = 4
    if sys.platform == "win32":
        x1 = 0  # Value unknown
        x2 = 0  # Value unknown

class Controller:
    def __init__(self) -> None: ...
    @property
    def position(self) -> tuple[int, int]: ...
    @position.setter
    def position(self, position: tuple[int, int]) -> None: ...
    def scroll(self, dx: int, dy: int) -> None: ...
    def press(self, button: Button) -> None: ...
    def release(self, button: Button) -> None: ...
    def move(self, dx: int, dy: int) -> None: ...
    def click(self, button: Button, count: int = 1) -> None: ...
    def __enter__(self) -> Self: ...
    def __exit__(
        self, exc_type: type[BaseException] | None, exc_val: BaseException | None, exc_tb: TracebackType | None
    ) -> None: ...

class Listener(AbstractListener):
    if sys.platform == "win32":
        WM_LBUTTONDOWN: int
        WM_LBUTTONUP: int
        WM_MBUTTONDOWN: int
        WM_MBUTTONUP: int
        WM_MOUSEMOVE: int
        WM_MOUSEWHEEL: int
        WM_MOUSEHWHEEL: int
        WM_RBUTTONDOWN: int
        WM_RBUTTONUP: int
        WM_XBUTTONDOWN: int
        WM_XBUTTONUP: int

        MK_XBUTTON1: int
        MK_XBUTTON2: int

        XBUTTON1: int
        XBUTTON2: int

        CLICK_BUTTONS: dict[int, tuple[Button, bool]]
        X_BUTTONS: dict[int, dict[int, tuple[Button, bool]]]
        SCROLL_BUTTONS: dict[int, tuple[int, int]]

    def __init__(
        self,
        on_move: Callable[[int, int], bool | None] | None = None,
        on_click: Callable[[int, int, Button, bool], bool | None] | None = None,
        on_scroll: Callable[[int, int, int, int], bool | None] | None = None,
        suppress: bool = False,
        **kwargs: Any,
    ) -> None: ...
