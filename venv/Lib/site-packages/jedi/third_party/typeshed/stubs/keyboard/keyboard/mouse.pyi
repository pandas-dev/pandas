import sys
from collections.abc import Callable, Iterable
from ctypes import c_long
from typing import Literal, SupportsInt, TypeVar
from typing_extensions import TypeAlias

from ._generic import GenericListener as _GenericListener
from ._mouse_event import (
    DOUBLE as DOUBLE,
    DOWN as DOWN,
    LEFT as LEFT,
    MIDDLE as MIDDLE,
    RIGHT as RIGHT,
    UP as UP,
    X2 as X2,
    ButtonEvent as ButtonEvent,
    MoveEvent as MoveEvent,
    WheelEvent as WheelEvent,
    X as X,
    _MouseButton,
    _MouseEvent,
    _MouseEventType,
)

# mypy doesn't support PEP 646's TypeVarTuple yet: https://github.com/python/mypy/issues/12280
# _Ts = TypeVarTuple("_Ts")
_Ts: TypeAlias = tuple[object, ...]
_Callback: TypeAlias = Callable[[_MouseEvent], bool | None]
_C = TypeVar("_C", bound=_Callback)

class _MouseListener(_GenericListener):
    def init(self) -> None: ...
    def pre_process_event(  # type: ignore[override]  # Mouse specific events and return
        self, event: _MouseEvent
    ) -> Literal[True]: ...
    def listen(self) -> None: ...

def is_pressed(button: _MouseButton = "left") -> bool: ...
def press(button: _MouseButton = "left") -> None: ...
def release(button: _MouseButton = "left") -> None: ...
def click(button: _MouseButton = "left") -> None: ...
def double_click(button: _MouseButton = "left") -> None: ...
def right_click() -> None: ...
def wheel(delta: int = 1) -> None: ...
def move(x: SupportsInt, y: SupportsInt, absolute: bool = True, duration: float = 0) -> None: ...
def drag(start_x: int, start_y: int, end_x: int, end_y: int, absolute: bool = True, duration: float = 0) -> None: ...
def on_button(
    callback: Callable[..., None],
    args: _Ts = (),
    # Omitting default: Darwin has no x and x2
    buttons: list[_MouseButton] | tuple[_MouseButton, ...] | _MouseButton = ...,
    # Omitting default: Darwin and Linux don't have "double", yet the defaults includes it
    types: list[_MouseEventType] | tuple[_MouseEventType, ...] | _MouseEventType = ...,
) -> _Callback: ...
def on_click(callback: Callable[..., None], args: _Ts = ()) -> _Callback: ...
def on_double_click(callback: Callable[..., None], args: _Ts = ()) -> _Callback: ...
def on_right_click(callback: Callable[..., None], args: _Ts = ()) -> _Callback: ...
def on_middle_click(callback: Callable[..., None], args: _Ts = ()) -> _Callback: ...
def wait(
    button: _MouseButton = "left",
    # Omitting default: Darwin and Linux don't have "double", yet the defaults includes it
    target_types: tuple[_MouseEventType, ...] = ...,
) -> None: ...

if sys.platform == "win32":
    def get_position() -> tuple[c_long, c_long]: ...

else:
    def get_position() -> tuple[int, int]: ...

def hook(callback: _C) -> _C: ...
def unhook(callback: _Callback) -> None: ...
def unhook_all() -> None: ...
def record(button: _MouseButton = "right", target_types: tuple[_MouseEventType] = ("down",)) -> _MouseEvent: ...
def play(
    events: Iterable[_MouseEvent],
    speed_factor: float = 1.0,
    include_clicks: bool = True,
    include_moves: bool = True,
    include_wheel: bool = True,
) -> None: ...

replay = play
hold = press
