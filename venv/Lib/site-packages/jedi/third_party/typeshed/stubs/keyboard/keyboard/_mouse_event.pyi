import sys
from typing import Literal, NamedTuple
from typing_extensions import TypeAlias

_MouseEvent: TypeAlias = ButtonEvent | WheelEvent | MoveEvent  # noqa: Y047  # Used outside

LEFT: Literal["left"]
RIGHT: Literal["right"]
MIDDLE: Literal["middle"]
X: Literal["x"]
X2: Literal["x2"]

UP: Literal["up"]
DOWN: Literal["down"]
DOUBLE: Literal["double"]
WHEEL: Literal["wheel"]

VERTICAL: Literal["vertical"]
HORIZONTAL: Literal["horizontal"]

if sys.platform == "linux" or sys.platform == "win32":
    _MouseButton: TypeAlias = Literal["left", "right", "middle", "x", "x2"]
else:
    _MouseButton: TypeAlias = Literal["left", "right", "middle"]

if sys.platform == "win32":
    _MouseEventType: TypeAlias = Literal["up", "down", "double", "wheel"]
else:
    _MouseEventType: TypeAlias = Literal["up", "down"]

class ButtonEvent(NamedTuple):
    event_type: _MouseEventType
    button: _MouseButton
    time: float

class WheelEvent(NamedTuple):
    delta: int
    time: float

class MoveEvent(NamedTuple):
    x: int
    y: int
    time: float
