import contextlib
from _typeshed import ConvertibleToInt
from collections.abc import Callable, Iterable, Sequence
from datetime import datetime
from typing import Final, NamedTuple, SupportsIndex, SupportsInt, TypeVar
from typing_extensions import ParamSpec, TypeAlias

from pyscreeze import (
    center as center,
    locate as locate,
    locateAll as locateAll,
    locateAllOnScreen as locateAllOnScreen,
    locateCenterOnScreen as locateCenterOnScreen,
    locateOnScreen as locateOnScreen,
    locateOnWindow as locateOnWindow,
    pixel as pixel,
    pixelMatchesColor as pixelMatchesColor,
    screenshot as screenshot,
)

_P = ParamSpec("_P")
_R = TypeVar("_R")
# Explicitly mentioning str despite being in the ConvertibleToInt Alias because it has a different meaning (filename on screen)
# Specifying non-None Y arg when X is a string or sequence raises an error
# TODO: This could be better represented through overloads
_NormalizeableXArg: TypeAlias = str | ConvertibleToInt | Sequence[ConvertibleToInt]

# Constants
KEY_NAMES: list[str]
KEYBOARD_KEYS: list[str]
LEFT: Final = "left"
MIDDLE: Final = "middle"
RIGHT: Final = "right"
PRIMARY: Final = "primary"
SECONDARY: Final = "secondary"
G_LOG_SCREENSHOTS_FILENAMES: list[str]
# Implementation details
QWERTY: Final[str]
QWERTZ: Final[str]
MINIMUM_SLEEP: Final[float]

# These are meant to be overridable
LOG_SCREENSHOTS: bool
LOG_SCREENSHOTS_LIMIT: int | None
# https://pyautogui.readthedocs.io/en/latest/index.html#fail-safes
FAILSAFE: bool
PAUSE: float
DARWIN_CATCH_UP_TIME: float
FAILSAFE_POINTS: list[tuple[int, int]]
# https://pyautogui.readthedocs.io/en/latest/mouse.htmln#mouse-movement
MINIMUM_DURATION: float

class PyAutoGUIException(Exception): ...
class FailSafeException(PyAutoGUIException): ...
class ImageNotFoundException(PyAutoGUIException): ...

def raisePyAutoGUIImageNotFoundException(wrappedFunction: Callable[_P, _R]) -> Callable[_P, _R]: ...
def mouseInfo() -> None: ...
def useImageNotFoundException(value: bool | None = None) -> None: ...
def isShiftCharacter(character: str) -> bool: ...

class Point(NamedTuple):
    x: int
    y: int

class Size(NamedTuple):
    width: int
    height: int

def getPointOnLine(x1: float, y1: float, x2: float, y2: float, n: float) -> tuple[float, float]: ...
def linear(n: float) -> float: ...
def position(x: int | None = None, y: int | None = None) -> Point: ...
def size() -> Size: ...

resolution = size

def onScreen(x: _NormalizeableXArg | None, y: SupportsInt | None = None) -> bool: ...
def mouseDown(
    x: _NormalizeableXArg | None = None,
    y: SupportsInt | None = None,
    # Docstring says `button` can also be `int`, but `.lower()` is called unconditionally in `_normalizeButton()`
    button: str = "primary",
    duration: float = 0.0,
    tween: Callable[[float], float] = ...,
    logScreenshot: bool | None = None,
    _pause: bool = True,
) -> None: ...
def mouseUp(
    x: _NormalizeableXArg | None = None,
    y: SupportsInt | None = None,
    # Docstring says `button` can also be `int`, but `.lower()` is called unconditionally in `_normalizeButton()`
    button: str = "primary",
    duration: float = 0.0,
    tween: Callable[[float], float] = ...,
    logScreenshot: bool | None = None,
    _pause: bool = True,
) -> None: ...
def click(
    x: _NormalizeableXArg | None = None,
    y: SupportsInt | None = None,
    clicks: SupportsIndex = 1,
    interval: float = 0.0,
    # Docstring says `button` can also be `int`, but `.lower()` is called unconditionally in `_normalizeButton()`
    button: str = "primary",
    duration: float = 0.0,
    tween: Callable[[float], float] = ...,
    logScreenshot: bool | None = None,
    _pause: bool = True,
) -> None: ...
def leftClick(
    x: _NormalizeableXArg | None = None,
    y: SupportsInt | None = None,
    interval: float = 0.0,
    duration: float = 0.0,
    tween: Callable[[float], float] = ...,
    logScreenshot: bool | None = None,
    _pause: bool = True,
) -> None: ...
def rightClick(
    x: _NormalizeableXArg | None = None,
    y: SupportsInt | None = None,
    interval: float = 0.0,
    duration: float = 0.0,
    tween: Callable[[float], float] = ...,
    logScreenshot: bool | None = None,
    _pause: bool = True,
) -> None: ...
def middleClick(
    x: _NormalizeableXArg | None = None,
    y: SupportsInt | None = None,
    interval: float = 0.0,
    duration: float = 0.0,
    tween: Callable[[float], float] = ...,
    logScreenshot: bool | None = None,
    _pause: bool = True,
) -> None: ...
def doubleClick(
    x: _NormalizeableXArg | None = None,
    y: SupportsInt | None = None,
    interval: float = 0.0,
    # Docstring says `button` can also be `int`, but `.lower()` is called unconditionally in `_normalizeButton()`
    button: str = "left",
    duration: float = 0.0,
    tween: Callable[[float], float] = ...,
    logScreenshot: bool | None = None,
    _pause: bool = True,
) -> None: ...
def tripleClick(
    x: _NormalizeableXArg | None = None,
    y: SupportsInt | None = None,
    interval: float = 0.0,
    # Docstring says `button` can also be `int`, but `.lower()` is called unconditionally in `_normalizeButton()`
    button: str = "left",
    duration: float = 0.0,
    tween: Callable[[float], float] = ...,
    logScreenshot: bool | None = None,
    _pause: bool = True,
) -> None: ...
def scroll(
    clicks: float,
    x: _NormalizeableXArg | None = None,
    y: SupportsInt | None = None,
    logScreenshot: bool | None = None,
    _pause: bool = True,
) -> None: ...
def hscroll(
    clicks: float,
    x: _NormalizeableXArg | None = None,
    y: SupportsInt | None = None,
    logScreenshot: bool | None = None,
    _pause: bool = True,
) -> None: ...
def vscroll(
    clicks: float,
    x: _NormalizeableXArg | None = None,
    y: SupportsInt | None = None,
    logScreenshot: bool | None = None,
    _pause: bool = True,
) -> None: ...
def moveTo(
    x: _NormalizeableXArg | None = None,
    y: SupportsInt | None = None,
    duration: float = 0.0,
    tween: Callable[[float], float] = ...,
    logScreenshot: bool = False,
    _pause: bool = True,
) -> None: ...
def moveRel(
    xOffset: _NormalizeableXArg | None = None,
    yOffset: SupportsInt | None = None,
    duration: float = 0.0,
    tween: Callable[[float], float] = ...,
    logScreenshot: bool = False,
    _pause: bool = True,
) -> None: ...

move = moveRel

def dragTo(
    x: _NormalizeableXArg | None = None,
    y: SupportsInt | None = None,
    duration: float = 0.0,
    tween: Callable[[float], float] = ...,
    # Docstring says `button` can also be `int`, but `.lower()` is called unconditionally in `_normalizeButton()`
    button: str = "primary",
    logScreenshot: bool | None = None,
    _pause: bool = True,
    mouseDownUp: bool = True,
) -> None: ...
def dragRel(
    xOffset: _NormalizeableXArg | None = 0,
    yOffset: SupportsInt | None = 0,
    duration: float = 0.0,
    tween: Callable[[float], float] = ...,
    # Docstring says `button` can also be `int`, but `.lower()` is called unconditionally in `_normalizeButton()`
    button: str = "primary",
    logScreenshot: bool | None = None,
    _pause: bool = True,
    mouseDownUp: bool = True,
) -> None: ...

drag = dragRel

def isValidKey(key: str) -> bool: ...
def keyDown(key: str, logScreenshot: bool | None = None, _pause: bool = True) -> None: ...
def keyUp(key: str, logScreenshot: bool | None = None, _pause: bool = True) -> None: ...
def press(
    keys: str | Iterable[str],
    presses: SupportsIndex = 1,
    interval: float = 0.0,
    logScreenshot: bool | None = None,
    _pause: bool = True,
) -> None: ...
def hold(
    keys: str | Iterable[str], logScreenshot: bool | None = None, _pause: bool = True
) -> contextlib._GeneratorContextManager[None]: ...
def typewrite(
    message: str | Sequence[str], interval: float = 0.0, logScreenshot: bool | None = None, _pause: bool = True
) -> None: ...

write = typewrite

def hotkey(*args: str, logScreenshot: bool | None = None, interval: float = 0.0) -> None: ...

shortcut = hotkey

def failSafeCheck() -> None: ...
def displayMousePosition(xOffset: float = 0, yOffset: float = 0) -> None: ...
def sleep(seconds: float) -> None: ...
def countdown(seconds: SupportsIndex) -> None: ...
def run(commandStr: str, _ssCount: Sequence[int] | None = None) -> None: ...
def printInfo(dontPrint: bool = False) -> str: ...
def getInfo() -> tuple[str, str, str, str, Size, datetime]: ...
