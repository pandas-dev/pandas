import contextlib
import enum
import sys
from collections.abc import Callable, Generator, Iterable, Iterator
from typing import Any, ClassVar, cast
from typing_extensions import Self

from pynput._util import AbstractListener

class KeyCode:
    _PLATFORM_EXTENSIONS: ClassVar[Iterable[str]]  # undocumented
    vk: int | None
    char: str | None
    is_dead: bool | None
    combining: str | None
    def __init__(self, vk: str | None = None, char: str | None = None, is_dead: bool = False, **kwargs: str) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...
    def join(self, key: Self) -> Self: ...
    @classmethod
    def from_vk(cls, vk: int, **kwargs: Any) -> Self: ...
    @classmethod
    def from_char(cls, char: str, **kwargs: Any) -> Self: ...
    @classmethod
    def from_dead(cls, char: str, **kwargs: Any) -> Self: ...

class Key(enum.Enum):
    alt = cast(KeyCode, ...)
    alt_l = cast(KeyCode, ...)
    alt_r = cast(KeyCode, ...)
    alt_gr = cast(KeyCode, ...)
    backspace = cast(KeyCode, ...)
    caps_lock = cast(KeyCode, ...)
    cmd = cast(KeyCode, ...)
    cmd_l = cast(KeyCode, ...)
    cmd_r = cast(KeyCode, ...)
    ctrl = cast(KeyCode, ...)
    ctrl_l = cast(KeyCode, ...)
    ctrl_r = cast(KeyCode, ...)
    delete = cast(KeyCode, ...)
    down = cast(KeyCode, ...)
    end = cast(KeyCode, ...)
    enter = cast(KeyCode, ...)
    esc = cast(KeyCode, ...)
    f1 = cast(KeyCode, ...)
    f2 = cast(KeyCode, ...)
    f3 = cast(KeyCode, ...)
    f4 = cast(KeyCode, ...)
    f5 = cast(KeyCode, ...)
    f6 = cast(KeyCode, ...)
    f7 = cast(KeyCode, ...)
    f8 = cast(KeyCode, ...)
    f9 = cast(KeyCode, ...)
    f10 = cast(KeyCode, ...)
    f11 = cast(KeyCode, ...)
    f12 = cast(KeyCode, ...)
    f13 = cast(KeyCode, ...)
    f14 = cast(KeyCode, ...)
    f15 = cast(KeyCode, ...)
    f16 = cast(KeyCode, ...)
    f17 = cast(KeyCode, ...)
    f18 = cast(KeyCode, ...)
    f19 = cast(KeyCode, ...)
    f20 = cast(KeyCode, ...)
    if sys.platform == "win32":
        f21 = cast(KeyCode, ...)
        f22 = cast(KeyCode, ...)
        f23 = cast(KeyCode, ...)
        f24 = cast(KeyCode, ...)
    home = cast(KeyCode, ...)
    left = cast(KeyCode, ...)
    page_down = cast(KeyCode, ...)
    page_up = cast(KeyCode, ...)
    right = cast(KeyCode, ...)
    shift = cast(KeyCode, ...)
    shift_l = cast(KeyCode, ...)
    shift_r = cast(KeyCode, ...)
    space = cast(KeyCode, ...)
    tab = cast(KeyCode, ...)
    up = cast(KeyCode, ...)
    media_play_pause = cast(KeyCode, ...)
    media_stop = cast(KeyCode, ...)
    media_volume_mute = cast(KeyCode, ...)
    media_volume_down = cast(KeyCode, ...)
    media_volume_up = cast(KeyCode, ...)
    media_previous = cast(KeyCode, ...)
    media_next = cast(KeyCode, ...)
    if sys.platform == "darwin":
        media_eject = cast(KeyCode, ...)
    insert = cast(KeyCode, ...)
    menu = cast(KeyCode, ...)
    num_lock = cast(KeyCode, ...)
    pause = cast(KeyCode, ...)
    print_screen = cast(KeyCode, ...)
    scroll_lock = cast(KeyCode, ...)

class Controller:
    _KeyCode: ClassVar[type[KeyCode]]  # undocumented
    _Key: ClassVar[type[Key]]  # undocumented

    if sys.platform == "linux":
        CTRL_MASK: ClassVar[int]
        SHIFT_MASK: ClassVar[int]

    class InvalidKeyException(Exception): ...
    class InvalidCharacterException(Exception): ...

    def __init__(self) -> None: ...
    def press(self, key: str | Key | KeyCode) -> None: ...
    def release(self, key: str | Key | KeyCode) -> None: ...
    def tap(self, key: str | Key | KeyCode) -> None: ...
    def touch(self, key: str | Key | KeyCode, is_press: bool) -> None: ...
    @contextlib.contextmanager
    def pressed(self, *args: str | Key | KeyCode) -> Generator[None]: ...
    def type(self, string: str) -> None: ...
    @property
    def modifiers(self) -> contextlib.AbstractContextManager[Iterator[set[Key]]]: ...
    @property
    def alt_pressed(self) -> bool: ...
    @property
    def alt_gr_pressed(self) -> bool: ...
    @property
    def ctrl_pressed(self) -> bool: ...
    @property
    def shift_pressed(self) -> bool: ...

class Listener(AbstractListener):
    def __init__(
        self,
        on_press: Callable[[Key | KeyCode | None], None] | None = None,
        on_release: Callable[[Key | KeyCode | None], None] | None = None,
        suppress: bool = False,
        **kwargs: Any,
    ) -> None: ...
    def canonical(self, key: Key | KeyCode) -> Key | KeyCode: ...
