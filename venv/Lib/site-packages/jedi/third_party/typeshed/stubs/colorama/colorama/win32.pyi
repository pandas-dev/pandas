import sys
from collections.abc import Callable
from typing import Literal

STDOUT: Literal[-11]
STDERR: Literal[-12]
ENABLE_VIRTUAL_TERMINAL_PROCESSING: int

if sys.platform == "win32":
    from ctypes import LibraryLoader, Structure, WinDLL, wintypes

    windll: LibraryLoader[WinDLL]
    COORD = wintypes._COORD

    class CONSOLE_SCREEN_BUFFER_INFO(Structure):
        dwSize: COORD
        dwCursorPosition: COORD
        wAttributes: wintypes.WORD
        srWindow: wintypes.SMALL_RECT
        dwMaximumWindowSize: COORD

    def winapi_test() -> bool: ...
    def GetConsoleScreenBufferInfo(stream_id: int = -11) -> CONSOLE_SCREEN_BUFFER_INFO: ...
    def SetConsoleTextAttribute(stream_id: int, attrs: wintypes.WORD) -> wintypes.BOOL: ...
    def SetConsoleCursorPosition(stream_id: int, position: COORD, adjust: bool = True) -> wintypes.BOOL: ...
    def FillConsoleOutputCharacter(stream_id: int, char: str, length: int, start: COORD) -> int: ...
    def FillConsoleOutputAttribute(stream_id: int, attr: int, length: int, start: COORD) -> wintypes.BOOL: ...
    def SetConsoleTitle(title: str) -> wintypes.BOOL: ...
    def GetConsoleMode(handle: int) -> int: ...
    def SetConsoleMode(handle: int, mode: int) -> None: ...

else:
    windll: None
    SetConsoleTextAttribute: Callable[..., None]
    winapi_test: Callable[..., None]
