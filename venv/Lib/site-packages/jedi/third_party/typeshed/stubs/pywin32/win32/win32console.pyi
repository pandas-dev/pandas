from typing import Literal, NoReturn, overload

import _win32typing
from win32.lib.pywintypes import error as error

def GetConsoleProcessList() -> tuple[int, ...]: ...
def CreateConsoleScreenBuffer(
    DesiredAccess=..., ShareMode=..., SecurityAttributes: _win32typing.PySECURITY_ATTRIBUTES | None = ..., Flags=...
) -> _win32typing.PyConsoleScreenBuffer: ...
def GetConsoleDisplayMode(): ...
def AttachConsole(ProcessId) -> None: ...
def AllocConsole() -> None: ...
def FreeConsole() -> None: ...
def GetConsoleCP(): ...
def GetConsoleOutputCP(): ...
def SetConsoleCP(CodePageId) -> None: ...
def SetConsoleOutputCP(CodePageID) -> None: ...
def GetConsoleSelectionInfo(): ...
def AddConsoleAlias(Source, Target, ExeName) -> None: ...
def GetConsoleAliases(ExeName: str) -> str: ...
def GetConsoleAliasExes(): ...
def GetConsoleWindow(): ...
def GetNumberOfConsoleFonts(): ...
def SetConsoleTitle(ConsoleTitle: str) -> None: ...
def GetConsoleTitle(): ...
@overload
def GenerateConsoleCtrlEvent(CtrlEvent: Literal[1], ProcessGroupId: Literal[0] = 0) -> NoReturn: ...
@overload
def GenerateConsoleCtrlEvent(CtrlEvent: Literal[0, 1], ProcessGroupId: int) -> None: ...
def GetStdHandle(StdHandle: int) -> _win32typing.PyConsoleScreenBuffer: ...

ATTACH_PARENT_PROCESS: int
BACKGROUND_BLUE: int
BACKGROUND_GREEN: int
BACKGROUND_INTENSITY: int
BACKGROUND_RED: int
COMMON_LVB_GRID_HORIZONTAL: int
COMMON_LVB_GRID_LVERTICAL: int
COMMON_LVB_GRID_RVERTICAL: int
COMMON_LVB_LEADING_BYTE: int
COMMON_LVB_REVERSE_VIDEO: int
COMMON_LVB_TRAILING_BYTE: int
COMMON_LVB_UNDERSCORE: int
CONSOLE_FULLSCREEN: int
CONSOLE_FULLSCREEN_HARDWARE: int
CONSOLE_FULLSCREEN_MODE: int
CONSOLE_MOUSE_DOWN: int
CONSOLE_MOUSE_SELECTION: int
CONSOLE_NO_SELECTION: int
CONSOLE_SELECTION_IN_PROGRESS: int
CONSOLE_SELECTION_NOT_EMPTY: int
CONSOLE_TEXTMODE_BUFFER: int
CONSOLE_WINDOWED_MODE: int
CTRL_BREAK_EVENT: int
CTRL_C_EVENT: int
ENABLE_ECHO_INPUT: int
ENABLE_LINE_INPUT: int
ENABLE_MOUSE_INPUT: int
ENABLE_PROCESSED_INPUT: int
ENABLE_PROCESSED_OUTPUT: int
ENABLE_WINDOW_INPUT: int
ENABLE_WRAP_AT_EOL_OUTPUT: int
FOCUS_EVENT: int
FOREGROUND_BLUE: int
FOREGROUND_GREEN: int
FOREGROUND_INTENSITY: int
FOREGROUND_RED: int
KEY_EVENT: int
LOCALE_USER_DEFAULT: int
MENU_EVENT: int
MOUSE_EVENT: int
PyCOORDType = _win32typing.PyCOORD
PyConsoleScreenBufferType = _win32typing.PyConsoleScreenBuffer
PyINPUT_RECORDType = _win32typing.PyINPUT_RECORD
PySMALL_RECTType = _win32typing.PySMALL_RECT
STD_ERROR_HANDLE: int
STD_INPUT_HANDLE: int
STD_OUTPUT_HANDLE: int
WINDOW_BUFFER_SIZE_EVENT: int
