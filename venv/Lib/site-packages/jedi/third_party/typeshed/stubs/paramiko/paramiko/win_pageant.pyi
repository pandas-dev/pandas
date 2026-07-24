import ctypes
import sys
from _typeshed import Incomplete
from typing import Literal
from typing_extensions import TypeAlias

if sys.platform == "win32":
    win32con_WM_COPYDATA: int
    def can_talk_to_agent() -> bool: ...

    ULONG_PTR: TypeAlias = ctypes.c_uint64 | ctypes.c_uint32

    class COPYDATASTRUCT(ctypes.Structure):
        num_data: Incomplete
        data_size: Incomplete
        data_loc: Incomplete

    class PageantConnection:
        def __init__(self) -> None: ...
        def send(self, data: bytes) -> None: ...
        def recv(self, n: int) -> Literal[""] | bytes: ...
        def close(self) -> None: ...
