from _typeshed import Unused
from collections.abc import Callable, Sequence, Sized
from typing import Any, Final, Literal, TypeVar

from Xlib.display import Display
from Xlib.protocol import display, rq
from Xlib.xobject import resource

_T = TypeVar("_T")
_S = TypeVar("_S", bound=Sized)

extname: Final = "RECORD"

FromServerTime: Final = 0x01
FromClientTime: Final = 0x02
FromClientSequence: Final = 0x04

CurrentClients: Final = 1
FutureClients: Final = 2
AllClients: Final = 3

FromServer: Final = 0
FromClient: Final = 1
ClientStarted: Final = 2
ClientDied: Final = 3
StartOfData: Final = 4
EndOfData: Final = 5
Record_Range8: rq.Struct
Record_Range16: rq.Struct
Record_ExtRange: rq.Struct
Record_Range: rq.Struct
Record_ClientInfo: rq.Struct

class RawField(rq.ValueField):
    structcode: None
    def pack_value(self, val: _S) -> tuple[_S, int, None]: ...  # type: ignore[override]
    def parse_binary_value(self, data: _T, display: Unused, length: Unused, format: Unused) -> tuple[_T, Literal[""]]: ...  # type: ignore[override]  # See: https://github.com/python-xlib/python-xlib/pull/249

class GetVersion(rq.ReplyRequest): ...

def get_version(self: Display | resource.Resource, major: int, minor: int) -> GetVersion: ...

class CreateContext(rq.Request): ...

def create_context(
    self: Display | resource.Resource,
    datum_flags: int,
    clients: Sequence[int],
    ranges: Sequence[
        tuple[
            tuple[int, int],
            tuple[int, int],
            tuple[int, int],
            tuple[int, int],
            tuple[int, int],
            tuple[int, int],
            tuple[int, int],
            bool,
            bool,
        ]
    ],
) -> int: ...

class RegisterClients(rq.Request): ...

def register_clients(
    self: Display | resource.Resource,
    context: int,
    element_header: int,
    clients: int,
    ranges: Sequence[
        tuple[
            tuple[int, int],
            tuple[int, int],
            tuple[int, int],
            tuple[int, int],
            tuple[int, int],
            tuple[int, int],
            tuple[int, int],
            bool,
            bool,
        ]
    ],
) -> None: ...

class UnregisterClients(rq.Request): ...

def unregister_clients(self: Display | resource.Resource, context: int, clients: Sequence[int]) -> None: ...

class GetContext(rq.ReplyRequest): ...

def get_context(self: Display | resource.Resource, context: int) -> GetContext: ...

class EnableContext(rq.ReplyRequest):
    def __init__(
        self,
        callback: Callable[[rq.DictWrapper | dict[str, Any]], Any],
        display: display.Display,
        defer: bool = False,
        *args: object | bool,
        **keys: object | bool,
    ) -> None: ...

def enable_context(
    self: Display | resource.Resource, context: int, callback: Callable[[rq.DictWrapper | dict[str, Any]], Any]
) -> None: ...

class DisableContext(rq.Request): ...

def disable_context(self: Display | resource.Resource, context: int) -> None: ...

class FreeContext(rq.Request): ...

def free_context(self: Display | resource.Resource, context: int) -> None: ...
def init(disp: Display, info: Unused) -> None: ...
