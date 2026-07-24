from typing import Final

from Xlib._typing import ErrorHandler
from Xlib.display import Display
from Xlib.protocol import request, rq
from Xlib.xobject import drawable

extname: Final = "MIT-SCREEN-SAVER"
NotifyMask: Final = 1
CycleMask: Final = 2
StateOff: Final = 0
StateOn: Final = 1
StateCycle: Final = 2
KindBlanked: Final = 0
KindInternal: Final = 1
KindExternal: Final = 2

class QueryVersion(rq.ReplyRequest): ...

def query_version(self: drawable.Drawable) -> QueryVersion: ...

class QueryInfo(rq.ReplyRequest): ...

def query_info(self: drawable.Drawable) -> QueryInfo: ...

class SelectInput(rq.Request): ...

def select_input(self: drawable.Drawable, mask: int) -> SelectInput: ...

class SetAttributes(rq.Request): ...

def set_attributes(
    self: drawable.Drawable,
    x: int,
    y: int,
    width: int,
    height: int,
    border_width: int,
    window_class: int = 0,
    depth: int = 0,
    visual: int = 0,
    onerror: ErrorHandler[object] | None = None,
    **keys: object,
) -> SetAttributes: ...

class UnsetAttributes(rq.Request): ...

def unset_attributes(self: drawable.Drawable, onerror: ErrorHandler[object] | None = None) -> UnsetAttributes: ...

class Notify(rq.Event): ...

def init(disp: Display, info: request.QueryExtension) -> None: ...
