from Xlib._typing import ErrorHandler
from Xlib.protocol.structs import _RGB3IntIterable
from Xlib.xobject import resource

class Cursor(resource.Resource):
    __cursor__ = resource.Resource.__resource__
    def free(self, onerror: ErrorHandler[object] | None = None) -> None: ...
    def recolor(
        self, foreground: _RGB3IntIterable, background: _RGB3IntIterable, onerror: ErrorHandler[object] | None = None
    ) -> None: ...
