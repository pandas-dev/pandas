from collections.abc import Sequence

from Xlib._typing import ErrorHandler
from Xlib.protocol import request
from Xlib.protocol.structs import _RGB3IntIterable
from Xlib.xobject import cursor, resource

class Fontable(resource.Resource):
    __fontable__ = resource.Resource.__resource__
    def query(self) -> request.QueryFont: ...
    def query_text_extents(self, string: str) -> request.QueryTextExtents: ...

class GC(Fontable):
    __gc__ = resource.Resource.__resource__
    def change(self, onerror: ErrorHandler[object] | None = None, **keys: object) -> None: ...
    def copy(self, src_gc: int, mask: int, onerror: ErrorHandler[object] | None = None) -> None: ...
    def set_dashes(self, offset: int, dashes: Sequence[int], onerror: ErrorHandler[object] | None = None) -> None: ...
    def set_clip_rectangles(
        self,
        x_origin: int,
        y_origin: int,
        rectangles: Sequence[dict[str, int]],
        ordering: int,
        onerror: ErrorHandler[object] | None = None,
    ) -> None: ...
    def free(self, onerror: ErrorHandler[object] | None = None) -> None: ...

class Font(Fontable):
    __font__ = resource.Resource.__resource__
    def close(self, onerror: ErrorHandler[object] | None = None) -> None: ...
    def create_glyph_cursor(
        self, mask: Font, source_char: int, mask_char: int, foreground: _RGB3IntIterable, background: _RGB3IntIterable
    ) -> cursor.Cursor: ...
