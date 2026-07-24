from _typeshed import StrPath
from _typeshed.wsgi import WSGIApplication
from collections.abc import Iterator
from typing import IO, Any

from webob.dec import wsgify
from webob.request import Request
from webob.response import Response

__all__ = ["FileApp", "DirectoryApp"]

BLOCK_SIZE: int

class FileApp:
    filename: StrPath
    kw: dict[str, Any]
    def __init__(self, filename: StrPath, **kw: Any) -> None: ...
    @wsgify
    def __call__(self, req: Request) -> WSGIApplication: ...

class FileIter:
    file: IO[bytes]
    def __init__(self, file: IO[bytes]) -> None: ...
    def app_iter_range(
        self, seek: int | None = None, limit: int | None = None, block_size: int | None = None
    ) -> Iterator[bytes]: ...
    __iter__ = app_iter_range

class DirectoryApp:
    path: StrPath
    index_page: str
    hide_index_with_redirect: bool
    fileapp_kw: dict[str, Any]
    def __init__(
        self, path: StrPath, index_page: str = "index.html", hide_index_with_redirect: bool = False, **kw: Any
    ) -> None: ...
    def make_fileapp(self, path: StrPath) -> FileApp: ...
    @wsgify
    def __call__(self, req: Request) -> Response | FileApp: ...
    def index(self, req: Request, path: StrPath) -> Response | FileApp: ...
