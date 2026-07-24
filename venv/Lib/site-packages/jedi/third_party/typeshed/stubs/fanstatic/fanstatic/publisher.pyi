from _typeshed import StrOrBytesPath
from _typeshed.wsgi import StartResponse, WSGIApplication, WSGIEnvironment
from collections.abc import Iterable
from typing import IO, Any, Literal

from fanstatic.core import Library
from fanstatic.registry import LibraryRegistry
from webob import Request, Response
from webob.dec import wsgify
from webob.static import DirectoryApp, FileApp

MINUTE_IN_SECONDS: Literal[60]
HOUR_IN_SECONDS: Literal[3600]
DAY_IN_SECONDS: Literal[86400]
YEAR_IN_SECONDS: int
FOREVER: int

class BundleApp(FileApp):
    filenames: list[str]
    def __init__(self, rootpath: str, bundle: IO[bytes], filenames: Iterable[StrOrBytesPath]) -> None: ...
    @wsgify
    def __call__(self, req: Request) -> Response: ...

class LibraryPublisher(DirectoryApp):
    ignores: list[str]
    library: Library
    cached_apps: dict[str, FileApp]
    def __init__(self, library: Library) -> None: ...
    @wsgify
    def __call__(self, req: Request) -> Response: ...

class Publisher:
    registry: LibraryRegistry
    directory_publishers: dict[str, LibraryPublisher]
    def __init__(self, registry: LibraryRegistry) -> None: ...
    @wsgify
    def __call__(self, request: Request) -> Response: ...

class Delegator:
    app: WSGIApplication
    publisher: Publisher
    publisher_signature: str
    trigger: str
    def __init__(self, app: WSGIApplication, publisher: Publisher, publisher_signature: str = "fanstatic") -> None: ...
    def is_resource(self, request: Request) -> bool: ...
    def __call__(self, environ: WSGIEnvironment, start_response: StartResponse) -> Iterable[bytes]: ...

def make_publisher(global_config: Any) -> Publisher: ...
