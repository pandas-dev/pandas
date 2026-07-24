from _typeshed import Unused
from collections.abc import Sequence
from typing import Final

from Xlib.display import Display
from Xlib.protocol import rq
from Xlib.xobject import resource

RES_MAJOR_VERSION: Final = 1
RES_MINOR_VERSION: Final = 2
extname: Final = "X-Resource"
ResQueryVersion: Final = 0
ResQueryClients: Final = 1
ResQueryClientResources: Final = 2
ResQueryClientPixmapBytes: Final = 3
ResQueryClientIds: Final = 4
ResQueryResourceBytes: Final = 5

class QueryVersion(rq.ReplyRequest): ...

def query_version(self: Display | resource.Resource, client_major: int = 1, client_minor: int = 2) -> QueryVersion: ...

Client: rq.Struct

class QueryClients(rq.ReplyRequest): ...

def query_clients(self: Display | resource.Resource) -> QueryClients: ...

Type: rq.Struct

class QueryClientResources(rq.ReplyRequest): ...

def query_client_resources(self: Display | resource.Resource, client: int) -> QueryClientResources: ...

class QueryClientPixmapBytes(rq.ReplyRequest): ...

def query_client_pixmap_bytes(self: Display | resource.Resource, client: int) -> QueryClientPixmapBytes: ...

class SizeOf(rq.LengthOf):
    item_size: int
    def __init__(self, name: str | list[str] | tuple[str, ...], size: int, item_size: int) -> None: ...
    def parse_value(self, length: int, display: Unused) -> int: ...  # type: ignore[override]

ClientXIDMask: Final = 0x1
LocalClientPIDMask: Final = 0x2
ClientIdSpec: rq.Struct
ClientIdValue: rq.Struct

class QueryClientIds(rq.ReplyRequest): ...

def query_client_ids(self: Display | resource.Resource, specs: Sequence[tuple[int, int]]) -> QueryClientIds: ...

ResourceIdSpec: rq.Struct
ResourceSizeSpec: rq.Struct
ResourceSizeValue: rq.Struct

class QueryResourceBytes(rq.ReplyRequest): ...

def query_resource_bytes(
    self: Display | resource.Resource, client: int, specs: Sequence[tuple[int, int]]
) -> QueryResourceBytes: ...
def init(disp: Display, info: Unused) -> None: ...
