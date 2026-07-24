import gzip
from _typeshed import ReadableBuffer
from typing import Final, Protocol, type_check_only
from xmlrpc.client import ExpatParser, Unmarshaller

@type_check_only
class _Readable(Protocol):
    def read(self, size: int | None = -1) -> bytes: ...

__origin__: Final = "xmlrpc.client"
MAX_DATA: Final = 31457280

def defused_gzip_decode(data: ReadableBuffer, limit: int | None = None) -> bytes: ...

class DefusedGzipDecodedResponse(gzip.GzipFile):
    limit: int
    readlength: int | None
    def __init__(self, response: _Readable, limit: int | None = None) -> None: ...
    def read(self, n: int) -> bytes: ...  # type: ignore[override]
    def close(self) -> None: ...

class DefusedExpatParser(ExpatParser):
    forbid_dtd: bool
    forbid_entities: bool
    forbid_external: bool
    def __init__(
        self, target: Unmarshaller, forbid_dtd: bool = False, forbid_entities: bool = True, forbid_external: bool = True
    ) -> None: ...
    def defused_start_doctype_decl(self, name: str, sysid: str | None, pubid: str | None, has_internal_subset: bool) -> None: ...
    def defused_entity_decl(
        self,
        name: str,
        is_parameter_entity: bool,
        value: str | None,
        base: str | None,
        sysid: str,
        pubid: str | None,
        notation_name: str | None,
    ) -> None: ...
    def defused_unparsed_entity_decl(
        self, name: str, base: str | None, sysid: str, pubid: str | None, notation_name: str
    ) -> None: ...
    def defused_external_entity_ref_handler(
        self, context: str, base: str | None, sysid: str | None, pubid: str | None
    ) -> None: ...

def monkey_patch() -> None: ...
def unmonkey_patch() -> None: ...
