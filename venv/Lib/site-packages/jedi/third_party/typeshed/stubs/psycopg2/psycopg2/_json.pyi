from collections.abc import Callable
from typing import Any
from typing_extensions import Self

from psycopg2._psycopg import _type, connection, cursor

JSON_OID: int
JSONARRAY_OID: int
JSONB_OID: int
JSONBARRAY_OID: int

class Json:
    adapted: Any
    def __init__(self, adapted: Any, dumps: Callable[..., str] | None = None) -> None: ...
    def __conform__(self, proto) -> Self | None: ...
    def dumps(self, obj: Any) -> str: ...
    def prepare(self, conn: connection | None) -> None: ...
    def getquoted(self) -> bytes: ...

def register_json(
    conn_or_curs: connection | cursor | None = None,
    globally: bool = False,
    loads: Callable[..., Any] | None = None,
    oid: int | None = None,
    array_oid: int | None = None,
    name: str = "json",
) -> tuple[_type, _type | None]: ...
def register_default_json(
    conn_or_curs: connection | cursor | None = None, globally: bool = False, loads: Callable[..., Any] | None = None
) -> tuple[_type, _type | None]: ...
def register_default_jsonb(
    conn_or_curs: connection | cursor | None = None, globally: bool = False, loads: Callable[..., Any] | None = None
) -> tuple[_type, _type | None]: ...
