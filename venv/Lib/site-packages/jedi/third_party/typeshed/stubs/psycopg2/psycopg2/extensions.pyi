from _typeshed import Unused
from collections.abc import Callable, Iterable
from typing import Any, TypeVar, overload

from psycopg2._psycopg import (
    BINARYARRAY as BINARYARRAY,
    BOOLEAN as BOOLEAN,
    BOOLEANARRAY as BOOLEANARRAY,
    BYTES as BYTES,
    BYTESARRAY as BYTESARRAY,
    DATE as DATE,
    DATEARRAY as DATEARRAY,
    DATETIMEARRAY as DATETIMEARRAY,
    DECIMAL as DECIMAL,
    DECIMALARRAY as DECIMALARRAY,
    FLOAT as FLOAT,
    FLOATARRAY as FLOATARRAY,
    INTEGER as INTEGER,
    INTEGERARRAY as INTEGERARRAY,
    INTERVAL as INTERVAL,
    INTERVALARRAY as INTERVALARRAY,
    LONGINTEGER as LONGINTEGER,
    LONGINTEGERARRAY as LONGINTEGERARRAY,
    PYDATE as PYDATE,
    PYDATEARRAY as PYDATEARRAY,
    PYDATETIME as PYDATETIME,
    PYDATETIMEARRAY as PYDATETIMEARRAY,
    PYDATETIMETZ as PYDATETIMETZ,
    PYDATETIMETZARRAY as PYDATETIMETZARRAY,
    PYINTERVAL as PYINTERVAL,
    PYINTERVALARRAY as PYINTERVALARRAY,
    PYTIME as PYTIME,
    PYTIMEARRAY as PYTIMEARRAY,
    ROWIDARRAY as ROWIDARRAY,
    STRINGARRAY as STRINGARRAY,
    TIME as TIME,
    TIMEARRAY as TIMEARRAY,
    UNICODE as UNICODE,
    UNICODEARRAY as UNICODEARRAY,
    AsIs as AsIs,
    Binary as Binary,
    Boolean as Boolean,
    Column as Column,
    ConnectionInfo as ConnectionInfo,
    DateFromPy as DateFromPy,
    Diagnostics as Diagnostics,
    Float as Float,
    Int as Int,
    IntervalFromPy as IntervalFromPy,
    ISQLQuote as ISQLQuote,
    Notify as Notify,
    QueryCanceledError as QueryCanceledError,
    QuotedString as QuotedString,
    TimeFromPy as TimeFromPy,
    TimestampFromPy as TimestampFromPy,
    TransactionRollbackError as TransactionRollbackError,
    Xid as Xid,
    _ISQLQuoteProto,
    _type,
    adapt as adapt,
    adapters as adapters,
    binary_types as binary_types,
    connection as connection,
    cursor as cursor,
    encodings as encodings,
    encrypt_password as encrypt_password,
    get_wait_callback as get_wait_callback,
    libpq_version as libpq_version,
    lobject as lobject,
    new_array_type as new_array_type,
    new_type as new_type,
    parse_dsn as parse_dsn,
    quote_ident as quote_ident,
    register_type as register_type,
    set_wait_callback as set_wait_callback,
    string_types as string_types,
)

ISOLATION_LEVEL_AUTOCOMMIT: int
ISOLATION_LEVEL_READ_UNCOMMITTED: int
ISOLATION_LEVEL_READ_COMMITTED: int
ISOLATION_LEVEL_REPEATABLE_READ: int
ISOLATION_LEVEL_SERIALIZABLE: int
ISOLATION_LEVEL_DEFAULT: Any
STATUS_SETUP: int
STATUS_READY: int
STATUS_BEGIN: int
STATUS_SYNC: int
STATUS_ASYNC: int
STATUS_PREPARED: int
STATUS_IN_TRANSACTION: int
POLL_OK: int
POLL_READ: int
POLL_WRITE: int
POLL_ERROR: int
TRANSACTION_STATUS_IDLE: int
TRANSACTION_STATUS_ACTIVE: int
TRANSACTION_STATUS_INTRANS: int
TRANSACTION_STATUS_INERROR: int
TRANSACTION_STATUS_UNKNOWN: int

_T = TypeVar("_T")

def register_adapter(typ: type[_T], callable: Callable[[_T], _ISQLQuoteProto]) -> None: ...

class SQL_IN:
    def __init__(self, seq: Iterable[object]) -> None: ...
    def prepare(self, conn: connection | None) -> None: ...
    def getquoted(self) -> bytes: ...

class NoneAdapter:
    def __init__(self, obj: Unused) -> None: ...
    def getquoted(self, _null: bytes = b"NULL") -> bytes: ...

@overload
def make_dsn(dsn: bytes) -> bytes: ...  # type: ignore[overload-overlap]
@overload
def make_dsn(dsn: None = None) -> str: ...
@overload
def make_dsn(dsn: str | bytes | None = None, **kwargs: Any) -> str: ...

JSON: _type
JSONARRAY: _type | None
JSONB: _type
JSONBARRAY: _type | None
