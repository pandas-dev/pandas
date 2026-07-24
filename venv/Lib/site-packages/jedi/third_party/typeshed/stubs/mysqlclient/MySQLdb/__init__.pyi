from _typeshed import Incomplete

from MySQLdb import connections as connections, constants as constants, converters as converters, cursors as cursors
from MySQLdb._mysql import (
    DatabaseError as DatabaseError,
    DataError as DataError,
    Error as Error,
    IntegrityError as IntegrityError,
    InterfaceError as InterfaceError,
    InternalError as InternalError,
    MySQLError as MySQLError,
    NotSupportedError as NotSupportedError,
    OperationalError as OperationalError,
    ProgrammingError as ProgrammingError,
    Warning as Warning,
    debug as debug,
    get_client_info as get_client_info,
    string_literal as string_literal,
)
from MySQLdb.connections import Connection as Connection
from MySQLdb.constants import FIELD_TYPE as FIELD_TYPE
from MySQLdb.release import version_info as version_info
from MySQLdb.times import (
    Date as Date,
    DateFromTicks as DateFromTicks,
    Time as Time,
    TimeFromTicks as TimeFromTicks,
    Timestamp as Timestamp,
    TimestampFromTicks as TimestampFromTicks,
)

threadsafety: int
apilevel: str
paramstyle: str

class DBAPISet(frozenset[Incomplete]):
    def __eq__(self, other): ...

STRING: Incomplete
BINARY: Incomplete
NUMBER: Incomplete
DATE: Incomplete
TIME: Incomplete
TIMESTAMP: Incomplete
DATETIME: Incomplete
ROWID: Incomplete

def Binary(x): ...
def Connect(*args, **kwargs) -> Connection: ...

connect = Connect

__all__ = [
    "BINARY",
    "Binary",
    "Connect",
    "Connection",
    "DATE",
    "Date",
    "Time",
    "Timestamp",
    "DateFromTicks",
    "TimeFromTicks",
    "TimestampFromTicks",
    "DataError",
    "DatabaseError",
    "Error",
    "FIELD_TYPE",
    "IntegrityError",
    "InterfaceError",
    "InternalError",
    "MySQLError",
    "NUMBER",
    "NotSupportedError",
    "DBAPISet",
    "OperationalError",
    "ProgrammingError",
    "ROWID",
    "STRING",
    "TIME",
    "TIMESTAMP",
    "Warning",
    "apilevel",
    "connect",
    "connections",
    "constants",
    "converters",
    "cursors",
    "debug",
    "get_client_info",
    "paramstyle",
    "string_literal",
    "threadsafety",
    "version_info",
]
