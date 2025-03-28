from typing import Final

from .connections import Connection as Connection
from .constants import FIELD_TYPE as FIELD_TYPE
from .converters import escape_dict as escape_dict, escape_sequence as escape_sequence, escape_string as escape_string
from .err import (
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
)
from .times import (
    Date as Date,
    DateFromTicks as DateFromTicks,
    Time as Time,
    TimeFromTicks as TimeFromTicks,
    Timestamp as Timestamp,
    TimestampFromTicks as TimestampFromTicks,
)

VERSION: Final[tuple[str | int, ...]]
VERSION_STRING: Final[str]

threadsafety: int
apilevel: str
paramstyle: str

class DBAPISet(frozenset[int]):
    def __ne__(self, other) -> bool: ...
    def __eq__(self, other) -> bool: ...
    def __hash__(self) -> int: ...

STRING: DBAPISet
BINARY: DBAPISet
NUMBER: DBAPISet
DATE: DBAPISet
TIME: DBAPISet
TIMESTAMP: DBAPISet
DATETIME: DBAPISet
ROWID: DBAPISet

def Binary(x) -> bytes: ...
def get_client_info() -> str: ...

__version__: str
version_info: tuple[int, int, int, str, int]
NULL: str

# pymysql/__init__.py says "Connect = connect = Connection = connections.Connection"
Connect = Connection
connect = Connection

def thread_safe() -> bool: ...
def install_as_MySQLdb() -> None: ...
