from _typeshed import Incomplete
from re import Pattern
from types import TracebackType
from typing import Any
from typing_extensions import LiteralString, Self, TypeAlias

from . import _mysql, cursors
from ._exceptions import (
    DatabaseError as DatabaseError,
    DataError as DataError,
    Error as Error,
    IntegrityError as IntegrityError,
    InterfaceError as InterfaceError,
    InternalError as InternalError,
    NotSupportedError as NotSupportedError,
    OperationalError as OperationalError,
    ProgrammingError as ProgrammingError,
    Warning as Warning,
)

# Any kind of object that can be passed to Connection.literal().
# The allowed types depend on the defined encoders, but the following
# types are always allowed.
_Literal: TypeAlias = str | bytearray | bytes | tuple[_Literal, ...] | list[_Literal] | Any

re_numeric_part: Pattern[str]

def numeric_part(s): ...

class Connection(_mysql.connection):
    default_cursor: type[cursors.Cursor]
    cursorclass: type[cursors.BaseCursor]
    encoders: Incomplete
    encoding: str
    messages: Incomplete
    def __init__(self, *args, **kwargs) -> None: ...
    def __enter__(self) -> Self: ...
    def __exit__(
        self, exc_type: type[BaseException] | None, exc_value: BaseException | None, traceback: TracebackType | None
    ) -> None: ...
    def autocommit(self, on: bool) -> None: ...
    def cursor(self, cursorclass: type[cursors.BaseCursor] | None = None): ...
    def query(self, query) -> None: ...
    def literal(self, o: _Literal) -> bytes: ...
    def begin(self) -> None: ...
    def warning_count(self): ...
    def set_character_set(self, charset: LiteralString, collation: LiteralString | None = None) -> None: ...
    def set_sql_mode(self, sql_mode) -> None: ...
    def show_warnings(self): ...
    Warning: type[BaseException]
    Error: type[BaseException]
    InterfaceError: type[BaseException]
    DatabaseError: type[BaseException]
    DataError: type[BaseException]
    OperationalError: type[BaseException]
    IntegrityError: type[BaseException]
    InternalError: type[BaseException]
    ProgrammingError: type[BaseException]
    NotSupportedError: type[BaseException]
