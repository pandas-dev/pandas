from _typeshed import Incomplete
from typing import ClassVar, Literal

import _win32typing

def odbc(connectionString: str, /) -> _win32typing.connection: ...
def SQLDataSources(direction, /) -> tuple[Incomplete, Incomplete]: ...

DATE: str
NUMBER: str
RAW: str
SQL_FETCH_ABSOLUTE: int
SQL_FETCH_FIRST: int
SQL_FETCH_FIRST_SYSTEM: int
SQL_FETCH_FIRST_USER: int
SQL_FETCH_LAST: int
SQL_FETCH_NEXT: int
SQL_FETCH_PRIOR: int
SQL_FETCH_RELATIVE: int
STRING: str
TYPES: tuple[Literal["STRING"], Literal["RAW"], Literal["NUMBER"], Literal["DATE"]]

class error(Exception):
    __name__: ClassVar[str] = "odbcError"

# These all pretend to come from a module called "dbi", but that module doesn't exist
class dataError(Exception): ...
class integrityError(Exception): ...
class internalError(Exception): ...
class noError(Exception): ...
class opError(Exception): ...
class progError(Exception): ...
