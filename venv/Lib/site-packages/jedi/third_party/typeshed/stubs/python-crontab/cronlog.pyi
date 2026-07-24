from _typeshed import StrOrBytesPath
from codecs import StreamReaderWriter
from collections.abc import Generator, Iterator
from datetime import datetime
from types import TracebackType
from typing_extensions import Self

MATCHER: str

class LogReader:
    filename: StrOrBytesPath
    mass: int
    size: int
    read: int
    pipe: StreamReaderWriter | None
    def __init__(self, filename: StrOrBytesPath, mass: int = ...) -> None: ...
    def __enter__(self) -> Self: ...
    def __exit__(
        self, error_type: type[BaseException] | None, value: BaseException | None, traceback: TracebackType | None
    ) -> None: ...
    def __iter__(self) -> Iterator[str]: ...
    def readlines(self, until: int = ...) -> Generator[tuple[int, str]]: ...

def cron_date_to_datetime(cron_str: str) -> datetime: ...

class CronLog(LogReader):
    user: str | None
    def __init__(self, filename: StrOrBytesPath = ..., user: str | None = ...) -> None: ...
    def for_program(self, command: str) -> ProgramLog: ...
    def __iter__(self) -> dict[str, str | None]: ...  # type: ignore[override]

class ProgramLog:
    log: CronLog
    command: str
    def __init__(self, log: CronLog, command: str) -> None: ...
    def __iter__(self) -> dict[str, str | None]: ...
