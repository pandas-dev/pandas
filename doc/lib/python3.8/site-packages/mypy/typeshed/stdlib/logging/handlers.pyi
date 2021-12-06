import datetime
import http.client
import ssl
import sys
from _typeshed import StrPath
from collections.abc import Callable
from logging import FileHandler, Handler, LogRecord
from socket import SocketKind, SocketType
from typing import Any, ClassVar, Optional, Pattern, Union

if sys.version_info >= (3, 7):
    from queue import Queue, SimpleQueue
else:
    from queue import Queue

DEFAULT_TCP_LOGGING_PORT: int
DEFAULT_UDP_LOGGING_PORT: int
DEFAULT_HTTP_LOGGING_PORT: int
DEFAULT_SOAP_LOGGING_PORT: int
SYSLOG_UDP_PORT: int
SYSLOG_TCP_PORT: int

class WatchedFileHandler(FileHandler):
    dev: int  # undocumented
    ino: int  # undocumented
    if sys.version_info >= (3, 9):
        def __init__(
            self,
            filename: StrPath,
            mode: str = ...,
            encoding: Optional[str] = ...,
            delay: bool = ...,
            errors: Optional[str] = ...,
        ) -> None: ...
    else:
        def __init__(self, filename: StrPath, mode: str = ..., encoding: Optional[str] = ..., delay: bool = ...) -> None: ...
    def _statstream(self) -> None: ...  # undocumented
    def reopenIfNeeded(self) -> None: ...

class BaseRotatingHandler(FileHandler):
    namer: Optional[Callable[[str], str]]
    rotator: Optional[Callable[[str, str], None]]
    if sys.version_info >= (3, 9):
        def __init__(
            self, filename: StrPath, mode: str, encoding: Optional[str] = ..., delay: bool = ..., errors: Optional[str] = ...
        ) -> None: ...
    else:
        def __init__(self, filename: StrPath, mode: str, encoding: Optional[str] = ..., delay: bool = ...) -> None: ...
    def rotation_filename(self, default_name: str) -> str: ...
    def rotate(self, source: str, dest: str) -> None: ...

class RotatingFileHandler(BaseRotatingHandler):
    maxBytes: str  # undocumented
    backupCount: str  # undocumented
    if sys.version_info >= (3, 9):
        def __init__(
            self,
            filename: StrPath,
            mode: str = ...,
            maxBytes: int = ...,
            backupCount: int = ...,
            encoding: Optional[str] = ...,
            delay: bool = ...,
            errors: Optional[str] = ...,
        ) -> None: ...
    else:
        def __init__(
            self,
            filename: StrPath,
            mode: str = ...,
            maxBytes: int = ...,
            backupCount: int = ...,
            encoding: Optional[str] = ...,
            delay: bool = ...,
        ) -> None: ...
    def doRollover(self) -> None: ...
    def shouldRollover(self, record: LogRecord) -> int: ...  # undocumented

class TimedRotatingFileHandler(BaseRotatingHandler):
    when: str  # undocumented
    backupCount: str  # undocumented
    utc: bool  # undocumented
    atTime: Optional[datetime.datetime]  # undocumented
    interval: int  # undocumented
    suffix: str  # undocumented
    dayOfWeek: int  # undocumented
    rolloverAt: int  # undocumented
    extMatch: Pattern[str]  # undocumented
    if sys.version_info >= (3, 9):
        def __init__(
            self,
            filename: StrPath,
            when: str = ...,
            interval: int = ...,
            backupCount: int = ...,
            encoding: Optional[str] = ...,
            delay: bool = ...,
            utc: bool = ...,
            atTime: Optional[datetime.datetime] = ...,
            errors: Optional[str] = ...,
        ) -> None: ...
    else:
        def __init__(
            self,
            filename: StrPath,
            when: str = ...,
            interval: int = ...,
            backupCount: int = ...,
            encoding: Optional[str] = ...,
            delay: bool = ...,
            utc: bool = ...,
            atTime: Optional[datetime.datetime] = ...,
        ) -> None: ...
    def doRollover(self) -> None: ...
    def shouldRollover(self, record: LogRecord) -> int: ...  # undocumented
    def computeRollover(self, currentTime: int) -> int: ...  # undocumented
    def getFilesToDelete(self) -> list[str]: ...  # undocumented

class SocketHandler(Handler):
    host: str  # undocumented
    port: Optional[int]  # undocumented
    address: Union[tuple[str, int], str]  # undocumented
    sock: Optional[SocketType]  # undocumented
    closeOnError: bool  # undocumented
    retryTime: Optional[float]  # undocumented
    retryStart: float  # undocumented
    retryFactor: float  # undocumented
    retryMax: float  # undocumented
    def __init__(self, host: str, port: Optional[int]) -> None: ...
    def makeSocket(self, timeout: float = ...) -> SocketType: ...  # timeout is undocumented
    def makePickle(self, record: LogRecord) -> bytes: ...
    def send(self, s: bytes) -> None: ...
    def createSocket(self) -> None: ...

class DatagramHandler(SocketHandler):
    def makeSocket(self) -> SocketType: ...  # type: ignore

class SysLogHandler(Handler):
    LOG_EMERG: int
    LOG_ALERT: int
    LOG_CRIT: int
    LOG_ERR: int
    LOG_WARNING: int
    LOG_NOTICE: int
    LOG_INFO: int
    LOG_DEBUG: int

    LOG_KERN: int
    LOG_USER: int
    LOG_MAIL: int
    LOG_DAEMON: int
    LOG_AUTH: int
    LOG_SYSLOG: int
    LOG_LPR: int
    LOG_NEWS: int
    LOG_UUCP: int
    LOG_CRON: int
    LOG_AUTHPRIV: int
    LOG_FTP: int

    if sys.version_info >= (3, 9):
        LOG_NTP: int
        LOG_SECURITY: int
        LOG_CONSOLE: int
        LOG_SOLCRON: int

    LOG_LOCAL0: int
    LOG_LOCAL1: int
    LOG_LOCAL2: int
    LOG_LOCAL3: int
    LOG_LOCAL4: int
    LOG_LOCAL5: int
    LOG_LOCAL6: int
    LOG_LOCAL7: int
    unixsocket: bool  # undocumented
    socktype: SocketKind  # undocumented
    ident: str  # undocumented
    append_nul: bool  # undocumented
    facility: int  # undocumented
    priority_names: ClassVar[dict[str, int]]  # undocumented
    facility_names: ClassVar[dict[str, int]]  # undocumented
    priority_map: ClassVar[dict[str, str]]  # undocumented
    def __init__(
        self, address: Union[tuple[str, int], str] = ..., facility: int = ..., socktype: Optional[SocketKind] = ...
    ) -> None: ...
    def encodePriority(self, facility: Union[int, str], priority: Union[int, str]) -> int: ...
    def mapPriority(self, levelName: str) -> str: ...

class NTEventLogHandler(Handler):
    def __init__(self, appname: str, dllname: Optional[str] = ..., logtype: str = ...) -> None: ...
    def getEventCategory(self, record: LogRecord) -> int: ...
    # TODO correct return value?
    def getEventType(self, record: LogRecord) -> int: ...
    def getMessageID(self, record: LogRecord) -> int: ...

class SMTPHandler(Handler):
    mailhost: str  # undocumented
    mailport: Optional[int]  # undocumented
    username: Optional[str]  # undocumented
    # password only exists as an attribute if passed credentials is a tuple or list
    password: str  # undocumented
    fromaddr: str  # undocumented
    toaddrs: list[str]  # undocumented
    subject: str  # undocumented
    secure: Union[tuple[()], tuple[str], tuple[str, str], None]  # undocumented
    timeout: float  # undocumented
    def __init__(
        self,
        mailhost: Union[str, tuple[str, int]],
        fromaddr: str,
        toaddrs: Union[str, list[str]],
        subject: str,
        credentials: Optional[tuple[str, str]] = ...,
        secure: Union[tuple[()], tuple[str], tuple[str, str], None] = ...,
        timeout: float = ...,
    ) -> None: ...
    def getSubject(self, record: LogRecord) -> str: ...

class BufferingHandler(Handler):
    capacity: int  # undocumented
    buffer: list[LogRecord]  # undocumented
    def __init__(self, capacity: int) -> None: ...
    def shouldFlush(self, record: LogRecord) -> bool: ...

class MemoryHandler(BufferingHandler):
    flushLevel: int  # undocumented
    target: Optional[Handler]  # undocumented
    flushOnClose: bool  # undocumented
    def __init__(
        self, capacity: int, flushLevel: int = ..., target: Optional[Handler] = ..., flushOnClose: bool = ...
    ) -> None: ...
    def setTarget(self, target: Optional[Handler]) -> None: ...

class HTTPHandler(Handler):
    host: str  # undocumented
    url: str  # undocumented
    method: str  # undocumented
    secure: bool  # undocumented
    credentials: Optional[tuple[str, str]]  # undocumented
    context: Optional[ssl.SSLContext]  # undocumented
    def __init__(
        self,
        host: str,
        url: str,
        method: str = ...,
        secure: bool = ...,
        credentials: Optional[tuple[str, str]] = ...,
        context: Optional[ssl.SSLContext] = ...,
    ) -> None: ...
    def mapLogRecord(self, record: LogRecord) -> dict[str, Any]: ...
    if sys.version_info >= (3, 9):
        def getConnection(self, host: str, secure: bool) -> http.client.HTTPConnection: ...  # undocumented

class QueueHandler(Handler):
    if sys.version_info >= (3, 7):
        queue: Union[SimpleQueue[Any], Queue[Any]]  # undocumented
        def __init__(self, queue: Union[SimpleQueue[Any], Queue[Any]]) -> None: ...
    else:
        queue: Queue[Any]  # undocumented
        def __init__(self, queue: Queue[Any]) -> None: ...
    def prepare(self, record: LogRecord) -> Any: ...
    def enqueue(self, record: LogRecord) -> None: ...

class QueueListener:
    handlers: tuple[Handler]  # undocumented
    respect_handler_level: bool  # undocumented
    if sys.version_info >= (3, 7):
        queue: Union[SimpleQueue[Any], Queue[Any]]  # undocumented
        def __init__(
            self, queue: Union[SimpleQueue[Any], Queue[Any]], *handlers: Handler, respect_handler_level: bool = ...
        ) -> None: ...
    else:
        queue: Queue[Any]  # undocumented
        def __init__(self, queue: Queue[Any], *handlers: Handler, respect_handler_level: bool = ...) -> None: ...
    def dequeue(self, block: bool) -> LogRecord: ...
    def prepare(self, record: LogRecord) -> Any: ...
    def start(self) -> None: ...
    def stop(self) -> None: ...
    def enqueue_sentinel(self) -> None: ...
    def handle(self, record: LogRecord) -> None: ...
