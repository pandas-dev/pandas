from typing_extensions import Buffer

from gunicorn.http import Message

class ParseException(Exception): ...

class NoMoreData(IOError):
    buf: Buffer

    def __init__(self, buf: Buffer | None = None) -> None: ...

class ConfigurationProblem(ParseException):
    info: str
    code: int

    def __init__(self, info: str) -> None: ...

class InvalidRequestLine(ParseException):
    req: str
    code: int

    def __init__(self, req: str) -> None: ...

class InvalidRequestMethod(ParseException):
    method: str

    def __init__(self, method: str) -> None: ...

class ExpectationFailed(ParseException):
    expect: str

    def __init__(self, expect: str) -> None: ...

class InvalidHTTPVersion(ParseException):
    version: str | tuple[int, int]

    def __init__(self, version: str | tuple[int, int]) -> None: ...

class InvalidHeader(ParseException):
    hdr: str
    req: Message | None

    def __init__(self, hdr: str, req: Message | None = None) -> None: ...

class ObsoleteFolding(ParseException):
    hdr: str

    def __init__(self, hdr: str) -> None: ...

class InvalidHeaderName(ParseException):
    hdr: str

    def __init__(self, hdr: str) -> None: ...

class UnsupportedTransferCoding(ParseException):
    hdr: str
    code: int

    def __init__(self, hdr: str) -> None: ...

class InvalidChunkSize(IOError):
    data: bytes

    def __init__(self, data: bytes) -> None: ...

class ChunkMissingTerminator(IOError):
    term: bytes

    def __init__(self, term: bytes) -> None: ...

class InvalidChunkExtension(IOError):
    reason: str
    def __init__(self, reason: str) -> None: ...

class LimitRequestLine(ParseException):
    size: int
    max_size: int | None

    def __init__(self, size: int, max_size: int | None = None) -> None: ...

class LimitRequestHeaders(ParseException):
    msg: str

    def __init__(self, msg: str) -> None: ...

class InvalidProxyLine(ParseException):
    line: str
    code: int

    def __init__(self, line: str) -> None: ...

class InvalidProxyHeader(ParseException):
    msg: str
    code: int

    def __init__(self, msg: str) -> None: ...

class ForbiddenProxyRequest(ParseException):
    host: str
    code: int

    def __init__(self, host: str) -> None: ...

class InvalidSchemeHeaders(ParseException): ...
