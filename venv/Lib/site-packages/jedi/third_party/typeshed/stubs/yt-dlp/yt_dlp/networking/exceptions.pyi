from ..utils import YoutubeDLError
from .common import RequestHandler, Response

class RequestError(YoutubeDLError):
    handler: RequestHandler | None
    cause: Exception | str | None
    def __init__(
        self, msg: str | None = None, cause: Exception | str | None = None, handler: RequestHandler | None = None
    ) -> None: ...

class UnsupportedRequest(RequestError): ...

class NoSupportingHandlers(RequestError):
    unsupported_errors: list[UnsupportedRequest]
    unexpected_errors: list[UnsupportedRequest]
    def __init__(self, unsupported_errors: list[UnsupportedRequest], unexpected_errors: list[Exception]) -> None: ...

class TransportError(RequestError): ...

class HTTPError(RequestError):
    response: Response
    status: int
    reason: str | None
    redirect_loop: bool

    def __init__(self, response: Response, redirect_loop: bool = False) -> None: ...
    def close(self) -> None: ...

class IncompleteRead(TransportError):
    def __init__(
        self,
        partial: int,
        expected: int | None = None,
        *,
        cause: Exception | str | None = None,
        handler: RequestHandler | None = None,
    ) -> None: ...

class SSLError(TransportError): ...
class CertificateVerifyError(SSLError): ...
class ProxyError(TransportError): ...

network_exceptions: tuple[type[RequestError], ...]
