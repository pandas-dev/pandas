from _typeshed.wsgi import StartResponse, WSGIEnvironment
from collections.abc import Iterable
from http.client import HTTPConnection, HTTPMessage, HTTPSConnection
from typing import ClassVar

__all__ = ["send_request_app", "SendRequest"]

class SendRequest:
    def __init__(self, HTTPConnection: type[HTTPConnection] = ..., HTTPSConnection: type[HTTPSConnection] = ...) -> None: ...
    def __call__(self, environ: WSGIEnvironment, start_response: StartResponse) -> Iterable[bytes]: ...
    filtered_headers: ClassVar[tuple[str, ...]]
    def parse_headers(self, message: HTTPMessage) -> list[tuple[str, str]]: ...

send_request_app: SendRequest
