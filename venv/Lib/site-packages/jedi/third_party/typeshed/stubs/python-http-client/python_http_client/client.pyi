from email.message import Message
from http.client import HTTPResponse
from typing import Any, Final

class Response:
    def __init__(self, response: HTTPResponse) -> None: ...
    @property
    def status_code(self) -> int: ...
    @property
    def body(self) -> bytes: ...
    @property
    def headers(self) -> Message: ...
    @property
    def to_dict(self) -> dict[str, Any] | None: ...  # dict of response from API if body is not empty

class Client:
    methods: Final[set[str]]
    host: str
    request_headers: dict[str, str]
    append_slash: bool
    timeout: int
    def __init__(
        self,
        host: str,
        request_headers: dict[str, str] | None = None,
        version: int | None = None,
        url_path: list[str] | None = None,
        append_slash: bool = False,
        timeout: int | None = None,
    ) -> None: ...
    def _(self, name: str) -> Client: ...
    def __getattr__(self, name: str) -> Client | Response: ...
