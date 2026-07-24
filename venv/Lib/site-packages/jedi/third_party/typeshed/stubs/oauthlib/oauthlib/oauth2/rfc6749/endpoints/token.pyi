from _typeshed import Incomplete
from logging import Logger
from typing import Literal

from oauthlib.common import Request, _HTTPMethod

from .base import BaseEndpoint

log: Logger

class TokenEndpoint(BaseEndpoint):
    valid_request_methods: tuple[Literal["POST"]]
    def __init__(self, default_grant_type: str, default_token_type: str, grant_types: dict[str, Incomplete]) -> None: ...
    @property
    def grant_types(self) -> dict[str, Incomplete]: ...
    @property
    def default_grant_type(self) -> str: ...
    @property
    def default_grant_type_handler(self): ...
    @property
    def default_token_type(self) -> str: ...
    def create_token_response(
        self,
        uri: str,
        http_method: _HTTPMethod = "POST",
        body: str | dict[str, str] | list[tuple[str, str]] | None = None,
        headers: dict[str, str] | None = None,
        credentials=None,
        grant_type_for_scope=None,
        claims=None,
    ): ...
    def validate_token_request(self, request: Request) -> None: ...
