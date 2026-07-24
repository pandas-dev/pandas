from logging import Logger
from typing import Literal

from oauthlib.common import Request, _HTTPMethod

from ..request_validator import RequestValidator
from .base import BaseEndpoint

log: Logger

class IntrospectEndpoint(BaseEndpoint):
    valid_token_types: tuple[Literal["access_token"], Literal["refresh_token"]]
    valid_request_methods: tuple[Literal["POST"]]
    request_validator: RequestValidator
    supported_token_types: tuple[str, ...]
    def __init__(self, request_validator: RequestValidator, supported_token_types: tuple[str, ...] | None = None) -> None: ...
    def create_introspect_response(
        self, uri: str, http_method: _HTTPMethod = "POST", body: str | None = None, headers: dict[str, str] | None = None
    ) -> tuple[dict[str, str], str, int]: ...
    def validate_introspect_request(self, request: Request) -> None: ...
