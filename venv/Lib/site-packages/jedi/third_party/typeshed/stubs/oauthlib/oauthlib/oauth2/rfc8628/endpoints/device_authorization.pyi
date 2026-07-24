from collections.abc import Callable
from logging import Logger

from oauthlib.common import Request, _HTTPMethod
from oauthlib.oauth2.rfc6749.endpoints.base import BaseEndpoint
from oauthlib.openid.connect.core.request_validator import RequestValidator

log: Logger

class DeviceAuthorizationEndpoint(BaseEndpoint):
    request_validator: RequestValidator
    user_code_generator: Callable[[None], str] | None
    def __init__(
        self,
        request_validator: RequestValidator,
        verification_uri: str,
        expires_in: int = 1800,
        interval: int | None = None,
        verification_uri_complete: str | None = None,
        user_code_generator: Callable[[None], str] | None = None,
    ) -> None: ...
    @property
    def interval(self) -> int | None: ...
    @property
    def expires_in(self) -> int: ...
    @property
    def verification_uri(self) -> str: ...
    def verification_uri_complete(self, user_code: str) -> str | None: ...
    def validate_device_authorization_request(self, request: Request) -> None: ...
    def create_device_authorization_response(
        self, uri: str, http_method: _HTTPMethod = "POST", body: str | None = None, headers: dict[str, str] | None = None
    ) -> tuple[dict[str, str], dict[str, str | int], int]: ...
