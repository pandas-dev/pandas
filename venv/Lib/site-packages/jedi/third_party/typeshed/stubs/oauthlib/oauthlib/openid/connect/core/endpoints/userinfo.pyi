from collections.abc import Mapping
from logging import Logger

from oauthlib.common import Request, _HTTPMethod
from oauthlib.oauth2.rfc6749.endpoints.base import BaseEndpoint
from oauthlib.oauth2.rfc6749.request_validator import RequestValidator as OAuth2RequestValidator
from oauthlib.oauth2.rfc6749.tokens import BearerToken

log: Logger

class UserInfoEndpoint(BaseEndpoint):
    bearer: BearerToken
    request_validator: OAuth2RequestValidator
    def __init__(self, request_validator: OAuth2RequestValidator) -> None: ...
    def create_userinfo_response(
        self,
        uri: str,
        http_method: _HTTPMethod = "GET",
        body: str | dict[str, str] | list[tuple[str, str]] | None = None,
        headers: Mapping[str, str] | None = None,
    ) -> tuple[dict[str, str], str, int]: ...
    def validate_userinfo_request(self, request: Request) -> None: ...
