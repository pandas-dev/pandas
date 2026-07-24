from _typeshed import Incomplete
from collections.abc import Generator
from typing import NoReturn
from typing_extensions import TypeAlias

from authlib.oauth2.auth import ClientAuth, TokenAuth
from authlib.oauth2.client import OAuth2Client as _OAuth2Client

from ..base_client import OAuthError

__all__ = ["OAuth2Auth", "OAuth2ClientAuth", "AsyncOAuth2Client", "OAuth2Client"]

_Response: TypeAlias = Incomplete  # actual type is httpx.Response
_Request: TypeAlias = Incomplete  # actual type is httpx.Request

# Inherits from httpx.Auth
class OAuth2Auth(TokenAuth):
    requires_request_body: bool
    def auth_flow(self, request: _Request) -> Generator[_Request, _Response]: ...

# Inherits from httpx.Auth
class OAuth2ClientAuth(ClientAuth):
    requires_request_body: bool
    def auth_flow(self, request: _Request) -> Generator[_Request, _Response]: ...

# Inherits from httpx.AsyncClient
class AsyncOAuth2Client(_OAuth2Client):
    SESSION_REQUEST_PARAMS: list[str]
    client_auth_class = OAuth2ClientAuth
    token_auth_class = OAuth2Auth
    oauth_error_class = OAuthError  # type: ignore[assignment]
    def __init__(
        self,
        client_id=None,
        client_secret=None,
        token_endpoint_auth_method=None,
        revocation_endpoint_auth_method=None,
        scope=None,
        redirect_uri=None,
        token=None,
        token_placement="header",
        update_token=None,
        leeway=60,
        **kwargs,
    ) -> None: ...
    async def request(self, method, url, withhold_token: bool = False, auth=..., **kwargs): ...
    async def stream(self, method, url, withhold_token: bool = False, auth=..., **kwargs) -> Generator[Incomplete]: ...
    async def ensure_active_token(self, token): ...  # type: ignore[override]

# Inherits from httpx.Client
class OAuth2Client(_OAuth2Client):
    SESSION_REQUEST_PARAMS: list[str]
    client_auth_class = OAuth2ClientAuth
    token_auth_class = OAuth2Auth
    oauth_error_class = OAuthError  # type: ignore[assignment]
    def __init__(
        self,
        client_id=None,
        client_secret=None,
        token_endpoint_auth_method=None,
        revocation_endpoint_auth_method=None,
        scope=None,
        redirect_uri=None,
        token=None,
        token_placement="header",
        update_token=None,
        **kwargs,
    ) -> None: ...
    @staticmethod
    def handle_error(error_type: str | None, error_description: str | None) -> NoReturn: ...
    def request(self, method, url, withhold_token: bool = False, auth=..., **kwargs): ...
    def stream(self, method, url, withhold_token: bool = False, auth=..., **kwargs): ...
