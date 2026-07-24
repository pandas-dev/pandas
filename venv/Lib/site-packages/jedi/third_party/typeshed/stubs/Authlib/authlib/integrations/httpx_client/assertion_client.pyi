from _typeshed import Incomplete

from authlib.oauth2.rfc7521 import AssertionClient as _AssertionClient

from ..base_client import OAuthError
from .oauth2_client import OAuth2Auth

__all__ = ["AsyncAssertionClient"]

# Inherits from httpx.AsyncClient
class AsyncAssertionClient(_AssertionClient):
    token_auth_class = OAuth2Auth
    oauth_error_class = OAuthError  # type: ignore[assignment]
    JWT_BEARER_GRANT_TYPE: Incomplete
    ASSERTION_METHODS: Incomplete
    DEFAULT_GRANT_TYPE: Incomplete
    def __init__(
        self,
        token_endpoint,
        issuer,
        subject,
        audience=None,
        grant_type=None,
        claims=None,
        token_placement="header",
        scope=None,
        **kwargs,
    ) -> None: ...
    async def request(self, method, url, withhold_token=False, auth=..., **kwargs): ...

# Inherits from httpx.Client
class AssertionClient(_AssertionClient):
    token_auth_class = OAuth2Auth
    oauth_error_class = OAuthError  # type: ignore[assignment]
    JWT_BEARER_GRANT_TYPE: Incomplete
    ASSERTION_METHODS: Incomplete
    DEFAULT_GRANT_TYPE: Incomplete
    def __init__(
        self,
        token_endpoint,
        issuer,
        subject,
        audience=None,
        grant_type=None,
        claims=None,
        token_placement="header",
        scope=None,
        **kwargs,
    ) -> None: ...
    def request(self, method, url, withhold_token=False, auth=..., **kwargs): ...
