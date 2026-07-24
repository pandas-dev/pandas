from _typeshed import Incomplete

from authlib.oauth2.auth import ClientAuth, TokenAuth
from authlib.oauth2.client import OAuth2Client

from ..base_client import OAuthError

__all__ = ["OAuth2Session", "OAuth2Auth"]

# Inherits from requests.auth.AuthBase
class OAuth2Auth(TokenAuth):
    def ensure_active_token(self) -> None: ...
    def __call__(self, req): ...

# Inherits from requests.auth.AuthBase
class OAuth2ClientAuth(ClientAuth):
    def __call__(self, req): ...

# Inherits from requests.Session
class OAuth2Session(OAuth2Client):
    client_auth_class = OAuth2ClientAuth
    token_auth_class = OAuth2Auth
    oauth_error_class = OAuthError  # type: ignore[assignment]
    SESSION_REQUEST_PARAMS: tuple[str, ...]  # type: ignore[assignment]
    default_timeout: Incomplete
    def __init__(
        self,
        client_id=None,
        client_secret=None,
        token_endpoint_auth_method=None,
        revocation_endpoint_auth_method=None,
        scope=None,
        state=None,
        redirect_uri=None,
        token=None,
        token_placement="header",
        update_token=None,
        leeway=60,
        default_timeout=None,
        **kwargs,
    ) -> None: ...
    def fetch_access_token(self, url=None, **kwargs): ...
    def request(self, method, url, withhold_token=False, auth=None, **kwargs): ...
