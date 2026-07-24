from _typeshed import Incomplete

from authlib.oauth2.rfc7521 import AssertionClient

from .oauth2_session import OAuth2Auth

class AssertionAuth(OAuth2Auth):
    def ensure_active_token(self): ...

# Inherits from requests.Session
class AssertionSession(AssertionClient):
    token_auth_class = AssertionAuth
    JWT_BEARER_GRANT_TYPE: Incomplete
    ASSERTION_METHODS: Incomplete
    DEFAULT_GRANT_TYPE: Incomplete
    default_timeout: Incomplete
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
        default_timeout=None,
        leeway=60,
        **kwargs,
    ) -> None: ...
    def request(self, method, url, withhold_token=False, auth=None, **kwargs): ...
