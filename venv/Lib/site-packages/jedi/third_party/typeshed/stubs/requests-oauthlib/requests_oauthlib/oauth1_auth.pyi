from logging import Logger
from typing import Any

from oauthlib.oauth1 import Client
from requests.auth import AuthBase

CONTENT_TYPE_FORM_URLENCODED: str
CONTENT_TYPE_MULTI_PART: str
log: Logger

class OAuth1(AuthBase):
    client_class: type[Client]
    client: Client
    force_include_body: bool
    def __init__(
        self,
        client_key,
        client_secret=None,
        resource_owner_key=None,
        resource_owner_secret=None,
        callback_uri=None,
        signature_method="HMAC-SHA1",
        signature_type="AUTH_HEADER",
        rsa_key=None,
        verifier=None,
        decoding: str | None = "utf-8",
        client_class: type[Client] | None = None,
        force_include_body: bool = False,
        *,
        realm=None,
        encoding: str = "utf-8",
        nonce=None,
        timestamp=None,
        **kwargs: Any,  # passed to client_class's __init__
    ) -> None: ...
