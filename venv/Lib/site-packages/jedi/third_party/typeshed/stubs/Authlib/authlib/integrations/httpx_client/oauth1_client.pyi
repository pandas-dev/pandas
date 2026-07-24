from _typeshed import Incomplete
from collections.abc import Generator
from typing import NoReturn
from typing_extensions import TypeAlias

from authlib.oauth1 import ClientAuth
from authlib.oauth1.client import OAuth1Client as _OAuth1Client

_Response: TypeAlias = Incomplete  # actual type is httpx.Response
_Request: TypeAlias = Incomplete  # actual type is httpx.Request

# Inherits from httpx.Auth
class OAuth1Auth(ClientAuth):
    requires_request_body: bool
    def auth_flow(self, request: _Request) -> Generator[_Request, _Response]: ...

# Inherits from httpx.AsyncClient
class AsyncOAuth1Client(_OAuth1Client):
    auth_class = OAuth1Auth
    def __init__(
        self,
        client_id,
        client_secret=None,
        token=None,
        token_secret=None,
        redirect_uri=None,
        rsa_key=None,
        verifier=None,
        signature_method=...,
        signature_type=...,
        force_include_body=False,
        **kwargs,
    ) -> None: ...
    async def fetch_access_token(self, url, verifier=None, **kwargs): ...
    @staticmethod
    def handle_error(error_type: str | None, error_description: str | None) -> NoReturn: ...

# Inherits from httpx.Client
class OAuth1Client(_OAuth1Client):
    auth_class = OAuth1Auth
    def __init__(
        self,
        client_id,
        client_secret=None,
        token=None,
        token_secret=None,
        redirect_uri=None,
        rsa_key=None,
        verifier=None,
        signature_method=...,
        signature_type=...,
        force_include_body=False,
        **kwargs,
    ) -> None: ...
    @staticmethod
    def handle_error(error_type: str | None, error_description: str | None) -> NoReturn: ...
