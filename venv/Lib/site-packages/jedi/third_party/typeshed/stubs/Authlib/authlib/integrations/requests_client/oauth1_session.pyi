from typing import NoReturn

from authlib.oauth1 import ClientAuth
from authlib.oauth1.client import OAuth1Client

# Inherits from requests.auth.AuthBase
class OAuth1Auth(ClientAuth):
    def __call__(self, req): ...

# Inherits from requests.Session
class OAuth1Session(OAuth1Client):
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
    def rebuild_auth(self, prepared_request, response) -> None: ...
    @staticmethod
    def handle_error(error_type: str | None, error_description: str | None) -> NoReturn: ...
