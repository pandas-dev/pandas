from _typeshed import Incomplete
from typing import NoReturn

from authlib.oauth2.rfc7009 import RevocationEndpoint

class JWTRevocationEndpoint(RevocationEndpoint):
    issuer: Incomplete
    def __init__(self, issuer, server=None, *args, **kwargs) -> None: ...
    def authenticate_token(self, request, client) -> NoReturn: ...
    def get_jwks(self): ...
