from logging import Logger

from authlib.oauth2.rfc6749 import BaseGrant, TokenEndpointMixin

log: Logger

class ClientCredentialsGrant(BaseGrant, TokenEndpointMixin):
    GRANT_TYPE: str
    def validate_token_request(self) -> None: ...
    def create_token_response(self): ...
