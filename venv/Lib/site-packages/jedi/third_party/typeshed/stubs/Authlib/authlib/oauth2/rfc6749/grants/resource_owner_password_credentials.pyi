from logging import Logger

from authlib.oauth2.rfc6749 import BaseGrant, TokenEndpointMixin

log: Logger

class ResourceOwnerPasswordCredentialsGrant(BaseGrant, TokenEndpointMixin):
    GRANT_TYPE: str
    def validate_token_request(self) -> None: ...
    def create_token_response(self): ...
    def authenticate_user(self, username, password): ...
