from logging import Logger

from oauthlib.common import Request

from ..tokens import TokenBase
from .base import GrantTypeBase

log: Logger

class ResourceOwnerPasswordCredentialsGrant(GrantTypeBase):
    def create_token_response(self, request: Request, token_handler: TokenBase) -> tuple[dict[str, str], str, int | None]: ...
    def validate_token_request(self, request: Request) -> None: ...
