from logging import Logger

from .base import BaseEndpoint as BaseEndpoint

log: Logger

class AccessTokenEndpoint(BaseEndpoint):
    def create_access_token(self, request, credentials): ...
    def create_access_token_response(self, uri, http_method: str = "GET", body=None, headers=None, credentials=None): ...
    def validate_access_token_request(self, request): ...
