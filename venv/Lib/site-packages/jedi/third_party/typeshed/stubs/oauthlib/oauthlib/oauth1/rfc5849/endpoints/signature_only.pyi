from logging import Logger

from .base import BaseEndpoint as BaseEndpoint

log: Logger

class SignatureOnlyEndpoint(BaseEndpoint):
    def validate_request(self, uri, http_method: str = "GET", body=None, headers=None): ...
