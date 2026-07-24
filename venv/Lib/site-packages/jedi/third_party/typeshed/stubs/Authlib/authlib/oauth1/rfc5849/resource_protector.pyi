from authlib.oauth1.rfc5849.base_server import BaseServer

from .wrapper import OAuth1Request

class ResourceProtector(BaseServer):
    def validate_request(self, method, uri, body, headers) -> OAuth1Request: ...
    def get_token_credential(self, request): ...
