from oauthlib.common import Request
from oauthlib.oauth2 import RequestValidator as OAuth2RequestValidator

class RequestValidator(OAuth2RequestValidator):
    def client_authentication_required(self, request: Request, *args, **kwargs) -> bool: ...
