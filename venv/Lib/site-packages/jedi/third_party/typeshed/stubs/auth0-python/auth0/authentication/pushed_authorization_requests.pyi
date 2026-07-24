from .base import AuthenticationBase

class PushedAuthorizationRequests(AuthenticationBase):
    def pushed_authorization_request(self, response_type: str, redirect_uri: str, **kwargs): ...
