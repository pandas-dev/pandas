from .base import AuthenticationBase

class Social(AuthenticationBase):
    def login(self, access_token: str, connection: str, scope: str = "openid"): ...
