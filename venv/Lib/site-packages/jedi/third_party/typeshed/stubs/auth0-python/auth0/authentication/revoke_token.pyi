from .base import AuthenticationBase

class RevokeToken(AuthenticationBase):
    def revoke_refresh_token(self, token: str): ...
