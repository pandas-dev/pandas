from .base import AuthenticationBase

class Delegated(AuthenticationBase):
    def get_token(
        self,
        target: str,
        api_type: str,
        grant_type: str,
        id_token: str | None = None,
        refresh_token: str | None = None,
        scope: str = "openid",
    ): ...
