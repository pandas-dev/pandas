from _typeshed import Incomplete

from .base import AuthenticationBase

class BackChannelLogin(AuthenticationBase):
    def back_channel_login(
        self,
        binding_message: str,
        login_hint: str,
        scope: str,
        authorization_details: str | list[dict[str, Incomplete]] | None = None,
        **kwargs,
    ): ...
