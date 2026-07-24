from _typeshed import Incomplete

from ..base_client import FrameworkIntegration

# actual type is django.dispatch.Signal
token_update: Incomplete

class DjangoIntegration(FrameworkIntegration):
    def update_token(self, token, refresh_token=None, access_token=None) -> None: ...
    @staticmethod
    def load_config(oauth, name, params): ...
