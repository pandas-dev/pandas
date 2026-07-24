from _typeshed import Incomplete

from ..base_client import FrameworkIntegration

token_update: Incomplete

class FlaskIntegration(FrameworkIntegration):
    def update_token(self, token, refresh_token=None, access_token=None) -> None: ...
    @staticmethod
    def load_config(oauth, name, params) -> dict[Incomplete, Incomplete]: ...
