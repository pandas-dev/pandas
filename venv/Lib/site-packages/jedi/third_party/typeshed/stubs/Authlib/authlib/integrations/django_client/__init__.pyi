from ..base_client import BaseOAuth, OAuthError as OAuthError
from .apps import DjangoOAuth1App as DjangoOAuth1App, DjangoOAuth2App as DjangoOAuth2App
from .integration import DjangoIntegration as DjangoIntegration, token_update as token_update

class OAuth(BaseOAuth):
    oauth1_client_cls = DjangoOAuth1App
    oauth2_client_cls = DjangoOAuth2App
    framework_integration_cls = DjangoIntegration

__all__ = ["OAuth", "DjangoOAuth1App", "DjangoOAuth2App", "DjangoIntegration", "token_update", "OAuthError"]
