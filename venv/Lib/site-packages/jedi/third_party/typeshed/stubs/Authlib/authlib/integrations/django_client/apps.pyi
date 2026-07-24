from ..base_client import BaseApp, OAuth1Mixin, OAuth2Mixin, OpenIDMixin
from ..requests_client import OAuth1Session, OAuth2Session

class DjangoAppMixin:
    def save_authorize_data(self, request, **kwargs) -> None: ...
    def authorize_redirect(self, request, redirect_uri=None, **kwargs): ...

class DjangoOAuth1App(DjangoAppMixin, OAuth1Mixin, BaseApp):
    client_cls = OAuth1Session
    def authorize_access_token(self, request, **kwargs): ...

class DjangoOAuth2App(DjangoAppMixin, OAuth2Mixin, OpenIDMixin, BaseApp):
    client_cls = OAuth2Session
    def authorize_access_token(self, request, **kwargs): ...
