from oauthlib.oauth2.rfc6749.tokens import OAuth2Token

from .base import Client

class MobileApplicationClient(Client):
    response_type: str
    def prepare_request_uri(
        self,
        uri,
        redirect_uri: str | None = None,
        scope: str | set[object] | tuple[object] | list[object] | None = None,
        state: str | None = None,
        **kwargs,
    ) -> str: ...
    token: OAuth2Token
    def parse_request_uri_response(
        self, uri: str, state: str | None = None, scope: str | set[object] | tuple[object] | list[object] | None = None
    ) -> OAuth2Token: ...
