from _typeshed import Incomplete
from logging import Logger

log: Logger

class BaseApp:
    client_cls: Incomplete
    OAUTH_APP_CONFIG: Incomplete
    def request(self, method, url, token=None, **kwargs): ...
    def get(self, url, **kwargs): ...
    def post(self, url, **kwargs): ...
    def patch(self, url, **kwargs): ...
    def put(self, url, **kwargs): ...
    def delete(self, url, **kwargs): ...

class _RequestMixin: ...

class OAuth1Base:
    client_cls: Incomplete
    framework: Incomplete
    name: Incomplete
    client_id: Incomplete
    client_secret: Incomplete
    request_token_url: Incomplete
    request_token_params: Incomplete
    access_token_url: Incomplete
    access_token_params: Incomplete
    authorize_url: Incomplete
    authorize_params: Incomplete
    api_base_url: Incomplete
    client_kwargs: Incomplete
    def __init__(
        self,
        framework,
        name=None,
        fetch_token=None,
        client_id=None,
        client_secret=None,
        request_token_url=None,
        request_token_params=None,
        access_token_url=None,
        access_token_params=None,
        authorize_url=None,
        authorize_params=None,
        api_base_url=None,
        client_kwargs=None,
        user_agent=None,
        **kwargs,
    ) -> None: ...

class OAuth1Mixin(_RequestMixin, OAuth1Base):
    def request(self, method, url, token=None, **kwargs): ...
    def create_authorization_url(self, redirect_uri=None, **kwargs) -> dict[Incomplete, Incomplete]: ...
    def fetch_access_token(self, request_token=None, **kwargs): ...

class OAuth2Base:
    client_cls: Incomplete
    framework: Incomplete
    name: Incomplete
    client_id: Incomplete
    client_secret: Incomplete
    access_token_url: Incomplete
    access_token_params: Incomplete
    authorize_url: Incomplete
    authorize_params: Incomplete
    api_base_url: Incomplete
    client_kwargs: Incomplete
    compliance_fix: Incomplete
    client_auth_methods: Incomplete
    server_metadata: Incomplete
    def __init__(
        self,
        framework,
        name=None,
        fetch_token=None,
        update_token=None,
        client_id=None,
        client_secret=None,
        access_token_url=None,
        access_token_params=None,
        authorize_url=None,
        authorize_params=None,
        api_base_url=None,
        client_kwargs=None,
        server_metadata_url=None,
        compliance_fix=None,
        client_auth_methods=None,
        user_agent=None,
        **kwargs,
    ) -> None: ...

class OAuth2Mixin(_RequestMixin, OAuth2Base):
    def request(self, method, url, token=None, **kwargs): ...
    def load_server_metadata(self) -> dict[Incomplete, Incomplete]: ...
    def create_authorization_url(self, redirect_uri=None, **kwargs) -> dict[Incomplete, Incomplete]: ...
    def fetch_access_token(self, redirect_uri=None, **kwargs): ...
