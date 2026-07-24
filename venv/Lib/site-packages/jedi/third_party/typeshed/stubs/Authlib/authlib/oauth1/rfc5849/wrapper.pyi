from _typeshed import Incomplete

class OAuth1Request:
    method: Incomplete
    uri: Incomplete
    body: Incomplete
    headers: Incomplete
    client: Incomplete | None
    credential: Incomplete | None
    user: Incomplete | None
    query: Incomplete
    query_params: Incomplete
    body_params: Incomplete
    auth_params: Incomplete
    realm: Incomplete
    signature_type: str | None
    oauth_params: Incomplete
    params: list[Incomplete]
    def __init__(self, method, uri, body=None, headers=None) -> None: ...
    @property
    def client_id(self): ...
    @property
    def client_secret(self): ...
    @property
    def rsa_public_key(self): ...
    @property
    def timestamp(self): ...
    @property
    def redirect_uri(self): ...
    @property
    def signature(self): ...
    @property
    def signature_method(self): ...
    @property
    def token(self): ...
    @property
    def token_secret(self): ...
