from _typeshed import Incomplete
from collections.abc import Callable
from typing import Final

CONTENT_TYPE_FORM_URLENCODED: Final = "application/x-www-form-urlencoded"
CONTENT_TYPE_MULTI_PART: Final = "multipart/form-data"

class ClientAuth:
    SIGNATURE_METHODS: dict[str, Callable[..., str]]
    @classmethod
    def register_signature_method(cls, name: str, sign: Callable[..., str]) -> None: ...
    client_id: Incomplete
    client_secret: Incomplete
    token: Incomplete
    token_secret: Incomplete
    redirect_uri: Incomplete
    signature_method: Incomplete
    signature_type: Incomplete
    rsa_key: Incomplete
    verifier: Incomplete
    realm: Incomplete
    force_include_body: Incomplete
    def __init__(
        self,
        client_id,
        client_secret=None,
        token=None,
        token_secret=None,
        redirect_uri=None,
        rsa_key=None,
        verifier=None,
        signature_method="HMAC-SHA1",
        signature_type="HEADER",
        realm=None,
        force_include_body: bool = False,
    ) -> None: ...
    def get_oauth_signature(self, method, uri, headers, body) -> str: ...
    def get_oauth_params(self, nonce, timestamp) -> list[tuple[str, Incomplete]]: ...
    def sign(self, method, uri, headers, body) -> tuple[Incomplete, Incomplete, Incomplete]: ...
    def prepare(self, method, uri, headers, body) -> tuple[Incomplete, ...]: ...

def generate_nonce() -> str: ...
def generate_timestamp() -> str: ...
