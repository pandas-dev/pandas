from _typeshed import Incomplete
from typing import Literal

from authlib.common.errors import AuthlibHTTPError

def invalid_error_characters(text: str) -> list[str]: ...

class OAuth2Error(AuthlibHTTPError):
    state: Incomplete
    redirect_uri: Incomplete
    redirect_fragment: Incomplete
    def __init__(
        self,
        description: str | None = None,
        uri=None,
        status_code=None,
        state=None,
        redirect_uri=None,
        redirect_fragment: bool = False,
        error=None,
    ) -> None: ...
    def get_body(self) -> list[tuple[Literal["error", "error_description", "error_uri"], str | None]]: ...
    def __call__(self, uri: str | None = None): ...
