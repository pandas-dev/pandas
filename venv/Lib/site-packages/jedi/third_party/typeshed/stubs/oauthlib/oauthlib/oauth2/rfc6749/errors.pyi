from _typeshed import Incomplete
from typing import NoReturn

from oauthlib.common import Request

class OAuth2Error(Exception):
    error: str | None
    status_code: int
    description: str
    uri: str | None
    state: str | None
    redirect_uri: str | None
    client_id: str | None
    scopes: Incomplete | None
    response_type: str | None
    response_mode: str | None
    grant_type: str | None
    def __init__(
        self,
        description: str | None = None,
        uri: str | None = None,
        state: str | None = None,
        status_code: int | None = None,
        request: Request | None = None,
    ) -> None: ...
    def in_uri(self, uri: str) -> str: ...
    @property
    def twotuples(self) -> list[tuple[str, str | None]]: ...
    @property
    def urlencoded(self) -> str: ...
    @property
    def json(self) -> str: ...
    @property
    def headers(self) -> dict[str, str]: ...

class TokenExpiredError(OAuth2Error):
    error: str

class InsecureTransportError(OAuth2Error):
    error: str
    description: str

class MismatchingStateError(OAuth2Error):
    error: str
    description: str

class MissingCodeError(OAuth2Error):
    error: str

class MissingTokenError(OAuth2Error):
    error: str

class MissingTokenTypeError(OAuth2Error):
    error: str

class FatalClientError(OAuth2Error): ...

class InvalidRequestFatalError(FatalClientError):
    error: str

class InvalidRedirectURIError(InvalidRequestFatalError):
    description: str

class MissingRedirectURIError(InvalidRequestFatalError):
    description: str

class MismatchingRedirectURIError(InvalidRequestFatalError):
    description: str

class InvalidClientIdError(InvalidRequestFatalError):
    description: str

class MissingClientIdError(InvalidRequestFatalError):
    description: str

class InvalidRequestError(OAuth2Error):
    error: str

class MissingResponseTypeError(InvalidRequestError):
    description: str

class MissingCodeChallengeError(InvalidRequestError):
    description: str

class MissingCodeVerifierError(InvalidRequestError):
    description: str

class AccessDeniedError(OAuth2Error):
    error: str

class UnsupportedResponseTypeError(OAuth2Error):
    error: str

class UnsupportedCodeChallengeMethodError(InvalidRequestError):
    description: str

class InvalidScopeError(OAuth2Error):
    error: str

class ServerError(OAuth2Error):
    error: str

class TemporarilyUnavailableError(OAuth2Error):
    error: str

class InvalidClientError(FatalClientError):
    error: str
    status_code: int

class InvalidGrantError(OAuth2Error):
    error: str
    status_code: int

class UnauthorizedClientError(OAuth2Error):
    error: str

class UnsupportedGrantTypeError(OAuth2Error):
    error: str

class UnsupportedTokenTypeError(OAuth2Error):
    error: str

class InvalidTokenError(OAuth2Error):
    error: str
    status_code: int
    description: str

class InsufficientScopeError(OAuth2Error):
    error: str
    status_code: int
    description: str

class ConsentRequired(OAuth2Error):
    error: str

class LoginRequired(OAuth2Error):
    error: str

class CustomOAuth2Error(OAuth2Error):
    def __init__(
        self,
        error: str,
        description: str | None = None,
        uri: str | None = None,
        state: str | None = None,
        status_code: int | None = None,
        request: Request | None = None,
    ) -> None: ...

def raise_from_error(error: str, params: dict[str, Incomplete] | None = None) -> NoReturn: ...
