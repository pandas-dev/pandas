from oauthlib.oauth2.rfc6749.errors import FatalClientError, OAuth2Error

class FatalOpenIDClientError(FatalClientError): ...
class OpenIDClientError(OAuth2Error): ...

class InteractionRequired(OpenIDClientError):
    error: str
    status_code: int

class LoginRequired(OpenIDClientError):
    error: str
    status_code: int

class AccountSelectionRequired(OpenIDClientError):
    error: str

class ConsentRequired(OpenIDClientError):
    error: str
    status_code: int

class InvalidRequestURI(OpenIDClientError):
    error: str
    description: str

class InvalidRequestObject(OpenIDClientError):
    error: str
    description: str

class RequestNotSupported(OpenIDClientError):
    error: str
    description: str

class RequestURINotSupported(OpenIDClientError):
    error: str
    description: str

class RegistrationNotSupported(OpenIDClientError):
    error: str
    description: str

class InvalidTokenError(OAuth2Error):
    error: str
    status_code: int
    description: str

class InsufficientScopeError(OAuth2Error):
    error: str
    status_code: int
    description: str

def raise_from_error(error: object, params: dict[str, str] | None = None) -> None: ...
