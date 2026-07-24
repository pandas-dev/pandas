from authlib.common.errors import AuthlibBaseError

class OAuthError(AuthlibBaseError):
    error: str

class MissingRequestTokenError(OAuthError):
    error: str

class MissingTokenError(OAuthError):
    error: str

class TokenExpiredError(OAuthError):
    error: str

class InvalidTokenError(OAuthError):
    error: str

class UnsupportedTokenTypeError(OAuthError):
    error: str

class MismatchingStateError(OAuthError):
    error: str
    description: str
