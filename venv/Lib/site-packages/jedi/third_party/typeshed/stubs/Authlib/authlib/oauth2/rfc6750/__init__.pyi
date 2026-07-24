from .errors import InsufficientScopeError as InsufficientScopeError, InvalidTokenError as InvalidTokenError
from .parameters import add_bearer_token as add_bearer_token
from .token import BearerTokenGenerator as BearerTokenGenerator
from .validator import BearerTokenValidator as BearerTokenValidator

__all__ = [
    "InvalidTokenError",
    "InsufficientScopeError",
    "add_bearer_token",
    "BearerToken",
    "BearerTokenGenerator",
    "BearerTokenValidator",
]

BearerToken = BearerTokenGenerator
