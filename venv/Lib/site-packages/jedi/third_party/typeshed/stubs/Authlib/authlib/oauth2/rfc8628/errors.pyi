from authlib.oauth2 import OAuth2Error

class AuthorizationPendingError(OAuth2Error):
    error: str

class SlowDownError(OAuth2Error):
    error: str

class ExpiredTokenError(OAuth2Error):
    error: str
