from authlib.oauth2 import OAuth2Error

class InteractionRequiredError(OAuth2Error):
    error: str

class LoginRequiredError(OAuth2Error):
    error: str

class AccountSelectionRequiredError(OAuth2Error):
    error: str

class ConsentRequiredError(OAuth2Error):
    error: str

class InvalidRequestURIError(OAuth2Error):
    error: str

class InvalidRequestObjectError(OAuth2Error):
    error: str

class RequestNotSupportedError(OAuth2Error):
    error: str

class RequestURINotSupportedError(OAuth2Error):
    error: str

class RegistrationNotSupportedError(OAuth2Error):
    error: str
