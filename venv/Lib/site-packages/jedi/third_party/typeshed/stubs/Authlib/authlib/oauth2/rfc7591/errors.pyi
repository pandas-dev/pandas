from authlib.oauth2 import OAuth2Error

class InvalidRedirectURIError(OAuth2Error):
    error: str

class InvalidClientMetadataError(OAuth2Error):
    error: str

class InvalidSoftwareStatementError(OAuth2Error):
    error: str

class UnapprovedSoftwareStatementError(OAuth2Error):
    error: str
