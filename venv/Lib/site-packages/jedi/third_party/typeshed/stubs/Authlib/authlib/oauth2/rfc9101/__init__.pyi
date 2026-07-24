from .authorization_server import JWTAuthenticationRequest as JWTAuthenticationRequest
from .discovery import AuthorizationServerMetadata as AuthorizationServerMetadata
from .registration import ClientMetadataClaims as ClientMetadataClaims

__all__ = ["AuthorizationServerMetadata", "JWTAuthenticationRequest", "ClientMetadataClaims"]
