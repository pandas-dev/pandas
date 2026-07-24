from . import (
    AccessTokenEndpoint as AccessTokenEndpoint,
    AuthorizationEndpoint as AuthorizationEndpoint,
    RequestTokenEndpoint as RequestTokenEndpoint,
    ResourceEndpoint as ResourceEndpoint,
)

class WebApplicationServer(RequestTokenEndpoint, AuthorizationEndpoint, AccessTokenEndpoint, ResourceEndpoint):
    def __init__(self, request_validator) -> None: ...
