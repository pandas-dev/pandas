from .authorization import AuthorizationEndpoint as AuthorizationEndpoint
from .introspect import IntrospectEndpoint as IntrospectEndpoint
from .metadata import MetadataEndpoint as MetadataEndpoint
from .pre_configured import (
    BackendApplicationServer as BackendApplicationServer,
    LegacyApplicationServer as LegacyApplicationServer,
    MobileApplicationServer as MobileApplicationServer,
    Server as Server,
    WebApplicationServer as WebApplicationServer,
)
from .resource import ResourceEndpoint as ResourceEndpoint
from .revocation import RevocationEndpoint as RevocationEndpoint
from .token import TokenEndpoint as TokenEndpoint
