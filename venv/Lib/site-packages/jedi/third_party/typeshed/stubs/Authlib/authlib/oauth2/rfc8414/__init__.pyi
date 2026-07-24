from .models import AuthorizationServerMetadata as AuthorizationServerMetadata
from .well_known import get_well_known_url as get_well_known_url

__all__ = ["AuthorizationServerMetadata", "get_well_known_url"]
