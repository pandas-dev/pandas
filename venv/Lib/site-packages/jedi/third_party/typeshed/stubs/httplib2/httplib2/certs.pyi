from collections.abc import Callable
from typing import Final

certifi_available: bool
certifi_where: Callable[[], str] | None
custom_ca_locater_available: bool
custom_ca_locater_where: Callable[..., str] | None
BUILTIN_CA_CERTS: Final[str]

def where() -> str: ...
