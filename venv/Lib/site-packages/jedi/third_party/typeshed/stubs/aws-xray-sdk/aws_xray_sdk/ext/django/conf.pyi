from _typeshed import Incomplete
from typing import Final

DEFAULTS: dict[str, str | bool | tuple[Incomplete] | list[Incomplete] | None]
XRAY_NAMESPACE: Final = "XRAY_RECORDER"
SUPPORTED_ENV_VARS: tuple[str, ...]

class XRaySettings:
    defaults: dict[str, str | bool | tuple[Incomplete] | list[Incomplete] | None]
    def __init__(self, user_settings=None) -> None: ...
    @property
    def user_settings(self): ...
    def __getattr__(self, attr): ...

settings: XRaySettings

def reload_settings(*, settings: str | None = None, value=None) -> None: ...
