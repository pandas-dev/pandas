from typing import Final
from typing_extensions import deprecated

__all__ = ["available", "enable", "disable", "enabled"]

available: Final = True
enabled: Final = True

@deprecated("Function `enable` is deprecated and no longer has any effect. Speedups are always available.")
def enable() -> None: ...
@deprecated("Function `disable` is deprecated and no longer has any effect. Speedups are always available.")
def disable() -> None: ...
