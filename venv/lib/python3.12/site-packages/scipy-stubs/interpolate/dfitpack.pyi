# This module is not meant for public use and will be removed in SciPy v2.0.0.

from typing import Any, Final, Protocol, type_check_only
from typing_extensions import deprecated

__all__ = ["spalde", "splder", "splev", "splint", "sproot"]

###

@type_check_only
class _DeprecatedFortranFunction(Protocol):
    __name__: str

    @deprecated(
        "The `scipy.interpolate.dfitpack` namespace is deprecated and will be removed in SciPy 2.0.0. "
        "Please use the `scipy.interpolate` namespace instead."
    )
    def __call__(self, /, *args: object, **kwds: object) -> Any: ...

###

spalde: Final[_DeprecatedFortranFunction] = ...
splder: Final[_DeprecatedFortranFunction] = ...
splev: Final[_DeprecatedFortranFunction] = ...
splint: Final[_DeprecatedFortranFunction] = ...
sproot: Final[_DeprecatedFortranFunction] = ...
