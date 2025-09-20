from typing import Final, TypedDict, type_check_only
from typing_extensions import CapsuleType, ReadOnly

from scipy._lib._ccallback import LowLevelCallable as LowLevelCallable

@type_check_only
class _CApiDict(TypedDict):
    _F_integrand: ReadOnly[CapsuleType]
    _F_integrand1: ReadOnly[CapsuleType]
    _F_integrand2: ReadOnly[CapsuleType]
    _F_integrand3: ReadOnly[CapsuleType]
    _F_integrand4: ReadOnly[CapsuleType]
    _set_action: ReadOnly[CapsuleType]

###

__pyx_capi__: Final[_CApiDict] = ...

nan: Final[float] = ...

def _ellipsoid(h2: float, k2: float, n: int, p: int, s: float) -> float: ...
def _ellipsoid_norm(h2: float, k2: float, n: int, p: int) -> float: ...
