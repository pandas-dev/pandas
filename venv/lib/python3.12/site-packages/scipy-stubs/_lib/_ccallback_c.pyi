from collections.abc import Callable
from typing import Any, Final, SupportsFloat, TypedDict, type_check_only
from typing_extensions import CapsuleType, ReadOnly, TypeIs

@type_check_only
class _CApiDict(TypedDict):
    plus1_cython: ReadOnly[CapsuleType]
    plus1b_cython: ReadOnly[CapsuleType]
    plus1bc_cython: ReadOnly[CapsuleType]
    sine: ReadOnly[CapsuleType]

###

__pyx_capi__: Final[_CApiDict] = ...  # undocumented
__test__: Final[dict[Any, Any]] = ...  # undocumented

def check_capsule(item: object) -> TypeIs[CapsuleType]: ...  # undocumented
def get_capsule_signature(capsule_obj: CapsuleType) -> str: ...  # undocumented
def get_raw_capsule(
    func_obj: CapsuleType | int, name_obj: str, context_obj: CapsuleType | int
) -> CapsuleType: ...  # undocumented

# namespace pollution
idx: int = 1  # undocumented
sig: tuple[bytes, int] = ...  # undocumented
sigs: list[tuple[bytes, int]] = ...  # undocumented

def test_call_cython(callback_obj: Callable[[float], SupportsFloat], value: float) -> float: ...  # undocumented
