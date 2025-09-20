import sys
from typing import Any, Literal, TypeAlias

if sys.version_info >= (3, 13):
    from typing import TypeAliasType, TypeVar
else:
    from typing_extensions import TypeAliasType, TypeVar


__all__ = [
    "AtLeast0D",
    "AtLeast1D",
    "AtLeast2D",
    "AtLeast3D",
    "AtMost0D",
    "AtMost1D",
    "AtMost2D",
    "AtMost3D",
    "NDim",
    "NDim0",
]


def __dir__() -> list[str]:
    return __all__


###
# undocumented aliases for internal use

Shape: TypeAlias = tuple[int, ...]
AnyShape: TypeAlias = tuple[Any, ...]


###

AxT = TypeVar("AxT", int, Any, default=int)

###
# Shape types with at least N dimensions. They're fully static by default, but can be
# turned "gradual" by passing `Any` as type argument.
AtLeast0D = TypeAliasType("AtLeast0D", tuple[AxT, ...], type_params=(AxT,))
AtLeast1D = TypeAliasType("AtLeast1D", tuple[int, *tuple[AxT, ...]], type_params=(AxT,))
AtLeast2D = TypeAliasType(
    "AtLeast2D",
    tuple[int, int, *tuple[AxT, ...]],
    type_params=(AxT,),
)
AtLeast3D = TypeAliasType(
    "AtLeast3D",
    tuple[int, int, int, *tuple[AxT, ...]],
    type_params=(AxT,),
)

###
# Shape types with at most N dimensions. Unlike the above, these are not gradual due to
# limitations in the Python type system.
AtMost0D = TypeAliasType("AtMost0D", tuple[()])
AtMost1D = TypeAliasType("AtMost1D", tuple[int] | AtMost0D)
AtMost2D = TypeAliasType("AtMost2D", tuple[int, int] | AtMost1D)
AtMost3D = TypeAliasType("AtMost3D", tuple[int, int, int] | AtMost2D)

###
# Integer literal types for the number of dimensions. Note that before `numpy>=2` the
# maximum number of dimensions was 32, but this was increased to 64 in `numpy>=2`.

# NOTE: on `numpy<2` this was at most 32
_1_64: TypeAlias = Literal[  # noqa: PYI042
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
    17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
    33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
    49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64,
]  # fmt: skip
_0_64: TypeAlias = Literal[0, _1_64]  # noqa: PYI042

NDim0 = TypeAliasType("NDim0", _1_64)
"""Integer literal between 1 and 64, inclusive."""

NDim = TypeAliasType("NDim", _0_64)
"""Integer literal between 0 and 64, inclusive."""
