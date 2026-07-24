from typing import (
    TYPE_CHECKING,
    Any,
    Final,
    Literal,
    Protocol,
    TypeAlias,
    TypeVar,
)

if TYPE_CHECKING:
    import numpy as np
    from numpy.dtypes import StringDType
    from numpy.polynomial._polybase import ABCPolyBase


__version__: Final = "20260602.2.4"

__all__ = (
    "NUMPY_GE_2_0",
    "NUMPY_GE_2_1",
    "NUMPY_GE_2_2",
    "NUMPY_GE_2_3",
    "NUMPY_GE_2_4",
    "NUMPY_GE_2_5",
    "ABCPolyBase",
    "CanArray",
    "LiteralFalse",
    "LiteralTrue",
    "StringDType",
    "__version__",
    "_check_version",
)


def __dir__() -> tuple[str, ...]:
    return __all__


__ALL_SET = frozenset(__all__)


def __getattr__(name: str, /) -> object:
    if name == "StringDType":
        from numpy.dtypes import StringDType

        return StringDType

    if name == "ABCPolyBase":
        from numpy.polynomial._polybase import ABCPolyBase

        return ABCPolyBase

    if name in __ALL_SET and name in globals():
        return globals()[name]

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


NUMPY_GE_2_0: Final = True
NUMPY_GE_2_1: Final = True
NUMPY_GE_2_2: Final = True
NUMPY_GE_2_3: Final = True
NUMPY_GE_2_4: Final = True
NUMPY_GE_2_5: Final = False


LiteralTrue: TypeAlias = "Literal[True] | np.bool[Literal[True]]"
LiteralFalse: TypeAlias = "Literal[False] | np.bool[Literal[False]]"


ShapeT_co = TypeVar("ShapeT_co", bound=tuple[int, ...], covariant=True)
DTypeT_co = TypeVar("DTypeT_co", bound="np.dtype[Any]", covariant=True)


class CanArray(Protocol[ShapeT_co, DTypeT_co]):
    def __array__(self, /) -> "np.ndarray[ShapeT_co, DTypeT_co]": ...


def _check_version() -> bool:
    """Check if the `numpy-typing-compat` version is compatible with `numpy`."""
    import numpy as np

    np_version = tuple(map(int, np.__version__.split(".", 2)[:2]))
    return (2, 4) <= np_version < (2, 5)
