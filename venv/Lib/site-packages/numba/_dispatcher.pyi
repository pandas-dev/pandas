# Although this module is not public API, `numba.core.Dispatcher` depends on
# `numba._dispatcher`, so it needs to be typed. The `compute_fingerprint` and
# `typeof_init` functions are annotated for the sake of completeness.

import array as _array
import mmap as _mmap

import numpy as np

from collections.abc import Callable
from typing import Any, type_check_only

from typing_extensions import Generic, ParamSpec, Protocol, TypeAlias, TypeVar

_SignatureT_co = TypeVar(
    "_SignatureT_co",
    bound=Callable[..., object],
    default=Callable[..., Any],
    covariant=True,
)
_ParamsT = ParamSpec("_ParamsT")
_ReturnT = TypeVar("_ReturnT")


@type_check_only
class _HasValue(Protocol):
    @property
    def value(self, /) -> object: ...


_Fingerprintable: TypeAlias = (
    None
    | bool
    | int
    | float
    | complex
    | bytes
    | bytearray
    | memoryview
    | _array.array
    | _mmap.mmap
    | _HasValue
    | np.generic
    | np.dtype
    | np.ndarray
    | tuple["_Fingerprintable", ...]
    | list["_Fingerprintable"]
    | set["_Fingerprintable"]
)

# NOTE: `Dispatcher` is not a generic type at runtime.
# undocumented
class Dispatcher(Generic[_SignatureT_co]):
    def __init__(
        self,
        tmaddrobj: int,
        argct: int,
        fold_args: bool,
        argnames: tuple[str, ...],
        defargs: tuple[object, ...],
        can_fallback: bool,
        has_stararg: bool = False,
        exact_match_required: bool = False,
        /,
    ) -> None: ...
    def __call__(
        self: Dispatcher[Callable[_ParamsT, _ReturnT]],
        /,
        *args: _ParamsT.args,
        **kwargs: _ParamsT.kwargs,
    ) -> _ReturnT: ...

# undocumented
def compute_fingerprint(val: _Fingerprintable, /) -> bytes: ...

# undocumented
def typeof_init(
    omittedarg_type: type[_HasValue],
    typecode_dict: dict[str, int],
    /,
) -> None: ...
