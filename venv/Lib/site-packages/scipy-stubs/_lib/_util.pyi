import inspect
import multiprocessing.pool as mpp
import sys
import types
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from contextlib import _GeneratorContextManager
from typing import Any, Concatenate, Final, Generic, NamedTuple, Never, Self, SupportsIndex, overload, override
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy._typing import ExitMixin

_T_co = TypeVar("_T_co", default=Any, covariant=True)
_T_contra = TypeVar("_T_contra", default=Never, contravariant=True)

###

copy_if_needed: Final = None

# NOTE: These aliases are implictly exported at runtime (I don't like this).
# NOTE: We can't use unions here because of a bug in mypy/stubtest
type IntNumber = Any
type DecimalNumber = Any
type _RNG = Any
type SeedType = Any

# see https://github.com/scipy/scipy/pull/25225
GeneratorType = TypeVar("GeneratorType", bound=_RNG)  # noqa: PYI001

if sys.version_info >= (3, 14):
    def wrapped_inspect_signature(callable: Callable[..., object]) -> inspect.Signature: ...

else:
    wrapped_inspect_signature = inspect.signature

class AxisError(ValueError, IndexError):
    __slots__ = "_msg", "axis", "ndim"

    _msg: Final[str | None]
    axis: Final[int | None]
    ndim: Final[onp.NDim | None]
    @overload
    def __init__(self, /, axis: str, ndim: None = None, msg_prefix: None = None) -> None: ...
    @overload
    def __init__(self, /, axis: int, ndim: onp.NDim, msg_prefix: str | None = None) -> None: ...

class FullArgSpec(NamedTuple):
    args: list[str]
    varargs: str | None
    varkw: str | None
    defaults: tuple[Any, ...] | None
    kwonlyargs: list[str]
    kwonlydefaults: dict[str, Any] | None
    annotations: dict[str, Any]

class _FunctionWrapper(Generic[_T_contra, _T_co]):
    f: Callable[Concatenate[_T_contra, ...], _T_co]
    args: tuple[Any, ...]
    @overload
    def __init__(self, /, f: Callable[[_T_contra], _T_co], args: tuple[()]) -> None: ...
    @overload
    def __init__(self, /, f: Callable[Concatenate[_T_contra, ...], _T_co], args: tuple[object, ...]) -> None: ...
    def __call__(self, /, x: _T_contra) -> _T_co: ...

class MapWrapper(ExitMixin):
    pool: int | mpp.Pool | None

    def __init__[VT, RT](self, /, pool: Callable[[Callable[[VT], RT], Iterable[VT]], Iterable[RT]] | int = 1) -> None: ...
    def __call__[VT, RT](self, /, func: Callable[[VT], RT], iterable: Iterable[VT]) -> Iterable[RT]: ...
    def __enter__(self, /) -> Self: ...
    def terminate(self, /) -> None: ...
    def join(self, /) -> None: ...
    def close(self, /) -> None: ...

# NOTE: At runtime this is a subclass of `dict`. But since that would require us to make this an invariant type,
# we instead pretend this "just" a `Mapping`, so that it can be covariant in the value type.
class _RichResult(Mapping[str, _T_co]):
    def __getattr__(self, name: str, /) -> _T_co: ...

    # These abstract methods of `Mapping` must be overridden even in stubs. And that's exactly why mixing protocols
    # and `abc` is not a good idea, but alas...
    @override
    def __getitem__(self, name: str, /) -> _T_co: ...
    @override
    def __iter__(self, /) -> Iterator[str]: ...
    @override
    def __len__(self) -> int: ...

#
def float_factorial(n: SupportsIndex) -> float: ...  # will be `np.inf` if `n >= 171`

#
def getfullargspec_no_self(func: Callable[..., object]) -> FullArgSpec: ...

#
@overload
def check_random_state[AnyRngT: (np.random.RandomState, np.random.Generator)](seed: AnyRngT) -> AnyRngT: ...
@overload
def check_random_state(seed: onp.ToJustInt | types.ModuleType | None) -> np.random.RandomState: ...

#
@overload
def rng_integers(
    gen: onp.random.RNG | None,
    low: onp.ToInt,
    high: onp.ToInt | None = None,
    size: tuple[()] | None = None,
    dtype: onp.AnyIntegerDType = "int64",
    endpoint: bool = False,
) -> npc.integer: ...
@overload
def rng_integers(
    gen: onp.random.RNG | None,
    low: onp.ToInt | onp.ToIntND,
    high: onp.ToInt | onp.ToIntND | None = None,
    size: SupportsIndex | Sequence[SupportsIndex] | None = None,
    dtype: onp.AnyIntegerDType = "int64",
    endpoint: bool = False,
) -> npc.integer | onp.ArrayND[npc.integer]: ...

#
def ignore_warns(expected_warning: type[Warning], *, match: str | None = None) -> _GeneratorContextManager[None]: ...

#
@overload
def normalize_axis_index(axis: int, ndim: onp.NDim) -> onp.NDim: ...
@overload
def normalize_axis_index[AxisT: npc.integer](axis: int | AxisT, ndim: AxisT) -> AxisT: ...
@overload
def normalize_axis_index[AxisT: npc.integer](axis: AxisT, ndim: onp.NDim | AxisT) -> AxisT: ...

#
def broadcastable(shape_a: tuple[int, ...], shape_b: tuple[int, ...]) -> bool: ...
