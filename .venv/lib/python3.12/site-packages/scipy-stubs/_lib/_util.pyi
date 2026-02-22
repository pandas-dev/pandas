import inspect
import multiprocessing.pool as mpp
import sys
import types
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from contextlib import _GeneratorContextManager
from typing import Any, Concatenate, Final, Generic, Literal, NamedTuple, Never, Self, TypeAlias, overload
from typing_extensions import TypeVar, override

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy._typing import ExitMixin

_AnyRNGT = TypeVar("_AnyRNGT", np.random.RandomState, np.random.Generator)

_VT = TypeVar("_VT")
_RT = TypeVar("_RT")

_T_co = TypeVar("_T_co", default=Any, covariant=True)
_T_contra = TypeVar("_T_contra", default=Never, contravariant=True)

_AxisT = TypeVar("_AxisT", bound=npc.integer)

###

np_long: Final[type[np.int32 | np.int64]] = ...  # `np.long` on `numpy>=2`, else `np.int_`
np_ulong: Final[type[np.uint32 | np.uint64]] = ...  # `np.ulong` on `numpy>=2`, else `np.uint`
copy_if_needed: Final[Literal[False] | None] = ...  # `None` on `numpy>=2`, otherwise `False`

###

# mypy<=1.19.0 workaround, see https://github.com/python/mypy/pull/20392
if sys.version_info >= (3, 14):
    __conditional_annotations__: Final[set[int]] = ...

# NOTE: These aliases are implictly exported at runtime
IntNumber: TypeAlias = int | npc.integer
DecimalNumber: TypeAlias = float | npc.floating | npc.integer
_RNG: TypeAlias = np.random.Generator | np.random.RandomState
SeedType: TypeAlias = IntNumber | _RNG | None
GeneratorType = TypeVar("GeneratorType", bound=_RNG)  # noqa: PYI001  # oof

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

    def __init__(self, /, pool: Callable[[Callable[[_VT], _RT], Iterable[_VT]], Iterable[_RT]] | int = 1) -> None: ...
    def __call__(self, /, func: Callable[[_VT], _RT], iterable: Iterable[_VT]) -> Iterable[_RT]: ...
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
def float_factorial(n: op.CanIndex) -> float: ...  # will be `np.inf` if `n >= 171`

#
def getfullargspec_no_self(func: Callable[..., object]) -> FullArgSpec: ...

#
@overload
def check_random_state(seed: _AnyRNGT) -> _AnyRNGT: ...
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
    endpoint: op.CanBool = False,
) -> npc.integer: ...
@overload
def rng_integers(
    gen: onp.random.RNG | None,
    low: onp.ToInt | onp.ToIntND,
    high: onp.ToInt | onp.ToIntND | None = None,
    size: op.CanIndex | Sequence[op.CanIndex] | None = None,
    dtype: onp.AnyIntegerDType = "int64",
    endpoint: op.CanBool = False,
) -> npc.integer | onp.ArrayND[npc.integer]: ...

#
def ignore_warns(expected_warning: type[Warning], *, match: str | None = None) -> _GeneratorContextManager[None]: ...

#
@overload
def normalize_axis_index(axis: int, ndim: onp.NDim) -> onp.NDim: ...
@overload
def normalize_axis_index(axis: int | _AxisT, ndim: _AxisT) -> _AxisT: ...
@overload
def normalize_axis_index(axis: _AxisT, ndim: onp.NDim | _AxisT) -> _AxisT: ...

#
@overload
def np_vecdot(x1: onp.ToIntStrict1D, x2: onp.ToIntStrict1D, /, *, axis: op.CanIndex = -1) -> npc.integer: ...
@overload
def np_vecdot(x1: onp.ToFloatStrict1D, x2: onp.ToJustFloatStrict1D, /, *, axis: op.CanIndex = -1) -> npc.floating: ...
@overload
def np_vecdot(x1: onp.ToJustFloatStrict1D, x2: onp.ToFloatStrict1D, /, *, axis: op.CanIndex = -1) -> npc.floating: ...
@overload
def np_vecdot(x1: onp.ToComplexStrict1D, x2: onp.ToJustComplexStrict1D, /, *, axis: op.CanIndex = -1) -> npc.complexfloating: ...
@overload
def np_vecdot(x1: onp.ToJustComplexStrict1D, x2: onp.ToComplexStrict1D, /, *, axis: op.CanIndex = -1) -> npc.complexfloating: ...

#
def broadcastable(shape_a: tuple[int, ...], shape_b: tuple[int, ...]) -> bool: ...
