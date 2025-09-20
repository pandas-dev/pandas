import sys
import types
from collections.abc import Callable, Iterable, Mapping, Sequence
from typing import Any, Literal as L, Protocol, TypeAlias  # noqa: N817

if sys.version_info >= (3, 13):
    from typing import TypeVar, runtime_checkable
else:
    from typing_extensions import TypeVar, runtime_checkable

import numpy as np

from ._shape import AnyShape
from optype._utils import set_module

__all__ = ["CanArrayFunction", "CanArrayUFunc", "UFunc"]


def __dir__() -> list[str]:
    return __all__


###


_AnyFunc: TypeAlias = Callable[..., Any]
_AnyArray: TypeAlias = np.ndarray[AnyShape, np.dtype[Any]]

_FT_co = TypeVar("_FT_co", bound=_AnyFunc, default=_AnyFunc, covariant=True)
_NInT_co = TypeVar("_NInT_co", bound=int, default=int, covariant=True)
_NoutT_co = TypeVar("_NoutT_co", bound=int, default=int, covariant=True)
_SigT_co = TypeVar("_SigT_co", bound=str | None, default=str | None, covariant=True)
_IdT_co = TypeVar(
    "_IdT_co",
    bound=complex | bytes | str | None,
    default=float | None,
    covariant=True,
)


###


@runtime_checkable
@set_module("optype.numpy")
class UFunc(Protocol[_FT_co, _NInT_co, _NoutT_co, _SigT_co, _IdT_co]):
    """
    A generic interface for `numpy.ufunc` "universal function" instances,
    e.g. `numpy.exp`, `numpy.add`, `numpy.frexp`, `numpy.divmod`.

    This also includes gufunc's (generalized universion functions), which
    have a specified `signature`, and aren't necessarily element-wise
    functions (which "regular" ufuncs are).
    At the moment (`numpy>=2.2,<2.4`), the only `GUFuncs` in the public numpy API
    are `matmul`, `matvec`, `vecdot`, and `vecmat`, and all four have `nin == 2`
    and `nout == 1`.
    """

    @property
    def __call__(self, /) -> _FT_co: ...

    # The number of positional-only parameters, within numpy this is
    # either 1 or 2, but e.g. `scipy.special.pro_rad2_cv` has 5.
    @property
    def nin(self, /) -> _NInT_co: ...
    # The number of output values, within numpy this is either 1 or 2,
    # but e.g. `scipy.special.ellipj` has 4.
    @property
    def nout(self, /) -> _NoutT_co: ...
    # A string i.f.f. this is a gufunc (generalized ufunc).
    @property
    def signature(self, /) -> _SigT_co: ...

    # If `signature is None and nin == 2 and nout == 1`, this *may* be set
    # to a python scalar s.t. `self(x, identity) == x` for all possible
    # `x`.
    # Within numpy==2.0.0, this is only the case for `multiply` (`1`),
    # `logaddexp` (`-inf`), `logaddexp2` (`-inf`), `logical_and` (`True`),
    # and `bitwise_and` (`-1`).
    # Note that the `complex` return annotation implicitly includes
    # `bool | int | float` (these are its supertypes).
    @property
    def identity(self, /) -> _IdT_co: ...
    # Within numpy this is always `nin + nout`, since each output value
    # comes with a corresponding (optional) `out` parameter.
    @property
    def nargs(self, /) -> int: ...
    # Equivalent to `len(types)`, within numpy this is at most 24, but for
    # 3rd party ufuncs it could be more.
    @property
    def ntypes(self, /) -> int: ...
    # A list of strings (`LiteralString` can't be used for compatibility
    # reasons), with signatures in terms of `numpy.dtype.char`, that match
    # `r'(\w{nin})->(\w{nout})'`.
    # For instance, `np.frexp` has `['e->ei', 'f->fi', 'd->di', 'g->gi']`.
    # Note that the `len` of each `types` elements is `nin + nout + 2`.
    # Also note that the elements aren't necessarily unique, because the
    # available data types are system dependent.
    @property
    def types(self, /) -> Sequence[str]: ...

    # NOTE: The following *methods* were incorrectly typed prior to NumPy 2.1:
    # https://github.com/numpy/numpy/pull/26847

    @property  # raises `ValueError` if `nout != 1 or bool(signature)`
    def at(self, /) -> Callable[..., None] | None: ...
    @property  # # raises `ValueError` if `nin != 2 or nout != 1 or bool(signature)`
    def outer(self, /) -> _AnyFunc | None: ...
    @property  # raises `ValueError` if `nin != 2 or nout != 1 or bool(signature)`
    def reduce(self, /) -> _AnyFunc | None: ...
    @property  # raises `ValueError` if `nin != 2 or nout != 1 or bool(signature)`
    def reduceat(self, /) -> Callable[..., _AnyArray] | None: ...
    @property  # raises `ValueError` if `nin != 2 or nout != 1 or bool(signature)`
    def accumulate(self, /) -> Callable[..., _AnyArray] | None: ...


_UFT_contra = TypeVar("_UFT_contra", bound=UFunc, default=np.ufunc, contravariant=True)
_T_co = TypeVar("_T_co", default=Any, covariant=True)

_MethodCommon: TypeAlias = L["__call__", "reduce", "reduceat", "accumulate", "outer"]


@runtime_checkable
@set_module("optype.numpy")
class CanArrayUFunc(Protocol[_UFT_contra, _T_co]):
    """
    Interface for ufunc operands.

    See Also:
        - https://numpy.org/devdocs/reference/arrays.classes.html
    """

    def __array_ufunc__(
        self,
        ufunc: _UFT_contra,
        method: L[_MethodCommon],
        /,
        *args: Any,
        **kwargs: Any,
    ) -> _T_co: ...


_FT_contra = TypeVar("_FT_contra", bound=_AnyFunc, default=_AnyFunc, contravariant=True)


@runtime_checkable
@set_module("optype.numpy")
class CanArrayFunction(Protocol[_FT_contra, _T_co]):
    def __array_function__(
        self,
        /,
        func: _FT_contra,
        # although this could be tighter, this ensures numpy.typing compat
        types: Iterable[type["CanArrayFunction"]],
        # ParamSpec can only be used on *args and **kwargs for some reason...
        args: tuple[Any, ...],
        kwargs: Mapping[str, Any],
    ) -> types.NotImplementedType | _T_co: ...
