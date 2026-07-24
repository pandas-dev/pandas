from collections.abc import Callable, Sequence
from typing import Any, Protocol, type_check_only

import numpy as np
from llvmlite.ir import IRBuilder, Value
from typing_extensions import Concatenate, Generic, Never, TypeVar

from numba.core.base import BaseContext
from numba.core.codegen import CodeLibrary
from numba.core.compiler import CompileResult
from numba.core.types import Type
from numba.core.typing.templates import Signature
from numba.np import npyimpl

_AtT_co = TypeVar(
    "_AtT_co",
    bound=Callable[Concatenate[Never, Never, ...], None],
    default=Callable[Concatenate[Any, Any, ...], None],
    covariant=True,
)
_ReduceT_co = TypeVar(
    "_ReduceT_co",
    bound=Callable[Concatenate[Never, ...], object],
    default=Callable[Concatenate[Any, ...], Any],
    covariant=True,
)
_ReduceAtT_co = TypeVar(
    "_ReduceAtT_co",
    bound=Callable[Concatenate[Never, Never, ...], object],
    default=Callable[Concatenate[Any, Any, ...], np.ndarray],
    covariant=True,
)
_AccumulateT_co = TypeVar(
    "_AccumulateT_co",
    bound=Callable[Concatenate[Never, ...], object],
    default=Callable[Concatenate[Any, ...], np.ndarray],
    covariant=True,
)
_OuterT_co = TypeVar(
    "_OuterT_co",
    bound=Callable[Concatenate[Never, Never, ...], object],
    default=Callable[Concatenate[Any, Any, ...], Any],
    covariant=True,
)
_IdentityT_co = TypeVar("_IdentityT_co", default=Any, covariant=True)
_SignatureT_co = TypeVar(
    "_SignatureT_co",
    bound=str | None,
    default=str | None,
    covariant=True,
)

@type_check_only
class _UFuncInterface(
    Protocol[
        _AtT_co,
        _ReduceT_co,
        _ReduceAtT_co,
        _AccumulateT_co,
        _OuterT_co,
        _IdentityT_co,
        _SignatureT_co,
    ]
):
    @property
    def at(self) -> _AtT_co: ...
    @property
    def reduce(self) -> _ReduceT_co: ...
    @property
    def reduceat(self) -> _ReduceAtT_co: ...
    @property
    def accumulate(self) -> _AccumulateT_co: ...
    @property
    def outer(self) -> _OuterT_co: ...
    @property
    def identity(self) -> _IdentityT_co: ...
    @property
    def signature(self) -> _SignatureT_co: ...
    @property
    def nin(self) -> int: ...
    @property
    def nout(self) -> int: ...
    @property
    def nargs(self) -> int: ...
    @property
    def ntypes(self) -> int: ...
    @property
    def types(self) -> list[str]: ...

@type_check_only
class _UFuncLowerable(Protocol):
    @property
    def __name__(self) -> str: ...
    @property
    def nin(self) -> int: ...
    @property
    def nout(self) -> int: ...
    def __call__(self, /, *args: Any, **kwargs: Any) -> Any: ...

# TODO: narrow upper bound?
_UFuncLowerableT_co = TypeVar(
    "_UFuncLowerableT_co",
    bound=_UFuncLowerable,
    default=_UFuncLowerable,
    covariant=True,
)
_KernelT_co = TypeVar(
    "_KernelT_co",
    bound=npyimpl._Kernel,
    default=npyimpl._Kernel,
    covariant=True,
)

###

class UfuncLowererBase(Generic[_UFuncLowerableT_co, _KernelT_co]):
    ufunc: _UFuncLowerableT_co
    kernel: type[_KernelT_co]
    libs: list[CodeLibrary]

    def __init__(
        self,
        ufunc: _UFuncLowerableT_co,
        make_kernel_fn: Callable[[_UFuncLowerableT_co], type[_KernelT_co]],
        make_ufunc_kernel_fn: Callable[
            [
                BaseContext,
                IRBuilder,
                Signature,
                Sequence[Value],
                _UFuncLowerableT_co,
                type[_KernelT_co],
            ],
            Value,
        ],
    ) -> None: ...
    def __call__(
        self,
        context: BaseContext,
        builder: IRBuilder,
        sig: Signature,
        args: Sequence[Value],
    ) -> Value: ...

class UfuncBase(
    Generic[
        _AtT_co,
        _ReduceT_co,
        _ReduceAtT_co,
        _AccumulateT_co,
        _OuterT_co,
        _IdentityT_co,
        _SignatureT_co,
    ]
):
    ufunc: _UFuncInterface[
        _AtT_co,
        _ReduceT_co,
        _ReduceAtT_co,
        _AccumulateT_co,
        _OuterT_co,
        _IdentityT_co,
        _SignatureT_co,
    ]

    # abstract-ish method
    match_signature: Callable[[tuple[Type, ...], Signature], bool]

    @property
    def nin(self) -> int: ...
    @property
    def nout(self) -> int: ...
    @property
    def nargs(self) -> int: ...
    @property
    def ntypes(self) -> int: ...
    @property
    def types(self) -> list[str]: ...
    @property
    def identity(self) -> _IdentityT_co: ...
    @property
    def signature(self) -> _SignatureT_co: ...

    #
    @property
    def accumulate(self) -> _AccumulateT_co: ...
    @property
    def at(self) -> _AtT_co: ...
    @property
    def outer(self) -> _OuterT_co: ...
    @property
    def reduce(self) -> _ReduceT_co: ...
    @property
    def reduceat(self) -> _ReduceAtT_co: ...

    #
    def disable_compile(self) -> None: ...
    def find_ewise_function(
        self, ewise_types: tuple[Type, ...]
    ) -> tuple[Signature, CompileResult] | tuple[None, None]: ...

