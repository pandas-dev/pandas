# pyright: reportIncompleteStub = false
from typing import (
    Any,
    Callable,
    Literal,
    overload,
)

import numba

from pandas._typing import F

def __getattr__(name: str) -> Any: ...  # incomplete
@overload
def jit(signature_or_function: F) -> F: ...
@overload
def jit(
    signature_or_function: str
    | list[str]
    | numba.core.types.abstract.Type
    | list[numba.core.types.abstract.Type] = ...,
    locals: dict = ...,  # TODO: Mapping of local variable names to Numba types
    cache: bool = ...,
    pipeline_class: numba.compiler.CompilerBase = ...,
    boundscheck: bool | None = ...,
    *,
    nopython: bool = ...,
    forceobj: bool = ...,
    looplift: bool = ...,
    error_model: Literal["python", "numpy"] = ...,
    inline: Literal["never", "always"] | Callable = ...,
    # TODO: If a callable is provided it will be called with the call expression
    # node that is requesting inlining, the caller's IR and callee's IR as
    # arguments, it is expected to return Truthy as to whether to inline.
    target: Literal["cpu", "gpu", "npyufunc", "cuda"] = ...,  # deprecated
    nogil: bool = ...,
    parallel: bool = ...,
) -> Callable[[F], F]: ...

njit = jit
generated_jit = jit
