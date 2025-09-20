# This module is not meant for public use and will be removed in SciPy v2.0.0.
from types import ModuleType
from typing_extensions import deprecated

__all__ = ["bicg", "bicgstab", "cg", "cgs", "gcrotmk", "gmres", "lgmres", "lsmr", "lsqr", "minres", "qmr", "test", "tfqmr"]

test: ModuleType

@deprecated("will be removed in SciPy v2.0.0")
def cgs(
    A: object,
    b: object,
    x0: object = ...,
    *,
    rtol: object = ...,
    atol: object = ...,
    maxiter: object = ...,
    M: object = ...,
    callback: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def gmres(
    A: object,
    b: object,
    x0: object = ...,
    *,
    rtol: object = ...,
    atol: object = ...,
    restart: object = ...,
    maxiter: object = ...,
    M: object = ...,
    callback: object = ...,
    callback_type: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def qmr(
    A: object,
    b: object,
    x0: object = ...,
    *,
    rtol: object = ...,
    atol: object = ...,
    maxiter: object = ...,
    M1: object = ...,
    M2: object = ...,
    callback: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def bicg(
    A: object,
    b: object,
    x0: object = ...,
    *,
    rtol: object = ...,
    atol: object = ...,
    maxiter: object = ...,
    M: object = ...,
    callback: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def bicgstab(
    A: object,
    b: object,
    x0: object = ...,
    *,
    rtol: object = ...,
    atol: object = ...,
    maxiter: object = ...,
    M: object = ...,
    callback: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def cg(
    A: object,
    b: object,
    x0: object = ...,
    *,
    rtol: object = ...,
    atol: object = ...,
    maxiter: object = ...,
    M: object = ...,
    callback: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tfqmr(
    A: object,
    b: object,
    x0: object = ...,
    *,
    rtol: object = ...,
    atol: object = ...,
    maxiter: object = ...,
    M: object = ...,
    callback: object = ...,
    show: object = ...,
) -> object: ...
def minres(
    A: object,
    b: object,
    x0: object = ...,
    *,
    rtol: object = ...,
    shift: object = ...,
    maxiter: object = ...,
    M: object = ...,
    callback: object = ...,
    show: object = ...,
    check: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lsqr(
    A: object,
    b: object,
    damp: object = ...,
    atol: object = ...,
    btol: object = ...,
    conlim: object = ...,
    iter_lim: object = ...,
    show: object = ...,
    calc_var: object = ...,
    x0: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def gcrotmk(
    A: object,
    b: object,
    x0: object = ...,
    *,
    rtol: object = ...,
    atol: object = ...,
    maxiter: object = ...,
    M: object = ...,
    callback: object = ...,
    m: object = ...,
    k: object = ...,
    CU: object = ...,
    discard_C: object = ...,
    truncate: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lgmres(
    A: object,
    b: object,
    x0: object = ...,
    *,
    rtol: object = ...,
    atol: object = ...,
    maxiter: object = ...,
    M: object = ...,
    callback: object = ...,
    inner_m: object = ...,
    outer_k: object = ...,
    outer_v: object = ...,
    store_outer_Av: object = ...,
    prepend_outer_v: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lsmr(
    A: object,
    b: object,
    damp: object = ...,
    atol: object = ...,
    btol: object = ...,
    conlim: object = ...,
    maxiter: object = ...,
    show: object = ...,
    x0: object = ...,
) -> object: ...
