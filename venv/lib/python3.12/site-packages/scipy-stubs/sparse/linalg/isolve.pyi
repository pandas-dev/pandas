# This module is not meant for public use and will be removed in SciPy v2.0.0.
from types import ModuleType
from typing_extensions import deprecated

__all__ = ["bicg", "bicgstab", "cg", "cgs", "gcrotmk", "gmres", "lgmres", "lsmr", "lsqr", "minres", "qmr", "test", "tfqmr"]

test: ModuleType

@deprecated("will be removed in SciPy v2.0.0")
def cgs(
    A: object,
    b: object,
    x0: object = None,
    *,
    rtol: object = 1e-5,
    atol: object = 0.0,
    maxiter: object = None,
    M: object = None,
    callback: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def gmres(
    A: object,
    b: object,
    x0: object = None,
    *,
    rtol: object = 1e-5,
    atol: object = 0.0,
    restart: object = None,
    maxiter: object = None,
    M: object = None,
    callback: object = None,
    callback_type: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def qmr(
    A: object,
    b: object,
    x0: object = None,
    *,
    rtol: object = 1e-5,
    atol: object = 0.0,
    maxiter: object = None,
    M1: object = None,
    M2: object = None,
    callback: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def bicg(
    A: object,
    b: object,
    x0: object = None,
    *,
    rtol: object = 1e-5,
    atol: object = 0.0,
    maxiter: object = None,
    M: object = None,
    callback: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def bicgstab(
    A: object,
    b: object,
    x0: object = None,
    *,
    rtol: object = 1e-5,
    atol: object = 0.0,
    maxiter: object = None,
    M: object = None,
    callback: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def cg(
    A: object,
    b: object,
    x0: object = None,
    *,
    rtol: object = 1e-5,
    atol: object = 0.0,
    maxiter: object = None,
    M: object = None,
    callback: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tfqmr(
    A: object,
    b: object,
    x0: object = None,
    *,
    rtol: object = 1e-5,
    atol: object = 0.0,
    maxiter: object = None,
    M: object = None,
    callback: object = None,
    show: object = False,
) -> object: ...
def minres(
    A: object,
    b: object,
    x0: object = None,
    *,
    rtol: object = 1e-5,
    shift: object = 0.0,
    maxiter: object = None,
    M: object = None,
    callback: object = None,
    show: object = False,
    check: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lsqr(
    A: object,
    b: object,
    damp: object = 0.0,
    atol: object = 1e-6,
    btol: object = 1e-6,
    conlim: object = ...,
    iter_lim: object = None,
    show: object = False,
    calc_var: object = False,
    x0: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def gcrotmk(
    A: object,
    b: object,
    x0: object = None,
    *,
    rtol: object = 1e-5,
    atol: object = 0.0,
    maxiter: object = 1000,
    M: object = None,
    callback: object = None,
    m: object = 20,
    k: object = None,
    CU: object = None,
    discard_C: object = False,
    truncate: object = "oldest",
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lgmres(
    A: object,
    b: object,
    x0: object = None,
    *,
    rtol: object = 1e-5,
    atol: object = 0.0,
    maxiter: object = 1000,
    M: object = None,
    callback: object = None,
    inner_m: object = 30,
    outer_k: object = 3,
    outer_v: object = None,
    store_outer_Av: object = True,
    prepend_outer_v: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lsmr(
    A: object,
    b: object,
    damp: object = 0.0,
    atol: object = 1e-6,
    btol: object = 1e-6,
    conlim: object = ...,
    maxiter: object = None,
    show: object = False,
    x0: object = None,
) -> object: ...
