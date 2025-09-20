from _typeshed import Incomplete
from collections.abc import Iterable, Iterator, Sequence
from typing import Literal as L, Protocol, SupportsIndex, TypeVar, overload, type_check_only

import numpy as np
import numpy.typing as npt
import optype.numpy as onp
import optype.numpy.compat as npc

from ._fblas import (
    caxpy as caxpy,
    ccopy as ccopy,
    cdotc as cdotc,
    cdotu as cdotu,
    cgbmv as cgbmv,
    cgemm as cgemm,
    cgemv as cgemv,
    cgerc as cgerc,
    cgeru as cgeru,
    chbmv as chbmv,
    chemm as chemm,
    chemv as chemv,
    cher as cher,
    cher2 as cher2,
    cher2k as cher2k,
    cherk as cherk,
    chpmv as chpmv,
    chpr as chpr,
    chpr2 as chpr2,
    crotg as crotg,
    cscal as cscal,
    cspmv as cspmv,
    cspr as cspr,
    csrot as csrot,
    csscal as csscal,
    cswap as cswap,
    csymm as csymm,
    csyr as csyr,
    csyr2k as csyr2k,
    csyrk as csyrk,
    ctbmv as ctbmv,
    ctbsv as ctbsv,
    ctpmv as ctpmv,
    ctpsv as ctpsv,
    ctrmm as ctrmm,
    ctrmv as ctrmv,
    ctrsm as ctrsm,
    ctrsv as ctrsv,
    dasum as dasum,
    daxpy as daxpy,
    dcopy as dcopy,
    ddot as ddot,
    dgbmv as dgbmv,
    dgemm as dgemm,
    dgemv as dgemv,
    dger as dger,
    dnrm2 as dnrm2,
    drot as drot,
    drotg as drotg,
    drotm as drotm,
    drotmg as drotmg,
    dsbmv as dsbmv,
    dscal as dscal,
    dspmv as dspmv,
    dspr as dspr,
    dspr2 as dspr2,
    dswap as dswap,
    dsymm as dsymm,
    dsymv as dsymv,
    dsyr as dsyr,
    dsyr2 as dsyr2,
    dsyr2k as dsyr2k,
    dsyrk as dsyrk,
    dtbmv as dtbmv,
    dtbsv as dtbsv,
    dtpmv as dtpmv,
    dtpsv as dtpsv,
    dtrmm as dtrmm,
    dtrmv as dtrmv,
    dtrsm as dtrsm,
    dtrsv as dtrsv,
    dzasum as dzasum,
    dznrm2 as dznrm2,
    icamax as icamax,
    idamax as idamax,
    isamax as isamax,
    izamax as izamax,
    sasum as sasum,
    saxpy as saxpy,
    scasum as scasum,
    scnrm2 as scnrm2,
    scopy as scopy,
    sdot as sdot,
    sgbmv as sgbmv,
    sgemm as sgemm,
    sgemv as sgemv,
    sger as sger,
    snrm2 as snrm2,
    srot as srot,
    srotg as srotg,
    srotm as srotm,
    srotmg as srotmg,
    ssbmv as ssbmv,
    sscal as sscal,
    sspmv as sspmv,
    sspr as sspr,
    sspr2 as sspr2,
    sswap as sswap,
    ssymm as ssymm,
    ssymv as ssymv,
    ssyr as ssyr,
    ssyr2 as ssyr2,
    ssyr2k as ssyr2k,
    ssyrk as ssyrk,
    stbmv as stbmv,
    stbsv as stbsv,
    stpmv as stpmv,
    stpsv as stpsv,
    strmm as strmm,
    strmv as strmv,
    strsm as strsm,
    strsv as strsv,
    zaxpy as zaxpy,
    zcopy as zcopy,
    zdotc as zdotc,
    zdotu as zdotu,
    zdrot as zdrot,
    zdscal as zdscal,
    zgbmv as zgbmv,
    zgemm as zgemm,
    zgemv as zgemv,
    zgerc as zgerc,
    zgeru as zgeru,
    zhbmv as zhbmv,
    zhemm as zhemm,
    zhemv as zhemv,
    zher as zher,
    zher2 as zher2,
    zher2k as zher2k,
    zherk as zherk,
    zhpmv as zhpmv,
    zhpr as zhpr,
    zhpr2 as zhpr2,
    zrotg as zrotg,
    zscal as zscal,
    zspmv as zspmv,
    zspr as zspr,
    zswap as zswap,
    zsymm as zsymm,
    zsyr as zsyr,
    zsyr2k as zsyr2k,
    zsyrk as zsyrk,
    ztbmv as ztbmv,
    ztbsv as ztbsv,
    ztpmv as ztpmv,
    ztpsv as ztpsv,
    ztrmm as ztrmm,
    ztrmv as ztrmv,
    ztrsm as ztrsm,
    ztrsv as ztrsv,
)

__all__ = ["find_best_blas_type", "get_blas_funcs"]

###

_VT_co = TypeVar("_VT_co", covariant=True)

# A slightly modified variant from https://github.com/python/typing/issues/256#issuecomment-1442633430
# This works because `str.__contains__` does not accept object (either in typeshed or at runtime)
@type_check_only
class _SequenceNotStr(Protocol[_VT_co]):
    @overload
    def __getitem__(self, index: SupportsIndex, /) -> _VT_co: ...
    @overload
    def __getitem__(self, index: slice, /) -> Sequence[_VT_co]: ...
    def __iter__(self, /) -> Iterator[_VT_co]: ...
    def __reversed__(self, /) -> Iterator[_VT_co]: ...
    def __contains__(self, value: object, /) -> bool: ...  # <-- the trick
    def __len__(self, /) -> int: ...
    def index(self, value: object, start: int = 0, stop: int = ..., /) -> int: ...
    def count(self, value: object, /) -> int: ...

# NOTE: used in `lapack.pyi`
@type_check_only
class _FortranFunction(Protocol):
    @property
    def dtype(self, /) -> np.dtype[Incomplete]: ...
    @property
    def int_dtype(self, /) -> np.dtype[npc.integer]: ...
    @property
    def module_name(self, /) -> str: ...
    @property
    def prefix(self, /) -> str: ...
    @property
    def typecode(self, /) -> str: ...
    def __call__(self, /, *args: object, **kwargs: object) -> Incomplete: ...

###

# see `scipy.linalg.blas._type_conv`
def find_best_blas_type(
    arrays: Sequence[onp.ArrayND] = (), dtype: npt.DTypeLike | None = None
) -> (
    tuple[L["s"], np.dtype[np.float32], bool]
    | tuple[L["f"], np.dtype[np.float64], bool]
    | tuple[L["c"], np.dtype[np.complex64], bool]
    | tuple[L["z"], np.dtype[np.complex128], bool]
): ...

#
@overload
def get_blas_funcs(
    names: str, arrays: Sequence[onp.ArrayND] = (), dtype: npt.DTypeLike | None = None, ilp64: L["preferred"] | bool = False
) -> _FortranFunction: ...
@overload
def get_blas_funcs(
    names: _SequenceNotStr[str],
    arrays: Sequence[onp.ArrayND] = (),
    dtype: npt.DTypeLike | None = None,
    ilp64: L["preferred"] | bool = False,
) -> list[_FortranFunction]: ...
@overload
def get_blas_funcs(
    names: Iterable[str],
    arrays: Sequence[onp.ArrayND] = (),
    dtype: npt.DTypeLike | None = None,
    ilp64: L["preferred"] | bool = False,
) -> list[_FortranFunction] | _FortranFunction: ...
