# pyright: reportDeprecated=false

from collections.abc import Callable, Iterable
from typing import Any, Concatenate, Final, Literal, TypeAlias, TypedDict, overload, type_check_only
from typing_extensions import deprecated

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["ODR", "Data", "Model", "OdrError", "OdrStop", "OdrWarning", "Output", "RealData", "odr", "odr_error", "odr_stop"]

# Technically this can be `np.floating | np.integer | np.bool`, but that'd
# force users to explicitly narrow the type in many places. In practice,
# `float64` is used most often, so we simplify this using the "any trick",
# even though it's technically type-unsafe this way.
_ToFloatScalar: TypeAlias = np.float64 | Any

_Float1D: TypeAlias = onp.Array1D[np.float64]
_Float2D: TypeAlias = onp.Array2D[np.float64]
_FloatND: TypeAlias = onp.ArrayND[np.float64]
_FCN: TypeAlias = Callable[Concatenate[_Float1D, _FloatND, ...], onp.ArrayND[_ToFloatScalar]]

_01: TypeAlias = Literal[0, 1]  # noqa: PYI042
_012: TypeAlias = Literal[0, 1, 2]  # noqa: PYI042
_0123: TypeAlias = Literal[0, 1, 2, 3]  # noqa: PYI042

# return value of `__odrpack.odr` with `full_output=False`
_RawOutput: TypeAlias = tuple[
    _Float1D,  # beta
    _Float1D,  # sd_beta
    _Float2D,  # cov_beta
]
_RawOutputFull: TypeAlias = tuple[
    _Float1D,  # beta
    _Float1D,  # sd_beta
    _Float2D,  # cov_beta
    _FullOutput,
]

@type_check_only
class _FullOutput(TypedDict):
    delta: _FloatND
    eps: _FloatND
    xplus: _FloatND
    y: _FloatND
    res_var: float
    sum_square: float
    sum_square_delta: float
    sum_square_eps: float
    inc_condnum: float
    rel_error: float
    work: _Float1D
    work_ind: dict[str, int]
    iwork: onp.Array1D[np.int32 | np.int64]
    info: int

###

@deprecated("`scipy.odr` is deprecated and will be removed in SciPy 1.19.0.")
class OdrWarning(UserWarning): ...

@deprecated("`scipy.odr` is deprecated and will be removed in SciPy 1.19.0.")
class OdrError(Exception): ...

@deprecated("`scipy.odr` is deprecated and will be removed in SciPy 1.19.0.")
class OdrStop(Exception): ...

odr_error = OdrError
odr_stop = OdrStop

@deprecated("`scipy.odr` is deprecated and will be removed in SciPy 1.19.0.")
class Data:
    x: Final[onp.ArrayND[_ToFloatScalar]]
    y: Final[onp.ArrayND[_ToFloatScalar] | _ToFloatScalar | None]
    we: onp.ArrayND[_ToFloatScalar] | _ToFloatScalar | None
    wd: onp.ArrayND[_ToFloatScalar] | _ToFloatScalar | None
    fix: Final[onp.ArrayND[npc.integer] | None]
    meta: Final[dict[str, Any]]

    def __init__(
        self,
        /,
        x: onp.ToFloatND,
        y: onp.ToFloat | onp.ToFloatND | None = None,
        we: onp.ToFloat | onp.ToFloatND | None = None,
        wd: onp.ToFloat | onp.ToFloatND | None = None,
        fix: onp.ToIntND | None = None,
        meta: dict[str, Any] | None = None,
    ) -> None: ...
    def set_meta(self, /, **kwds: object) -> None: ...

@deprecated("`scipy.odr` is deprecated and will be removed in SciPy 1.19.0.")
class RealData(Data):
    sx: Final[onp.ArrayND[_ToFloatScalar] | None]
    sy: Final[onp.ArrayND[_ToFloatScalar] | None]
    covx: Final[onp.ArrayND[_ToFloatScalar] | None]
    covy: Final[onp.ArrayND[_ToFloatScalar] | None]

    # readonly attributes that are dynamically computed in `RealData.__getattr__`
    wd: onp.ArrayND[_ToFloatScalar] | None  # pyrefly: ignore[bad-override]
    we: onp.ArrayND[_ToFloatScalar] | None  # pyrefly: ignore[bad-override]

    @overload
    def __init__(
        self,
        /,
        x: onp.ToFloatND,
        y: onp.ToFloat | onp.ToFloatND | None = None,
        sx: onp.ToFloatND | None = None,
        sy: onp.ToFloat | onp.ToFloatND | None = None,
        covx: None = None,
        covy: None = None,
        fix: onp.ToIntND | None = None,
        meta: dict[str, Any] | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        /,
        x: onp.ToFloatND,
        y: onp.ToFloat | onp.ToFloatND | None,
        sx: None,
        sy: onp.ToFloat | onp.ToFloatND | None,
        covx: onp.ToFloatND,
        covy: None = None,
        fix: onp.ToIntND | None = None,
        meta: dict[str, Any] | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        /,
        x: onp.ToFloatND,
        y: onp.ToFloat | onp.ToFloatND | None,
        sx: onp.ToFloatND | None,
        sy: None,
        covx: None,
        covy: onp.ToFloatND,
        fix: onp.ToIntND | None = None,
        meta: dict[str, Any] | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        /,
        x: onp.ToFloatND,
        y: onp.ToFloat | onp.ToFloatND | None,
        sx: None,
        sy: None,
        covx: onp.ToFloatND,
        covy: onp.ToFloatND,
        fix: onp.ToIntND | None = None,
        meta: dict[str, Any] | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        /,
        x: onp.ToFloatND,
        y: onp.ToFloat | onp.ToFloatND | None = None,
        sx: None = None,
        sy: onp.ToFloat | onp.ToFloatND | None = None,
        *,
        covx: onp.ToFloatND,
        covy: None = None,
        fix: onp.ToIntND | None = None,
        meta: dict[str, Any] | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        /,
        x: onp.ToFloatND,
        y: onp.ToFloat | onp.ToFloatND | None = None,
        sx: onp.ToFloatND | None = None,
        sy: None = None,
        *,
        covx: None = None,
        covy: onp.ToFloatND,
        fix: onp.ToIntND | None = None,
        meta: dict[str, Any] | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        /,
        x: onp.ToFloatND,
        y: onp.ToFloat | onp.ToFloatND | None = None,
        sx: None = None,
        sy: None = None,
        *,
        covx: onp.ToFloatND,
        covy: onp.ToFloatND,
        fix: onp.ToIntND | None = None,
        meta: dict[str, Any] | None = None,
    ) -> None: ...

@deprecated("`scipy.odr` is deprecated and will be removed in SciPy 1.19.0.")
class Model:
    fcn: Final[_FCN]
    fjacb: Final[_FCN]
    fjacd: Final[_FCN]
    extra_args: Final[tuple[Any, ...]]
    implicit: Final[bool | Literal[0, 1]]
    meta: Final[dict[str, Any]]

    def __init__(
        self,
        /,
        fcn: _FCN,
        fjacb: _FCN | None = None,
        fjacd: _FCN | None = None,
        extra_args: Iterable[object] | None = None,
        estimate: onp.ToFloat1D | None = None,
        implicit: bool | Literal[0, 1] = 0,
        meta: dict[str, Any] | None = None,
    ) -> None: ...
    def set_meta(self, /, **kwds: object) -> None: ...

@deprecated("`scipy.odr` is deprecated and will be removed in SciPy 1.19.0.")
class Output:
    beta: Final[onp.Array1D[np.float64]]
    sd_beta: Final[onp.Array1D[np.float64]]
    cov_beta: Final[onp.Array2D[np.float64]]

    # the following attributes are only available if `full_output=True` was used
    delta: Final[onp.Array1D[np.float64] | onp.Array2D[np.float64]]
    eps: Final[onp.Array1D[np.float64] | onp.Array2D[np.float64]]
    xplus: Final[onp.Array1D[np.float64] | onp.Array2D[np.float64]]
    y: Final[onp.Array1D[np.float64] | onp.Array2D[np.float64]]
    res_var: Final[float]
    sum_square: Final[float]
    sum_square_delta: Final[float]
    sum_square_eps: Final[float]
    inv_condnum: Final[float]
    rel_error: Final[float]
    work: Final[onp.Array1D[np.float64]]
    work_ind: Final[dict[str, int]]
    info: Final[int]
    stopreason: Final[list[str]]

    def __init__(self, /, output: _RawOutput | _RawOutputFull) -> None: ...
    def pprint(self, /) -> None: ...

@deprecated("`scipy.odr` is deprecated and will be removed in SciPy 1.19.0.")
class ODR:
    data: Final[Data]
    model: Final[Model]
    output: Output | None

    beta0: Final[onp.Array1D[_ToFloatScalar]]
    delta0: Final[onp.Array1D[_ToFloatScalar] | None]
    ifixx: Final[onp.Array1D[np.int32] | None]
    ifixb: Final[onp.Array1D[np.int32] | None]
    errfile: Final[str | None]
    rptfile: Final[str | None]
    ndigit: Final[int | None]
    taufac: Final[float | None]
    sstol: Final[float | None]
    partol: Final[float | None]
    stpb: Final[onp.Array1D[_ToFloatScalar] | None]
    stpd: Final[onp.Array1D[_ToFloatScalar] | None]
    sclb: Final[onp.Array1D[_ToFloatScalar] | None]
    scld: Final[onp.Array1D[_ToFloatScalar] | None]

    job: int | None
    iprint: int | None
    maxit: int | None
    work: onp.Array1D[np.float64] | None
    iwork: onp.Array1D[np.int32 | np.int64] | None

    def __init__(
        self,
        /,
        data: Data,
        model: Model,
        beta0: onp.ToFloat1D | None = None,
        delta0: onp.ToFloat1D | None = None,
        ifixb: onp.ToInt1D | None = None,
        ifixx: onp.ToIntND | None = None,
        job: int | None = None,
        iprint: int | None = None,
        errfile: str | None = None,
        rptfile: str | None = None,
        ndigit: int | None = None,
        taufac: float | None = None,  # = 1
        sstol: float | None = None,  # = eps**(1/2)
        partol: float | None = None,  # = eps**(2/3) (explicit), = eps**(1/3) (implicit)
        maxit: int | None = None,  # = 10
        stpb: onp.ToFloat1D | None = None,
        stpd: onp.ToFloatND | None = None,
        sclb: onp.ToFloat1D | None = None,
        scld: onp.ToFloatND | None = None,
        work: onp.ArrayND[np.float64] | None = None,
        iwork: onp.ArrayND[np.int32 | np.int64] | None = None,
        overwrite: bool = False,
    ) -> None: ...
    def set_job(
        self,
        /,
        fit_type: _012 | None = None,
        deriv: _0123 | None = None,
        var_calc: _012 | None = None,
        del_init: _01 | None = None,
        restart: _01 | None = None,
    ) -> None: ...
    def set_iprint(
        self,
        /,
        init: _012 | None = None,
        so_init: _012 | None = None,
        iter: _012 | None = None,
        so_iter: _012 | None = None,
        iter_step: _012 | None = None,
        final: _012 | None = None,
        so_final: _012 | None = None,
    ) -> None: ...
    def run(self, /) -> Output: ...
    def restart(self, /, iter: int | None = None) -> Output: ...

@overload
@deprecated("`scipy.odr` is deprecated and will be removed in SciPy 1.19.0.")
def odr(
    fcn: _FCN,
    beta0: onp.ToFloat1D,
    y: onp.ToFloat | onp.ToFloatND,
    x: onp.ToFloatND,
    we: onp.ToFloat | onp.ToFloatND | None = None,
    wd: onp.ToFloat | onp.ToFloatND | None = None,
    fjacb: _FCN | None = None,
    fjacd: _FCN | None = None,
    extra_args: tuple[object, ...] | None = None,
    ifixx: onp.ToIntND | None = None,
    ifixb: onp.ToInt1D | None = None,
    job: int = 0,
    iprint: int = 0,
    errfile: str | None = None,
    rptfile: str | None = None,
    ndigit: int = 0,
    taufac: float = 0.0,
    sstol: float = -1.0,
    partol: float = -1.0,
    maxit: int = -1,
    stpb: onp.ToFloat1D | None = None,
    stpd: onp.ToFloatND | None = None,
    sclb: onp.ToFloat1D | None = None,
    scld: onp.ToFloatND | None = None,
    work: onp.ArrayND[np.float64] | None = None,
    iwork: onp.ArrayND[np.int32 | np.int64] | None = None,
    full_output: onp.ToFalse = 0,
) -> _RawOutput: ...
@overload
@deprecated("`scipy.odr` is deprecated and will be removed in SciPy 1.19.0.")
def odr(
    fcn: _FCN,
    beta0: onp.ToFloat1D,
    y: onp.ToFloat | onp.ToFloatND,
    x: onp.ToFloatND,
    we: onp.ToFloat | onp.ToFloatND | None = None,
    wd: onp.ToFloat | onp.ToFloatND | None = None,
    fjacb: _FCN | None = None,
    fjacd: _FCN | None = None,
    extra_args: tuple[object, ...] | None = None,
    ifixx: onp.ToIntND | None = None,
    ifixb: onp.ToInt1D | None = None,
    job: int = 0,
    iprint: int = 0,
    errfile: str | None = None,
    rptfile: str | None = None,
    ndigit: int = 0,
    taufac: float = 0.0,
    sstol: float = -1.0,
    partol: float = -1.0,
    maxit: int = -1,
    stpb: onp.ToFloat1D | None = None,
    stpd: onp.ToFloatND | None = None,
    sclb: onp.ToFloat1D | None = None,
    scld: onp.ToFloatND | None = None,
    work: onp.ArrayND[np.float64] | None = None,
    iwork: onp.ArrayND[np.int32 | np.int64] | None = None,
    *,
    full_output: onp.ToTrue,
) -> _RawOutputFull: ...
