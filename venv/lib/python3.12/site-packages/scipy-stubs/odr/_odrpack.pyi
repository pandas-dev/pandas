from collections.abc import Callable, Mapping
from typing import Concatenate, Final, Literal, TypeAlias, TypedDict, overload, type_check_only

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["ODR", "Data", "Model", "OdrError", "OdrStop", "OdrWarning", "Output", "RealData", "odr", "odr_error", "odr_stop"]

_ToIntScalar: TypeAlias = npc.integer | np.bool_
_ToFloatScalar: TypeAlias = npc.floating | _ToIntScalar

_Float1D: TypeAlias = onp.Array1D[np.float64]
_Float2D: TypeAlias = onp.Array2D[np.float64]
_FloatND: TypeAlias = onp.ArrayND[np.float64]
_FCN: TypeAlias = Callable[Concatenate[_Float1D, _FloatND, ...], onp.ArrayND[_ToFloatScalar]]

_01: TypeAlias = Literal[0, 1]  # noqa: PYI042
_012: TypeAlias = Literal[0, 1, 2]  # noqa: PYI042
_0123: TypeAlias = Literal[0, 1, 2, 3]  # noqa: PYI042

@type_check_only
class _FullOutput(TypedDict):
    delta: _Float1D
    eps: _Float1D
    xplus: _Float1D
    y: _Float1D
    res_var: float
    sum_square: float
    sum_square_delta: float
    sum_square_eps: float
    inc_condnum: float
    rel_error: float
    work: _Float1D
    work_ind: dict[str, int]
    iwork: onp.Array1D[np.int32]
    info: int

###

odr_error = OdrError
odr_stop = OdrStop

class OdrWarning(UserWarning): ...
class OdrError(Exception): ...
class OdrStop(Exception): ...

class Data:
    x: Final[onp.ArrayND[_ToFloatScalar]]
    y: Final[_ToFloatScalar | onp.ArrayND[_ToFloatScalar] | None]
    we: Final[_ToFloatScalar | onp.ArrayND[_ToFloatScalar] | None]
    wd: Final[_ToFloatScalar | onp.ArrayND[_ToFloatScalar] | None]
    fix: Final[onp.ArrayND[_ToIntScalar] | None]
    meta: Final[Mapping[str, object]]

    def __init__(
        self,
        /,
        x: onp.ToFloatND,
        y: onp.ToFloat | onp.ToFloatND | None = None,
        we: onp.ToFloat | onp.ToFloatND | None = None,
        wd: onp.ToFloat | onp.ToFloatND | None = None,
        fix: onp.ToIntND | None = None,
        meta: Mapping[str, object] | None = None,
    ) -> None: ...
    def set_meta(self, /, **kwds: object) -> None: ...

class RealData(Data):
    sx: Final[onp.ArrayND[_ToFloatScalar] | None]
    sy: Final[onp.ArrayND[_ToFloatScalar] | None]
    covx: Final[onp.ArrayND[_ToFloatScalar] | None]
    covy: Final[onp.ArrayND[_ToFloatScalar] | None]

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
        meta: Mapping[str, object] | None = None,
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
        meta: Mapping[str, object] | None = None,
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
        meta: Mapping[str, object] | None = None,
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
        meta: Mapping[str, object] | None = None,
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
        meta: Mapping[str, object] | None = None,
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
        meta: Mapping[str, object] | None = None,
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
        meta: Mapping[str, object] | None = None,
    ) -> None: ...

class Model:
    fcn: Final[_FCN]
    fjacb: Final[_FCN]
    fjacd: Final[_FCN]
    extra_args: Final[tuple[object, ...]]
    covx: Final[onp.ArrayND[_ToFloatScalar] | None]
    implicit: Final[onp.ToBool]
    meta: Final[Mapping[str, object]]

    def __init__(
        self,
        /,
        fcn: _FCN,
        fjacb: _FCN | None = None,
        fjacd: _FCN | None = None,
        extra_args: tuple[object, ...] | None = None,
        estimate: onp.ToFloat1D | None = None,
        implicit: onp.ToBool = 0,
        meta: Mapping[str, object] | None = None,
    ) -> None: ...
    def set_meta(self, /, **kwds: object) -> None: ...

class Output:
    beta: Final[onp.Array1D[_ToFloatScalar]]
    sd_beta: Final[onp.Array1D[_ToFloatScalar]]
    cov_beta: Final[onp.Array1D[_ToFloatScalar]]
    stopreason: Final[list[str]]

    def __init__(self, /, output: onp.ArrayND[_ToFloatScalar]) -> None: ...
    def pprint(self, /) -> None: ...

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
) -> tuple[_Float1D, _Float1D, _Float2D]: ...
@overload
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
) -> tuple[_Float1D, _Float1D, _Float2D, _FullOutput]: ...
