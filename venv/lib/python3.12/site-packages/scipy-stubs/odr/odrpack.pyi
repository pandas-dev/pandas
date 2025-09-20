# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing import Any, Final
from typing_extensions import deprecated

from ._odrpack import ODR as _ODR, Data as _Data, Model as _Model, Output as _Output, RealData as _RealData

__all__ = ["ODR", "Data", "Model", "OdrError", "OdrStop", "OdrWarning", "Output", "RealData", "odr", "odr_error", "odr_stop"]

__MESSAGE: Final = "will be removed in SciPy v2.0.0"

@deprecated(__MESSAGE)
class OdrWarning(UserWarning): ...

@deprecated(__MESSAGE)
class OdrError(Exception): ...

@deprecated(__MESSAGE)
class OdrStop(Exception): ...

@deprecated(__MESSAGE)
class Data(_Data): ...

@deprecated(__MESSAGE)
class RealData(_RealData): ...

@deprecated(__MESSAGE)
class Model(_Model): ...

@deprecated(__MESSAGE)
class Output(_Output): ...

@deprecated(__MESSAGE)
class ODR(_ODR): ...

@deprecated(__MESSAGE)
def odr(
    fcn: object,
    beta0: object,
    y: object,
    x: object,
    we: object | None = None,
    wd: object | None = None,
    fjacb: object | None = None,
    fjacd: object | None = None,
    extra_args: tuple[object, ...] | None = None,
    ifixx: object | None = None,
    ifixb: object | None = None,
    job: int = 0,
    iprint: int = 0,
    errfile: str | None = None,
    rptfile: str | None = None,
    ndigit: int = 0,
    taufac: float = 0.0,
    sstol: float = -1.0,
    partol: float = -1.0,
    maxit: int = -1,
    stpb: object | None = None,
    stpd: object | None = None,
    sclb: object | None = None,
    scld: object | None = None,
    work: object | None = None,
    iwork: object | None = None,
    full_output: int = 0,
) -> Any: ...

odr_error = OdrError  # pyright: ignore[reportDeprecated]
odr_stop = OdrStop  # pyright: ignore[reportDeprecated]
