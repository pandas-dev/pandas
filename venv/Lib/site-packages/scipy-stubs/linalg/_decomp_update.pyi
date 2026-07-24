from typing import Literal, overload

import numpy as np
import optype.numpy as onp

__all__ = ["qr_delete", "qr_insert", "qr_update"]

###

type _FloatND = onp.ArrayND[np.float32 | np.float64]
type _FloatQR = tuple[_FloatND, _FloatND]

type _ComplexND = onp.ArrayND[np.complex64 | np.complex128]
type _ComplexQR = _FloatQR | tuple[_ComplexND, _ComplexND]

type _Which = Literal["row", "col"]

###

@overload
def qr_delete(
    Q: onp.ToFloatND,
    R: onp.ToFloatND,
    k: onp.ToJustInt,
    p: onp.ToJustInt = 1,
    which: _Which = "row",
    overwrite_qr: bool = False,
    check_finite: bool = True,
) -> _FloatQR: ...
@overload
def qr_delete(
    Q: onp.ToComplexND,
    R: onp.ToComplexND,
    k: onp.ToJustInt,
    p: onp.ToJustInt = 1,
    which: _Which = "row",
    overwrite_qr: bool = False,
    check_finite: bool = True,
) -> _ComplexQR: ...

#
@overload
def qr_insert(
    Q: onp.ToFloatND,
    R: onp.ToFloatND,
    u: onp.ToFloatND,
    k: onp.ToJustInt,
    which: _Which = "row",
    rcond: onp.ToFloat | None = None,
    overwrite_qru: bool = False,
    check_finite: bool = True,
) -> _FloatQR: ...
@overload
def qr_insert(
    Q: onp.ToComplexND,
    R: onp.ToComplexND,
    u: onp.ToComplexND,
    k: onp.ToJustInt,
    which: _Which = "row",
    rcond: onp.ToFloat | None = None,
    overwrite_qru: bool = False,
    check_finite: bool = True,
) -> _ComplexQR: ...

#
@overload
def qr_update(
    Q: onp.ToFloatND,
    R: onp.ToFloatND,
    u: onp.ToFloatND,
    v: onp.ToFloatND,
    overwrite_qruv: bool = False,
    check_finite: bool = True,
) -> _FloatQR: ...
@overload
def qr_update(
    Q: onp.ToComplexND,
    R: onp.ToComplexND,
    u: onp.ToComplexND,
    v: onp.ToComplexND,
    overwrite_qruv: bool = False,
    check_finite: bool = True,
) -> _ComplexQR: ...
