from typing import Literal, TypeAlias, overload

import numpy as np
import optype.numpy as onp

__all__ = ["qr_delete", "qr_insert", "qr_update"]

_FloatND: TypeAlias = onp.ArrayND[np.float32 | np.float64]
_FloatQR: TypeAlias = tuple[_FloatND, _FloatND]

_ComplexND: TypeAlias = onp.ArrayND[np.complex64 | np.complex128]
_ComplexQR: TypeAlias = _FloatQR | tuple[_ComplexND, _ComplexND]

_Which: TypeAlias = Literal["row", "col"]

###

@overload
def qr_delete(
    Q: onp.ToFloatND,
    R: onp.ToFloatND,
    k: onp.ToJustInt,
    p: onp.ToJustInt = 1,
    which: _Which = "row",
    overwrite_qr: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> _FloatQR: ...
@overload
def qr_delete(
    Q: onp.ToComplexND,
    R: onp.ToComplexND,
    k: onp.ToJustInt,
    p: onp.ToJustInt = 1,
    which: _Which = "row",
    overwrite_qr: onp.ToBool = False,
    check_finite: onp.ToBool = True,
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
    overwrite_qru: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> _FloatQR: ...
@overload
def qr_insert(
    Q: onp.ToComplexND,
    R: onp.ToComplexND,
    u: onp.ToComplexND,
    k: onp.ToJustInt,
    which: _Which = "row",
    rcond: onp.ToFloat | None = None,
    overwrite_qru: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> _ComplexQR: ...

#
@overload
def qr_update(
    Q: onp.ToFloatND,
    R: onp.ToFloatND,
    u: onp.ToFloatND,
    v: onp.ToFloatND,
    overwrite_qruv: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> _FloatQR: ...
@overload
def qr_update(
    Q: onp.ToComplexND,
    R: onp.ToComplexND,
    u: onp.ToComplexND,
    v: onp.ToComplexND,
    overwrite_qruv: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> _ComplexQR: ...
