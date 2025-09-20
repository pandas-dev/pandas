from typing import TypeAlias, overload

import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["cho_factor", "cho_solve", "cho_solve_banded", "cholesky", "cholesky_banded"]

_Float2D: TypeAlias = onp.Array2D[npc.floating]
_FloatND: TypeAlias = onp.ArrayND[npc.floating]
_Complex2D: TypeAlias = onp.Array2D[npc.inexact]
_ComplexND: TypeAlias = onp.ArrayND[npc.inexact]

###

@overload
def cholesky(
    a: onp.ToFloatStrict2D, lower: onp.ToBool = False, overwrite_a: onp.ToBool = False, check_finite: onp.ToBool = True
) -> _Float2D: ...
@overload
def cholesky(
    a: onp.ToFloatND, lower: onp.ToBool = False, overwrite_a: onp.ToBool = False, check_finite: onp.ToBool = True
) -> _FloatND: ...
@overload
def cholesky(
    a: onp.ToComplexStrict2D, lower: onp.ToBool = False, overwrite_a: onp.ToBool = False, check_finite: onp.ToBool = True
) -> _Complex2D: ...
@overload
def cholesky(
    a: onp.ToComplexND, lower: onp.ToBool = False, overwrite_a: onp.ToBool = False, check_finite: onp.ToBool = True
) -> _ComplexND: ...

#
@overload
def cho_factor(
    a: onp.ToFloatND, lower: onp.ToBool = False, overwrite_a: onp.ToBool = False, check_finite: onp.ToBool = True
) -> tuple[_FloatND, bool]: ...
@overload
def cho_factor(
    a: onp.ToComplexND, lower: onp.ToBool = False, overwrite_a: onp.ToBool = False, check_finite: onp.ToBool = True
) -> tuple[_ComplexND, bool]: ...

#
@overload
def cho_solve(
    c_and_lower: tuple[onp.ToFloatStrict2D, onp.ToBool],
    b: onp.ToFloatStrict1D,
    overwrite_b: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> _Float2D: ...
@overload
def cho_solve(
    c_and_lower: tuple[onp.ToFloatND, onp.ToBool],
    b: onp.ToFloatND,
    overwrite_b: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> _FloatND: ...
@overload
def cho_solve(
    c_and_lower: tuple[onp.ToComplexStrict2D, onp.ToBool],
    b: onp.ToComplexStrict1D,
    overwrite_b: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> _Complex2D: ...
@overload
def cho_solve(
    c_and_lower: tuple[onp.ToComplexND, onp.ToBool],
    b: onp.ToComplexND,
    overwrite_b: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> _ComplexND: ...

#
@overload
def cholesky_banded(
    ab: onp.ToFloatStrict2D, overwrite_ab: onp.ToBool = False, lower: onp.ToBool = False, check_finite: onp.ToBool = True
) -> _Float2D: ...
@overload
def cholesky_banded(
    ab: onp.ToFloatND, overwrite_ab: onp.ToBool = False, lower: onp.ToBool = False, check_finite: onp.ToBool = True
) -> _FloatND: ...
@overload
def cholesky_banded(
    ab: onp.ToComplexStrict2D, overwrite_ab: onp.ToBool = False, lower: onp.ToBool = False, check_finite: onp.ToBool = True
) -> _Complex2D: ...
@overload
def cholesky_banded(
    ab: onp.ToComplexND, overwrite_ab: onp.ToBool = False, lower: onp.ToBool = False, check_finite: onp.ToBool = True
) -> _ComplexND: ...

#
@overload
def cho_solve_banded(
    cb_and_lower: tuple[onp.ToFloatStrict2D, onp.ToBool],
    b: onp.ToComplexStrict1D,
    overwrite_b: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> _Complex2D: ...
@overload
def cho_solve_banded(
    cb_and_lower: tuple[onp.ToFloatND, onp.ToBool],
    b: onp.ToComplexND,
    overwrite_b: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> _ComplexND: ...
@overload
def cho_solve_banded(
    cb_and_lower: tuple[onp.ToComplexStrict2D, onp.ToBool],
    b: onp.ToComplexStrict1D,
    overwrite_b: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> _Complex2D: ...
@overload
def cho_solve_banded(
    cb_and_lower: tuple[onp.ToComplexND, onp.ToBool],
    b: onp.ToComplexND,
    overwrite_b: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> _ComplexND: ...
