from typing import Final, Literal, TypedDict, overload, type_check_only

import numpy as np
import optype.numpy as onp

from ._optimize import OptimizeResult as _OptimizeResult

@type_check_only
class _CommonOptions(TypedDict, total=False):
    maximize: onp.ToBool
    rng: onp.random.ToRNG | None
    partial_match: onp.ToInt2D | None

@type_check_only
class _FAQOptions(_CommonOptions, TypedDict, total=False):
    P0: onp.ToFloat2D | Literal["barycenter", "randomized"]
    shuffle: onp.ToBool
    maxiter: onp.ToJustInt
    tol: onp.ToFloat

@type_check_only
class _2OptOptions(_CommonOptions, TypedDict, total=False):
    partial_guess: onp.ToInt2D | None

###

QUADRATIC_ASSIGNMENT_METHODS: Final = ["faq", "2opt"]

class OptimizeResult(_OptimizeResult):
    col_ind: onp.Array1D[np.intp]
    fun: float | np.float64
    nit: int

@overload
def quadratic_assignment(
    A: onp.ToFloat2D, B: onp.ToFloat2D, method: Literal["faq"] = "faq", options: _FAQOptions | None = None
) -> OptimizeResult: ...
@overload
def quadratic_assignment(
    A: onp.ToFloat2D, B: onp.ToFloat2D, method: Literal["2opt"], options: _2OptOptions | None = None
) -> OptimizeResult: ...
