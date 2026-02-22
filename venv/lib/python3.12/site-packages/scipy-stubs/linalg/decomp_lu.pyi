# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

from ._misc import LinAlgWarning as _LinAlgWarning

__all__ = ["LinAlgWarning", "get_lapack_funcs", "lu", "lu_factor", "lu_solve"]

@deprecated("will be removed in SciPy v2.0.0")
class LinAlgWarning(_LinAlgWarning): ...

@deprecated("will be removed in SciPy v2.0.0")
def get_lapack_funcs(names: object, arrays: object = (), dtype: object = None, ilp64: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lu(
    a: object, permute_l: object = False, overwrite_a: object = False, check_finite: object = True, p_indices: object = False
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lu_factor(a: object, overwrite_a: object = False, check_finite: object = True) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lu_solve(
    lu_and_piv: object, b: object, trans: object = 0, overwrite_b: object = False, check_finite: object = True
) -> object: ...
