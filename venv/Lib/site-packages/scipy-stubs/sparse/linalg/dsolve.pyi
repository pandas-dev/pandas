# This module is not meant for public use and will be removed in SciPy v2.0.0.
from types import ModuleType
from typing import final
from typing_extensions import deprecated

from . import _dsolve

__all__ = ["MatrixRankWarning", "SuperLU", "factorized", "spilu", "splu", "spsolve", "spsolve_triangular", "test", "use_solver"]

test: ModuleType

@deprecated("will be removed in SciPy v2.0.0")
class MatrixRankWarning(_dsolve.MatrixRankWarning): ...

@final
@deprecated("will be removed in SciPy v2.0.0")
# pyrefly: ignore [invalid-inheritance]
class SuperLU(_dsolve.SuperLU): ...  # type: ignore[misc]  # pyright: ignore[reportGeneralTypeIssues]  # ty: ignore[subclass-of-final-class]

@deprecated("will be removed in SciPy v2.0.0")
def use_solver(**kwargs: object) -> None: ...
@deprecated("will be removed in SciPy v2.0.0")
def spsolve(A: object, b: object, permc_spec: object = None, use_umfpack: object = True) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def splu(
    A: object,
    permc_spec: object = None,
    diag_pivot_thresh: object = None,
    relax: object = None,
    panel_size: object = None,
    options: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def spilu(
    A: object,
    drop_tol: object = None,
    fill_factor: object = None,
    drop_rule: object = None,
    permc_spec: object = None,
    diag_pivot_thresh: object = None,
    relax: object = None,
    panel_size: object = None,
    options: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def factorized(A: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def spsolve_triangular(
    A: object,
    b: object,
    lower: object = True,
    overwrite_A: object = False,
    overwrite_b: object = False,
    unit_diagonal: object = False,
) -> object: ...
