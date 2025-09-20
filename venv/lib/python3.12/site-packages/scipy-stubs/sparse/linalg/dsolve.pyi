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
class SuperLU(_dsolve.SuperLU): ...  # type: ignore[misc]  # pyright: ignore[reportGeneralTypeIssues]

@deprecated("will be removed in SciPy v2.0.0")
def use_solver(**kwargs: object) -> None: ...
@deprecated("will be removed in SciPy v2.0.0")
def spsolve(A: object, b: object, permc_spec: object = ..., use_umfpack: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def splu(
    A: object,
    permc_spec: object = ...,
    diag_pivot_thresh: object = ...,
    relax: object = ...,
    panel_size: object = ...,
    options: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def spilu(
    A: object,
    drop_tol: object = ...,
    fill_factor: object = ...,
    drop_rule: object = ...,
    permc_spec: object = ...,
    diag_pivot_thresh: object = ...,
    relax: object = ...,
    panel_size: object = ...,
    options: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def factorized(A: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def spsolve_triangular(
    A: object, b: object, lower: object = ..., overwrite_A: object = ..., overwrite_b: object = ..., unit_diagonal: object = ...
) -> object: ...
