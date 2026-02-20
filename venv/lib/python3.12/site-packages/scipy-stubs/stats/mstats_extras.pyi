# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

__all__ = [
    "compare_medians_ms",
    "hdmedian",
    "hdquantiles",
    "hdquantiles_sd",
    "idealfourths",
    "median_cihs",
    "mjci",
    "mquantiles_cimj",
    "rsh",
    "trimmed_mean_ci",
]

@deprecated("will be removed in SciPy v2.0.0")
def hdquantiles(data: object, prob: object = (0.25, 0.5, 0.75), axis: object = None, var: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def hdmedian(data: object, axis: int = -1, var: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def hdquantiles_sd(data: object, prob: object = (0.25, 0.5, 0.75), axis: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def trimmed_mean_ci(
    data: object, limits: object = (0.2, 0.2), inclusive: object = (True, True), alpha: object = 0.05, axis: object = None
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def mjci(data: object, prob: object = (0.25, 0.5, 0.75), axis: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def mquantiles_cimj(data: object, prob: object = (0.25, 0.5, 0.75), alpha: object = 0.05, axis: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def median_cihs(data: object, alpha: object = 0.05, axis: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def compare_medians_ms(group_1: object, group_2: object, axis: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def idealfourths(data: object, axis: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def rsh(data: object, points: object = None) -> object: ...
