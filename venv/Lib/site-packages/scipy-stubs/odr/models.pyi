# This module is not meant for public use and will be removed in SciPy v2.0.0.
from ._models import exponential, multilinear, polynomial, quadratic, unilinear
from ._odrpack import Model

__all__ = ["Model", "exponential", "multilinear", "polynomial", "quadratic", "unilinear"]
