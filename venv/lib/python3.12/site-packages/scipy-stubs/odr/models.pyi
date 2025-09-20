# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

from ._models import exponential, multilinear, quadratic, unilinear
from ._odrpack import Model as _Model

__all__ = ["Model", "exponential", "multilinear", "polynomial", "quadratic", "unilinear"]

@deprecated("will be removed in SciPy v2.0.0")
class Model(_Model): ...

@deprecated("will be removed in SciPy v2.0.0")
def polynomial(order: object) -> Model: ...  # pyright: ignore[reportDeprecated]
