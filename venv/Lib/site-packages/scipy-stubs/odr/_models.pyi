# pyright: reportDeprecated=false

from typing import Final, type_check_only
from typing_extensions import deprecated

import optype.numpy as onp

from ._odrpack import Model

__all__ = ["Model", "exponential", "multilinear", "polynomial", "quadratic", "unilinear"]

@type_check_only
class _NamedModel(Model):
    name: Final[str]
    equ: Final[str]
    TeXequ: Final[str]

@type_check_only
class _SimpleModel(_NamedModel):
    def __init__(self, /) -> None: ...

###

@deprecated("`scipy.odr` is deprecated and will be removed in SciPy 1.19.0.")
class _MultilinearModel(_SimpleModel): ...

@deprecated("`scipy.odr` is deprecated and will be removed in SciPy 1.19.0.")
class _ExponentialModel(_SimpleModel): ...

@deprecated("`scipy.odr` is deprecated and will be removed in SciPy 1.19.0.")
class _UnilinearModel(_SimpleModel): ...

@deprecated("`scipy.odr` is deprecated and will be removed in SciPy 1.19.0.")
class _QuadraticModel(_SimpleModel): ...

@deprecated("`scipy.odr` is deprecated and will be removed in SciPy 1.19.0.")
def polynomial(order: onp.ToInt | onp.ToInt1D) -> _NamedModel: ...

multilinear: Final[_MultilinearModel] = ...
exponential: Final[_ExponentialModel] = ...
unilinear: Final[_UnilinearModel] = ...
quadratic: Final[_QuadraticModel] = ...
