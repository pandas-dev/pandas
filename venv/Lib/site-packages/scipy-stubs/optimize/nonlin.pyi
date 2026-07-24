# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

from ._nonlin import BroydenFirst as _BroydenFirst, InverseJacobian as _InverseJacobian, KrylovJacobian as _KrylovJacobian

__all__ = [
    "BroydenFirst",
    "InverseJacobian",
    "KrylovJacobian",
    "anderson",
    "broyden1",
    "broyden2",
    "diagbroyden",
    "excitingmixing",
    "linearmixing",
    "newton_krylov",
]

@deprecated("will be removed in SciPy v2.0.0")
class BroydenFirst(_BroydenFirst): ...

@deprecated("will be removed in SciPy v2.0.0")
class InverseJacobian(_InverseJacobian): ...

@deprecated("will be removed in SciPy v2.0.0")
class KrylovJacobian(_KrylovJacobian): ...

# Deprecated
broyden1: object
broyden2: object
anderson: object
linearmixing: object
diagbroyden: object
excitingmixing: object
newton_krylov: object
