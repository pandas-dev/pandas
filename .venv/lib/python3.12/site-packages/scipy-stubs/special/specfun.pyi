# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

__all__ = ["lqmn", "pbdv"]

# originally defined in scipy/special/_specfun.pyx
@deprecated("will be removed in SciPy v2.0.0")
def pbdv(v: object, x: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lqmn(m: object, n: object, z: object) -> object: ...
