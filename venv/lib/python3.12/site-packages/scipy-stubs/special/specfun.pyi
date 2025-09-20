# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

__all__ = ["clpmn", "lpmn", "lpn", "lqmn", "pbdv"]

# originally defined in scipy/special/_specfun.pyx
@deprecated("will be removed in SciPy v2.0.0")
def pbdv(v: object, x: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def clpmn(m: object, n: object, z: object, type: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lpmn(m: object, n: object, z: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lpn(n: object, z: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lqmn(m: object, n: object, z: object) -> object: ...
