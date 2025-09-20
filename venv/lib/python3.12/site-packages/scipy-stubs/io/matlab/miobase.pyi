# This module is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

__all__ = ["MatReadError", "MatReadWarning", "MatWriteError"]

@deprecated("will be removed in SciPy v2.0.0")
class MatReadError(Exception): ...

@deprecated("will be removed in SciPy v2.0.0")
class MatWriteError(Exception): ...

@deprecated("will be removed in SciPy v2.0.0")
class MatReadWarning(UserWarning): ...
