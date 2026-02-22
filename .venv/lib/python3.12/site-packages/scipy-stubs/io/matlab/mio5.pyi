# This module is not meant for public use and will be removed in SciPy v2.0.0.

from typing import Self
from typing_extensions import deprecated

import numpy as np

__all__ = ["MatReadError", "MatReadWarning", "MatWriteError", "MatlabFunction", "MatlabObject", "mat_struct", "varmats_from_mat"]

@deprecated("will be removed in SciPy v2.0.0")
class MatReadError(Exception): ...

@deprecated("will be removed in SciPy v2.0.0")
class MatWriteError(Exception): ...

@deprecated("will be removed in SciPy v2.0.0")
class MatReadWarning(UserWarning): ...

@deprecated("will be removed in SciPy v2.0.0")
class MatlabObject(np.ndarray[tuple[int, ...], np.dtype[np.generic]]):
    def __new__(cls, input_array: object, classname: object = None) -> Self: ...

@deprecated("will be removed in SciPy v2.0.0")
class MatlabFunction(MatlabObject):  # pyright: ignore[reportDeprecated]  # ty: ignore[deprecated]
    def __new__(cls, input_array: object) -> Self: ...

@deprecated("will be removed in SciPy v2.0.0")
class mat_struct: ...

@deprecated("will be removed in SciPy v2.0.0")
def varmats_from_mat(file_obj: object) -> object: ...
