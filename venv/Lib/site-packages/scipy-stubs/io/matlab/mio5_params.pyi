# This module is not meant for public use and will be removed in SciPy v2.0.0.

from typing import Self
from typing_extensions import deprecated

__all__ = ["MatlabFunction", "MatlabObject", "MatlabOpaque", "mat_struct"]

@deprecated("will be removed in SciPy v2.0.0")
class MatlabFunction:
    def __new__(cls, input_array: object) -> Self: ...

@deprecated("will be removed in SciPy v2.0.0")
class MatlabObject:
    def __new__(cls, input_array: object, classname: object = None) -> Self: ...
    def __array_finalize__(self, /, obj: Self) -> None: ...

@deprecated("will be removed in SciPy v2.0.0")
class MatlabOpaque:
    def __new__(cls, input_array: object) -> Self: ...

@deprecated("will be removed in SciPy v2.0.0")
class mat_struct: ...
