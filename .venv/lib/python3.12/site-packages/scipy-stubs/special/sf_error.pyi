# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

__all__ = ["SpecialFunctionError", "SpecialFunctionWarning"]

@deprecated("will be removed in SciPy v2.0.0")
class SpecialFunctionWarning(Warning): ...

@deprecated("will be removed in SciPy v2.0.0")
class SpecialFunctionError(Exception): ...
