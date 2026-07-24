import ctypes as ct
import sys
from ctypes import _CData as CType  # noqa: PLC2701
from typing import Any, TypeAliasType

if sys.version_info >= (3, 13):
    from typing import TypeVar
else:
    from typing_extensions import TypeVar

__all__ = "CScalar", "CType"

_T = TypeVar("_T", default=Any)
CScalar = TypeAliasType("CScalar", ct._SimpleCData[_T], type_params=(_T,))  # noqa: SLF001, UP040
