import ctypes as ct
from ctypes import _CData as CType  # noqa: PLC2701
from typing import Any
from typing_extensions import TypeAliasType, TypeVar

__all__ = "CScalar", "CType"

_T = TypeVar("_T", default=Any)
CScalar = TypeAliasType("CScalar", ct._SimpleCData[_T], type_params=(_T,))  # noqa: SLF001
