import ctypes as ct
from _ctypes import CFuncPtr as _CFuncPtr
from types import ModuleType
from typing import ClassVar, Generic, Literal, NoReturn, Protocol, Self, TypeAlias, final, overload, type_check_only
from typing_extensions import CapsuleType as PyCapsule, TypeVar, TypeVarTuple, Unpack, override

# some quick interfaces for the relevant `cffi` types

@type_check_only
@final
class _CFFIBackendType(Protocol):
    cname: str
    kind: Literal[
        "primitive",
        "bool",
        "int",
        "float",
        "char",
        "byte",
        "pointer",
        "charp",
        "bytep",
        "voidp",
        "generic",
        "struct",
        "union",
        "enum",
        "anonymous",
        "typedef",
        "function",
    ]

_CTs = TypeVarTuple("_CTs", default=Unpack[tuple[_CFFIType, ...]])
_CT_co = TypeVar("_CT_co", covariant=True, bound=_CFFIType, default=_CFFIType)

@type_check_only
class _CFFIType(Protocol):
    is_array_type: ClassVar[bool]
    is_raw_function: ClassVar[bool]

    def is_integer_type(self, /) -> bool: ...
    def has_c_name(self, /) -> bool: ...
    def get_c_name(self, /, replace_with: str = "", context: str = "a C file", quals: int = 0) -> str: ...
    def get_cached_btype(self, /, ffi: object, finishlist: list[object], can_delay: bool = False) -> _CFFIBackendType: ...

    # virtual
    def build_backend_type(self, /, ffi: object, finishlist: list[object]) -> _CFFIBackendType: ...
    @property
    def c_name_with_marker(self, /) -> str: ...

@type_check_only
@final
class _CFFIVoid(_CFFIType, Protocol):
    is_array_type: ClassVar[bool] = False
    is_raw_function: ClassVar[bool] = False

    def __init__(self, /) -> None: ...

@type_check_only
class _CFFIFunc(_CFFIType, Protocol[_CT_co, *_CTs]):
    is_array_type: ClassVar[bool] = False

    @property
    def args(self, /) -> tuple[*_CTs]: ...
    @property
    def result(self, /) -> _CT_co: ...
    @property
    def ellipsis(self, /) -> bool: ...
    @property
    def abi(self, /) -> int | str | None: ...
    def __init__(self, /, args: tuple[*_CTs], result: _CT_co, ellipsis: bool, abi: int | None = None) -> None: ...

@type_check_only
@final
class _CFFIFuncPtr(_CFFIFunc[_CT_co, *_CTs], Protocol[_CT_co, *_CTs]):
    is_raw_function: ClassVar[bool] = False

    def as_raw_function(self, /) -> _CFFIFunc[_CT_co, *_CTs]: ...

@type_check_only
@final
class _CFFIPointerType(_CFFIType, Protocol[_CT_co]):
    is_array_type: ClassVar[bool] = False
    is_raw_function: ClassVar[bool] = False

    @property
    def totype(self, /) -> _CT_co: ...
    @property
    def quals(self, /) -> int: ...
    def __init__(self, /, totype: _CT_co, quals: int = 0) -> None: ...

_CFFIVoidP: TypeAlias = _CFFIPointerType[_CFFIVoid]

# helper aliases

_Function: TypeAlias = PyCapsule | PyCFuncPtr | _CFFIFuncPtr | CData
_UserData: TypeAlias = PyCapsule | ct.c_void_p | _CFFIVoidP

_FuncT_co = TypeVar("_FuncT_co", bound=_Function, covariant=True, default=_Function)
_DataT = TypeVar("_DataT", bound=_UserData | None)
_DataT_co = TypeVar("_DataT_co", bound=_UserData | None, covariant=True, default=None)

ffi: Literal[False] | None

# public api

@final
class PyCFuncPtr(_CFuncPtr): ...

@final
class CData: ...

class LowLevelCallable(tuple[PyCapsule, _FuncT_co, _DataT_co], Generic[_FuncT_co, _DataT_co]):
    __slots__: ClassVar[tuple[()]] = ()

    @property
    def function(self, /) -> _FuncT_co: ...
    @property
    def user_data(self, /) -> _DataT_co: ...
    @property
    def signature(self, /) -> str: ...
    @overload
    def __new__(cls, function: Self, user_data: _DataT_co | None = None, signature: str | None = None) -> Self: ...
    @overload
    def __new__(cls, function: _FuncT_co, user_data: _DataT_co, signature: str | None = None) -> Self: ...
    @classmethod
    @overload
    def from_cython(
        cls, module: ModuleType, name: str, user_data: None = None, signature: str | None = None
    ) -> LowLevelCallable[PyCapsule, None]: ...
    @classmethod
    @overload
    def from_cython(
        cls, module: ModuleType, name: str, user_data: _DataT, signature: str | None = None
    ) -> LowLevelCallable[PyCapsule, _DataT]: ...

    # NOTE: `__getitem__` will always raise a `ValueError`
    @override
    def __getitem__(self, idx: object, /) -> NoReturn: ...
