from _typeshed import Incomplete
from collections.abc import Callable
from ctypes import _CFunctionType
import functools
from typing import Any, Generic, TypeVar

from _cffi_backend import _CDataBase as CData
from typing_extensions import ParamSpec

from .compiler import Compiler
from .types import Type

_SignatureT_co = TypeVar(
    "_SignatureT_co",
    bound=Callable[..., object],
    default=Callable[..., Any],
    covariant=True,
)
_ParamsT = ParamSpec("_ParamsT")
_ReturnT = TypeVar("_ReturnT")

class CFunc(Generic[_SignatureT_co]):
    __name__: str
    __qualname__: str
    __wrapped__: _SignatureT_co

    def __init__(
        self,
        /,
        pyfunc: _SignatureT_co,
        sig: tuple[tuple[Type, ...], Type | None],
        locals: dict[str, Any],
        options: dict[str, Incomplete],  # TODO
        pipeline_class: type[Compiler] = Compiler,
    ) -> None: ...
    def __call__(
        self: CFunc[Callable[_ParamsT, _ReturnT]],
        /,
        *args: _ParamsT.args,
        **kwargs: _ParamsT.kwargs,
    ) -> _ReturnT: ...

    #
    def enable_caching(self) -> None: ...
    def compile(self) -> None: ...
    def inspect_llvm(self) -> str: ...

    #
    @property
    def native_name(self) -> str: ...
    @property
    def address(self) -> int | None: ...
    @property
    def cache_hits(self) -> int: ...
    @functools.cached_property
    def cffi(self) -> CData: ...
    @functools.cached_property
    def ctypes(self) -> _CFunctionType: ...
