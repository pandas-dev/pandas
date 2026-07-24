import inspect
from builtins import dict as _dict  # alias to avoid conflicts with attribute name
from collections.abc import Callable, Generator, Iterator
from contextlib import _GeneratorContextManager
from inspect import Signature, getfullargspec as getfullargspec, iscoroutinefunction as iscoroutinefunction
from re import Pattern
from typing import Any, Final, Literal, TypeVar
from typing_extensions import ParamSpec

_C = TypeVar("_C", bound=Callable[..., Any])
_Func = TypeVar("_Func", bound=Callable[..., Any])
_T = TypeVar("_T")
_P = ParamSpec("_P")

DEF: Final[Pattern[str]]
POS: Final[Literal[inspect._ParameterKind.POSITIONAL_OR_KEYWORD]]
EMPTY: Final[type[inspect._empty]]

class FunctionMaker:
    args: list[str]
    varargs: str | None
    varkw: str | None
    defaults: tuple[Any, ...] | None
    kwonlyargs: list[str]
    kwonlydefaults: _dict[str, Any] | None
    shortsignature: str | None
    name: str
    doc: str | None
    module: str | None
    annotations: _dict[str, Any]
    signature: str
    dict: _dict[str, Any]
    def __init__(
        self,
        func: Callable[..., Any] | None = ...,
        name: str | None = ...,
        signature: str | None = ...,
        defaults: tuple[Any, ...] | None = ...,
        doc: str | None = ...,
        module: str | None = ...,
        funcdict: _dict[str, Any] | None = ...,
    ) -> None: ...
    def update(self, func: Any, **kw: Any) -> None: ...
    def make(
        self, src_templ: str, evaldict: _dict[str, Any] | None = ..., addsource: bool = ..., **attrs: Any
    ) -> Callable[..., Any]: ...
    @classmethod
    def create(
        cls,
        obj: Any,
        body: str,
        evaldict: _dict[str, Any],
        defaults: tuple[Any, ...] | None = ...,
        doc: str | None = ...,
        module: str | None = ...,
        addsource: bool = ...,
        **attrs: Any,
    ) -> Callable[..., Any]: ...

def fix(args: tuple[Any, ...], kwargs: dict[str, Any], sig: Signature) -> tuple[tuple[Any, ...], dict[str, Any]]: ...
def decorate(func: _Func, caller: Callable[..., Any], extras: tuple[Any, ...] = ..., kwsyntax: bool = False) -> _Func: ...
def decoratorx(caller: Callable[..., Any]) -> Callable[..., Any]: ...
def decorator(
    caller: Callable[..., Any], _func: Callable[..., Any] | None = None, kwsyntax: bool = False
) -> Callable[[Callable[..., Any]], Callable[..., Any]]: ...

class ContextManager(_GeneratorContextManager[_T]):
    def __init__(self, g: Callable[..., Generator[_T]], *a: Any, **k: Any) -> None: ...
    def __call__(self, func: _C) -> _C: ...

def contextmanager(func: Callable[_P, Iterator[_T]]) -> Callable[_P, ContextManager[_T]]: ...
def append(a: type, vancestors: list[type]) -> None: ...
def dispatch_on(*dispatch_args: Any) -> Callable[[Callable[..., Any]], Callable[..., Any]]: ...
