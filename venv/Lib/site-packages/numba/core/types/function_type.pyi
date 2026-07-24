import abc
from collections.abc import Collection, Sequence
from typing import (
    Any,
    ClassVar,
    Final,
    Iterable,
    Never,
    TypeAlias,
    runtime_checkable,
    type_check_only,
)

from typing_extensions import Protocol, Self, TypeVar, override

from numba.core.compiler import CompileResult
from numba.core.dispatcher import Dispatcher
from numba.core.typing.context import Context
from numba.core.typing.templates import Signature

from .abstract import Type

__all__ = [
    "FunctionType",
    "UndefinedFunctionType",
    "FunctionPrototype",
    "WrapperAddressProtocol",
    "CompileResultWAP",
]

_T = TypeVar("_T")
_T_co = TypeVar("_T_co", covariant=True)

@type_check_only
class _CanUnliteral(Protocol[_T_co]):
    def __unliteral__(self, /) -> _T_co: ...

@type_check_only
class _HasLiteralType(Protocol[_T_co]):
    @property
    def literal_type(self, /) -> _T_co: ...

# represents the return type of `types.unliteral` for literals
_LiteralLike: TypeAlias = _CanUnliteral[_T] | _HasLiteralType[_T]

###

class FunctionType(Type):
    cconv: ClassVar[None] = None

    nargs: Final[int]
    signature: Final[Signature]
    ftype: Final[FunctionPrototype]
    _key: str

    @override
    def __init__(self, signature: Signature | _LiteralLike[Signature]) -> None: ...

    #
    @property
    @override
    def key(self) -> str: ...
    @property  # type: ignore[misc]
    @override
    def name(self) -> str: ...

    #
    def is_precise(self) -> bool: ...
    def get_precise(self) -> FunctionType: ...
    def dump(self, tab: str = "") -> None: ...
    def get_call_type(
        self,
        context: Context,
        args: Sequence[Type],
        kws: dict[Never, Never] | None,  # empty dict or None
    ) -> Signature: ...
    def check_signature(self, other_sig: Signature) -> bool: ...
    def unify(self, context: Context, other: Type) -> Self | None: ...

class UndefinedFunctionType(FunctionType):
    dispatchers: Final[Collection[Dispatcher]]

    @override
    def __init__(
        self,
        nargs: int,
        dispatchers: Collection[Dispatcher],
    ) -> None: ...

class FunctionPrototype(Type):
    cconv: ClassVar[None] = None
    rtype: Type
    atypes: tuple[Type, ...]

    @override
    def __init__(self, rtype: Type, atypes: Iterable[Type]) -> None: ...

    #
    @property
    @override
    def key(self) -> str: ...

@runtime_checkable
class WrapperAddressProtocol(Protocol):
    @abc.abstractmethod
    def __wrapper_address__(self) -> int: ...
    @abc.abstractmethod
    def signature(self) -> Signature: ...

class CompileResultWAP(WrapperAddressProtocol):
    def __init__(self, cres: CompileResult) -> None: ...
    def dump(self, tab: str = "") -> None: ...

    #
    @override
    def __wrapper_address__(self) -> int: ...
    @override
    def signature(self) -> Signature: ...

    #
    def __call__(self, *args: Any, **kwargs: Any) -> Any: ...
