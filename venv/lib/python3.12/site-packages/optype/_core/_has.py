import sys
import types
from collections.abc import Callable, Iterable, Mapping
from typing import Any, ClassVar, LiteralString, TypeAlias

if sys.version_info >= (3, 13):
    from typing import (
        ParamSpec,
        Protocol,
        TypeVar,
        TypeVarTuple,
        override,
        runtime_checkable,
    )
else:
    from typing_extensions import (
        ParamSpec,
        Protocol,
        TypeVar,
        TypeVarTuple,
        override,
        runtime_checkable,
    )

__all__ = [
    "HasAnnotations",
    "HasClass",
    "HasCode",
    "HasDict",
    "HasDoc",
    "HasFunc",
    "HasMatchArgs",
    "HasModule",
    "HasName",
    "HasNames",
    "HasQualname",
    "HasSelf",
    "HasSlots",
    "HasTypeParams",
    "HasWrapped",
]


def __dir__() -> list[str]:
    return __all__


###

_TypeParams: TypeAlias = tuple[TypeVar | ParamSpec | TypeVarTuple, ...]

_TypeT = TypeVar("_TypeT", bound=type)
_ObjectT_co = TypeVar("_ObjectT_co", default=object, covariant=True)
_FuncT_co = TypeVar("_FuncT_co", bound=Callable[..., object], covariant=True)
_TypeParamsT = TypeVar("_TypeParamsT", bound=_TypeParams, default=_TypeParams)

__AnyMapping: TypeAlias = "Mapping[str, object]"
__AnyDict: TypeAlias = dict[str, Any]
_DictT = TypeVar("_DictT", bound=__AnyMapping, default=__AnyDict)
_DictT_co = TypeVar("_DictT_co", bound=__AnyMapping, default=__AnyDict, covariant=True)

_NameT = TypeVar("_NameT", bound=str, default=str)
_QualNameT = TypeVar("_QualNameT", bound=str, default=_NameT)
_StrT_co = TypeVar("_StrT_co", bound=str, default=str, covariant=True)


###


@runtime_checkable
class HasMatchArgs(Protocol):
    __match_args__: ClassVar[tuple[LiteralString, ...] | list[LiteralString]]


@runtime_checkable
class HasSlots(Protocol):
    __slots__: ClassVar[LiteralString | Iterable[LiteralString]]


@runtime_checkable
class HasDict(Protocol[_DictT]):  # type: ignore[misc]
    # the typeshed annotations for `builtins.object.__dict__` too narrow
    __dict__: _DictT  # type: ignore[assignment]  # pyright: ignore[reportIncompatibleVariableOverride]


@runtime_checkable
class HasClass(Protocol[_TypeT]):
    """
    Can be seen as the **invariant** inverse of `type[T]`, i.e. `HasClass[type[T]]`
    represents (but is not equivalent to) `T`. However, `HasClass` is stricter, and
    rejects types that aren't fully static, such as `int | str`.

    It works best to annotate input types within `.pyi` stubs. If we, for example, have
    `def typeof(obj, /): return type(obj)` at runtime, we can stub it like this:

    ```pyi
    def typeof[TypeT: type](obj: HasClass[TypeT], /) -> TypeT: ...
    ```

    It behaves the same as `type` if we use fully static types:

    ```pyi
    just_int: int

    type(just_int)  # type[int]
    typeof(just_int)  # type[int]
    ```

    but when we have a dynamic type, the `HasClass` will reject it, but `type` will
    accept it:

    ```pyi
    int_or_str: int | str

    type(int_or_str)  # type[int | str]
    typeof(int_or_str)  # Error: misc (mypy), reportArgumentType (pyright)
    ```

    The inferred `type[int | str]` type by `type` does not represent a concrete type
    that can exist at runtime. The (stricter) `typeof` function therefore doesn't
    accept it, and both mypy and pyright will report it as an error. So `typeof` can
    be seen as more "realistic", or "less abstract", in that way, and it better reflects
    the possible runtime outcomes of `typeof`, than that `type` does of itself.
    """

    @property  # type: ignore[override]  # mypy bug
    @override
    def __class__(self) -> _TypeT: ...
    @__class__.setter
    @override
    def __class__(self, __class__: _TypeT, /) -> None: ...


@runtime_checkable
class HasModule(Protocol[_StrT_co]):
    __module__: _StrT_co


@runtime_checkable
class HasName(Protocol[_NameT]):
    __name__: _NameT


@runtime_checkable
class HasQualname(Protocol[_NameT]):  # pyright: ignore[reportInvalidTypeVarUse]
    __qualname__: _NameT


@runtime_checkable
class HasNames(  # pyright: ignore[reportIncompatibleVariableOverride, reportInvalidTypeVarUse]
    HasName[_NameT],
    HasQualname[_QualNameT],
    Protocol[_NameT, _QualNameT],
): ...


# docs and type hints


@runtime_checkable
class HasDoc(Protocol[_StrT_co]):
    # note that docstrings are stripped if ran with e.g. `python -OO`
    __doc__: _StrT_co | None


@runtime_checkable
class HasAnnotations(Protocol[_DictT_co]):  # pyright: ignore[reportInvalidTypeVarUse]
    __annotations__: _DictT_co  # type: ignore[assignment]  # pyright: ignore[reportIncompatibleVariableOverride]


@runtime_checkable
class HasTypeParams(Protocol[_TypeParamsT]):
    __type_params__: _TypeParamsT


# functions and methods


@runtime_checkable
class HasFunc(Protocol[_FuncT_co]):
    @property
    def __func__(self, /) -> _FuncT_co: ...


@runtime_checkable
class HasWrapped(Protocol[_FuncT_co]):
    @property
    def __wrapped__(self, /) -> _FuncT_co: ...


@runtime_checkable
class HasSelf(Protocol[_ObjectT_co]):
    @property
    def __self__(self, /) -> _ObjectT_co: ...


@runtime_checkable
class HasCode(Protocol):
    __code__: types.CodeType
