import sys
import types
from collections.abc import Iterable
from typing import Literal, TypeAlias, cast, get_args as _get_args, overload

if sys.version_info >= (3, 13):
    from typing import TypeAliasType, TypeIs, is_protocol
else:
    from typing_extensions import TypeAliasType, TypeIs, is_protocol


from ._core import _can as _c
from .types import (
    AnnotatedAlias,
    GenericType,
    LiteralAlias,
    UnionAlias,
    WrappedFinalType,
)

__all__ = (
    "get_args",
    "get_protocol_members",
    "get_protocols",
    "is_final",
    "is_generic_alias",
    "is_iterable",
    "is_protocol",
    "is_runtime_protocol",
    "is_union_type",
)


def __dir__() -> tuple[str, ...]:
    return __all__


###


_ToIter: TypeAlias = Iterable[object] | _c.CanGetitem[int, object]


def is_iterable(obj: object, /) -> TypeIs[_ToIter]:
    """
    Check whether the object can be iterated over, i.e. if it can be used in
    a `for` loop, or if it can be passed to `builtins.iter`.

    Note:
        Contrary to popular *belief*, this isn't limited to objects that
        implement `__iter___`, as suggested by the name of
        `Iterable`.

        Sequence-like objects that implement `__getitem__` for consecutive
        `int` keys that start at `0` (or raise `IndexEeror` if out of bounds),
        can also be used in a `for` loop.
        In fact, all builtins that accept "iterables" also accept these
        sequence-likes at runtime.

    See also:
        - [`optype.types.Iterable`][optype.types.Iterable]
    """
    if isinstance(obj, type):
        # workaround the false-positive bug in `@runtime_checkable` for types
        return False

    if isinstance(obj, _c.CanIter):
        return True
    if isinstance(obj, _c.CanGetitem):
        # check if obj is a sequence-like
        if isinstance(obj, _c.CanLen):
            return True

        # not all sequence-likes implement __len__, e.g. `ctypes.pointer`
        try:
            obj[0]  # pyright: ignore[reportArgumentType]
        except (IndexError, StopIteration):
            pass
        except (KeyError, ValueError, TypeError):
            return False

        return True

    return False


@overload
def is_final(final_cls_or_method: WrappedFinalType, /) -> Literal[True]: ...
@overload
def is_final(cls: object, /) -> bool: ...
def is_final(arg: object, /) -> bool:
    """
    Check if the type, method, classmethod, staticmethod, or property, is
    decorated with `@typing.final` or `@typing_extensions.final`.

    IMPORTANT: A final `@property` won't be recognized unless `@final` is
    applied before the `@property` decorator, i.e. directly on the method.

    Do this:

    ```python
    class Europe:
        @property
        @final
        def countdown(self): ...
    ```

    Don't do this:

    ```python
    class Europe:
        @final
        @property
        def countdown(self): ...
    ```

    The reason for this is that `builtins.property.__final__` cannot be
    written to, so `@final` won't be able to set `__final__ = True` on it.
    Only the getter method of a `@property` is checked.

    NOTE: Accessing a `classmethod` or `staticmethod` of a class can't be done
    through regular attribute access, but should be done with
    `inspect.getattr_static`.
    But unlike `@property`, it doesn't matter in which order you use the
    `@classmethod` and `@staticmethod` decorators and `@final`.
    """
    if getattr(arg, "__final__", None) is True:
        return True

    if isinstance(arg, type):
        try:
            _ = type("_", (arg,), {})
        except TypeError:
            return True

    if (wrapped := getattr(arg, "__wrapped__", arg)) is not arg:
        return is_final(wrapped)

    return isinstance(arg, property) and arg.fget is not None and is_final(arg.fget)


def _get_alias(tp: type | object, /) -> type | object:
    seen: set[TypeAliasType] = set()
    for _ in range(sys.getrecursionlimit()):
        if isinstance(tp, TypeAliasType):
            seen.add(tp)
            tp = cast("type | object", tp.__value__)
            if tp in seen:
                raise RecursionError("type alias of itself")
            continue
        if isinstance(tp, AnnotatedAlias):
            assert len(tp.__args__) == 1
            tp = tp.__args__[0]
            continue
        return tp

    raise RecursionError


def is_runtime_protocol(type_expr: type | object, /) -> bool:
    """
    Check if `type_expr` is a `typing[_extensions].Protocol` that's decorated
    with `typing[_extensions].runtime_checkable`.
    """
    if isinstance(type_expr, AnnotatedAlias | TypeAliasType):
        type_expr = _get_alias(type_expr)
    return (
        bool(getattr(type_expr, "_is_runtime_protocol", False))
        and isinstance(type_expr, type)
        and is_protocol(type_expr)
    )


def is_union_type(type_expr: type | object, /) -> TypeIs[types.UnionType | UnionAlias]:
    if isinstance(type_expr, AnnotatedAlias | TypeAliasType):
        type_expr = _get_alias(type_expr)
    return isinstance(type_expr, types.UnionType | UnionAlias)


def is_generic_alias(
    type_expr: type | object,
    /,
) -> TypeIs[GenericType | types.GenericAlias]:
    if isinstance(type_expr, AnnotatedAlias | TypeAliasType):
        type_expr = _get_alias(type_expr)
    return (
        isinstance(type_expr, GenericType | types.GenericAlias)
        and not isinstance(type_expr, types.UnionType | UnionAlias)
    )  # fmt: skip


def get_args(tp: type | object, /) -> tuple[type | object, ...]:
    """
    A less broken implementation of `typing[_extensions].get_args()` that

    - unpacks `Annotated` and `TypeAliasType`,
    - recursively flattens unions / nested `Literal`s, and
    - raises `TypeError` if `tp` if isn't a generic type (alias).

    Raises:
        NotImplementedError: If `tp` is a string.
        TypeError: If `tp` is not a generic type or type-alias.
    """
    if isinstance(tp, str):
        raise NotImplementedError("str")

    should_raise = True
    if isinstance(tp, AnnotatedAlias | TypeAliasType):
        tp = _get_alias(tp)
        should_raise = False

    if isinstance(tp, types.UnionType | UnionAlias | LiteralAlias):
        arg: object
        args: list[object] = []
        for arg in tp.__args__:
            if isinstance(arg, TypeAliasType | AnnotatedAlias):
                arg = _get_alias(arg)  # noqa: PLW2901
            if isinstance(arg, types.UnionType | UnionAlias | LiteralAlias):
                args.extend(get_args(arg))
            else:
                args.append(arg)
        return tuple(args)

    if hasattr(tp, "__origin__") and hasattr(tp, "__args__"):
        return _get_args(tp)

    if should_raise and not isinstance(tp, type):
        raise TypeError(repr(tp))

    return ()


def get_protocol_members(cls: type, /) -> frozenset[str]:
    """
    A variant of `typing[_extensions].get_protocol_members()` that

    - doesn't hide `__annotations__`, `__dict__`, `__init__` or `__new__`,
    - doesn't add a `__hash__` if there's an `__eq__` method, and
    - doesn't include methods of base types from different modules.

    Raises:
        TypeError: If `cls` is not a protocol type.
    """
    if not is_protocol(cls):
        msg = f"{cls!r} is not a protocol"
        raise TypeError(msg)

    annotations, module = cls.__annotations__, cls.__module__

    member_dict: types.MappingProxyType[str, object] = vars(cls)
    members = annotations.keys() | {
        name
        for name, v in member_dict.items()
        if (
            callable(f := getattr(v, "__func__", v))
            or (isinstance(v, property) and (f := v.fget) is not None)
        )
        and f.__module__ == module
    }

    # this hack here is plagiarized from the (often incorrect)
    # `typing_extensions.get_protocol_members`.
    # Maybe the `typing.get_protocol_member`s` that's coming in 3.13 will
    # won't be as broken. I have little hope though...
    if (protocol_attrs := getattr(cls, "__protocol_attrs__", None)) is not None:
        members |= protocol_attrs
    else:
        # python <3.12
        from optype._inspect import _get_protocol_attrs  # noqa: PLC0415

        members |= _get_protocol_attrs(cls)

    # sometimes __protocol_attrs__ hallicunates some non-existing dunders.
    # the `getattr_static` avoids potential descriptor magic
    from inspect import getattr_static  # noqa: PLC0415

    members = {
        member
        for member in members
        if member in annotations
        or getattr_static(cls, member) is not None
    }  # fmt: skip

    # also include any of the parents
    for supercls in cls.mro()[1:]:
        if is_protocol(supercls):
            members |= get_protocol_members(supercls)

    return frozenset(members)


def get_protocols(
    module: types.ModuleType,
    /,
    private: bool = False,
) -> frozenset[type]:
    """Return the protocol types within the given module."""
    members: list[str] | tuple[str, ...]
    if private:
        members = dir(module)
    elif hasattr(module, "__all__"):
        members = cast("list[str] | tuple[str, ...]", module.__all__)
    else:
        members = [k for k in dir(module) if not k.startswith("_")]

    return frozenset({
        cls for name in members
        if is_protocol(cls := cast("type", getattr(module, name)))
    })  # fmt: skip
