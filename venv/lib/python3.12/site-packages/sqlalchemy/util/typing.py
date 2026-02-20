# util/typing.py
# Copyright (C) 2022-2026 the SQLAlchemy authors and contributors
# <see AUTHORS file>
#
# This module is part of SQLAlchemy and is released under
# the MIT License: https://www.opensource.org/licenses/mit-license.php
# mypy: allow-untyped-defs, allow-untyped-calls

from __future__ import annotations

import builtins
from collections import deque
import collections.abc as collections_abc
import re
import sys
import typing
from typing import Any
from typing import Callable
from typing import Dict
from typing import ForwardRef
from typing import Generic
from typing import Iterable
from typing import Mapping
from typing import NewType
from typing import NoReturn
from typing import Optional
from typing import overload
from typing import Set
from typing import Tuple
from typing import Type
from typing import TYPE_CHECKING
from typing import TypeVar
from typing import Union

import typing_extensions

from . import compat

if True:  # zimports removes the tailing comments
    from typing_extensions import Annotated as Annotated  # 3.8
    from typing_extensions import Concatenate as Concatenate  # 3.10
    from typing_extensions import (
        dataclass_transform as dataclass_transform,  # 3.11,
    )
    from typing_extensions import Final as Final  # 3.8
    from typing_extensions import final as final  # 3.8
    from typing_extensions import get_args as get_args  # 3.10
    from typing_extensions import get_origin as get_origin  # 3.10
    from typing_extensions import Literal as Literal  # 3.8
    from typing_extensions import NotRequired as NotRequired  # 3.11
    from typing_extensions import ParamSpec as ParamSpec  # 3.10
    from typing_extensions import Protocol as Protocol  # 3.8
    from typing_extensions import SupportsIndex as SupportsIndex  # 3.8
    from typing_extensions import TypeAlias as TypeAlias  # 3.10
    from typing_extensions import TypedDict as TypedDict  # 3.8
    from typing_extensions import TypeGuard as TypeGuard  # 3.10
    from typing_extensions import Self as Self  # 3.11
    from typing_extensions import TypeAliasType as TypeAliasType  # 3.12
    from typing_extensions import Never as Never  # 3.11
    from typing_extensions import LiteralString as LiteralString  # 3.11

_T = TypeVar("_T", bound=Any)
_KT = TypeVar("_KT")
_KT_co = TypeVar("_KT_co", covariant=True)
_KT_contra = TypeVar("_KT_contra", contravariant=True)
_VT = TypeVar("_VT")
_VT_co = TypeVar("_VT_co", covariant=True)

if compat.py310:
    # why they took until py310 to put this in stdlib is beyond me,
    # I've been wanting it since py27
    from types import NoneType as NoneType
else:
    NoneType = type(None)  # type: ignore


def is_fwd_none(typ: Any) -> bool:
    return isinstance(typ, ForwardRef) and typ.__forward_arg__ == "None"


_AnnotationScanType = Union[
    Type[Any], str, ForwardRef, NewType, TypeAliasType, "GenericProtocol[Any]"
]


class ArgsTypeProtocol(Protocol):
    """protocol for types that have ``__args__``

    there's no public interface for this AFAIK

    """

    __args__: Tuple[_AnnotationScanType, ...]


class GenericProtocol(Protocol[_T]):
    """protocol for generic types.

    this since Python.typing _GenericAlias is private

    """

    __args__: Tuple[_AnnotationScanType, ...]
    __origin__: Type[_T]

    # Python's builtin _GenericAlias has this method, however builtins like
    # list, dict, etc. do not, even though they have ``__origin__`` and
    # ``__args__``
    #
    # def copy_with(self, params: Tuple[_AnnotationScanType, ...]) -> Type[_T]:
    #     ...


# copied from TypeShed, required in order to implement
# MutableMapping.update()
class SupportsKeysAndGetItem(Protocol[_KT, _VT_co]):
    def keys(self) -> Iterable[_KT]: ...

    def __getitem__(self, __k: _KT) -> _VT_co: ...


# work around https://github.com/microsoft/pyright/issues/3025
_LiteralStar = Literal["*"]


def de_stringify_annotation(
    cls: Type[Any],
    annotation: _AnnotationScanType,
    originating_module: str,
    locals_: Mapping[str, Any],
    *,
    str_cleanup_fn: Optional[Callable[[str, str], str]] = None,
    include_generic: bool = False,
    _already_seen: Optional[Set[Any]] = None,
) -> Type[Any]:
    """Resolve annotations that may be string based into real objects.

    This is particularly important if a module defines "from __future__ import
    annotations", as everything inside of __annotations__ is a string. We want
    to at least have generic containers like ``Mapped``, ``Union``, ``List``,
    etc.

    """
    # looked at typing.get_type_hints(), looked at pydantic.  We need much
    # less here, and we here try to not use any private typing internals
    # or construct ForwardRef objects which is documented as something
    # that should be avoided.

    original_annotation = annotation

    if is_fwd_ref(annotation):
        annotation = annotation.__forward_arg__

    if isinstance(annotation, str):
        if str_cleanup_fn:
            annotation = str_cleanup_fn(annotation, originating_module)

        annotation = eval_expression(
            annotation, originating_module, locals_=locals_, in_class=cls
        )

    if (
        include_generic
        and is_generic(annotation)
        and not is_literal(annotation)
    ):
        if _already_seen is None:
            _already_seen = set()

        if annotation in _already_seen:
            # only occurs recursively.  outermost return type
            # will always be Type.
            # the element here will be either ForwardRef or
            # Optional[ForwardRef]
            return original_annotation  # type: ignore
        else:
            _already_seen.add(annotation)

        elements = tuple(
            de_stringify_annotation(
                cls,
                elem,
                originating_module,
                locals_,
                str_cleanup_fn=str_cleanup_fn,
                include_generic=include_generic,
                _already_seen=_already_seen,
            )
            for elem in annotation.__args__
        )

        return _copy_generic_annotation_with(annotation, elements)

    return annotation  # type: ignore


def fixup_container_fwd_refs(
    type_: _AnnotationScanType,
) -> _AnnotationScanType:
    """Correct dict['x', 'y'] into dict[ForwardRef('x'), ForwardRef('y')]
    and similar for list, set

    """

    if (
        is_generic(type_)
        and get_origin(type_)
        in (
            dict,
            set,
            list,
            collections_abc.MutableSet,
            collections_abc.MutableMapping,
            collections_abc.MutableSequence,
            collections_abc.Mapping,
            collections_abc.Sequence,
        )
        # fight, kick and scream to struggle to tell the difference between
        # dict[] and typing.Dict[] which DO NOT compare the same and DO NOT
        # behave the same yet there is NO WAY to distinguish between which type
        # it is using public attributes
        and not re.match(
            "typing.(?:Dict|List|Set|.*Mapping|.*Sequence|.*Set)", repr(type_)
        )
    ):
        # compat with py3.10 and earlier
        return get_origin(type_).__class_getitem__(  # type: ignore
            tuple(
                [
                    ForwardRef(elem) if isinstance(elem, str) else elem
                    for elem in get_args(type_)
                ]
            )
        )
    return type_


def _copy_generic_annotation_with(
    annotation: GenericProtocol[_T], elements: Tuple[_AnnotationScanType, ...]
) -> Type[_T]:
    if hasattr(annotation, "copy_with"):
        # List, Dict, etc. real generics
        return annotation.copy_with(elements)  # type: ignore
    else:
        # Python builtins list, dict, etc.
        return annotation.__origin__[elements]  # type: ignore


def eval_expression(
    expression: str,
    module_name: str,
    *,
    locals_: Optional[Mapping[str, Any]] = None,
    in_class: Optional[Type[Any]] = None,
) -> Any:
    try:
        base_globals: Dict[str, Any] = sys.modules[module_name].__dict__
    except KeyError as ke:
        raise NameError(
            f"Module {module_name} isn't present in sys.modules; can't "
            f"evaluate expression {expression}"
        ) from ke

    try:
        if in_class is not None:
            cls_namespace = dict(in_class.__dict__)
            cls_namespace.setdefault(in_class.__name__, in_class)

            # see #10899.  We want the locals/globals to take precedence
            # over the class namespace in this context, even though this
            # is not the usual way variables would resolve.
            cls_namespace.update(base_globals)

            annotation = eval(expression, cls_namespace, locals_)
        else:
            annotation = eval(expression, base_globals, locals_)
    except Exception as err:
        raise NameError(
            f"Could not de-stringify annotation {expression!r}"
        ) from err
    else:
        return annotation


def eval_name_only(
    name: str,
    module_name: str,
    *,
    locals_: Optional[Mapping[str, Any]] = None,
) -> Any:
    if "." in name:
        return eval_expression(name, module_name, locals_=locals_)

    try:
        base_globals: Dict[str, Any] = sys.modules[module_name].__dict__
    except KeyError as ke:
        raise NameError(
            f"Module {module_name} isn't present in sys.modules; can't "
            f"resolve name {name}"
        ) from ke

    # name only, just look in globals.  eval() works perfectly fine here,
    # however we are seeking to have this be faster, as this occurs for
    # every Mapper[] keyword, etc. depending on configuration
    try:
        return base_globals[name]
    except KeyError as ke:
        # check in builtins as well to handle `list`, `set` or `dict`, etc.
        try:
            return builtins.__dict__[name]
        except KeyError:
            pass

        raise NameError(
            f"Could not locate name {name} in module {module_name}"
        ) from ke


def resolve_name_to_real_class_name(name: str, module_name: str) -> str:
    try:
        obj = eval_name_only(name, module_name)
    except NameError:
        return name
    else:
        return getattr(obj, "__name__", name)


def is_pep593(type_: Optional[Any]) -> bool:
    return type_ is not None and get_origin(type_) in _type_tuples.Annotated


def is_non_string_iterable(obj: Any) -> TypeGuard[Iterable[Any]]:
    return isinstance(obj, collections_abc.Iterable) and not isinstance(
        obj, (str, bytes)
    )


def is_literal(type_: Any) -> bool:
    return get_origin(type_) in _type_tuples.Literal


def is_newtype(type_: Optional[_AnnotationScanType]) -> TypeGuard[NewType]:
    return hasattr(type_, "__supertype__")

    # doesn't work in 3.8, 3.7 as it passes a closure, not an
    # object instance
    # isinstance(type, type_instances.NewType)


def is_generic(type_: _AnnotationScanType) -> TypeGuard[GenericProtocol[Any]]:
    return hasattr(type_, "__args__") and hasattr(type_, "__origin__")


def is_pep695(type_: _AnnotationScanType) -> TypeGuard[TypeAliasType]:
    # NOTE: a generic TAT does not instance check as TypeAliasType outside of
    # python 3.10. For sqlalchemy use cases it's fine to consider it a TAT
    # though.
    # NOTE: things seems to work also without this additional check
    if is_generic(type_):
        return is_pep695(type_.__origin__)
    return isinstance(type_, _type_instances.TypeAliasType)


def flatten_newtype(type_: NewType) -> Type[Any]:
    super_type = type_.__supertype__
    while is_newtype(super_type):
        super_type = super_type.__supertype__
    return super_type  # type: ignore[return-value]


def pep695_values(type_: _AnnotationScanType) -> Set[Any]:
    """Extracts the value from a TypeAliasType, recursively exploring unions
    and inner TypeAliasType to flatten them into a single set.

    Forward references are not evaluated, so no recursive exploration happens
    into them.
    """
    _seen = set()

    def recursive_value(inner_type):
        if inner_type in _seen:
            # recursion are not supported (at least it's flagged as
            # an error by pyright). Just avoid infinite loop
            return inner_type
        _seen.add(inner_type)
        if not is_pep695(inner_type):
            return inner_type
        value = inner_type.__value__
        if not is_union(value):
            return value
        return [recursive_value(t) for t in value.__args__]

    res = recursive_value(type_)
    if isinstance(res, list):
        types = set()
        stack = deque(res)
        while stack:
            t = stack.popleft()
            if isinstance(t, list):
                stack.extend(t)
            else:
                types.add(None if t is NoneType or is_fwd_none(t) else t)
        return types
    else:
        return {res}


def is_fwd_ref(
    type_: _AnnotationScanType,
    check_generic: bool = False,
    check_for_plain_string: bool = False,
) -> TypeGuard[ForwardRef]:
    if check_for_plain_string and isinstance(type_, str):
        return True
    elif isinstance(type_, _type_instances.ForwardRef):
        return True
    elif check_generic and is_generic(type_):
        return any(
            is_fwd_ref(
                arg, True, check_for_plain_string=check_for_plain_string
            )
            for arg in type_.__args__
        )
    else:
        return False


@overload
def de_optionalize_union_types(type_: str) -> str: ...


@overload
def de_optionalize_union_types(type_: Type[Any]) -> Type[Any]: ...


@overload
def de_optionalize_union_types(
    type_: _AnnotationScanType,
) -> _AnnotationScanType: ...


def de_optionalize_union_types(
    type_: _AnnotationScanType,
) -> _AnnotationScanType:
    """Given a type, filter out ``Union`` types that include ``NoneType``
    to not include the ``NoneType``.

    Contains extra logic to work on non-flattened unions, unions that contain
    ``None`` (seen in py38, 37)

    """

    if is_fwd_ref(type_):
        return _de_optionalize_fwd_ref_union_types(type_, False)

    elif is_union(type_) and includes_none(type_):
        if compat.py39:
            typ = set(type_.__args__)
        else:
            # py38, 37 - unions are not automatically flattened, can contain
            # None rather than NoneType
            stack_of_unions = deque([type_])
            typ = set()
            while stack_of_unions:
                u_typ = stack_of_unions.popleft()
                for elem in u_typ.__args__:
                    if is_union(elem):
                        stack_of_unions.append(elem)
                    else:
                        typ.add(elem)

            typ.discard(None)  # type: ignore

        typ = {t for t in typ if t is not NoneType and not is_fwd_none(t)}

        return make_union_type(*typ)

    else:
        return type_


@overload
def _de_optionalize_fwd_ref_union_types(
    type_: ForwardRef, return_has_none: Literal[True]
) -> bool: ...


@overload
def _de_optionalize_fwd_ref_union_types(
    type_: ForwardRef, return_has_none: Literal[False]
) -> _AnnotationScanType: ...


def _de_optionalize_fwd_ref_union_types(
    type_: ForwardRef, return_has_none: bool
) -> Union[_AnnotationScanType, bool]:
    """return the non-optional type for Optional[], Union[None, ...], x|None,
    etc. without de-stringifying forward refs.

    unfortunately this seems to require lots of hardcoded heuristics

    """

    annotation = type_.__forward_arg__

    mm = re.match(r"^(.+?)\[(.+)\]$", annotation)
    if mm:
        g1 = mm.group(1).split(".")[-1]
        if g1 == "Optional":
            return True if return_has_none else ForwardRef(mm.group(2))
        elif g1 == "Union":
            if "[" in mm.group(2):
                # cases like "Union[Dict[str, int], int, None]"
                elements: list[str] = []
                current: list[str] = []
                ignore_comma = 0
                for char in mm.group(2):
                    if char == "[":
                        ignore_comma += 1
                    elif char == "]":
                        ignore_comma -= 1
                    elif ignore_comma == 0 and char == ",":
                        elements.append("".join(current).strip())
                        current.clear()
                        continue
                    current.append(char)
            else:
                elements = re.split(r",\s*", mm.group(2))
            parts = [ForwardRef(elem) for elem in elements if elem != "None"]
            if return_has_none:
                return len(elements) != len(parts)
            else:
                return make_union_type(*parts) if parts else Never  # type: ignore[return-value] # noqa: E501
        else:
            return False if return_has_none else type_

    pipe_tokens = re.split(r"\s*\|\s*", annotation)
    has_none = "None" in pipe_tokens
    if return_has_none:
        return has_none
    if has_none:
        anno_str = "|".join(p for p in pipe_tokens if p != "None")
        return ForwardRef(anno_str) if anno_str else Never  # type: ignore[return-value] # noqa: E501

    return type_


def make_union_type(*types: _AnnotationScanType) -> Type[Any]:
    """Make a Union type."""

    return Union[types]  # type: ignore


def includes_none(type_: Any) -> bool:
    """Returns if the type annotation ``type_`` allows ``None``.

    This function supports:
    * forward refs
    * unions
    * pep593 - Annotated
    * pep695 - TypeAliasType (does not support looking into
    fw reference of other pep695)
    * NewType
    * plain types like ``int``, ``None``, etc
    """
    if is_fwd_ref(type_):
        return _de_optionalize_fwd_ref_union_types(type_, True)
    if is_union(type_):
        return any(includes_none(t) for t in get_args(type_))
    if is_pep593(type_):
        return includes_none(get_args(type_)[0])
    if is_pep695(type_):
        return any(includes_none(t) for t in pep695_values(type_))
    if is_newtype(type_):
        return includes_none(type_.__supertype__)
    try:
        return type_ in (NoneType, None) or is_fwd_none(type_)
    except TypeError:
        # if type_ is Column, mapped_column(), etc. the use of "in"
        # resolves to ``__eq__()`` which then gives us an expression object
        # that can't resolve to boolean.  just catch it all via exception
        return False


def is_a_type(type_: Any) -> bool:
    return (
        isinstance(type_, type)
        or get_origin(type_) is not None
        or getattr(type_, "__module__", None)
        in ("typing", "typing_extensions")
        or type(type_).__mro__[0].__module__ in ("typing", "typing_extensions")
    )


def is_union(type_: Any) -> TypeGuard[ArgsTypeProtocol]:
    return is_origin_of(type_, "Union", "UnionType")


def is_origin_of_cls(
    type_: Any, class_obj: Union[Tuple[Type[Any], ...], Type[Any]]
) -> bool:
    """return True if the given type has an __origin__ that shares a base
    with the given class"""

    origin = get_origin(type_)
    if origin is None:
        return False

    return isinstance(origin, type) and issubclass(origin, class_obj)


def is_origin_of(
    type_: Any, *names: str, module: Optional[str] = None
) -> bool:
    """return True if the given type has an __origin__ with the given name
    and optional module."""

    origin = get_origin(type_)
    if origin is None:
        return False

    return _get_type_name(origin) in names and (
        module is None or origin.__module__.startswith(module)
    )


def _get_type_name(type_: Type[Any]) -> str:
    if compat.py310:
        return type_.__name__
    else:
        typ_name = getattr(type_, "__name__", None)
        if typ_name is None:
            typ_name = getattr(type_, "_name", None)

        return typ_name  # type: ignore


class DescriptorProto(Protocol):
    def __get__(self, instance: object, owner: Any) -> Any: ...

    def __set__(self, instance: Any, value: Any) -> None: ...

    def __delete__(self, instance: Any) -> None: ...


_DESC = TypeVar("_DESC", bound=DescriptorProto)


class DescriptorReference(Generic[_DESC]):
    """a descriptor that refers to a descriptor.

    used for cases where we need to have an instance variable referring to an
    object that is itself a descriptor, which typically confuses typing tools
    as they don't know when they should use ``__get__`` or not when referring
    to the descriptor assignment as an instance variable. See
    sqlalchemy.orm.interfaces.PropComparator.prop

    """

    if TYPE_CHECKING:

        def __get__(self, instance: object, owner: Any) -> _DESC: ...

        def __set__(self, instance: Any, value: _DESC) -> None: ...

        def __delete__(self, instance: Any) -> None: ...


_DESC_co = TypeVar("_DESC_co", bound=DescriptorProto, covariant=True)


class RODescriptorReference(Generic[_DESC_co]):
    """a descriptor that refers to a descriptor.

    same as :class:`.DescriptorReference` but is read-only, so that subclasses
    can define a subtype as the generically contained element

    """

    if TYPE_CHECKING:

        def __get__(self, instance: object, owner: Any) -> _DESC_co: ...

        def __set__(self, instance: Any, value: Any) -> NoReturn: ...

        def __delete__(self, instance: Any) -> NoReturn: ...


_FN = TypeVar("_FN", bound=Optional[Callable[..., Any]])


class CallableReference(Generic[_FN]):
    """a descriptor that refers to a callable.

    works around mypy's limitation of not allowing callables assigned
    as instance variables


    """

    if TYPE_CHECKING:

        def __get__(self, instance: object, owner: Any) -> _FN: ...

        def __set__(self, instance: Any, value: _FN) -> None: ...

        def __delete__(self, instance: Any) -> None: ...


class _TypingInstances:
    def __getattr__(self, key: str) -> tuple[type, ...]:
        types = tuple(
            {
                t
                for t in [
                    getattr(typing, key, None),
                    getattr(typing_extensions, key, None),
                ]
                if t is not None
            }
        )
        if not types:
            raise AttributeError(key)
        self.__dict__[key] = types
        return types


_type_tuples = _TypingInstances()
if TYPE_CHECKING:
    _type_instances = typing_extensions
else:
    _type_instances = _type_tuples

LITERAL_TYPES = _type_tuples.Literal
