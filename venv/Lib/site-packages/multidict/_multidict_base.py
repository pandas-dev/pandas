import sys
from collections.abc import (
    Container,
    ItemsView,
    Iterable,
    KeysView,
    Mapping,
    Set,
    ValuesView,
)
from typing import Literal, Union

if sys.version_info >= (3, 10):
    from types import NotImplementedType
else:
    from typing import Any as NotImplementedType

if sys.version_info >= (3, 11):
    from typing import assert_never
else:
    from typing_extensions import assert_never


def _abc_itemsview_register(view_cls: type[object]) -> None:
    ItemsView.register(view_cls)


def _abc_keysview_register(view_cls: type[object]) -> None:
    KeysView.register(view_cls)


def _abc_valuesview_register(view_cls: type[object]) -> None:
    ValuesView.register(view_cls)


def _viewbaseset_richcmp(
    view: set[object], other: object, op: Literal[0, 1, 2, 3, 4, 5]
) -> Union[bool, NotImplementedType]:
    if op == 0:  # <
        if not isinstance(other, Set):
            return NotImplemented  # type: ignore[no-any-return]
        return len(view) < len(other) and view <= other
    elif op == 1:  # <=
        if not isinstance(other, Set):
            return NotImplemented  # type: ignore[no-any-return]
        if len(view) > len(other):
            return False
        for elem in view:
            if elem not in other:
                return False
        return True
    elif op == 2:  # ==
        if not isinstance(other, Set):
            return NotImplemented  # type: ignore[no-any-return]
        return len(view) == len(other) and view <= other
    elif op == 3:  # !=
        return not view == other
    elif op == 4:  # >
        if not isinstance(other, Set):
            return NotImplemented  # type: ignore[no-any-return]
        return len(view) > len(other) and view >= other
    elif op == 5:  # >=
        if not isinstance(other, Set):
            return NotImplemented  # type: ignore[no-any-return]
        if len(view) < len(other):
            return False
        for elem in other:
            if elem not in view:
                return False
        return True
    else:  # pragma: no cover
        assert_never(op)


def _viewbaseset_and(
    view: set[object], other: object
) -> Union[set[object], NotImplementedType]:
    if not isinstance(other, Iterable):
        return NotImplemented  # type: ignore[no-any-return]
    if isinstance(view, Set):
        view = set(iter(view))
    if isinstance(other, Set):
        other = set(iter(other))
    if not isinstance(other, Set):
        other = set(iter(other))
    return view & other


def _viewbaseset_or(
    view: set[object], other: object
) -> Union[set[object], NotImplementedType]:
    if not isinstance(other, Iterable):
        return NotImplemented  # type: ignore[no-any-return]
    if isinstance(view, Set):
        view = set(iter(view))
    if isinstance(other, Set):
        other = set(iter(other))
    if not isinstance(other, Set):
        other = set(iter(other))
    return view | other


def _viewbaseset_sub(
    view: set[object], other: object
) -> Union[set[object], NotImplementedType]:
    if not isinstance(other, Iterable):
        return NotImplemented  # type: ignore[no-any-return]
    if isinstance(view, Set):
        view = set(iter(view))
    if isinstance(other, Set):
        other = set(iter(other))
    if not isinstance(other, Set):
        other = set(iter(other))
    return view - other


def _viewbaseset_xor(
    view: set[object], other: object
) -> Union[set[object], NotImplementedType]:
    if not isinstance(other, Iterable):
        return NotImplemented  # type: ignore[no-any-return]
    if isinstance(view, Set):
        view = set(iter(view))
    if isinstance(other, Set):
        other = set(iter(other))
    if not isinstance(other, Set):
        other = set(iter(other))
    return view ^ other


def _itemsview_isdisjoint(view: Container[object], other: Iterable[object]) -> bool:
    "Return True if two sets have a null intersection."
    for v in other:
        if v in view:
            return False
    return True


def _itemsview_repr(view: Iterable[tuple[object, object]]) -> str:
    lst = []
    for k, v in view:
        lst.append("{!r}: {!r}".format(k, v))
    body = ", ".join(lst)
    return "{}({})".format(view.__class__.__name__, body)


def _keysview_isdisjoint(view: Container[object], other: Iterable[object]) -> bool:
    "Return True if two sets have a null intersection."
    for k in other:
        if k in view:
            return False
    return True


def _keysview_repr(view: Iterable[object]) -> str:
    lst = []
    for k in view:
        lst.append("{!r}".format(k))
    body = ", ".join(lst)
    return "{}({})".format(view.__class__.__name__, body)


def _valuesview_repr(view: Iterable[object]) -> str:
    lst = []
    for v in view:
        lst.append("{!r}".format(v))
    body = ", ".join(lst)
    return "{}({})".format(view.__class__.__name__, body)


def _mdrepr(md: Mapping[object, object]) -> str:
    lst = []
    for k, v in md.items():
        lst.append("'{}': {!r}".format(k, v))
    body = ", ".join(lst)
    return "<{}({})>".format(md.__class__.__name__, body)
