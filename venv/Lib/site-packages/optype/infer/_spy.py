# ruff: noqa: TD002, TD003, PYI034

"""Recording proxy objects that trace the operations performed on them."""

from collections.abc import Generator, Iterator
from contextvars import ContextVar
from typing import Any, NamedTuple, final, override


class _Fork(BaseException): ...


_fork: ContextVar[Iterator[bool] | None] = ContextVar("_fork", default=None)


def _decide() -> bool:
    if (plan := _fork.get()) is None:
        return True
    if (value := next(plan, None)) is None:
        raise _Fork
    return value


class _TraceItem(NamedTuple):
    attr: str
    args: tuple[Any, ...]
    kwargs: dict[str, Any]
    return_: Any


class _Spy:
    __optype_trace__: list[_TraceItem]

    def __optype_trace_add__[OutT](  # noqa: PLW3201
        self,
        attr: str,
        args: tuple[Any, ...],
        kwargs: dict[str, Any],
        out: OutT,
    ) -> OutT:
        self.__optype_trace__.append(_TraceItem(attr, args, kwargs, out))
        return out

    def __init__(self, /, *_args: object, **_kwargs: object) -> None:
        self.__optype_trace__ = []


@final
class _SpyStr(str, _Spy):
    __slots__ = ()  # pyrefly:ignore[implicit-any-attribute]


@final
class _SpyBytes(bytes, _Spy):
    __slots__ = ()  # pyrefly:ignore[implicit-any-attribute]


class _SpyObject(_Spy):  # noqa: PLR0904
    __optype_element__: "_SpyObject | None" = None
    __optype_iterator__: bool = False

    ###

    def __getattr__(self, attr: str, /) -> "_SpyObject | None":
        # TODO: specialize for known special dunder attrs, e.g. `__name__: str`
        if attr.startswith("__optype"):
            return None

        return self.__optype_trace_add__("__getattr__", (attr,), {}, _SpyObject())

    @override
    def __setattr__(self, attr: str, value: object, /) -> None:
        if attr.startswith("__optype"):
            return super().__setattr__(attr, value)

        return self.__optype_trace_add__("__setattr__", (attr, value), {}, None)

    @override
    def __delattr__(self, attr: str, /) -> None:
        self.__optype_trace_add__("__delattr__", (attr,), {}, None)

    @override
    def __dir__(self, /) -> "_SpyObject":
        # TODO: maybe return specialized `_SpyObject & Iterable[str]`
        return self.__optype_trace_add__("__dir__", (), {}, _SpyObject())

    # TODO: __class__ -> _SpyType

    ###

    @override
    def __repr__(self, /) -> _SpyStr:
        return self.__optype_trace_add__("__repr__", (), {}, _SpyStr())

    @override
    def __str__(self, /) -> _SpyStr:
        return self.__optype_trace_add__("__str__", (), {}, _SpyStr())

    @override
    def __format__(self, format_spec: str, /) -> _SpyStr:
        return self.__optype_trace_add__("__format__", (format_spec,), {}, _SpyStr())

    def __bytes__(self, /) -> _SpyBytes:
        return self.__optype_trace_add__("__bytes__", (), {}, _SpyBytes())

    ###

    @override
    def __eq__(self, other: object, /) -> "_SpyObject":  # type:ignore[override] # pyright:ignore[reportIncompatibleMethodOverride] # ty:ignore[invalid-method-override]
        return self.__optype_trace_add__("__eq__", (other,), {}, _SpyObject())

    @override
    def __ne__(self, other: object, /) -> "_SpyObject":  # type:ignore[override] # pyright:ignore[reportIncompatibleMethodOverride] # ty:ignore[invalid-method-override]
        return self.__optype_trace_add__("__ne__", (other,), {}, _SpyObject())

    def __lt__(self, other: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__lt__", (other,), {}, _SpyObject())

    def __le__(self, other: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__le__", (other,), {}, _SpyObject())

    def __gt__(self, other: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__gt__", (other,), {}, _SpyObject())

    def __ge__(self, other: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__ge__", (other,), {}, _SpyObject())

    ###

    @override
    def __hash__(self, /) -> int:
        return self.__optype_trace_add__("__hash__", (), {}, super().__hash__())

    def __bool__(self, /) -> bool:
        return self.__optype_trace_add__("__bool__", (), {}, _decide())

    ###

    def __get__(self, instance: object, owner: type | None = None) -> "_SpyObject":
        return self.__optype_trace_add__("__get__", (instance, owner), {}, _SpyObject())

    def __set__(self, instance: object, value: object, /) -> None:
        self.__optype_trace_add__("__set__", (instance, value), {}, None)

    def __delete__(self, instance: object, /) -> None:
        self.__optype_trace_add__("__delete__", (instance,), {}, None)

    # TODO: __objclass__ -> _SpyType

    def __set_name__(self, owner: type, name: str, /) -> None:
        self.__optype_trace_add__("__set_name__", (owner, name), {}, None)

    ###

    # TODO: __mro_entries__
    # TODO: __instancecheck__
    # TODO: __subclasscheck__

    ###

    def __call__(self, /, *args: object, **kwargs: object) -> "_SpyObject":
        return self.__optype_trace_add__("__call__", args, kwargs, _SpyObject())

    ###

    def __len__(self, /) -> int:
        # TODO: over-recorded as required when probed optionally (e.g. `list()` via
        # `length_hint`); detect optional protocols by forking value-vs-raise
        return self.__optype_trace_add__("__len__", (), {}, _decide())

    # no need for `__length_hint__`

    def __getitem__(self, key: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__getitem__", (key,), {}, _SpyObject())

    def __setitem__(self, key: object, value: object, /) -> None:
        return self.__optype_trace_add__("__setitem__", (key, value), {}, None)

    def __delitem__(self, key: object, /) -> None:
        return self.__optype_trace_add__("__delitem__", (key,), {}, None)

    # no need for `__missing__`

    def __iter__(self, /) -> "_SpyObject":
        if self.__optype_iterator__:
            return self  # an iterator is its own iterable (idempotent `iter()`)
        return self.__optype_trace_add__("__iter__", (), {}, _iterator_of(self))

    def __reversed__(self, /) -> "_SpyObject":
        return self.__optype_trace_add__("__reversed__", (), {}, _iterator_of(self))

    def __contains__(self, item: object, /) -> bool:
        return self.__optype_trace_add__("__contains__", (item,), {}, _decide())

    # return `Any` instead of `_SpyObject` to avoid an LSP error for `__dir__`
    def __next__(self, /) -> Any:
        if any(item.attr == "__next__" for item in self.__optype_trace__):
            raise StopIteration
        return self.__optype_trace_add__("__next__", (), {}, _element_of(self))

    ###

    def __add__(self, rhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__add__", (rhs,), {}, _SpyObject())

    def __sub__(self, rhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__sub__", (rhs,), {}, _SpyObject())

    def __mul__(self, rhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__mul__", (rhs,), {}, _SpyObject())

    def __matmul__(self, rhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__matmul__", (rhs,), {}, _SpyObject())

    def __truediv__(self, rhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__truediv__", (rhs,), {}, _SpyObject())

    def __floordiv__(self, rhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__floordiv__", (rhs,), {}, _SpyObject())

    def __mod__(self, rhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__mod__", (rhs,), {}, _SpyObject())

    def __divmod__(self, rhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__divmod__", (rhs,), {}, _SpyObject())

    def __pow__(self, rhs: object, /, *args: object) -> "_SpyObject":
        return self.__optype_trace_add__("__pow__", (rhs, *args), {}, _SpyObject())

    def __lshift__(self, rhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__lshift__", (rhs,), {}, _SpyObject())

    def __rshift__(self, rhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__rshift__", (rhs,), {}, _SpyObject())

    def __and__(self, rhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__and__", (rhs,), {}, _SpyObject())

    def __xor__(self, rhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__xor__", (rhs,), {}, _SpyObject())

    def __or__(self, rhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__or__", (rhs,), {}, _SpyObject())

    ###

    def __radd__(self, lhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__radd__", (lhs,), {}, _SpyObject())

    def __rsub__(self, lhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__rsub__", (lhs,), {}, _SpyObject())

    def __rmul__(self, lhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__rmul__", (lhs,), {}, _SpyObject())

    def __rmatmul__(self, lhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__rmatmul__", (lhs,), {}, _SpyObject())

    def __rtruediv__(self, lhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__rtruediv__", (lhs,), {}, _SpyObject())

    def __rfloordiv__(self, lhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__rfloordiv__", (lhs,), {}, _SpyObject())

    def __rmod__(self, lhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__rmod__", (lhs,), {}, _SpyObject())

    def __rdivmod__(self, lhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__rdivmod__", (lhs,), {}, _SpyObject())

    def __rpow__(self, lhs: object, /, *args: object) -> "_SpyObject":
        return self.__optype_trace_add__("__rpow__", (lhs, *args), {}, _SpyObject())

    def __rlshift__(self, lhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__rlshift__", (lhs,), {}, _SpyObject())

    def __rrshift__(self, lhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__rrshift__", (lhs,), {}, _SpyObject())

    def __rand__(self, lhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__rand__", (lhs,), {}, _SpyObject())

    def __rxor__(self, lhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__rxor__", (lhs,), {}, _SpyObject())

    def __ror__(self, lhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__ror__", (lhs,), {}, _SpyObject())

    ###

    def __iadd__(self, rhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__iadd__", (rhs,), {}, _SpyObject())

    def __isub__(self, rhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__isub__", (rhs,), {}, _SpyObject())

    def __imul__(self, rhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__imul__", (rhs,), {}, _SpyObject())

    def __imatmul__(self, rhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__imatmul__", (rhs,), {}, _SpyObject())

    def __itruediv__(self, rhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__itruediv__", (rhs,), {}, _SpyObject())

    def __ifloordiv__(self, rhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__ifloordiv__", (rhs,), {}, _SpyObject())

    def __imod__(self, rhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__imod__", (rhs,), {}, _SpyObject())

    def __ipow__(self, rhs: object, /, *args: object) -> "_SpyObject":
        return self.__optype_trace_add__("__ipow__", (rhs, *args), {}, _SpyObject())

    def __ilshift__(self, rhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__ilshift__", (rhs,), {}, _SpyObject())

    def __irshift__(self, rhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__irshift__", (rhs,), {}, _SpyObject())

    def __iand__(self, rhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__iand__", (rhs,), {}, _SpyObject())

    def __ixor__(self, rhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__ixor__", (rhs,), {}, _SpyObject())

    def __ior__(self, rhs: object, /) -> "_SpyObject":
        return self.__optype_trace_add__("__ior__", (rhs,), {}, _SpyObject())

    ###

    def __neg__(self, /) -> "_SpyObject":
        return self.__optype_trace_add__("__neg__", (), {}, _SpyObject())

    def __pos__(self, /) -> "_SpyObject":
        return self.__optype_trace_add__("__pos__", (), {}, _SpyObject())

    def __abs__(self, /) -> "_SpyObject":
        return self.__optype_trace_add__("__abs__", (), {}, _SpyObject())

    def __invert__(self, /) -> "_SpyObject":
        return self.__optype_trace_add__("__invert__", (), {}, _SpyObject())

    ###

    def __complex__(self, /) -> complex:
        return self.__optype_trace_add__("__complex__", (), {}, 0j)

    def __float__(self, /) -> float:
        return self.__optype_trace_add__("__float__", (), {}, 0.0)

    def __int__(self, /) -> int:
        return self.__optype_trace_add__("__int__", (), {}, int(_decide()))

    def __index__(self, /) -> int:
        return self.__optype_trace_add__("__index__", (), {}, int(_decide()))

    ###

    def __round__(self, /, *args: object) -> "_SpyObject":
        return self.__optype_trace_add__("__round__", args, {}, _SpyObject())

    def __trunc__(self, /) -> "_SpyObject":
        return self.__optype_trace_add__("__trunc__", (), {}, _SpyObject())

    def __floor__(self, /) -> "_SpyObject":
        return self.__optype_trace_add__("__floor__", (), {}, _SpyObject())

    def __ceil__(self, /) -> "_SpyObject":
        return self.__optype_trace_add__("__ceil__", (), {}, _SpyObject())

    ###

    def __enter__(self, /) -> "_SpyObject":
        return self.__optype_trace_add__("__enter__", (), {}, _SpyObject())

    def __exit__(self, /, *args: object) -> None:
        # TODO: maybe fork and return falsy/truthy in case of exception??
        return self.__optype_trace_add__("__exit__", args, {}, None)

    ###

    def __buffer__(self, flags: int, /) -> memoryview:
        return self.__optype_trace_add__("__buffer__", (flags,), {}, memoryview(b""))

    def __release_buffer__(self, buffer: memoryview, /) -> None:
        return self.__optype_trace_add__("__releasebuffer__", (buffer,), {}, None)

    ###

    def __await__(self, /) -> Generator[Any, None, "_SpyObject"]:
        out = _SpyObject()

        def spy_generator() -> Generator[Any, None, "_SpyObject"]:
            yield from ()
            return out  # noqa: B901

        return self.__optype_trace_add__("__await__", (), {}, spy_generator())

    def __aiter__(self, /) -> "_SpyObject":
        return self.__optype_trace_add__("__aiter__", (), {}, _SpyObject())

    def __anext__(self, /) -> "_SpyObject":
        return self.__optype_trace_add__("__anext__", (), {}, _SpyObject())

    def __aenter__(self, /) -> "_SpyObject":
        return self.__optype_trace_add__("__aenter__", (), {}, _SpyObject())

    def __aexit__(self, /, *args: object) -> None:
        return self.__optype_trace_add__("__aexit__", args, {}, None)


# Free functions, not methods: a method would be an unrecorded hole in the proxy.
def _element_of(spy: _SpyObject) -> _SpyObject:
    if (element := spy.__optype_element__) is None:
        element = _SpyObject()
        spy.__optype_element__ = element  # ty:ignore[invalid-assignment]
    return element


def _iterator_of(spy: _SpyObject) -> _SpyObject:
    iterator = _SpyObject()
    iterator.__optype_element__ = _element_of(spy)  # ty:ignore[invalid-assignment]
    iterator.__optype_iterator__ = True
    return iterator
