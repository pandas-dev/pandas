import enum
import sys
from array import array
from collections.abc import (
    Callable,
    ItemsView,
    Iterable,
    Iterator,
    KeysView,
    Mapping,
    ValuesView,
)
from typing import (
    TYPE_CHECKING,
    Generic,
    NoReturn,
    TypeVar,
    Union,
    cast,
    overload,
)

from ._abc import MDArg, MultiMapping, MutableMultiMapping, SupportsKeys

if sys.version_info >= (3, 11):
    from typing import Self
else:
    from typing_extensions import Self


class istr(str):
    """Case insensitive str."""

    __is_istr__ = True


_V = TypeVar("_V")
_T = TypeVar("_T")

_SENTINEL = enum.Enum("_SENTINEL", "sentinel")
sentinel = _SENTINEL.sentinel

_version = array("Q", [0])


class _Impl(Generic[_V]):
    __slots__ = ("_items", "_version")

    def __init__(self) -> None:
        self._items: list[tuple[str, str, _V]] = []
        self.incr_version()

    def incr_version(self) -> None:
        global _version
        v = _version
        v[0] += 1
        self._version = v[0]

    if sys.implementation.name != "pypy":

        def __sizeof__(self) -> int:
            return object.__sizeof__(self) + sys.getsizeof(self._items)


class _Iter(Generic[_T]):
    __slots__ = ("_size", "_iter")

    def __init__(self, size: int, iterator: Iterator[_T]):
        self._size = size
        self._iter = iterator

    def __iter__(self) -> Self:
        return self

    def __next__(self) -> _T:
        return next(self._iter)

    def __length_hint__(self) -> int:
        return self._size


class _ViewBase(Generic[_V]):
    def __init__(self, impl: _Impl[_V]):
        self._impl = impl

    def __len__(self) -> int:
        return len(self._impl._items)


class _ItemsView(_ViewBase[_V], ItemsView[str, _V]):
    def __contains__(self, item: object) -> bool:
        if not isinstance(item, (tuple, list)) or len(item) != 2:
            return False
        for i, k, v in self._impl._items:
            if item[0] == k and item[1] == v:
                return True
        return False

    def __iter__(self) -> _Iter[tuple[str, _V]]:
        return _Iter(len(self), self._iter(self._impl._version))

    def _iter(self, version: int) -> Iterator[tuple[str, _V]]:
        for i, k, v in self._impl._items:
            if version != self._impl._version:
                raise RuntimeError("Dictionary changed during iteration")
            yield k, v

    def __repr__(self) -> str:
        lst = []
        for item in self._impl._items:
            lst.append("{!r}: {!r}".format(item[1], item[2]))
        body = ", ".join(lst)
        return "{}({})".format(self.__class__.__name__, body)


class _ValuesView(_ViewBase[_V], ValuesView[_V]):
    def __contains__(self, value: object) -> bool:
        for item in self._impl._items:
            if item[2] == value:
                return True
        return False

    def __iter__(self) -> _Iter[_V]:
        return _Iter(len(self), self._iter(self._impl._version))

    def _iter(self, version: int) -> Iterator[_V]:
        for item in self._impl._items:
            if version != self._impl._version:
                raise RuntimeError("Dictionary changed during iteration")
            yield item[2]

    def __repr__(self) -> str:
        lst = []
        for item in self._impl._items:
            lst.append("{!r}".format(item[2]))
        body = ", ".join(lst)
        return "{}({})".format(self.__class__.__name__, body)


class _KeysView(_ViewBase[_V], KeysView[str]):
    def __contains__(self, key: object) -> bool:
        for item in self._impl._items:
            if item[1] == key:
                return True
        return False

    def __iter__(self) -> _Iter[str]:
        return _Iter(len(self), self._iter(self._impl._version))

    def _iter(self, version: int) -> Iterator[str]:
        for item in self._impl._items:
            if version != self._impl._version:
                raise RuntimeError("Dictionary changed during iteration")
            yield item[1]

    def __repr__(self) -> str:
        lst = []
        for item in self._impl._items:
            lst.append("{!r}".format(item[1]))
        body = ", ".join(lst)
        return "{}({})".format(self.__class__.__name__, body)


class _Base(MultiMapping[_V]):
    _impl: _Impl[_V]

    def _title(self, key: str) -> str:
        return key

    @overload
    def getall(self, key: str) -> list[_V]: ...
    @overload
    def getall(self, key: str, default: _T) -> Union[list[_V], _T]: ...
    def getall(
        self, key: str, default: Union[_T, _SENTINEL] = sentinel
    ) -> Union[list[_V], _T]:
        """Return a list of all values matching the key."""
        identity = self._title(key)
        res = [v for i, k, v in self._impl._items if i == identity]
        if res:
            return res
        if not res and default is not sentinel:
            return default
        raise KeyError("Key not found: %r" % key)

    @overload
    def getone(self, key: str) -> _V: ...
    @overload
    def getone(self, key: str, default: _T) -> Union[_V, _T]: ...
    def getone(
        self, key: str, default: Union[_T, _SENTINEL] = sentinel
    ) -> Union[_V, _T]:
        """Get first value matching the key.

        Raises KeyError if the key is not found and no default is provided.
        """
        identity = self._title(key)
        for i, k, v in self._impl._items:
            if i == identity:
                return v
        if default is not sentinel:
            return default
        raise KeyError("Key not found: %r" % key)

    # Mapping interface #

    def __getitem__(self, key: str) -> _V:
        return self.getone(key)

    @overload
    def get(self, key: str, /) -> Union[_V, None]: ...
    @overload
    def get(self, key: str, /, default: _T) -> Union[_V, _T]: ...
    def get(self, key: str, default: Union[_T, None] = None) -> Union[_V, _T, None]:
        """Get first value matching the key.

        If the key is not found, returns the default (or None if no default is provided)
        """
        return self.getone(key, default)

    def __iter__(self) -> Iterator[str]:
        return iter(self.keys())

    def __len__(self) -> int:
        return len(self._impl._items)

    def keys(self) -> KeysView[str]:
        """Return a new view of the dictionary's keys."""
        return _KeysView(self._impl)

    def items(self) -> ItemsView[str, _V]:
        """Return a new view of the dictionary's items *(key, value) pairs)."""
        return _ItemsView(self._impl)

    def values(self) -> _ValuesView[_V]:
        """Return a new view of the dictionary's values."""
        return _ValuesView(self._impl)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Mapping):
            return NotImplemented
        if isinstance(other, _Base):
            lft = self._impl._items
            rht = other._impl._items
            if len(lft) != len(rht):
                return False
            for (i1, k2, v1), (i2, k2, v2) in zip(lft, rht):
                if i1 != i2 or v1 != v2:
                    return False
            return True
        if len(self._impl._items) != len(other):
            return False
        for k, v in self.items():
            nv = other.get(k, sentinel)
            if v != nv:
                return False
        return True

    def __contains__(self, key: object) -> bool:
        if not isinstance(key, str):
            return False
        identity = self._title(key)
        for i, k, v in self._impl._items:
            if i == identity:
                return True
        return False

    def __repr__(self) -> str:
        body = ", ".join("'{}': {!r}".format(k, v) for k, v in self.items())
        return "<{}({})>".format(self.__class__.__name__, body)


class MultiDict(_Base[_V], MutableMultiMapping[_V]):
    """Dictionary with the support for duplicate keys."""

    def __init__(self, arg: MDArg[_V] = None, /, **kwargs: _V):
        self._impl = _Impl()

        self._extend(arg, kwargs, self.__class__.__name__, self._extend_items)

    if sys.implementation.name != "pypy":

        def __sizeof__(self) -> int:
            return object.__sizeof__(self) + sys.getsizeof(self._impl)

    def __reduce__(self) -> tuple[type[Self], tuple[list[tuple[str, _V]]]]:
        return (self.__class__, (list(self.items()),))

    def _title(self, key: str) -> str:
        return key

    def _key(self, key: str) -> str:
        if isinstance(key, str):
            return key
        else:
            raise TypeError("MultiDict keys should be either str or subclasses of str")

    def add(self, key: str, value: _V) -> None:
        identity = self._title(key)
        self._impl._items.append((identity, self._key(key), value))
        self._impl.incr_version()

    def copy(self) -> Self:
        """Return a copy of itself."""
        cls = self.__class__
        return cls(self.items())

    __copy__ = copy

    def extend(self, arg: MDArg[_V] = None, /, **kwargs: _V) -> None:
        """Extend current MultiDict with more values.

        This method must be used instead of update.
        """
        self._extend(arg, kwargs, "extend", self._extend_items)

    def _extend(
        self,
        arg: MDArg[_V],
        kwargs: Mapping[str, _V],
        name: str,
        method: Callable[[list[tuple[str, str, _V]]], None],
    ) -> None:
        if arg:
            if isinstance(arg, (MultiDict, MultiDictProxy)) and not kwargs:
                items = arg._impl._items
            else:
                if hasattr(arg, "keys"):
                    arg = cast(SupportsKeys[_V], arg)
                    arg = [(k, arg[k]) for k in arg.keys()]
                if kwargs:
                    arg = list(arg)
                    arg.extend(list(kwargs.items()))
                items = []
                for item in arg:
                    if not len(item) == 2:
                        raise TypeError(
                            "{} takes either dict or list of (key, value) "
                            "tuples".format(name)
                        )
                    items.append((self._title(item[0]), self._key(item[0]), item[1]))

            method(items)
        else:
            method(
                [
                    (self._title(key), self._key(key), value)
                    for key, value in kwargs.items()
                ]
            )

    def _extend_items(self, items: Iterable[tuple[str, str, _V]]) -> None:
        for identity, key, value in items:
            self.add(key, value)

    def clear(self) -> None:
        """Remove all items from MultiDict."""
        self._impl._items.clear()
        self._impl.incr_version()

    # Mapping interface #

    def __setitem__(self, key: str, value: _V) -> None:
        self._replace(key, value)

    def __delitem__(self, key: str) -> None:
        identity = self._title(key)
        items = self._impl._items
        found = False
        for i in range(len(items) - 1, -1, -1):
            if items[i][0] == identity:
                del items[i]
                found = True
        if not found:
            raise KeyError(key)
        else:
            self._impl.incr_version()

    @overload
    def setdefault(
        self: "MultiDict[Union[_T, None]]", key: str, default: None = None
    ) -> Union[_T, None]: ...
    @overload
    def setdefault(self, key: str, default: _V) -> _V: ...
    def setdefault(self, key: str, default: Union[_V, None] = None) -> Union[_V, None]:  # type: ignore[misc]
        """Return value for key, set value to default if key is not present."""
        identity = self._title(key)
        for i, k, v in self._impl._items:
            if i == identity:
                return v
        self.add(key, default)  # type: ignore[arg-type]
        return default

    @overload
    def popone(self, key: str) -> _V: ...
    @overload
    def popone(self, key: str, default: _T) -> Union[_V, _T]: ...
    def popone(
        self, key: str, default: Union[_T, _SENTINEL] = sentinel
    ) -> Union[_V, _T]:
        """Remove specified key and return the corresponding value.

        If key is not found, d is returned if given, otherwise
        KeyError is raised.

        """
        identity = self._title(key)
        for i in range(len(self._impl._items)):
            if self._impl._items[i][0] == identity:
                value = self._impl._items[i][2]
                del self._impl._items[i]
                self._impl.incr_version()
                return value
        if default is sentinel:
            raise KeyError(key)
        else:
            return default

    # Type checking will inherit signature for pop() if we don't confuse it here.
    if not TYPE_CHECKING:
        pop = popone

    @overload
    def popall(self, key: str) -> list[_V]: ...
    @overload
    def popall(self, key: str, default: _T) -> Union[list[_V], _T]: ...
    def popall(
        self, key: str, default: Union[_T, _SENTINEL] = sentinel
    ) -> Union[list[_V], _T]:
        """Remove all occurrences of key and return the list of corresponding
        values.

        If key is not found, default is returned if given, otherwise
        KeyError is raised.

        """
        found = False
        identity = self._title(key)
        ret = []
        for i in range(len(self._impl._items) - 1, -1, -1):
            item = self._impl._items[i]
            if item[0] == identity:
                ret.append(item[2])
                del self._impl._items[i]
                self._impl.incr_version()
                found = True
        if not found:
            if default is sentinel:
                raise KeyError(key)
            else:
                return default
        else:
            ret.reverse()
            return ret

    def popitem(self) -> tuple[str, _V]:
        """Remove and return an arbitrary (key, value) pair."""
        if self._impl._items:
            i = self._impl._items.pop(0)
            self._impl.incr_version()
            return i[1], i[2]
        else:
            raise KeyError("empty multidict")

    def update(self, arg: MDArg[_V] = None, /, **kwargs: _V) -> None:
        """Update the dictionary from *other*, overwriting existing keys."""
        self._extend(arg, kwargs, "update", self._update_items)

    def _update_items(self, items: list[tuple[str, str, _V]]) -> None:
        if not items:
            return
        used_keys: dict[str, int] = {}
        for identity, key, value in items:
            start = used_keys.get(identity, 0)
            for i in range(start, len(self._impl._items)):
                item = self._impl._items[i]
                if item[0] == identity:
                    used_keys[identity] = i + 1
                    self._impl._items[i] = (identity, key, value)
                    break
            else:
                self._impl._items.append((identity, key, value))
                used_keys[identity] = len(self._impl._items)

        # drop tails
        i = 0
        while i < len(self._impl._items):
            item = self._impl._items[i]
            identity = item[0]
            pos = used_keys.get(identity)
            if pos is None:
                i += 1
                continue
            if i >= pos:
                del self._impl._items[i]
            else:
                i += 1

        self._impl.incr_version()

    def _replace(self, key: str, value: _V) -> None:
        key = self._key(key)
        identity = self._title(key)
        items = self._impl._items

        for i in range(len(items)):
            item = items[i]
            if item[0] == identity:
                items[i] = (identity, key, value)
                # i points to last found item
                rgt = i
                self._impl.incr_version()
                break
        else:
            self._impl._items.append((identity, key, value))
            self._impl.incr_version()
            return

        # remove all tail items
        # Mypy bug: https://github.com/python/mypy/issues/14209
        i = rgt + 1  # type: ignore[possibly-undefined]
        while i < len(items):
            item = items[i]
            if item[0] == identity:
                del items[i]
            else:
                i += 1


class CIMultiDict(MultiDict[_V]):
    """Dictionary with the support for duplicate case-insensitive keys."""

    def _title(self, key: str) -> str:
        return key.title()


class MultiDictProxy(_Base[_V]):
    """Read-only proxy for MultiDict instance."""

    def __init__(self, arg: Union[MultiDict[_V], "MultiDictProxy[_V]"]):
        if not isinstance(arg, (MultiDict, MultiDictProxy)):
            raise TypeError(
                "ctor requires MultiDict or MultiDictProxy instance"
                ", not {}".format(type(arg))
            )

        self._impl = arg._impl

    def __reduce__(self) -> NoReturn:
        raise TypeError("can't pickle {} objects".format(self.__class__.__name__))

    def copy(self) -> MultiDict[_V]:
        """Return a copy of itself."""
        return MultiDict(self.items())


class CIMultiDictProxy(MultiDictProxy[_V]):
    """Read-only proxy for CIMultiDict instance."""

    def __init__(self, arg: Union[MultiDict[_V], MultiDictProxy[_V]]):
        if not isinstance(arg, (CIMultiDict, CIMultiDictProxy)):
            raise TypeError(
                "ctor requires CIMultiDict or CIMultiDictProxy instance"
                ", not {}".format(type(arg))
            )

        self._impl = arg._impl

    def _title(self, key: str) -> str:
        return key.title()

    def copy(self) -> CIMultiDict[_V]:
        """Return a copy of itself."""
        return CIMultiDict(self.items())


def getversion(md: Union[MultiDict[object], MultiDictProxy[object]]) -> int:
    if not isinstance(md, _Base):
        raise TypeError("Parameter should be multidict or proxy")
    return md._impl._version
