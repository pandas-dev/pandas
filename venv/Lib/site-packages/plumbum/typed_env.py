from __future__ import annotations

__lazy_modules__ = {"inspect"}

import abc
import inspect
import os
import typing
from collections.abc import Callable, Iterator, MutableMapping
from typing import Any, Final, Generic, TypeVar

if typing.TYPE_CHECKING:
    from ._compat.typing import Self


class _NoDefault:
    __slots__ = ()


NO_DEFAULT: Final[_NoDefault] = _NoDefault()


K = TypeVar("K", bound=str)
V = TypeVar("V")
T = TypeVar("T")


# must not inherit from AttributeError, so not to mess with python's attribute-lookup flow
class EnvironmentVariableError(KeyError):
    pass


class _BaseVar(Generic[V], metaclass=abc.ABCMeta):
    __slots__ = ("default", "name", "names")

    def __init__(
        self,
        name: str | tuple[str, ...] | list[str],
        default: V | _NoDefault = NO_DEFAULT,
    ):
        self.names: tuple[str, ...] = (
            tuple(name) if isinstance(name, (tuple, list)) else (name,)
        )
        self.name = self.names[0]
        self.default = default

    @abc.abstractmethod
    def convert(self, value: str) -> V:
        pass

    @typing.overload
    def __get__(self, instance: TypedEnv, owner: type[TypedEnv]) -> V: ...

    @typing.overload
    def __get__(self, instance: None, owner: type[TypedEnv]) -> Self: ...

    def __get__(self, instance: TypedEnv | None, owner: type[TypedEnv]) -> V | Self:
        if not instance:
            return self
        try:
            return self.convert(instance._raw_get(*self.names))
        except EnvironmentVariableError:
            if isinstance(self.default, _NoDefault):
                raise
            return self.default

    def __set__(self, instance: TypedEnv, value: V) -> None:
        instance[self.name] = value


class Str(_BaseVar[str]):
    __slots__ = ()

    def convert(self, value: str) -> str:
        return value


class Bool(_BaseVar[bool]):
    """
    Converts 'yes|true|1|no|false|0' to the appropriate boolean value.
    Case-insensitive. Throws a ``ValueError`` for any other value.
    """

    __slots__ = ()

    def convert(self, value: str) -> bool:
        value = value.lower()
        if value not in {"yes", "no", "true", "false", "1", "0"}:
            raise ValueError(f"Unrecognized boolean value: {value!r}")
        return value in {"yes", "true", "1"}

    def __set__(self, instance: TypedEnv, value: bool) -> None:
        instance[self.name] = "yes" if value else "no"


class Int(_BaseVar[int]):
    __slots__ = ()

    def convert(self, value: str) -> int:
        return int(value)


class Float(_BaseVar[float]):
    __slots__ = ()

    def convert(self, value: str) -> float:
        return float(value)


class CSV(_BaseVar[list[str]]):
    __slots__ = ("separator", "type")
    """
    Comma-separated-strings get split using the ``separator`` (',' by default) into
    a list of objects of type ``type`` (``str`` by default).
    """

    def __init__(
        self,
        name: str,
        default: list[str] | _NoDefault = NO_DEFAULT,
        type: Callable[[str], str] = str,
        separator: str = ",",
    ):  # pylint:disable=redefined-builtin
        super().__init__(name, default=default)
        self.type = type
        self.separator = separator

    def __set__(self, instance: TypedEnv, value: list[str]) -> None:
        instance[self.name] = self.separator.join(map(str, value))

    def convert(self, value: str) -> list[str]:
        return [self.type(v.strip()) for v in value.split(self.separator)]


class TypedEnv(MutableMapping[str, str]):
    """
    This object can be used in 'exploratory' mode:

        nv = TypedEnv()
        print(nv.HOME)

    It can also be used as a parser and validator of environment variables:

    class MyEnv(TypedEnv):
        username = TypedEnv.Str("USER", default='anonymous')
        path = TypedEnv.CSV("PATH", separator=":")
        tmp = TypedEnv.Str("TMP TEMP".split())  # support 'fallback' var-names

    nv = MyEnv()

    print(nv.username)

    for p in nv.path:
        print(p)

    try:
        print(p.tmp)
    except EnvironmentVariableError:
        print("TMP/TEMP is not defined")
    else:
        assert False
    """

    __slots__ = ("_defined_keys", "_env")

    Str = Str
    Bool = Bool
    Int = Int
    Float = Float
    CSV = CSV

    def __init__(self, env: MutableMapping[str, str] | None = None):
        if env is None:
            env = os.environ
        self._env = env
        self._defined_keys = {
            k
            for (k, v) in inspect.getmembers(self.__class__)
            if isinstance(v, _BaseVar)
        }

    def __iter__(self) -> Iterator[str]:
        # Iterate the declared keys if any are defined (so that dict(env),
        # keys(), etc. go through the typed descriptors); otherwise fall back
        # to the raw environment keys. Unlike __dir__, this deliberately does
        # not include class members.
        if self._defined_keys:
            return iter(sorted(self._defined_keys))
        return iter(self._env)

    def __len__(self) -> int:
        return len(self._defined_keys) if self._defined_keys else len(self._env)

    def __delitem__(self, name: str) -> None:
        del self._env[name]

    def __setitem__(self, name: str, value: Any) -> None:
        self._env[name] = str(value)

    def _raw_get(self, *key_names: str) -> str:
        for key in key_names:
            value = self._env.get(key, NO_DEFAULT)
            if not isinstance(value, _NoDefault):
                return value
        raise EnvironmentVariableError(key_names[0])

    def __contains__(self, key: object) -> bool:
        # Stay consistent with __iter__/__len__: when descriptors are defined,
        # membership is over the declared keys; otherwise the raw environment.
        if self._defined_keys:
            return key in self._defined_keys
        return key in self._env

    def __getattr__(self, name: str) -> str:
        # if we're here then there was no descriptor defined
        try:
            return self._raw_get(name)
        except EnvironmentVariableError:
            raise AttributeError(
                f"{self.__class__} has no attribute {name!r}"
            ) from None

    def __getitem__(self, key: str) -> Any:
        return getattr(self, key)  # delegate through the descriptors

    @typing.overload
    def get(self, key: str) -> str | None: ...

    @typing.overload
    def get(self, key: str, default: str = ...) -> str: ...

    @typing.overload
    def get(self, key: str, default: T = ...) -> str | T: ...

    def get(self, key: str, default: str | None | T = None) -> str | None | T:
        return getattr(self, key, default)

    def __dir__(self) -> list[str]:
        if self._defined_keys:
            # return only defined
            return sorted(self._defined_keys)
        # return whatever is in the environment (for convenience)
        members = set(self._env.keys())
        members.update(dir(self.__class__))
        return sorted(members)


__all__ = [
    "CSV",
    "NO_DEFAULT",
    "Bool",
    "EnvironmentVariableError",
    "Float",
    "Int",
    "Str",
    "TypedEnv",
]


def __dir__() -> list[str]:
    return list(__all__)
