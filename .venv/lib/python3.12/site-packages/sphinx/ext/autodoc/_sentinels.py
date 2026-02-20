from __future__ import annotations

TYPE_CHECKING = False
if TYPE_CHECKING:
    from typing import Final, Literal, NoReturn, Self, _SpecialForm


class _Sentinel:
    """Create a unique sentinel object."""

    __slots__ = ('_name',)

    _name: str

    def __new__(cls, name: str, /) -> Self:
        sentinel = super().__new__(cls)
        object.__setattr__(sentinel, '_name', str(name))
        return sentinel

    def __repr__(self) -> str:
        return self._name

    def __setattr__(self, key: str, value: object) -> NoReturn:
        msg = f'{self._name} is immutable'
        raise AttributeError(msg)

    def __or__(self, other: object) -> _SpecialForm:
        from typing import Union

        return Union[self, other]  # NoQA: UP007  # ty: ignore[invalid-type-form, invalid-return-type]

    def __ror__(self, other: object) -> _SpecialForm:
        from typing import Union

        return Union[other, self]  # NoQA: UP007  # ty: ignore[invalid-type-form, invalid-return-type]

    def __getstate__(self) -> NoReturn:
        msg = f'Cannot pickle {self._name}'
        raise TypeError(msg)


class _All(_Sentinel):
    """A special value for :*-members: that matches to any member."""

    def __contains__(self, item: object) -> Literal[True]:
        return True

    def append(self, item: object) -> None:
        pass  # nothing


class _Empty(_Sentinel):
    """A special value for :exclude-members: that never matches to any member."""

    def __contains__(self, item: object) -> Literal[False]:
        return False


if TYPE_CHECKING:
    # For the sole purpose of satisfying the type checker.
    # fmt: off
    import enum
    class _AllTC(enum.Enum):
        ALL = enum.auto()

        def __contains__(self, item: object) -> Literal[True]: return True
        def __add__(self, other: object) -> Self: pass
    type ALL_T = Literal[_AllTC.ALL]
    ALL: Final[ALL_T] = _AllTC.ALL

    class _EmptyTC(enum.Enum):
        EMPTY = enum.auto()

        def __contains__(self, item: object) -> Literal[False]: return False
    type EMPTY_T = Literal[_EmptyTC.EMPTY]
    EMPTY: Final[EMPTY_T] = _EmptyTC.EMPTY

    class _SentinelTC(enum.Enum):
        INSTANCE_ATTR = enum.auto()
        RUNTIME_INSTANCE_ATTRIBUTE = enum.auto()
        SLOTS_ATTR = enum.auto()
        SUPPRESS = enum.auto()
        UNINITIALIZED_ATTR = enum.auto()
    type INSTANCE_ATTR_T = Literal[_SentinelTC.INSTANCE_ATTR]
    type RUNTIME_INSTANCE_ATTRIBUTE_T = Literal[
        _SentinelTC.RUNTIME_INSTANCE_ATTRIBUTE
    ]
    type SLOTS_ATTR_T = Literal[_SentinelTC.SLOTS_ATTR]
    type SUPPRESS_T = Literal[_SentinelTC.SUPPRESS]
    type UNINITIALIZED_ATTR_T = Literal[_SentinelTC.UNINITIALIZED_ATTR]
    INSTANCE_ATTR: Final[INSTANCE_ATTR_T] = _SentinelTC.INSTANCE_ATTR
    RUNTIME_INSTANCE_ATTRIBUTE: Final[RUNTIME_INSTANCE_ATTRIBUTE_T] = (
        _SentinelTC.RUNTIME_INSTANCE_ATTRIBUTE
    )
    SLOTS_ATTR: Final[SLOTS_ATTR_T] = _SentinelTC.SLOTS_ATTR
    SUPPRESS: Final[SUPPRESS_T] = _SentinelTC.SUPPRESS
    UNINITIALIZED_ATTR: Final[UNINITIALIZED_ATTR_T] = _SentinelTC.UNINITIALIZED_ATTR
    # fmt: on
else:
    ALL = _All('ALL')
    EMPTY = _Empty('EMPTY')
    INSTANCE_ATTR = _Sentinel('INSTANCE_ATTR')
    RUNTIME_INSTANCE_ATTRIBUTE = _Sentinel('RUNTIME_INSTANCE_ATTRIBUTE')
    SLOTS_ATTR = _Sentinel('SLOTS_ATTR')
    SUPPRESS = _Sentinel('SUPPRESS')
    UNINITIALIZED_ATTR = _Sentinel('UNINITIALIZED_ATTR')
