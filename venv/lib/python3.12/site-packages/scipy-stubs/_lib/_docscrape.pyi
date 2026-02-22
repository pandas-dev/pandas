from collections.abc import Callable, Iterable, Iterator, Mapping, MutableSequence, Sequence
from typing import ClassVar, Final, Generic, Literal, LiteralString, NamedTuple, TypeAlias
from typing_extensions import TypeVar, override

import optype as op

_FunctionT = TypeVar("_FunctionT", bound=Callable[..., object], default=Callable[..., object])
_MutableStrSequenceT = TypeVar("_MutableStrSequenceT", bound=MutableSequence[str])
_AnyStrT = TypeVar("_AnyStrT", str, LiteralString)

_SectionValue: TypeAlias = str | list[str] | Mapping[str, list[str]]

###

class ParseError(Exception): ...

class Parameter(NamedTuple):
    name: str
    type: str
    desc: str

class Reader:
    def __init__(self, /, data: str | list[str]) -> None: ...
    def __getitem__(self, n: op.CanIndex, /) -> str: ...
    def reset(self, /) -> None: ...
    def read(self, /) -> str: ...
    def seek_next_non_empty_line(self, /) -> None: ...
    def eof(self, /) -> bool: ...
    def read_to_condition(self, /, condition_func: Callable[[str], object]) -> list[str]: ...
    def read_to_next_empty_line(self, /) -> list[str]: ...
    def read_to_next_unindented_line(self, /) -> list[str]: ...
    def peek(self, /, n: int = 0) -> str: ...
    def is_empty(self, /) -> bool: ...

class NumpyDocString(Mapping[str, _SectionValue]):
    empty_description: ClassVar[str] = ".."
    sections: ClassVar[Mapping[str, _SectionValue]]

    def __init__(self, /, docstring: str, config: Mapping[str, object] | None = None) -> None: ...
    @override
    def __str__(self, /, func_role: str = "") -> str: ...
    @override
    def __len__(self, /) -> int: ...
    @override
    def __iter__(self, /) -> Iterator[str]: ...
    @override
    def __getitem__(self, key: str, /) -> _SectionValue: ...
    def __setitem__(self, key: str, val: _SectionValue, /) -> None: ...

class FunctionDoc(NumpyDocString, Generic[_FunctionT]):
    def __init__(
        self,
        /,
        func: _FunctionT,
        role: Literal["func", "meth"] = "func",
        doc: str | None = None,
        config: Mapping[str, object] | None = None,
    ) -> None: ...
    @override
    def __str__(self, /) -> str: ...  # type: ignore[override] # pyright: ignore[reportIncompatibleMethodOverride]  # ty: ignore[invalid-method-override]  # noqa: PYI029
    def get_func(self, /) -> _FunctionT: ...

class ClassDoc(NumpyDocString):
    extra_public_methods: ClassVar[Sequence[str]]
    show_inherited_members: Final[bool]

    @property
    def methods(self, /) -> list[str]: ...
    @property
    def properties(self, /) -> list[str]: ...
    def __init__(
        self,
        /,
        cls: type,
        doc: str | None = None,
        modulename: str = "",
        # NOTE: we can't set this to `type[FunctionDoc]` because of a mypy bug
        func_doc: type = ...,
        config: Mapping[str, object] | None = None,
    ) -> None: ...

class ObjDoc(NumpyDocString):
    def __init__(self, /, obj: object, doc: str | None = None, config: Mapping[str, object] | None = None) -> None: ...

def strip_blank_lines(l: _MutableStrSequenceT) -> _MutableStrSequenceT: ...
def dedent_lines(lines: Iterable[_AnyStrT]) -> _AnyStrT: ...
def get_doc_object(
    obj: object,
    what: str | None = None,
    doc: str | None = None,
    config: Mapping[str, object] | None = None,
    # NOTE: due to a mypy bug, these can't be made any more specific
    class_doc: type = ...,
    func_doc: type = ...,
    obj_doc: type = ...,
) -> NumpyDocString: ...
