from collections.abc import Callable, Mapping, Sequence
from typing import ClassVar, Final, Generic, Literal, NamedTuple, TypeAlias
from typing_extensions import TypeVar, override

import optype as op

_FT = TypeVar("_FT", bound=Callable[..., object], default=Callable[..., object])
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
    def __iter__(self, /) -> op.CanIterSelf[str]: ...
    @override
    def __getitem__(self, key: str, /) -> _SectionValue: ...
    def __setitem__(self, key: str, val: _SectionValue, /) -> None: ...

class FunctionDoc(NumpyDocString, Generic[_FT]):
    def __init__(
        self,
        /,
        func: _FT,
        role: Literal["func", "meth"] = "func",
        doc: str | None = None,
        config: Mapping[str, object] | None = None,
    ) -> None: ...
    @override
    def __str__(self, /) -> str: ...  # type: ignore[override]  # noqa: PYI029  # pyright: ignore[reportIncompatibleMethodOverride]
    def get_func(self, /) -> _FT: ...

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

def strip_blank_lines(l: list[str]) -> list[str]: ...
def dedent_lines(lines: op.CanIter[op.CanIterSelf[str]]) -> str: ...
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
