from collections.abc import Generator
from itertools import zip_longest as zip_longest
from types import FunctionType
from typing import TypeVar

_T = TypeVar("_T")

text_type = str
string_type = str

def with_str_method(cls: _T) -> _T: ...
def with_repr_method(cls: _T) -> _T: ...
def get_methods(cls: object) -> Generator[tuple[str, FunctionType]]: ...
