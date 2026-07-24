from collections.abc import Generator, Mapping, MutableMapping
from io import StringIO
from typing import ClassVar, Literal
from typing_extensions import TypeAlias

def b2u(string: bytes | str) -> str: ...
def u2b(string: str | bytes) -> bytes: ...

_Quotes: TypeAlias = Literal["'", '"']
_ContextType: TypeAlias = Literal["ARG", "ENV", "LABEL"]

class WordSplitter:
    SQUOTE: ClassVar[_Quotes]
    DQUOTE: ClassVar[_Quotes]
    stream: StringIO
    args: Mapping[str, str] | None
    envs: Mapping[str, str] | None
    quotes: _Quotes | None
    escaped: bool
    def __init__(self, s: str, args: Mapping[str, str] | None = None, envs: Mapping[str, str] | None = None) -> None: ...
    def dequote(self) -> str: ...
    def split(self, maxsplit: int | None = None, dequote: bool = True) -> Generator[str | None]: ...

def extract_key_values(
    env_replace: bool, args: Mapping[str, str], envs: Mapping[str, str], instruction_value: str
) -> list[tuple[str, str]]: ...
def get_key_val_dictionary(
    instruction_value: str,
    env_replace: bool = False,
    args: Mapping[str, str] | None = None,
    envs: Mapping[str, str] | None = None,
) -> dict[str, str]: ...

class Context:
    args: MutableMapping[str, str]
    envs: MutableMapping[str, str]
    labels: MutableMapping[str, str]
    line_args: Mapping[str, str]
    line_envs: Mapping[str, str]
    line_labels: Mapping[str, str]
    def __init__(
        self,
        args: MutableMapping[str, str] | None = None,
        envs: MutableMapping[str, str] | None = None,
        labels: MutableMapping[str, str] | None = None,
        line_args: Mapping[str, str] | None = None,
        line_envs: Mapping[str, str] | None = None,
        line_labels: Mapping[str, str] | None = None,
    ) -> None: ...
    def set_line_value(self, context_type: _ContextType, value: Mapping[str, str]) -> None: ...
    def get_line_value(self, context_type: _ContextType) -> Mapping[str, str]: ...
    def get_values(self, context_type: _ContextType) -> Mapping[str, str]: ...
