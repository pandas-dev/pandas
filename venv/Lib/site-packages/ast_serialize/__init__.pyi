from typing import TypedDict, type_check_only
from typing_extensions import NotRequired, TypeAlias

__all__ = ["parse"]

_TypeIgnores: TypeAlias = list[tuple[int, list[str]]]

@type_check_only
class ParseError(TypedDict):
    line: int
    column: int
    message: str
    blocker: NotRequired[bool]
    code: NotRequired[str]

@type_check_only
class _ASTData(TypedDict):
    is_partial_package: bool
    uses_template_strings: bool
    mypy_ignores: _TypeIgnores
    source_hash: str
    mypy_comments: list[tuple[int, str]]

def parse(
    fnam: str,
    source: str | bytes | None = None,
    skip_function_bodies: bool = False,
    python_version: tuple[int, int] | None = None,
    platform: str | None = None,
    always_true: list[str] | None = None,
    always_false: list[str] | None = None,
    cache_version: int = 0,
) -> tuple[bytes, list[ParseError], _TypeIgnores, bytes, _ASTData]:
    ...
