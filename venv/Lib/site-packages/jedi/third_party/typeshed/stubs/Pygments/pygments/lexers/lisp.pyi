from _typeshed import Incomplete
from collections.abc import Iterator
from typing import ClassVar

from ..lexer import RegexLexer
from ..token import _TokenType

__all__ = [
    "SchemeLexer",
    "CommonLispLexer",
    "HyLexer",
    "RacketLexer",
    "NewLispLexer",
    "EmacsLispLexer",
    "ShenLexer",
    "CPSALexer",
    "XtlangLexer",
    "FennelLexer",
]

class SchemeLexer(RegexLexer):
    valid_name: ClassVar[str]
    token_end: ClassVar[str]
    def get_tokens_unprocessed(self, text: str) -> Iterator[tuple[int, _TokenType, str]]: ...  # type: ignore[override]
    number_rules: ClassVar[dict[Incomplete, Incomplete]]
    def decimal_cb(self, match) -> Iterator[tuple[Incomplete, Incomplete, Incomplete]]: ...

class CommonLispLexer(RegexLexer):
    nonmacro: ClassVar[str]
    constituent: ClassVar[str]
    terminated: ClassVar[str]
    symbol: ClassVar[str]
    def get_tokens_unprocessed(self, text: str) -> Iterator[tuple[int, _TokenType, str]]: ...  # type: ignore[override]

class HyLexer(RegexLexer):
    special_forms: ClassVar[tuple[str, ...]]
    declarations: ClassVar[tuple[str, ...]]
    hy_builtins: ClassVar[tuple[str, ...]]
    hy_core: ClassVar[tuple[str, ...]]
    builtins: ClassVar[tuple[str, ...]]
    valid_name: ClassVar[str]

class RacketLexer(RegexLexer): ...

class NewLispLexer(RegexLexer):
    builtins: ClassVar[tuple[str, ...]]
    valid_name: ClassVar[str]

class EmacsLispLexer(RegexLexer):
    nonmacro: ClassVar[str]
    constituent: ClassVar[str]
    terminated: ClassVar[str]
    symbol: ClassVar[str]
    macros: ClassVar[set[str]]
    special_forms: ClassVar[set[str]]
    builtin_function: ClassVar[set[str]]
    builtin_function_highlighted: ClassVar[set[str]]
    lambda_list_keywords: ClassVar[set[str]]
    error_keywords: ClassVar[set[str]]
    def get_tokens_unprocessed(self, text: str) -> Iterator[tuple[int, _TokenType, str]]: ...  # type: ignore[override]

class ShenLexer(RegexLexer):
    DECLARATIONS: ClassVar[tuple[str, ...]]
    SPECIAL_FORMS: ClassVar[tuple[str, ...]]
    BUILTINS: ClassVar[tuple[str, ...]]
    BUILTINS_ANYWHERE: ClassVar[tuple[str, ...]]
    MAPPINGS: ClassVar[dict[str, Incomplete]]

    valid_symbol_chars: ClassVar[str]
    valid_name: ClassVar[str]
    symbol_name: ClassVar[str]
    variable: ClassVar[str]

    def get_tokens_unprocessed(self, text: str) -> Iterator[tuple[int, _TokenType, str]]: ...  # type: ignore[override]

class CPSALexer(RegexLexer):
    valid_name: ClassVar[str]

class XtlangLexer(RegexLexer):
    common_keywords: ClassVar[tuple[str, ...]]
    scheme_keywords: ClassVar[tuple[str, ...]]
    xtlang_bind_keywords: ClassVar[tuple[str, ...]]
    xtlang_keywords: ClassVar[tuple[str, ...]]
    common_functions: ClassVar[tuple[str, ...]]
    scheme_functions: ClassVar[tuple[str, ...]]
    xtlang_functions: ClassVar[tuple[str, ...]]

    valid_scheme_name: ClassVar[str]
    valid_xtlang_name: ClassVar[str]
    valid_xtlang_type: ClassVar[str]

class FennelLexer(RegexLexer):
    special_forms: ClassVar[tuple[str, ...]]
    declarations: ClassVar[tuple[str, ...]]
    builtins: ClassVar[tuple[str, ...]]
    valid_name: ClassVar[str]
