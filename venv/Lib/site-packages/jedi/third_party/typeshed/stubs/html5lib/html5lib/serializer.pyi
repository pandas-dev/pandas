from _typeshed import Incomplete
from collections.abc import Generator
from typing import Literal, overload

def htmlentityreplace_errors(exc: UnicodeError) -> tuple[str | bytes, int]: ...
@overload
def serialize(
    input,
    tree: Literal["dom", "genshi", "lxml", "etree"] = "etree",
    encoding: Literal[""] | None = None,
    *,
    quote_attr_values: Literal["legacy", "spec", "always"] = "legacy",
    quote_char: str = '"',
    use_best_quote_char: bool = ...,  # default value depends on whether quote_char was passed
    omit_optional_tags: bool = True,
    minimize_boolean_attributes: bool = True,
    use_trailing_solidus: bool = False,
    space_before_trailing_solidus: bool = True,
    escape_lt_in_attrs: bool = False,
    escape_rcdata: bool = False,
    resolve_entities: bool = True,
    alphabetical_attributes: bool = False,
    inject_meta_charset: bool = True,
    strip_whitespace: bool = False,
    sanitize: bool = False,
) -> str: ...
@overload
def serialize(
    input,
    tree: Literal["dom", "genshi", "lxml", "etree"] = "etree",
    encoding: str = ...,
    *,
    quote_attr_values: Literal["legacy", "spec", "always"] = "legacy",
    quote_char: str = '"',
    use_best_quote_char: bool = ...,  # default value depends on whether quote_char was passed
    omit_optional_tags: bool = True,
    minimize_boolean_attributes: bool = True,
    use_trailing_solidus: bool = False,
    space_before_trailing_solidus: bool = True,
    escape_lt_in_attrs: bool = False,
    escape_rcdata: bool = False,
    resolve_entities: bool = True,
    alphabetical_attributes: bool = False,
    inject_meta_charset: bool = True,
    strip_whitespace: bool = False,
    sanitize: bool = False,
) -> bytes: ...

class HTMLSerializer:
    quote_attr_values: Literal["legacy", "spec", "always"]
    quote_char: str
    use_best_quote_char: bool
    omit_optional_tags: bool
    minimize_boolean_attributes: bool
    use_trailing_solidus: bool
    space_before_trailing_solidus: bool
    escape_lt_in_attrs: bool
    escape_rcdata: bool
    resolve_entities: bool
    alphabetical_attributes: bool
    inject_meta_charset: bool
    strip_whitespace: bool
    sanitize: bool
    options: tuple[str, ...]
    errors: list[Incomplete]
    strict: bool
    def __init__(
        self,
        *,
        quote_attr_values: Literal["legacy", "spec", "always"] = "legacy",
        quote_char: str = '"',
        use_best_quote_char: bool = ...,  # default value depends on whether quote_char was passed
        omit_optional_tags: bool = True,
        minimize_boolean_attributes: bool = True,
        use_trailing_solidus: bool = False,
        space_before_trailing_solidus: bool = True,
        escape_lt_in_attrs: bool = False,
        escape_rcdata: bool = False,
        resolve_entities: bool = True,
        alphabetical_attributes: bool = False,
        inject_meta_charset: bool = True,
        strip_whitespace: bool = False,
        sanitize: bool = False,
    ) -> None: ...
    def encode(self, string: str) -> str | bytes: ...  # result depends on self.encoding
    def encodeStrict(self, string: str) -> str | bytes: ...  # result depends on self.encoding
    encoding: str | None
    @overload
    def serialize(self, treewalker, encoding: Literal[""] | None = None) -> Generator[str]: ...
    @overload
    def serialize(self, treewalker, encoding: str = ...) -> Generator[bytes]: ...
    @overload
    def render(self, treewalker, encoding: Literal[""] | None = None) -> str: ...
    @overload
    def render(self, treewalker, encoding: str = ...) -> bytes: ...
    def serializeError(self, data="XXX ERROR MESSAGE NEEDED") -> None: ...

class SerializeError(Exception): ...
