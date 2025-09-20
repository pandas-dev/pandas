import configparser
import sys
from collections.abc import Mapping, Sequence

class Parser(configparser.RawConfigParser):
    # __init__ signature was taken from RawConfigParser, but with no allow_no_value argument and all arguments are keyword-only
    if sys.version_info >= (3, 13):
        def __init__(
            self,
            *,
            defaults: Mapping[str, str] | None = None,
            dict_type: type[Mapping[str, str]] = ...,
            delimiters: Sequence[str] = ("=", ":"),
            comment_prefixes: Sequence[str] = ("#", ";"),
            inline_comment_prefixes: Sequence[str] | None = None,
            strict: bool = True,
            empty_lines_in_values: bool = True,
            default_section: str = "DEFAULT",
            interpolation: configparser.Interpolation | None = ...,
            converters: configparser._ConvertersMap = ...,
            allow_unnamed_section: bool = False,
        ) -> None: ...
    else:
        def __init__(
            self,
            *,
            defaults: Mapping[str, str] | None = None,
            dict_type: type[Mapping[str, str]] = ...,
            delimiters: Sequence[str] = ("=", ":"),
            comment_prefixes: Sequence[str] = ("#", ";"),
            inline_comment_prefixes: Sequence[str] | None = None,
            strict: bool = True,
            empty_lines_in_values: bool = True,
            default_section: str = "DEFAULT",
            interpolation: configparser.Interpolation | None = ...,
            converters: configparser._ConvertersMap = ...,
        ) -> None: ...

    def optionxform(self, key: str) -> str: ...
    def get(self, section: configparser._SectionName, option: str) -> str: ...  # type: ignore[override]
