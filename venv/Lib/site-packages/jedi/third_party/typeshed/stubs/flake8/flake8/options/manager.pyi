import argparse
from _typeshed import Incomplete
from collections.abc import Callable, Sequence
from enum import Enum
from logging import Logger
from typing import Any

from ..plugins.finder import Plugins

LOG: Logger

class _ARG(Enum):
    NO = 1

class Option:
    short_option_name: Incomplete
    long_option_name: Incomplete
    option_args: Incomplete
    action: Incomplete
    default: Incomplete
    type: Incomplete
    dest: Incomplete
    nargs: Incomplete
    const: Incomplete
    choices: Incomplete
    help: Incomplete
    metavar: Incomplete
    required: Incomplete
    option_kwargs: Incomplete
    parse_from_config: Incomplete
    comma_separated_list: Incomplete
    normalize_paths: Incomplete
    config_name: Incomplete
    def __init__(
        self,
        short_option_name: str | _ARG = ...,
        long_option_name: str | _ARG = ...,
        action: str | type[argparse.Action] | _ARG = ...,
        default: Any | _ARG = ...,
        type: Callable[..., Any] | _ARG = ...,
        dest: str | _ARG = ...,
        nargs: int | str | _ARG = ...,
        const: Any | _ARG = ...,
        choices: Sequence[Any] | _ARG = ...,
        help: str | _ARG = ...,
        metavar: str | _ARG = ...,
        required: bool | _ARG = ...,
        parse_from_config: bool = False,
        comma_separated_list: bool = False,
        normalize_paths: bool = False,
    ) -> None: ...
    @property
    def filtered_option_kwargs(self) -> dict[str, Any]: ...
    def normalize(self, value: Any, *normalize_args: str) -> Any: ...
    def to_argparse(self) -> tuple[list[str], dict[str, Any]]: ...

class OptionManager:
    formatter_names: Incomplete
    parser: Incomplete
    config_options_dict: Incomplete
    options: Incomplete
    extended_default_ignore: Incomplete
    extended_default_select: Incomplete
    def __init__(
        self, *, version: str, plugin_versions: str, parents: list[argparse.ArgumentParser], formatter_names: list[str]
    ) -> None: ...
    def register_plugins(self, plugins: Plugins) -> None: ...
    def add_option(self, *args: Any, **kwargs: Any) -> None: ...
    def extend_default_ignore(self, error_codes: Sequence[str]) -> None: ...
    def extend_default_select(self, error_codes: Sequence[str]) -> None: ...
    def parse_args(self, args: Sequence[str] | None = None, values: argparse.Namespace | None = None) -> argparse.Namespace: ...
