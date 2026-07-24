from collections.abc import Callable, MutableMapping, Sequence
from typing import Any, Final, Literal, overload
from typing_extensions import deprecated

import click

__all__ = ["DefaultGroup"]
__version__: Final[str]

class DefaultGroup(click.Group):
    ignore_unknown_options: bool
    default_cmd_name: str | None
    default_if_no_args: bool
    # type hints were taken from click lib
    def __init__(
        self,
        name: str | None = None,
        commands: MutableMapping[str, click.Command] | Sequence[click.Command] | None = None,
        *,
        ignore_unknown_options: Literal[True] | None = True,
        default: str | None = None,
        default_if_no_args: bool = False,
        invoke_without_command: bool = False,
        no_args_is_help: bool | None = None,
        subcommand_metavar: str | None = None,
        chain: bool = False,
        result_callback: Callable[..., Any] | None = None,  # Any is specified in click lib
        context_settings: MutableMapping[str, Any] | None = None,  # Any is specified in click lib
        callback: Callable[..., Any] | None = None,  # Any is specified in click lib
        params: list[click.Parameter] | None = None,
        help: str | None = None,
        epilog: str | None = None,
        short_help: str | None = None,
        options_metavar: str | None = "[OPTIONS]",
        add_help_option: bool = True,
        hidden: bool = False,
        deprecated: bool = False,
    ) -> None: ...
    def set_default_command(self, command: click.Command) -> None: ...
    def parse_args(self, ctx: click.Context, args: list[str]) -> list[str]: ...
    def get_command(self, ctx: click.Context, cmd_name: str) -> click.Command | None: ...
    def resolve_command(self, ctx: click.Context, args: list[str]) -> tuple[str | None, click.Command | None, list[str]]: ...
    def format_commands(self, ctx: click.Context, formatter: click.HelpFormatter) -> None: ...
    @overload
    def command(
        self,
        __func: Callable[..., Any],
        /,
        *,
        name: str | None = ...,
        cls: type[click.Command] | None = ...,
        default: Literal[False] = False,
    ) -> click.Command: ...
    @overload
    @deprecated("Use default param of `DefaultGroup` or `set_default_command()` instead")
    def command(
        self,
        __func: Callable[..., Any],
        /,
        *,
        name: str | None = ...,
        cls: type[click.Command] | None = ...,
        default: Literal[True],
    ) -> click.Command: ...
    @overload
    def command(
        self, *, name: str | None = ..., cls: type[click.Command] | None = ..., default: Literal[False] = False
    ) -> Callable[[Callable[..., Any]], click.Command]: ...
    @overload
    @deprecated("Use default param of `DefaultGroup` or `set_default_command()` instead")
    def command(
        self, *, name: str | None = ..., cls: type[click.Command] | None = ..., default: Literal[True]
    ) -> Callable[[Callable[..., Any]], click.Command]: ...
    @overload
    def command(self, *args: Any, **kwargs: Any) -> Callable[[Callable[..., Any]], click.Command] | click.Command: ...

class DefaultCommandFormatter:
    group: click.Group
    formatter: click.HelpFormatter
    mark: str
    def __init__(self, group: click.Group, formatter: click.HelpFormatter, mark: str = "*") -> None: ...
    def write_dl(self, rows: Sequence[tuple[str, str]], col_max: int = 30, col_spacing: int = -2) -> None: ...
    def __getattr__(self, attr: str) -> Any: ...  # attribute access is forwarded to click.HelpFormatter
