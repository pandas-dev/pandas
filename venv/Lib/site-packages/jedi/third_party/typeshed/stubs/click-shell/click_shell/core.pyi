from collections.abc import Callable
from logging import Logger
from typing import Any

import click

from ._cmd import ClickCmd

logger: Logger

def get_invoke(command: click.Command) -> Callable[[ClickCmd, str], bool]: ...
def get_help(command: click.Command) -> Callable[[ClickCmd], None]: ...
def get_complete(command: click.Command) -> Callable[[ClickCmd, str, str, int, int], list[str]]: ...

class ClickShell(ClickCmd):
    def add_command(self, cmd: click.Command, name: str) -> None: ...

def make_click_shell(
    ctx: click.Context,
    prompt: str | Callable[[], str] | Callable[[click.Context, str], str] | None = None,
    intro: str | None = None,
    hist_file: str | None = None,
) -> ClickShell: ...

class Shell(click.Group):
    def __init__(
        self,
        prompt: str | Callable[[], str] | Callable[[click.Context, str], str] | None = None,
        intro: str | None = None,
        hist_file: str | None = None,
        on_finished: Callable[[click.Context], None] | None = None,
        **attrs: Any,
    ) -> None: ...
    def add_command(self, cmd: click.Command, name: str | None = None) -> None: ...
    def invoke(self, ctx: click.Context) -> Any: ...
