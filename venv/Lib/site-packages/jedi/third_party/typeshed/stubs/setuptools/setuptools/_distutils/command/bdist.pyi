from _typeshed import Unused
from collections.abc import Callable
from typing import ClassVar
from typing_extensions import deprecated

from ..cmd import Command

def show_formats() -> None: ...

class ListCompat(dict[str, tuple[str, str]]):
    @deprecated("format_commands is now a dict. append is deprecated")
    def append(self, item: Unused) -> None: ...

class bdist(Command):
    description: ClassVar[str]
    user_options: ClassVar[list[tuple[str, str | None, str]]]
    boolean_options: ClassVar[list[str]]
    help_options: ClassVar[list[tuple[str, str | None, str, Callable[[], Unused]]]]
    no_format_option: ClassVar[tuple[str, ...]]
    default_format: ClassVar[dict[str, str]]
    format_commands: ClassVar[ListCompat]
    format_command = format_commands

    def initialize_options(self) -> None: ...
    def finalize_options(self) -> None: ...
    def run(self) -> None: ...
