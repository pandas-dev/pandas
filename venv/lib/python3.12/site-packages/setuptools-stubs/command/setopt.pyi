from _typeshed import Incomplete
from abc import abstractmethod
from typing import ClassVar

from .. import Command

__all__ = ["config_file", "edit_config", "option_base", "setopt"]

def config_file(kind: str = "local"): ...
def edit_config(filename, settings, dry_run: bool = False) -> None: ...

class option_base(Command):
    user_options: ClassVar[list[tuple[str, str, str]]]
    boolean_options: ClassVar[list[str]]
    global_config: Incomplete
    user_config: Incomplete
    filename: Incomplete
    def initialize_options(self) -> None: ...
    def finalize_options(self) -> None: ...
    @abstractmethod
    def run(self) -> None: ...

class setopt(option_base):
    description: str
    user_options: ClassVar[list[tuple[str, str, str]]]
    boolean_options: ClassVar[list[str]]
    command: Incomplete
    option: Incomplete
    set_value: Incomplete
    remove: Incomplete
    def initialize_options(self) -> None: ...
    def finalize_options(self) -> None: ...
    def run(self) -> None: ...
