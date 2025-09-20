from _typeshed import Incomplete
from typing import ClassVar

from .. import Command, namespaces

class install_egg_info(namespaces.Installer, Command):
    description: str
    user_options: ClassVar[list[tuple[str, str, str]]]
    install_dir: Incomplete
    def initialize_options(self) -> None: ...
    source: Incomplete
    target: str
    outputs: list[str]
    def finalize_options(self) -> None: ...
    def run(self) -> None: ...
    def get_outputs(self): ...
    def copytree(self): ...
