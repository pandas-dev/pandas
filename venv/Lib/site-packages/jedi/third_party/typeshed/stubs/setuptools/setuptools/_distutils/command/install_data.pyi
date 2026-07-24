from _typeshed import Incomplete
from typing import ClassVar

from ..cmd import Command

class install_data(Command):
    description: ClassVar[str]
    user_options: ClassVar[list[tuple[str, str | None, str]]]
    boolean_options: Incomplete
    install_dir: Incomplete
    outfiles: Incomplete
    root: Incomplete
    force: bool
    data_files: Incomplete
    warn_dir: bool
    def initialize_options(self) -> None: ...
    def finalize_options(self) -> None: ...
    def run(self) -> None: ...
    def get_inputs(self): ...
    def get_outputs(self): ...
