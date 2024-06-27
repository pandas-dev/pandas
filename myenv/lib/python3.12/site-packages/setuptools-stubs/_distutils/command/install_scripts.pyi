from _typeshed import Incomplete

from ..cmd import Command

class install_scripts(Command):
    description: str
    user_options: Incomplete
    boolean_options: Incomplete
    install_dir: Incomplete
    force: int
    build_dir: Incomplete
    skip_build: Incomplete
    def initialize_options(self) -> None: ...
    def finalize_options(self) -> None: ...
    outfiles: list[str]
    def run(self) -> None: ...
    def get_inputs(self): ...
    def get_outputs(self): ...
