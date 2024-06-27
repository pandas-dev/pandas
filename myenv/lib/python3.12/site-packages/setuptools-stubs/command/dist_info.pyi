from typing import Any, ClassVar

from .._distutils.cmd import Command

class dist_info(Command):
    description: str
    user_options: Any
    boolean_options: ClassVar[list[str]]
    negative_opt: ClassVar[dict[str, str]]
    egg_base: Any
    def initialize_options(self) -> None: ...
    def finalize_options(self) -> None: ...
    def run(self) -> None: ...
