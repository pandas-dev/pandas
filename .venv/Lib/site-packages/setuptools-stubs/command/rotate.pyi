from _typeshed import Incomplete
from typing import ClassVar

from .. import Command

class rotate(Command):
    description: str
    user_options: ClassVar[list[tuple[str, str, str]]]
    boolean_options: ClassVar[list[str]]
    match: Incomplete
    dist_dir: Incomplete
    keep: Incomplete
    def initialize_options(self) -> None: ...
    def finalize_options(self) -> None: ...
    def run(self) -> None: ...
