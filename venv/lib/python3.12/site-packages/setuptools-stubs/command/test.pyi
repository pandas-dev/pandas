from typing import ClassVar, NoReturn
from typing_extensions import deprecated

from .. import Command

@deprecated(
    """\
The test command is disabled and references to it are deprecated. \
Please remove any references to `setuptools.command.test` in all supported versions of the affected package.\
"""
)
class test(Command):
    description: ClassVar[str]
    user_options: ClassVar[list[tuple[str, str, str]]]
    def initialize_options(self) -> None: ...
    def finalize_options(self) -> None: ...
    def run(self) -> NoReturn: ...
