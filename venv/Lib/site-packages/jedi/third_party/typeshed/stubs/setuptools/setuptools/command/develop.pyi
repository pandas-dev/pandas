from _typeshed import Incomplete
from typing import ClassVar
from typing_extensions import deprecated

from setuptools import Command
from setuptools.warnings import SetuptoolsDeprecationWarning

class develop(Command):
    user_options: ClassVar[list[tuple[str, str | None, str]]]
    boolean_options: ClassVar[list[str]]
    install_dir: Incomplete
    no_deps: bool
    user: bool
    prefix: Incomplete
    index_url: Incomplete
    def run(self) -> None: ...
    @deprecated(
        "develop command is deprecated. Please avoid running `setup.py` and `develop`. "
        "Instead, use standards-based tools like pip or uv."
    )
    def initialize_options(self) -> None: ...
    def finalize_options(self) -> None: ...

class DevelopDeprecationWarning(SetuptoolsDeprecationWarning): ...
