from _typeshed import Incomplete
from typing import ClassVar

from pkg_resources import Distribution

from .. import namespaces
from .easy_install import easy_install

class develop(namespaces.DevelopInstaller, easy_install):
    description: str
    user_options: ClassVar[list[tuple[str, str | None, str]]]
    boolean_options: ClassVar[list[str]]
    command_consumes_arguments: bool
    multi_version: bool
    def run(self) -> None: ...  # type: ignore[override]
    uninstall: Incomplete
    egg_path: Incomplete
    setup_path: Incomplete
    always_copy_from: str
    def initialize_options(self) -> None: ...
    args: list[Incomplete]
    egg_link: str
    egg_base: Incomplete
    dist: Distribution
    def finalize_options(self) -> None: ...
    def install_for_development(self) -> None: ...
    def uninstall_link(self) -> None: ...
    def install_egg_scripts(self, dist): ...
    def install_wrapper_scripts(self, dist): ...

class VersionlessRequirement:
    def __init__(self, dist) -> None: ...
    def __getattr__(self, name: str): ...
    def as_requirement(self): ...
