from typing import Any

from .. import namespaces
from .easy_install import easy_install

class develop(namespaces.DevelopInstaller, easy_install):
    description: str
    user_options: Any
    boolean_options: Any
    command_consumes_arguments: bool
    multi_version: bool
    def run(self) -> None: ...  # type: ignore[override]
    uninstall: Any
    egg_path: Any
    setup_path: Any
    always_copy_from: str
    def initialize_options(self) -> None: ...
    args: Any
    egg_link: Any
    egg_base: Any
    dist: Any
    def finalize_options(self) -> None: ...
    def install_for_development(self) -> None: ...
    def uninstall_link(self) -> None: ...
    def install_egg_scripts(self, dist): ...
    def install_wrapper_scripts(self, dist): ...

class VersionlessRequirement:
    def __init__(self, dist) -> None: ...
    def __getattr__(self, name: str): ...
    def as_requirement(self): ...
