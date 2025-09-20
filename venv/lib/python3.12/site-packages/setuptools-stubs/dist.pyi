from _typeshed import Incomplete, StrPath
from collections.abc import Iterable, Iterator, MutableMapping
from importlib import metadata
from typing import Literal, TypeVar, overload

from . import Command, SetuptoolsDeprecationWarning
from ._distutils.cmd import Command as _Command
from ._distutils.dist import Distribution as _Distribution
from .command.alias import alias
from .command.bdist_egg import bdist_egg
from .command.bdist_rpm import bdist_rpm
from .command.bdist_wheel import bdist_wheel
from .command.build import build
from .command.build_clib import build_clib
from .command.build_ext import build_ext
from .command.build_py import build_py
from .command.develop import develop
from .command.dist_info import dist_info
from .command.easy_install import easy_install
from .command.editable_wheel import editable_wheel
from .command.egg_info import egg_info
from .command.install import install
from .command.install_egg_info import install_egg_info
from .command.install_lib import install_lib
from .command.install_scripts import install_scripts
from .command.rotate import rotate
from .command.saveopts import saveopts
from .command.sdist import sdist
from .command.setopt import setopt

_CommandT = TypeVar("_CommandT", bound=_Command)

__all__ = ["Distribution"]

class Distribution(_Distribution):
    include_package_data: bool | None
    exclude_package_data: dict[str, list[str]] | None
    src_root: str | None
    dependency_links: list[str]
    setup_requires: list[str]
    def __init__(self, attrs: MutableMapping[str, Incomplete] | None = None) -> None: ...
    def parse_config_files(self, filenames: Iterable[StrPath] | None = None, ignore_option_errors: bool = False) -> None: ...
    def fetch_build_eggs(self, requires: str | Iterable[str]) -> list[metadata.Distribution]: ...
    def get_egg_cache_dir(self) -> str: ...
    def fetch_build_egg(self, req): ...
    # NOTE: Commands that setuptools doesn't re-expose are considered deprecated (they must be imported from distutils directly)
    # So we're not listing them here. This list comes directly from the setuptools/command folder. Minus the test command.
    @overload  # type: ignore[override]
    def get_command_obj(self, command: Literal["alias"], create: Literal[1, True] = 1) -> alias: ...
    @overload
    def get_command_obj(self, command: Literal["bdist_egg"], create: Literal[1, True] = 1) -> bdist_egg: ...
    @overload
    def get_command_obj(self, command: Literal["bdist_rpm"], create: Literal[1, True] = 1) -> bdist_rpm: ...  # type: ignore[overload-overlap]
    @overload
    def get_command_obj(self, command: Literal["bdist_wheel"], create: Literal[1, True] = 1) -> bdist_wheel: ...
    @overload
    def get_command_obj(self, command: Literal["build"], create: Literal[1, True] = 1) -> build: ...  # type: ignore[overload-overlap]
    @overload
    def get_command_obj(self, command: Literal["build_clib"], create: Literal[1, True] = 1) -> build_clib: ...  # type: ignore[overload-overlap]
    @overload
    def get_command_obj(self, command: Literal["build_ext"], create: Literal[1, True] = 1) -> build_ext: ...  # type: ignore[overload-overlap]
    @overload
    def get_command_obj(self, command: Literal["build_py"], create: Literal[1, True] = 1) -> build_py: ...  # type: ignore[overload-overlap]
    @overload
    def get_command_obj(self, command: Literal["develop"], create: Literal[1, True] = 1) -> develop: ...
    @overload
    def get_command_obj(self, command: Literal["dist_info"], create: Literal[1, True] = 1) -> dist_info: ...  # type: ignore[overload-overlap]
    @overload
    def get_command_obj(self, command: Literal["easy_install"], create: Literal[1, True] = 1) -> easy_install: ...
    @overload
    def get_command_obj(self, command: Literal["editable_wheel"], create: Literal[1, True] = 1) -> editable_wheel: ...
    @overload
    def get_command_obj(self, command: Literal["egg_info"], create: Literal[1, True] = 1) -> egg_info: ...
    @overload
    def get_command_obj(self, command: Literal["install"], create: Literal[1, True] = 1) -> install: ...  # type: ignore[overload-overlap]
    @overload
    def get_command_obj(self, command: Literal["install_egg_info"], create: Literal[1, True] = 1) -> install_egg_info: ...
    @overload
    def get_command_obj(self, command: Literal["install_lib"], create: Literal[1, True] = 1) -> install_lib: ...  # type: ignore[overload-overlap]
    @overload
    def get_command_obj(self, command: Literal["install_scripts"], create: Literal[1, True] = 1) -> install_scripts: ...  # type: ignore[overload-overlap]
    @overload
    def get_command_obj(self, command: Literal["rotate"], create: Literal[1, True] = 1) -> rotate: ...
    @overload
    def get_command_obj(self, command: Literal["saveopts"], create: Literal[1, True] = 1) -> saveopts: ...
    @overload
    def get_command_obj(self, command: Literal["sdist"], create: Literal[1, True] = 1) -> sdist: ...  # type: ignore[overload-overlap]
    @overload
    def get_command_obj(self, command: Literal["setopt"], create: Literal[1, True] = 1) -> setopt: ...
    @overload
    def get_command_obj(self, command: str, create: Literal[1, True] = 1) -> Command: ...
    # Not replicating the overloads for "Command | None", user may use "isinstance"
    @overload
    def get_command_obj(self, command: str, create: Literal[0, False]) -> Command | None: ...
    @overload
    def get_command_class(self, command: Literal["alias"]) -> type[alias]: ...
    @overload
    def get_command_class(self, command: Literal["bdist_egg"]) -> type[bdist_egg]: ...
    @overload
    def get_command_class(self, command: Literal["bdist_rpm"]) -> type[bdist_rpm]: ...  # type: ignore[overload-overlap]
    @overload
    def get_command_class(self, command: Literal["bdist_wheel"]) -> type[bdist_wheel]: ...
    @overload
    def get_command_class(self, command: Literal["build"]) -> type[build]: ...  # type: ignore[overload-overlap]
    @overload
    def get_command_class(self, command: Literal["build_clib"]) -> type[build_clib]: ...  # type: ignore[overload-overlap]
    @overload
    def get_command_class(self, command: Literal["build_ext"]) -> type[build_ext]: ...  # type: ignore[overload-overlap]
    @overload
    def get_command_class(self, command: Literal["build_py"]) -> type[build_py]: ...  # type: ignore[overload-overlap]
    @overload
    def get_command_class(self, command: Literal["develop"]) -> type[develop]: ...
    @overload
    def get_command_class(self, command: Literal["dist_info"]) -> type[dist_info]: ...  # type: ignore[overload-overlap]
    @overload
    def get_command_class(self, command: Literal["easy_install"]) -> type[easy_install]: ...
    @overload
    def get_command_class(self, command: Literal["editable_wheel"]) -> type[editable_wheel]: ...
    @overload
    def get_command_class(self, command: Literal["egg_info"]) -> type[egg_info]: ...
    @overload
    def get_command_class(self, command: Literal["install"]) -> type[install]: ...  # type: ignore[overload-overlap]
    @overload
    def get_command_class(self, command: Literal["install_egg_info"]) -> type[install_egg_info]: ...
    @overload
    def get_command_class(self, command: Literal["install_lib"]) -> type[install_lib]: ...  # type: ignore[overload-overlap]
    @overload
    def get_command_class(self, command: Literal["install_scripts"]) -> type[install_scripts]: ...  # type: ignore[overload-overlap]
    @overload
    def get_command_class(self, command: Literal["rotate"]) -> type[rotate]: ...
    @overload
    def get_command_class(self, command: Literal["saveopts"]) -> type[saveopts]: ...
    @overload
    def get_command_class(self, command: Literal["sdist"]) -> type[sdist]: ...  # type: ignore[overload-overlap]
    @overload
    def get_command_class(self, command: Literal["setopt"]) -> type[setopt]: ...
    @overload
    def get_command_class(self, command: str) -> type[Command]: ...
    @overload  # type: ignore[override]
    def reinitialize_command(self, command: Literal["alias"], reinit_subcommands: bool = False) -> alias: ...
    @overload
    def reinitialize_command(self, command: Literal["bdist_egg"], reinit_subcommands: bool = False) -> bdist_egg: ...
    @overload
    def reinitialize_command(self, command: Literal["bdist_rpm"], reinit_subcommands: bool = False) -> bdist_rpm: ...  # type: ignore[overload-overlap]
    @overload
    def reinitialize_command(self, command: Literal["bdist_wheel"], reinit_subcommands: bool = False) -> bdist_wheel: ...
    @overload
    def reinitialize_command(self, command: Literal["build"], reinit_subcommands: bool = False) -> build: ...  # type: ignore[overload-overlap]
    @overload
    def reinitialize_command(self, command: Literal["build_clib"], reinit_subcommands: bool = False) -> build_clib: ...  # type: ignore[overload-overlap]
    @overload
    def reinitialize_command(self, command: Literal["build_ext"], reinit_subcommands: bool = False) -> build_ext: ...  # type: ignore[overload-overlap]
    @overload
    def reinitialize_command(self, command: Literal["build_py"], reinit_subcommands: bool = False) -> build_py: ...  # type: ignore[overload-overlap]
    @overload
    def reinitialize_command(self, command: Literal["develop"], reinit_subcommands: bool = False) -> develop: ...
    @overload
    def reinitialize_command(self, command: Literal["dist_info"], reinit_subcommands: bool = False) -> dist_info: ...  # type: ignore[overload-overlap]
    @overload
    def reinitialize_command(self, command: Literal["easy_install"], reinit_subcommands: bool = False) -> easy_install: ...
    @overload
    def reinitialize_command(self, command: Literal["editable_wheel"], reinit_subcommands: bool = False) -> editable_wheel: ...
    @overload
    def reinitialize_command(self, command: Literal["egg_info"], reinit_subcommands: bool = False) -> egg_info: ...
    @overload
    def reinitialize_command(self, command: Literal["install"], reinit_subcommands: bool = False) -> install: ...  # type: ignore[overload-overlap]
    @overload
    def reinitialize_command(
        self, command: Literal["install_egg_info"], reinit_subcommands: bool = False
    ) -> install_egg_info: ...
    @overload
    def reinitialize_command(self, command: Literal["install_lib"], reinit_subcommands: bool = False) -> install_lib: ...  # type: ignore[overload-overlap]
    @overload
    def reinitialize_command(self, command: Literal["install_scripts"], reinit_subcommands: bool = False) -> install_scripts: ...  # type: ignore[overload-overlap]
    @overload
    def reinitialize_command(self, command: Literal["rotate"], reinit_subcommands: bool = False) -> rotate: ...
    @overload
    def reinitialize_command(self, command: Literal["saveopts"], reinit_subcommands: bool = False) -> saveopts: ...
    @overload
    def reinitialize_command(self, command: Literal["sdist"], reinit_subcommands: bool = False) -> sdist: ...  # type: ignore[overload-overlap]
    @overload
    def reinitialize_command(self, command: Literal["setopt"], reinit_subcommands: bool = False) -> setopt: ...
    @overload
    def reinitialize_command(self, command: str, reinit_subcommands: bool = False) -> Command: ...
    @overload
    def reinitialize_command(self, command: _CommandT, reinit_subcommands: bool = False) -> _CommandT: ...
    def include(self, **attrs) -> None: ...
    def exclude_package(self, package: str) -> None: ...
    def has_contents_for(self, package: str) -> bool: ...
    def exclude(self, **attrs) -> None: ...
    def get_cmdline_options(self) -> dict[str, dict[str, str | None]]: ...
    def iter_distribution_names(self) -> Iterator[str]: ...
    def handle_display_options(self, option_order): ...

class DistDeprecationWarning(SetuptoolsDeprecationWarning): ...
