from _typeshed import Incomplete
from abc import abstractmethod
from collections.abc import Mapping, Sequence
from typing import Any, Literal, TypedDict, TypeVar, overload, type_check_only
from typing_extensions import NotRequired

from ._distutils.cmd import Command as _Command
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
from .depends import Require as Require
from .discovery import _Finder
from .dist import Distribution as Distribution
from .extension import Extension as Extension
from .warnings import SetuptoolsDeprecationWarning as SetuptoolsDeprecationWarning

_CommandT = TypeVar("_CommandT", bound=_Command)

__all__ = [
    "setup",
    "Distribution",
    "Command",
    "Extension",
    "Require",
    "SetuptoolsDeprecationWarning",
    "find_packages",
    "find_namespace_packages",
]

__version__: str

@type_check_only
class _BuildInfo(TypedDict):
    sources: list[str] | tuple[str, ...]
    obj_deps: NotRequired[dict[str, list[str] | tuple[str, ...]]]
    macros: NotRequired[list[tuple[str] | tuple[str, str | None]]]
    include_dirs: NotRequired[list[str]]
    cflags: NotRequired[list[str]]

find_packages = _Finder.find
find_namespace_packages = _Finder.find

def setup(
    *,
    name: str = ...,
    version: str = ...,
    description: str = ...,
    long_description: str = ...,
    long_description_content_type: str = ...,
    author: str = ...,
    author_email: str = ...,
    maintainer: str = ...,
    maintainer_email: str = ...,
    url: str = ...,
    download_url: str = ...,
    packages: list[str] = ...,
    py_modules: list[str] = ...,
    scripts: list[str] = ...,
    ext_modules: Sequence[Extension] = ...,
    classifiers: list[str] = ...,
    distclass: type[Distribution] = ...,
    script_name: str = ...,
    script_args: list[str] = ...,
    options: Mapping[str, Incomplete] = ...,
    license: str = ...,
    keywords: list[str] | str = ...,
    platforms: list[str] | str = ...,
    cmdclass: Mapping[str, type[_Command]] = ...,
    data_files: list[tuple[str, list[str]]] = ...,
    package_dir: Mapping[str, str] = ...,
    obsoletes: list[str] = ...,
    provides: list[str] = ...,
    requires: list[str] = ...,
    command_packages: list[str] = ...,
    command_options: Mapping[str, Mapping[str, tuple[Incomplete, Incomplete]]] = ...,
    package_data: Mapping[str, list[str]] = ...,
    include_package_data: bool = ...,
    # libraries for `Distribution` or `build_clib`, not `Extension`, `build_ext` or `CCompiler`
    libraries: list[tuple[str, _BuildInfo]] = ...,
    headers: list[str] = ...,
    ext_package: str = ...,
    include_dirs: list[str] = ...,
    password: str = ...,
    fullname: str = ...,
    # Custom Distributions could accept more params
    **attrs: Any,
) -> Distribution: ...

class Command(_Command):
    command_consumes_arguments: bool
    distribution: Distribution
    # Any: Dynamic command subclass attributes
    def __init__(self, dist: Distribution, **kw: Any) -> None: ...
    # Note: Commands that setuptools doesn't re-expose are considered deprecated (they must be imported from distutils directly)
    # So we're not listing them here. This list comes directly from the setuptools/command folder. Minus the test command.
    @overload  # type: ignore[override]
    def get_finalized_command(self, command: Literal["alias"], create: bool | Literal[0, 1] = 1) -> alias: ...
    @overload
    def get_finalized_command(self, command: Literal["bdist_egg"], create: bool | Literal[0, 1] = 1) -> bdist_egg: ...
    @overload
    def get_finalized_command(self, command: Literal["bdist_rpm"], create: bool | Literal[0, 1] = 1) -> bdist_rpm: ...  # type: ignore[overload-overlap]
    @overload
    def get_finalized_command(self, command: Literal["bdist_wheel"], create: bool | Literal[0, 1] = 1) -> bdist_wheel: ...
    @overload
    def get_finalized_command(self, command: Literal["build"], create: bool | Literal[0, 1] = 1) -> build: ...  # type: ignore[overload-overlap]
    @overload
    def get_finalized_command(self, command: Literal["build_clib"], create: bool | Literal[0, 1] = 1) -> build_clib: ...  # type: ignore[overload-overlap]
    @overload
    def get_finalized_command(self, command: Literal["build_ext"], create: bool | Literal[0, 1] = 1) -> build_ext: ...  # type: ignore[overload-overlap]
    @overload
    def get_finalized_command(self, command: Literal["build_py"], create: bool | Literal[0, 1] = 1) -> build_py: ...  # type: ignore[overload-overlap]
    @overload
    def get_finalized_command(self, command: Literal["develop"], create: bool | Literal[0, 1] = 1) -> develop: ...
    @overload
    def get_finalized_command(self, command: Literal["dist_info"], create: bool | Literal[0, 1] = 1) -> dist_info: ...  # type: ignore[overload-overlap]
    @overload
    def get_finalized_command(self, command: Literal["easy_install"], create: bool | Literal[0, 1] = 1) -> easy_install: ...
    @overload
    def get_finalized_command(self, command: Literal["editable_wheel"], create: bool | Literal[0, 1] = 1) -> editable_wheel: ...
    @overload
    def get_finalized_command(self, command: Literal["egg_info"], create: bool | Literal[0, 1] = 1) -> egg_info: ...
    @overload
    def get_finalized_command(self, command: Literal["install"], create: bool | Literal[0, 1] = 1) -> install: ...  # type: ignore[overload-overlap]
    @overload
    def get_finalized_command(
        self, command: Literal["install_egg_info"], create: bool | Literal[0, 1] = 1
    ) -> install_egg_info: ...
    @overload
    def get_finalized_command(self, command: Literal["install_lib"], create: bool | Literal[0, 1] = 1) -> install_lib: ...  # type: ignore[overload-overlap]
    @overload
    def get_finalized_command(self, command: Literal["install_scripts"], create: bool | Literal[0, 1] = 1) -> install_scripts: ...  # type: ignore[overload-overlap]
    @overload
    def get_finalized_command(self, command: Literal["rotate"], create: bool | Literal[0, 1] = 1) -> rotate: ...
    @overload
    def get_finalized_command(self, command: Literal["saveopts"], create: bool | Literal[0, 1] = 1) -> saveopts: ...
    @overload
    def get_finalized_command(self, command: Literal["sdist"], create: bool | Literal[0, 1] = 1) -> sdist: ...  # type: ignore[overload-overlap]
    @overload
    def get_finalized_command(self, command: Literal["setopt"], create: bool | Literal[0, 1] = 1) -> setopt: ...
    @overload
    def get_finalized_command(self, command: str, create: bool | Literal[0, 1] = 1) -> Command: ...
    @overload  # type: ignore[override] # Extra **kw param
    def reinitialize_command(self, command: Literal["alias"], reinit_subcommands: bool = False, **kw) -> alias: ...
    @overload
    def reinitialize_command(self, command: Literal["bdist_egg"], reinit_subcommands: bool = False, **kw) -> bdist_egg: ...
    @overload
    def reinitialize_command(self, command: Literal["bdist_rpm"], reinit_subcommands: bool = False, **kw) -> bdist_rpm: ...
    @overload
    def reinitialize_command(self, command: Literal["bdist_wheel"], reinit_subcommands: bool = False, **kw) -> bdist_wheel: ...
    @overload
    def reinitialize_command(self, command: Literal["build"], reinit_subcommands: bool = False, **kw) -> build: ...
    @overload
    def reinitialize_command(self, command: Literal["build_clib"], reinit_subcommands: bool = False, **kw) -> build_clib: ...
    @overload
    def reinitialize_command(self, command: Literal["build_ext"], reinit_subcommands: bool = False, **kw) -> build_ext: ...
    @overload
    def reinitialize_command(self, command: Literal["build_py"], reinit_subcommands: bool = False, **kw) -> build_py: ...
    @overload
    def reinitialize_command(self, command: Literal["develop"], reinit_subcommands: bool = False, **kw) -> develop: ...
    @overload
    def reinitialize_command(self, command: Literal["dist_info"], reinit_subcommands: bool = False, **kw) -> dist_info: ...
    @overload
    def reinitialize_command(self, command: Literal["easy_install"], reinit_subcommands: bool = False, **kw) -> easy_install: ...
    @overload
    def reinitialize_command(
        self, command: Literal["editable_wheel"], reinit_subcommands: bool = False, **kw
    ) -> editable_wheel: ...
    @overload
    def reinitialize_command(self, command: Literal["egg_info"], reinit_subcommands: bool = False, **kw) -> egg_info: ...
    @overload
    def reinitialize_command(self, command: Literal["install"], reinit_subcommands: bool = False, **kw) -> install: ...
    @overload
    def reinitialize_command(
        self, command: Literal["install_egg_info"], reinit_subcommands: bool = False, **kw
    ) -> install_egg_info: ...
    @overload
    def reinitialize_command(self, command: Literal["install_lib"], reinit_subcommands: bool = False, **kw) -> install_lib: ...
    @overload
    def reinitialize_command(
        self, command: Literal["install_scripts"], reinit_subcommands: bool = False, **kw
    ) -> install_scripts: ...
    @overload
    def reinitialize_command(self, command: Literal["rotate"], reinit_subcommands: bool = False, **kw) -> rotate: ...
    @overload
    def reinitialize_command(self, command: Literal["saveopts"], reinit_subcommands: bool = False, **kw) -> saveopts: ...
    @overload
    def reinitialize_command(self, command: Literal["sdist"], reinit_subcommands: bool = False, **kw) -> sdist: ...
    @overload
    def reinitialize_command(self, command: Literal["setopt"], reinit_subcommands: bool = False, **kw) -> setopt: ...
    @overload
    def reinitialize_command(self, command: str, reinit_subcommands: bool = False, **kw) -> Command: ...
    @overload
    def reinitialize_command(self, command: _CommandT, reinit_subcommands: bool = False, **kw) -> _CommandT: ...
    @abstractmethod
    def initialize_options(self) -> None: ...
    @abstractmethod
    def finalize_options(self) -> None: ...
    @abstractmethod
    def run(self) -> None: ...

class sic(str): ...
