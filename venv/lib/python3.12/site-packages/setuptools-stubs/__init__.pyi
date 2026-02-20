from _typeshed import Incomplete, StrPath
from abc import abstractmethod
from collections.abc import ItemsView, Iterable, Mapping, Sequence
from typing import Any, Literal, Protocol, TypedDict, TypeVar, overload, type_check_only
from typing_extensions import Never, NotRequired

from ._distutils.cmd import Command as _Command
from ._distutils.dist import Distribution as _Distribution
from ._distutils.extension import Extension as _Extension
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
_DistributionT = TypeVar("_DistributionT", bound=_Distribution, default=Distribution)
_KT = TypeVar("_KT")
_VT_co = TypeVar("_VT_co", covariant=True)

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

# We need any Command subclass to be valid
# Any: pyright would accept using covariance in __setitem__, but mypy won't let a dict be assignable to this protocol
# This is unsound, but it's a quirk of setuptools' internals
@type_check_only
class _DictLike(Protocol[_KT, _VT_co]):
    # See note about using _VT_co instead of Any
    def get(self, key: _KT, default: Any | None = None, /) -> _VT_co | None: ...
    def items(self) -> ItemsView[_KT, _VT_co]: ...
    def keys(self) -> Iterable[_KT]: ...
    def __getitem__(self, key: _KT, /) -> _VT_co: ...
    def __contains__(self, x: object, /) -> bool: ...

@type_check_only
class _MutableDictLike(_DictLike[_KT, _VT_co], Protocol):
    # See note about using _VT_co instead of Any
    def __setitem__(self, key: _KT, value: Any, /) -> None: ...
    def setdefault(self, key: _KT, default: Any, /) -> _VT_co: ...

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
    # Attributes from distutils.dist.DistributionMetadata.set_*
    # These take priority over attributes from distutils.dist.DistributionMetadata.__init__
    keywords: str | Iterable[str] = ...,
    platforms: str | Iterable[str] = ...,
    classifiers: str | Iterable[str] = ...,
    requires: Iterable[str] = ...,
    provides: Iterable[str] = ...,
    obsoletes: Iterable[str] = ...,
    # Attributes from distutils.dist.DistributionMetadata.__init__
    # These take priority over attributes from distutils.dist.Distribution.__init__
    name: str | None = None,
    version: str | None = None,
    author: str | None = None,
    author_email: str | None = None,
    maintainer: str | None = None,
    maintainer_email: str | None = None,
    url: str | None = None,
    license: str | None = None,
    description: str | None = None,
    long_description: str | None = None,
    download_url: str | None = None,
    # Attributes from distutils.dist.Distribution.__init__ (except self.metadata)
    # These take priority over attributes from distutils.dist.Distribution.display_option_names
    verbose: bool = True,
    help: bool = False,
    cmdclass: _MutableDictLike[str, type[_Command]] = {},
    command_packages: str | list[str] | None = None,
    script_name: StrPath | None = ...,  # default is actually set in distutils.core.setup
    script_args: list[str] | None = ...,  # default is actually set in distutils.core.setup
    command_options: _MutableDictLike[str, _DictLike[str, tuple[str, str]]] = {},
    packages: list[str] | None = None,
    package_dir: Mapping[str, str] | None = None,
    py_modules: list[str] | None = None,
    libraries: list[tuple[str, _BuildInfo]] | None = None,
    headers: list[str] | None = None,
    ext_modules: Sequence[_Extension] | None = None,
    ext_package: str | None = None,
    include_dirs: list[str] | None = None,
    extra_path: Never = ...,  # Deprecated
    scripts: list[str] | None = None,
    data_files: list[tuple[str, Sequence[str]]] | None = None,
    password: str = "",
    command_obj: _MutableDictLike[str, _Command] = {},
    have_run: _MutableDictLike[str, bool] = {},
    # kwargs used directly in distutils.dist.Distribution.__init__
    options: Mapping[str, Mapping[str, str]] | None = None,
    licence: Never = ...,  # Deprecated
    # Attributes from distutils.dist.Distribution.display_option_names
    # (this can more easily be copied from the `if TYPE_CHECKING` block)
    help_commands: bool = False,
    fullname: str | Literal[False] = False,
    contact: str | Literal[False] = False,
    contact_email: str | Literal[False] = False,
    # kwargs used directly in setuptools.dist.Distribution.__init__
    # and attributes from setuptools.dist.Distribution.__init__
    package_data: _DictLike[str, list[str]] = {},
    dist_files: list[tuple[str, str, str]] = [],
    include_package_data: bool | None = None,
    exclude_package_data: _DictLike[str, list[str]] | None = None,
    src_root: str | None = None,
    dependency_links: list[str] = [],
    setup_requires: list[str] = [],
    # From Distribution._DISTUTILS_UNSUPPORTED_METADATA set in Distribution._set_metadata_defaults
    long_description_content_type: str | None = None,
    project_urls: _DictLike[Incomplete, Incomplete] = {},
    provides_extras: _MutableDictLike[Incomplete, Incomplete] = {},
    license_expression: str | None = None,
    license_file: Never = ...,  # Deprecated
    license_files: Iterable[str] | None = None,
    install_requires: str | Iterable[str] = [],
    extras_require: _DictLike[Incomplete, Incomplete] = {},
    # kwargs used directly in distutils.core.setup
    distclass: type[_DistributionT] = Distribution,  # type: ignore[assignment] # noqa: Y011
    # Custom Distributions could accept more params
    **attrs: Any,
) -> _DistributionT: ...

class Command(_Command):
    command_consumes_arguments: bool
    distribution: Distribution
    dry_run: bool
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
