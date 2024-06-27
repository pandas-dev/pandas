from _typeshed import StrPath
from abc import abstractmethod
from collections.abc import Iterable, Mapping, Sequence
from typing import Any

from ._distutils.cmd import Command as _Command
from .depends import Require as Require
from .dist import Distribution as Distribution
from .extension import Extension as Extension
from .warnings import SetuptoolsDeprecationWarning as SetuptoolsDeprecationWarning

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

# Pytype fails with the following:
# find_packages = PackageFinder.find
# find_namespace_packages = PEP420PackageFinder.find
def find_packages(where: StrPath = ".", exclude: Iterable[str] = (), include: Iterable[str] = ("*",)) -> list[str]: ...
def find_namespace_packages(where: StrPath = ".", exclude: Iterable[str] = (), include: Iterable[str] = ("*",)) -> list[str]: ...
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
    options: Mapping[str, Any] = ...,
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
    command_options: Mapping[str, Mapping[str, tuple[Any, Any]]] = ...,
    package_data: Mapping[str, list[str]] = ...,
    include_package_data: bool = ...,
    libraries: list[str] = ...,
    headers: list[str] = ...,
    ext_package: str = ...,
    include_dirs: list[str] = ...,
    password: str = ...,
    fullname: str = ...,
    **attrs: Any,
) -> Distribution: ...

class Command(_Command):
    command_consumes_arguments: bool
    distribution: Distribution
    def __init__(self, dist: Distribution, **kw: Any) -> None: ...
    def ensure_string_list(self, option: str | list[str]) -> None: ...
    def reinitialize_command(self, command: _Command | str, reinit_subcommands: bool = False, **kw: Any) -> _Command: ...  # type: ignore[override]
    @abstractmethod
    def initialize_options(self) -> None: ...
    @abstractmethod
    def finalize_options(self) -> None: ...
    @abstractmethod
    def run(self) -> None: ...

class sic(str): ...
