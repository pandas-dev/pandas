from _typeshed import Incomplete, StrPath
from abc import abstractmethod
from collections.abc import Callable, Iterable
from typing import Any, ClassVar, Generic, TypeVar
from typing_extensions import TypeAlias

from .._distutils.dist import DistributionMetadata
from ..dist import Distribution
from . import expand

SingleCommandOptions: TypeAlias = dict[str, tuple[str, Any]]
AllCommandOptions: TypeAlias = dict[str, SingleCommandOptions]
Target = TypeVar("Target", Distribution, DistributionMetadata)  # noqa: Y001 # Exists at runtime

def read_configuration(
    filepath: StrPath, find_others: bool = False, ignore_option_errors: bool = False
) -> dict[Incomplete, Incomplete]: ...
def apply_configuration(dist: Distribution, filepath: StrPath) -> Distribution: ...
def configuration_to_dict(
    handlers: Iterable[ConfigHandler[Distribution] | ConfigHandler[DistributionMetadata]],
) -> dict[Incomplete, Incomplete]: ...
def parse_configuration(
    distribution: Distribution, command_options: AllCommandOptions, ignore_option_errors: bool = False
) -> tuple[ConfigMetadataHandler, ConfigOptionsHandler]: ...

class ConfigHandler(Generic[Target]):
    section_prefix: str
    aliases: ClassVar[dict[str, str]]
    ignore_option_errors: Incomplete
    target_obj: Target
    sections: dict[str, SingleCommandOptions]
    set_options: list[str]
    ensure_discovered: expand.EnsurePackagesDiscovered
    def __init__(
        self,
        target_obj: Target,
        options: AllCommandOptions,
        ignore_option_errors,
        ensure_discovered: expand.EnsurePackagesDiscovered,
    ) -> None: ...
    @property
    @abstractmethod
    def parsers(self) -> dict[str, Callable[..., Incomplete]]: ...
    def __setitem__(self, option_name, value): ...
    def parse_section(self, section_options) -> None: ...
    def parse(self) -> None: ...

class ConfigMetadataHandler(ConfigHandler[DistributionMetadata]):
    section_prefix: str
    aliases: ClassVar[dict[str, str]]
    strict_mode: bool
    package_dir: dict[Incomplete, Incomplete] | None
    root_dir: StrPath | None
    def __init__(
        self,
        target_obj: DistributionMetadata,
        options: AllCommandOptions,
        ignore_option_errors: bool,
        ensure_discovered: expand.EnsurePackagesDiscovered,
        package_dir: dict[Incomplete, Incomplete] | None = None,
        root_dir: StrPath | None = ".",
    ) -> None: ...
    @property
    def parsers(self) -> dict[str, Callable[..., Incomplete]]: ...

class ConfigOptionsHandler(ConfigHandler[Distribution]):
    section_prefix: str
    root_dir: str | None
    package_dir: dict[str, str]
    def __init__(
        self,
        target_obj: Distribution,
        options: AllCommandOptions,
        ignore_option_errors: bool,
        ensure_discovered: expand.EnsurePackagesDiscovered,
    ) -> None: ...
    @property
    def parsers(self) -> dict[str, Callable[..., Incomplete]]: ...
    def parse_section_packages__find(self, section_options): ...
    def parse_section_entry_points(self, section_options) -> None: ...
    def parse_section_package_data(self, section_options) -> None: ...
    def parse_section_exclude_package_data(self, section_options) -> None: ...
    def parse_section_extras_require(self, section_options) -> None: ...
    def parse_section_data_files(self, section_options) -> None: ...
