from _typeshed import Incomplete, StrPath
from typing import Generic, TypeVar

from .._distutils.dist import DistributionMetadata
from ..dist import Distribution
from . import expand

SingleCommandOptions: Incomplete
AllCommandOptions: Incomplete
Target = TypeVar("Target", bound=Distribution | DistributionMetadata)  # noqa: Y001 # Exists at runtime

def read_configuration(
    filepath: StrPath, find_others: bool = False, ignore_option_errors: bool = False
) -> dict[Incomplete, Incomplete]: ...
def apply_configuration(dist: Distribution, filepath: StrPath) -> Distribution: ...
def configuration_to_dict(
    handlers: tuple[ConfigHandler[Distribution | DistributionMetadata], ...]
) -> dict[Incomplete, Incomplete]: ...
def parse_configuration(
    distribution: Distribution, command_options: AllCommandOptions, ignore_option_errors: bool = False
) -> tuple[ConfigMetadataHandler, ConfigOptionsHandler]: ...

class ConfigHandler(Generic[Target]):
    section_prefix: str
    aliases: dict[str, str]
    ignore_option_errors: Incomplete
    target_obj: Incomplete
    sections: Incomplete
    set_options: Incomplete
    ensure_discovered: Incomplete
    def __init__(
        self,
        target_obj: Target,
        options: AllCommandOptions,
        ignore_option_errors,
        ensure_discovered: expand.EnsurePackagesDiscovered,
    ) -> None: ...
    @property
    def parsers(self) -> None: ...
    def __setitem__(self, option_name, value): ...
    def parse_section(self, section_options) -> None: ...
    def parse(self) -> None: ...

class ConfigMetadataHandler(ConfigHandler[DistributionMetadata]):
    section_prefix: str
    aliases: Incomplete
    strict_mode: bool
    package_dir: Incomplete
    root_dir: Incomplete
    def __init__(
        self,
        target_obj: DistributionMetadata,
        options: AllCommandOptions,
        ignore_option_errors: bool,
        ensure_discovered: expand.EnsurePackagesDiscovered,
        package_dir: dict[Incomplete, Incomplete] | None = None,
        root_dir: StrPath = ".",
    ) -> None: ...
    @property
    def parsers(self): ...

class ConfigOptionsHandler(ConfigHandler[Distribution]):
    section_prefix: str
    root_dir: Incomplete
    package_dir: Incomplete
    def __init__(
        self,
        target_obj: Distribution,
        options: AllCommandOptions,
        ignore_option_errors: bool,
        ensure_discovered: expand.EnsurePackagesDiscovered,
    ) -> None: ...
    @property
    def parsers(self): ...
    def parse_section_packages__find(self, section_options): ...
    def parse_section_entry_points(self, section_options) -> None: ...
    def parse_section_package_data(self, section_options) -> None: ...
    def parse_section_exclude_package_data(self, section_options) -> None: ...
    def parse_section_extras_require(self, section_options): ...
    def parse_section_data_files(self, section_options) -> None: ...
