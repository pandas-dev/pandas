import optparse
from _typeshed import Incomplete, StrPath
from collections.abc import Iterable, Mapping, Sequence
from configparser import RawConfigParser
from typing import Any, ClassVar, Final, Literal, Protocol, overload, type_check_only
from typing_extensions import deprecated

from docutils import SettingsSpec
from docutils.utils import DependencyList

__docformat__: Final = "reStructuredText"

@type_check_only
class _OptionValidator(Protocol):
    def __call__(
        self,
        setting: str,
        value: str | None,
        option_parser: OptionParser,
        /,
        config_parser: ConfigParser | None = None,
        config_section: str | None = None,
    ) -> Any: ...

@deprecated("Deprecated and will be removed with the switch to from optparse to argparse.")
def store_multiple(option: optparse.Option, opt: str, value, parser: OptionParser, *args: str, **kwargs) -> None: ...
@deprecated("Deprecated and will be removed with the switch to from optparse to argparse.")
def read_config_file(option: optparse.Option, opt: str, value, parser: OptionParser) -> None: ...
def validate_encoding(
    setting: str,
    value: str | None = None,
    option_parser: OptionParser | None = None,
    config_parser: ConfigParser | None = None,
    config_section: str | None = None,
) -> str: ...
def validate_encoding_error_handler(
    setting: str,
    value: str | None = None,
    option_parser: OptionParser | None = None,
    config_parser: ConfigParser | None = None,
    config_section: str | None = None,
) -> str: ...
def validate_encoding_and_error_handler(
    setting: str,
    value: str | None = None,
    option_parser: OptionParser | None = None,
    config_parser: ConfigParser | None = None,
    config_section: str | None = None,
) -> str: ...
def validate_boolean(
    setting: str | bool,
    value: str | None = None,
    option_parser: OptionParser | None = None,
    config_parser: ConfigParser | None = None,
    config_section: str | None = None,
) -> bool: ...
def validate_ternary(
    setting: str | bool,
    value: str | None = None,
    option_parser: OptionParser | None = None,
    config_parser: ConfigParser | None = None,
    config_section: str | None = None,
) -> str | bool | None: ...
def validate_nonnegative_int(
    setting: str | int,
    value: str | None = None,
    option_parser: OptionParser | None = None,
    config_parser: ConfigParser | None = None,
    config_section: str | None = None,
) -> int: ...
def validate_threshold(
    setting: str | int,
    value: str | None = None,
    option_parser: OptionParser | None = None,
    config_parser: ConfigParser | None = None,
    config_section: str | None = None,
) -> int: ...
def validate_colon_separated_string_list(
    setting: str | list[str],
    value: str | None = None,
    option_parser: OptionParser | None = None,
    config_parser: ConfigParser | None = None,
    config_section: str | None = None,
) -> list[str]: ...
def validate_comma_separated_list(
    setting: str | list[str],
    value: str | None = None,
    option_parser: OptionParser | None = None,
    config_parser: ConfigParser | None = None,
    config_section: str | None = None,
) -> list[str]: ...
def validate_math_output(
    setting: str,
    value: str | None = None,
    option_parser: OptionParser | None = None,
    config_parser: ConfigParser | None = None,
    config_section: str | None = None,
) -> tuple[()] | tuple[str, str]: ...
def validate_url_trailing_slash(
    setting: str,
    value: str | None = None,
    option_parser: OptionParser | None = None,
    config_parser: ConfigParser | None = None,
    config_section: str | None = None,
) -> str: ...
def validate_dependency_file(
    setting: str | None,
    value: str | None = None,
    option_parser: OptionParser | None = None,
    config_parser: ConfigParser | None = None,
    config_section: str | None = None,
) -> DependencyList: ...
def validate_strip_class(
    setting: str,
    value: str | None = None,
    option_parser: OptionParser | None = None,
    config_parser: ConfigParser | None = None,
    config_section: str | None = None,
) -> list[str]: ...
def validate_smartquotes_locales(
    setting: str | list[str | tuple[str, str]],
    value: str | None = None,
    option_parser: OptionParser | None = None,
    config_parser: ConfigParser | None = None,
    config_section: str | None = None,
) -> list[tuple[str, Sequence[str]]]: ...
def make_paths_absolute(
    pathdict: dict[str, list[StrPath] | StrPath], keys: tuple[str], base_path: StrPath | None = None
) -> None: ...
@deprecated("The `frontend.make_one_path_absolute` will be removed in Docutils 2.0 or later.")
def make_one_path_absolute(base_path: StrPath, path: StrPath) -> str: ...
def filter_settings_spec(settings_spec, *exclude, **replace) -> tuple[Any, ...]: ...
@deprecated("The `frontend.Values` class will be removed in Docutils 2.0 or later.")
class Values(optparse.Values):
    record_dependencies: DependencyList
    def __init__(self, defaults: dict[str, Any] | None = None) -> None: ...
    def update(self, other_dict: Values | Mapping[str, Incomplete], option_parser: OptionParser) -> None: ...
    def copy(self) -> Values: ...
    def setdefault(self, name: str, default): ...

@deprecated("The `frontend.Option` class will be removed in Docutils 2.0 or later.")
class Option(optparse.Option):
    ATTRS: list[str]
    validator: _OptionValidator
    overrides: str | None
    def __init__(self, *args: str | None, **kwargs) -> None: ...

@deprecated(
    "The `frontend.OptionParser` class will be replaced by a subclass of `argparse.ArgumentParser` in Docutils 2.0 or later."
)
class OptionParser(optparse.OptionParser, SettingsSpec):
    standard_config_files: ClassVar[list[str]]
    threshold_choices: ClassVar[tuple[str, ...]]
    thresholds: ClassVar[dict[str, int]]
    booleans: ClassVar[dict[str, bool]]
    default_error_encoding: ClassVar[str]
    default_error_encoding_error_handler: ClassVar[str]
    config_section: ClassVar[str]
    version_template: ClassVar[str]
    details: str
    lists: dict[str, Literal[True]]
    config_files: list[str]
    relative_path_settings: ClassVar[tuple[str, ...]]
    version: str
    components: tuple[SettingsSpec, ...]
    def __init__(
        self,
        components: Iterable[SettingsSpec | type[SettingsSpec]] = (),
        defaults: Mapping[str, Any] | None = None,
        read_config_files: bool | None = False,
        *args,
        **kwargs,
    ) -> None: ...
    def populate_from_components(self, components: Iterable[SettingsSpec]) -> None: ...
    @classmethod
    def get_standard_config_files(cls) -> Sequence[StrPath]: ...
    def get_standard_config_settings(self) -> Values: ...
    def get_config_file_settings(self, config_file: str) -> dict[str, Incomplete]: ...
    def check_values(self, values: Values, args: list[str]) -> Values: ...  # type: ignore[override]
    def check_args(self, args: list[str]) -> tuple[str | None, str | None]: ...
    def get_default_values(self) -> Values: ...
    def get_option_by_dest(self, dest: str) -> Option: ...

class ConfigParser(RawConfigParser):
    old_settings: ClassVar[dict[str, tuple[str, str]]]
    old_warning: ClassVar[str]
    not_utf8_error: ClassVar[str]
    @overload  # type: ignore[override]
    def read(self, filenames: str | Sequence[str]) -> list[str]: ...
    @overload
    @deprecated("The `option_parser` parameter is deprecated and will be removed in Docutils 0.24.")
    def read(self, filenames: str | Sequence[str], option_parser: OptionParser | None) -> list[str]: ...
    def handle_old_config(self, filename: str) -> None: ...
    def validate_settings(self, filename: str, option_parser: OptionParser) -> None: ...
    def optionxform(self, optionstr: str) -> str: ...

class ConfigDeprecationWarning(FutureWarning): ...

def get_default_settings(*components: SettingsSpec) -> Values: ...
