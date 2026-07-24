from _typeshed import Incomplete, StrPath
from typing import Final
from typing_extensions import deprecated

from docutils import SettingsSpec
from docutils.io import FileInput, Input, Output
from docutils.parsers import Parser
from docutils.readers import Reader
from docutils.utils import SystemMessage
from docutils.writers import Writer, _WriterParts

__docformat__: Final = "reStructuredText"

class Publisher:
    document: Incomplete | None
    reader: Reader[Incomplete]
    parser: Parser
    writer: Writer[Incomplete]
    source: Input[Incomplete]
    source_class: Incomplete
    destination: Output | None
    destination_class: Incomplete
    settings: dict[str, Incomplete]
    def __init__(
        self,
        reader: Reader[Incomplete] | None = None,
        parser: Parser | None = None,
        writer: Writer[Incomplete] | None = None,
        source: Input[Incomplete] | None = None,
        source_class=...,
        destination: Output | None = None,
        destination_class=...,
        settings: dict[str, Incomplete] | None = None,
    ) -> None: ...
    def set_reader(self, reader: str, parser: Parser | None = None, parser_name: str | None = None) -> None: ...
    def set_writer(self, writer_name: str) -> None: ...
    @deprecated("The `Publisher.set_components()` will be removed in Docutils 2.0.")
    def set_components(self, reader_name: str, parser_name: str, writer_name: str) -> None: ...
    def get_settings(
        self,
        usage: str | None = None,
        description: str | None = None,
        settings_spec: SettingsSpec | None = None,
        config_section: str | None = None,
        **defaults,
    ): ...
    def process_programmatic_settings(self, settings_spec, settings_overrides, config_section) -> None: ...
    def process_command_line(
        self,
        argv: list[str] | None = None,
        usage=None,
        description: str | None = None,
        settings_spec=None,
        config_section=None,
        **defaults,
    ) -> None: ...
    def set_io(self, source_path: StrPath | None = None, destination_path: StrPath | None = None) -> None: ...
    def set_source(self, source: str | None = None, source_path: StrPath | None = None) -> None: ...
    def set_destination(self, destination=None, destination_path: StrPath | None = None) -> None: ...
    def apply_transforms(self) -> None: ...
    def publish(
        self,
        argv: list[str] | None = None,
        usage: str | None = None,
        description: str | None = None,
        settings_spec=None,
        settings_overrides=None,
        config_section: str | None = None,
        enable_exit_status: bool = False,
    ): ...
    def debugging_dumps(self) -> None: ...
    def prompt(self) -> None: ...
    def report_Exception(self, error: BaseException) -> None: ...
    def report_SystemMessage(self, error: SystemMessage) -> None: ...
    def report_UnicodeError(self, error: UnicodeEncodeError) -> None: ...

default_usage: Final[str]
default_description: Final[str]

def publish_cmdline(
    reader: Reader[Incomplete] | None = None,
    reader_name: str | None = None,
    parser: Parser | None = None,
    parser_name: str | None = None,
    writer: Writer[Incomplete] | None = None,
    writer_name: str | None = None,
    settings=None,
    settings_spec=None,
    settings_overrides=None,
    config_section: str | None = None,
    enable_exit_status: bool = True,
    argv: list[str] | None = None,
    usage: str = "%prog [options] [<source> [<destination>]]",
    description: str = ...,
): ...
def publish_file(
    source=None,
    source_path: StrPath | None = None,
    destination=None,
    destination_path: StrPath | None = None,
    reader=None,
    reader_name: str | None = None,
    parser=None,
    parser_name: str | None = None,
    writer=None,
    writer_name: str | None = None,
    settings=None,
    settings_spec=None,
    settings_overrides=None,
    config_section: str | None = None,
    enable_exit_status: bool = False,
): ...
def publish_string(
    source,
    source_path: StrPath | None = None,
    destination_path: StrPath | None = None,
    reader=None,
    reader_name: str | None = None,
    parser=None,
    parser_name: str | None = None,
    writer=None,
    writer_name: str | None = None,
    settings=None,
    settings_spec=None,
    settings_overrides=None,
    config_section: str | None = None,
    enable_exit_status: bool = False,
): ...
def publish_parts(
    source,
    source_path: StrPath | None = None,
    source_class=...,
    destination_path: StrPath | None = None,
    reader=None,
    reader_name: str | None = None,
    parser=None,
    parser_name: str | None = None,
    writer=None,
    writer_name: str | None = None,
    settings=None,
    settings_spec=None,
    settings_overrides: dict[str, Incomplete] | None = None,
    config_section: str | None = None,
    enable_exit_status: bool = False,
) -> _WriterParts: ...
def publish_doctree(
    source,
    source_path: StrPath | None = None,
    source_class=...,
    reader=None,
    reader_name: str | None = None,
    parser=None,
    parser_name: str | None = None,
    settings=None,
    settings_spec=None,
    settings_overrides=None,
    config_section: str | None = None,
    enable_exit_status: bool = False,
): ...
def publish_from_doctree(
    document,
    destination_path: StrPath | None = None,
    writer=None,
    writer_name: str | None = None,
    settings=None,
    settings_spec=None,
    settings_overrides=None,
    config_section: str | None = None,
    enable_exit_status: bool = False,
): ...
@deprecated("The `publish_cmdline_to_binary()` is deprecated  by `publish_cmdline()` and will be removed in Docutils 0.24.")
def publish_cmdline_to_binary(
    reader=None,
    reader_name: str = "standalone",
    parser=None,
    parser_name: str = "restructuredtext",
    writer=None,
    writer_name: str = "pseudoxml",
    settings=None,
    settings_spec=None,
    settings_overrides=None,
    config_section: str | None = None,
    enable_exit_status: bool = True,
    argv: list[str] | None = None,
    usage: str = "%prog [options] [<source> [<destination>]]",
    description: str = ...,
    destination=None,
    destination_class=...,
): ...
def publish_programmatically(
    source_class: type[FileInput],
    source,
    source_path: StrPath | None,
    destination_class,
    destination,
    destination_path: StrPath | None,
    reader,
    reader_name: str,
    parser,
    parser_name: str,
    writer,
    writer_name: str,
    settings,
    settings_spec,
    settings_overrides,
    config_section: str,
    enable_exit_status: bool,
) -> tuple[str | bytes | None, Publisher]: ...
def rst2something(writer: str, documenttype: str, doc_path: str = "") -> None: ...
def rst2html() -> None: ...
def rst2html4() -> None: ...
def rst2html5() -> None: ...
def rst2latex() -> None: ...
def rst2man() -> None: ...
def rst2odt() -> None: ...
def rst2pseudoxml() -> None: ...
def rst2s5() -> None: ...
def rst2xetex() -> None: ...
def rst2xml() -> None: ...
