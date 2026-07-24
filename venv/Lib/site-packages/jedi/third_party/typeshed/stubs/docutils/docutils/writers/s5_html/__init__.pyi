import re
from _typeshed import StrPath
from typing import ClassVar, Final

from docutils import nodes
from docutils.writers import html4css1

__docformat__: Final = "reStructuredText"
themes_dir_path: Final[str]

def find_theme(name: StrPath) -> str: ...

class Writer(html4css1.Writer):
    settings_default_overrides: ClassVar[dict[str, int]]
    config_section_dependencies: ClassVar[tuple[str, ...]]
    translator_class: type[S5HTMLTranslator]

class S5HTMLTranslator(html4css1.HTMLTranslator):
    s5_stylesheet_template: ClassVar[str]
    disable_current_slide: ClassVar[str]
    layout_template: ClassVar[str]
    default_theme: ClassVar[str]
    base_theme_file: ClassVar[str]
    direct_theme_files: ClassVar[tuple[str, ...]]
    indirect_theme_files: ClassVar[tuple[str, ...]]
    required_theme_files: ClassVar[tuple[str, ...]]
    theme_file_path: str | None
    s5_footer: list[str]
    s5_header: list[str]
    section_count: int
    theme_files_copied: dict[str, bool]
    def __init__(self, document: nodes.document, /) -> None: ...
    def setup_theme(self) -> None: ...
    def copy_theme(self) -> None: ...
    files_to_skip_pattern: re.Pattern[str]
    def copy_file(self, name, source_dir, dest_dir): ...
    def depart_document(self, node: nodes.document) -> None: ...
    def depart_footer(self, node: nodes.footer) -> None: ...
    def depart_header(self, node: nodes.header) -> None: ...
    def visit_section(self, node: nodes.section) -> None: ...
    def visit_subtitle(self, node: nodes.subtitle) -> None: ...
    def visit_title(self, node: nodes.title) -> None: ...
