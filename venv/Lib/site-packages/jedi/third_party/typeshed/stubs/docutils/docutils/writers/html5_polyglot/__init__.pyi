from pathlib import Path
from typing import ClassVar, Final

from docutils.writers import _html_base

__docformat__: Final = "reStructuredText"

class Writer(_html_base.Writer):
    default_stylesheets: ClassVar[list[str]]
    default_stylesheet_dirs: ClassVar[list[str]]
    default_template: ClassVar[Path]
    translator_class: type[HTMLTranslator]

class HTMLTranslator(_html_base.HTMLTranslator):
    supported_block_tags: set[str]
    supported_inline_tags: set[str]
