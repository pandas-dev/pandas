from collections.abc import Iterable
from re import Pattern
from typing import Any, ClassVar

from markdown.core import Markdown
from markdown.extensions import Extension
from markdown.preprocessors import Preprocessor

class FencedCodeExtension(Extension):
    def __init__(self, **kwargs) -> None: ...

class FencedBlockPreprocessor(Preprocessor):
    FENCED_BLOCK_RE: ClassVar[Pattern[str]]
    checked_for_deps: bool
    codehilite_conf: dict[str, Any]
    use_attr_list: bool
    bool_options: list[str]
    def __init__(self, md: Markdown, config: dict[str, Any]) -> None: ...
    def handle_attrs(self, attrs: Iterable[tuple[str, str]]) -> tuple[str, list[str], dict[str, Any]]: ...

def makeExtension(**kwargs) -> FencedCodeExtension: ...
