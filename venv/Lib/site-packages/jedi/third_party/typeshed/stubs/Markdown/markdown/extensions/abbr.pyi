from re import Pattern
from typing import ClassVar
from typing_extensions import deprecated
from xml.etree.ElementTree import Element

from markdown.blockparser import BlockParser
from markdown.blockprocessors import BlockProcessor
from markdown.core import Markdown
from markdown.extensions import Extension
from markdown.inlinepatterns import InlineProcessor
from markdown.treeprocessors import Treeprocessor

class AbbrExtension(Extension):
    def reset(self) -> None: ...
    def reset_glossary(self) -> None: ...
    def load_glossary(self, dictionary: dict[str, str]) -> None: ...

class AbbrTreeprocessor(Treeprocessor):
    RE: Pattern[str] | None
    abbrs: dict[str, str]
    def __init__(self, md: Markdown | None = None, abbrs: dict[str, str] | None = None) -> None: ...
    def create_element(self, title: str, text: str, tail: str) -> Element: ...
    def iter_element(self, el: Element, parent: Element | None = None) -> None: ...

# Techinically it is the same type as `AbbrPreprocessor` just not deprecated.
class AbbrBlockprocessor(BlockProcessor):
    RE: ClassVar[Pattern[str]]
    abbrs: dict[str, str]
    def __init__(self, parser: BlockParser, abbrs: dict[str, str]) -> None: ...

@deprecated("This class will be removed in the future; use `AbbrTreeprocessor` instead.")
class AbbrPreprocessor(AbbrBlockprocessor): ...

@deprecated("This class will be removed in the future; use `AbbrTreeprocessor` instead.")
class AbbrInlineProcessor(InlineProcessor):
    title: str
    def __init__(self, pattern: str, title: str) -> None: ...

def makeExtension(**kwargs) -> AbbrExtension: ...
