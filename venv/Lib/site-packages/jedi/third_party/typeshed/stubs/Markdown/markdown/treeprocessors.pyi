from re import Pattern
from typing import ClassVar
from typing_extensions import TypeGuard
from xml.etree.ElementTree import Element

from markdown import util
from markdown.core import Markdown

def build_treeprocessors(md: Markdown, **kwargs) -> util.Registry[Treeprocessor]: ...
def isString(s: object) -> TypeGuard[str]: ...

class Treeprocessor(util.Processor):
    def run(self, root: Element) -> Element | None: ...

class InlineProcessor(Treeprocessor):
    inlinePatterns: util.Registry[InlineProcessor]
    ancestors: list[str]
    def __init__(self, md: Markdown) -> None: ...
    stashed_nodes: dict[str, Element | str]
    parent_map: dict[Element[str], Element[str]]
    def run(self, tree: Element, ancestors: list[str] | None = None) -> Element: ...

class PrettifyTreeprocessor(Treeprocessor): ...

class UnescapeTreeprocessor(Treeprocessor):
    RE: ClassVar[Pattern[str]]
    def unescape(self, text: str) -> str: ...
