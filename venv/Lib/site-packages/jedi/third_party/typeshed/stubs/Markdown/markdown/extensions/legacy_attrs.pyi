from re import Pattern
from xml.etree.ElementTree import Element

from markdown.extensions import Extension
from markdown.treeprocessors import Treeprocessor

ATTR_RE: Pattern[str]

class LegacyAttrs(Treeprocessor):
    def run(self, doc: Element) -> None: ...
    def handleAttributes(self, el: Element, txt: str) -> str: ...

class LegacyAttrExtension(Extension): ...

def makeExtension(**kwargs) -> LegacyAttrExtension: ...
