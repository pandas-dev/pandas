from re import Pattern

from markdown.blockprocessors import BlockProcessor, ListIndentProcessor
from markdown.extensions import Extension

class DefListProcessor(BlockProcessor):
    RE: Pattern[str]
    NO_INDENT_RE: Pattern[str]

class DefListIndentProcessor(ListIndentProcessor): ...
class DefListExtension(Extension): ...

def makeExtension(**kwargs) -> DefListExtension: ...
