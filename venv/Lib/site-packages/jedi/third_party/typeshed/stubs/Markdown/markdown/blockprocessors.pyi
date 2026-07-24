from logging import Logger
from re import Match, Pattern
from typing import Any, ClassVar
from xml.etree.ElementTree import Element

from markdown.blockparser import BlockParser
from markdown.core import Markdown

logger: Logger

def build_block_parser(md: Markdown, **kwargs: Any) -> BlockParser: ...

class BlockProcessor:
    parser: BlockParser
    tab_length: int
    def __init__(self, parser: BlockParser) -> None: ...
    def lastChild(self, parent: Element) -> Element | None: ...
    def detab(self, text: str, length: int | None = None) -> tuple[str, str]: ...
    def looseDetab(self, text: str, level: int = 1) -> str: ...
    def test(self, parent: Element, block: str) -> bool: ...
    def run(self, parent: Element, blocks: list[str]) -> bool | None: ...

class ListIndentProcessor(BlockProcessor):
    ITEM_TYPES: list[str]
    LIST_TYPES: list[str]
    INDENT_RE: Pattern[str]
    def __init__(self, parser: BlockParser) -> None: ...  # Note: This was done because the args are sent as-is.
    def create_item(self, parent: Element, block: str) -> None: ...
    def get_level(self, parent: Element, block: str) -> tuple[int, Element]: ...

class CodeBlockProcessor(BlockProcessor): ...

class BlockQuoteProcessor(BlockProcessor):
    RE: Pattern[str]
    def clean(self, line: str) -> str: ...

class OListProcessor(BlockProcessor):
    TAG: ClassVar[str]
    STARTSWITH: ClassVar[str]
    LAZY_OL: ClassVar[bool]
    SIBLING_TAGS: ClassVar[list[str]]
    RE: Pattern[str]
    CHILD_RE: Pattern[str]
    INDENT_RE: Pattern[str]
    def __init__(self, parser: BlockParser) -> None: ...
    def get_items(self, block: str) -> list[str]: ...

class UListProcessor(OListProcessor):
    def __init__(self, parser: BlockParser) -> None: ...

class HashHeaderProcessor(BlockProcessor):
    RE: ClassVar[Pattern[str]]

class SetextHeaderProcessor(BlockProcessor):
    RE: ClassVar[Pattern[str]]

class HRProcessor(BlockProcessor):
    RE: ClassVar[str]
    SEARCH_RE: ClassVar[Pattern[str]]
    match: Match[str]

class EmptyBlockProcessor(BlockProcessor): ...

class ReferenceProcessor(BlockProcessor):
    RE: ClassVar[Pattern[str]]

class ParagraphProcessor(BlockProcessor): ...
