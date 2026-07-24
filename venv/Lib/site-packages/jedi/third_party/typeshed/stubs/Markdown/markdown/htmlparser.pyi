import html.parser as htmlparser
import re
from _frozen_importlib import ModuleSpec
from collections.abc import Sequence

from markdown import Markdown

spec: ModuleSpec
commentclose: re.Pattern[str]
blank_line_re: re.Pattern[str]

class HTMLExtractor(htmlparser.HTMLParser):
    empty_tags: set[str]
    lineno_start_cache: list[int]
    md: Markdown
    def __init__(self, md: Markdown, *args, **kwargs): ...
    inraw: bool
    intail: bool
    stack: list[str]
    cleandoc: list[str]
    @property
    def line_offset(self) -> int: ...
    def at_line_start(self) -> bool: ...
    def get_endtag_text(self, tag: str) -> str: ...
    def handle_starttag(self, tag: str, attrs: Sequence[tuple[str, str]]) -> None: ...  # type: ignore[override]
    def handle_empty_tag(self, data: str, is_block: bool) -> None: ...
    def handle_decl(self, data: str) -> None: ...
    def parse_bogus_comment(self, i: int, report: int = 0) -> int: ...
    def get_starttag_text(self) -> str: ...
