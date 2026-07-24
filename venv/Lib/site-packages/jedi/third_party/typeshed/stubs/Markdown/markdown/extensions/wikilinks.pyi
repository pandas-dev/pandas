from typing import Any

from markdown.core import Markdown
from markdown.extensions import Extension
from markdown.inlinepatterns import InlineProcessor

def build_url(label: str, base: str, end: str) -> str: ...

class WikiLinkExtension(Extension):
    def __init__(self, **kwargs) -> None: ...
    md: Markdown

class WikiLinksInlineProcessor(InlineProcessor):
    config: dict[str, Any]
    def __init__(self, pattern: str, config: dict[str, Any]) -> None: ...

def makeExtension(**kwargs) -> WikiLinkExtension: ...
