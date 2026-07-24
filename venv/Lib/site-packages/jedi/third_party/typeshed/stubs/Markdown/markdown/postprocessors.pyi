import re
from typing import ClassVar
from typing_extensions import deprecated

from markdown.core import Markdown

from . import util

def build_postprocessors(md: Markdown, **kwargs) -> util.Registry[Postprocessor]: ...

class Postprocessor(util.Processor):
    def run(self, text: str) -> str: ...

class RawHtmlPostprocessor(Postprocessor):
    BLOCK_LEVEL_REGEX: ClassVar[re.Pattern[str]]
    def isblocklevel(self, html: str) -> bool: ...
    def stash_to_string(self, text: str) -> str: ...

class AndSubstitutePostprocessor(Postprocessor): ...

@deprecated(
    "This class is deprecated and will be removed in the future; "
    "use [`UnescapeTreeprocessor`][markdown.treeprocessors.UnescapeTreeprocessor] instead."
)
class UnescapePostprocessor(Postprocessor):
    RE: ClassVar[re.Pattern[str]]
    def unescape(self, m: re.Match[str]) -> str: ...
