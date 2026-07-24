from collections.abc import Sequence
from typing import Final

from docutils.parsers.rst import Directive

__docformat__: Final = "reStructuredText"

class Contents(Directive):
    backlinks_values: Sequence[str]
    def backlinks(arg): ...

class Sectnum(Directive): ...
class Header(Directive): ...
class Footer(Directive): ...
