from typing import ClassVar, Final, TypeVar

from docutils.parsers.rst import states
from docutils.readers import standalone

__docformat__: Final = "reStructuredText"

_S = TypeVar("_S", bound=str | bytes)

class Reader(standalone.Reader[_S]):
    settings_default_overrides: ClassVar[dict[str, int]]
    inliner_class: ClassVar[type[states.Inliner]]
