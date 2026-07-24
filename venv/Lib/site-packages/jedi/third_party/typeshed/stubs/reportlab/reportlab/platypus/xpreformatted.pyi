from reportlab.lib.styles import PropertySet
from reportlab.platypus.paragraph import Paragraph, ParaLines
from reportlab.platypus.paraparser import ParaFrag

class XPreformatted(Paragraph):
    def __init__(
        self,
        text: str,
        # NOTE: This should be a ParagraphStyle
        style: PropertySet,
        bulletText: str | None = None,
        frags: list[ParaFrag] | None = None,
        caseSensitive: int = 1,
        dedent: int = 0,
    ) -> None: ...
    def breakLinesCJK(self, width: float | list[float] | tuple[float, ...]) -> ParaLines | ParaFrag: ...

class PythonPreformatted(XPreformatted):
    formats: dict[str, tuple[str, str]]
    def __init__(
        self,
        text: str,
        # NOTE: This should be a ParagraphStyle
        style: PropertySet,
        bulletText: str | None = None,
        dedent: int = 0,
        frags: list[ParaFrag] | None = None,
    ) -> None: ...
    def escapeHtml(self, text: str) -> str: ...
    def fontify(self, code: str) -> str: ...

__all__ = ("XPreformatted", "PythonPreformatted")
