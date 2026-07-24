from collections.abc import Container
from typing import Final

ALLOWED_CSS_PROPERTIES: Final[frozenset[str]]
ALLOWED_SVG_PROPERTIES: Final[frozenset[str]]

class CSSSanitizer:
    allowed_css_properties: Container[str]
    allowed_svg_properties: Container[str]

    def __init__(self, allowed_css_properties: Container[str] = ..., allowed_svg_properties: Container[str] = ...) -> None: ...
    def sanitize_css(self, style: str) -> str: ...
