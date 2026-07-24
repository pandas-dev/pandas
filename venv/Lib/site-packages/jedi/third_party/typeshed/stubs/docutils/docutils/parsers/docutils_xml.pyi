import xml.etree.ElementTree as ET
from typing import ClassVar, Final

from docutils import nodes, parsers

__docformat__: Final = "reStructuredText"

class Parser(parsers.Parser):
    config_section_dependencies: ClassVar[tuple[str, ...]]
    settings_default_overrides: ClassVar[dict[str, bool]]

class Unknown(nodes.Special, nodes.Inline, nodes.Element): ...

def parse_element(inputstring: str, document: nodes.document | None = None) -> nodes.Element: ...
def element2node(element: ET.Element | None, document: nodes.document | None = None, unindent: bool = True) -> nodes.Element: ...
def append_text(node: nodes.Element, text: str | None, unindent: bool | None) -> None: ...
