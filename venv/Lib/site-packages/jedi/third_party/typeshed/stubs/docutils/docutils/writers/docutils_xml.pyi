from collections.abc import Callable
from typing import ClassVar, Final
from xml.sax.handler import ContentHandler
from xml.sax.xmlreader import Locator, XMLReader

import docutils
from docutils import nodes, writers
from docutils.frontend import Values

__docformat__: Final = "reStructuredText"

class RawXmlError(docutils.ApplicationError): ...

class Writer(writers.Writer[str]):
    settings_defaults: ClassVar[dict[str, str]]
    config_section: ClassVar[str]
    config_section_dependencies: ClassVar[tuple[str, ...]]
    translator_class: type[XMLTranslator]
    visitor: XMLTranslator

class XMLTranslator(nodes.GenericNodeVisitor):
    doctype: ClassVar[str]
    generator: ClassVar[str]
    xmlparser: ClassVar[XMLReader]
    warn: Callable[..., nodes.system_message]
    error: Callable[..., nodes.system_message]
    settings: Values
    indent: str
    newline: str
    level: int
    in_simple: int
    fixed_text: int
    output: list[str]
    the_handle: TestXml
    def __init__(self, document: nodes.document) -> None: ...
    simple_nodes: ClassVar[tuple[type[nodes.Element], ...]]
    def default_visit(self, node: nodes.Element) -> None: ...  # type: ignore[override]
    def default_departure(self, node: nodes.Element) -> None: ...  # type: ignore[override]
    def visit_Text(self, node: nodes.Text) -> None: ...
    def depart_Text(self, node: nodes.Text) -> None: ...
    def visit_raw(self, node: nodes.raw) -> None: ...

class TestXml(ContentHandler):
    locator: Locator
    def setDocumentLocator(self, locator: Locator) -> None: ...
