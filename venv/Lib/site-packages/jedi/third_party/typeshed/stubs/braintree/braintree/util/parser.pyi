from xml.dom.minidom import Document

from .generator import _XML

binary_type = bytes

class Parser:
    doc: Document
    def __init__(self, xml: str | bytes) -> None: ...
    def parse(self) -> _XML: ...
