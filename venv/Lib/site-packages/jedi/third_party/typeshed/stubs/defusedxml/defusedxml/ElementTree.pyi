from _typeshed import ReadableBuffer
from collections.abc import Sequence
from typing import Final
from xml.etree.ElementTree import (
    Element,
    ElementTree,
    ParseError as ParseError,
    XMLParser as _XMLParser,
    _FileRead,
    _IterParseIterator,
    _Target,
    tostring as tostring,
)

__origin__: Final = "xml.etree.ElementTree"

class DefusedXMLParser(_XMLParser):
    forbid_dtd: bool
    forbid_entities: bool
    forbid_external: bool
    def __init__(
        self,
        html: object = ...,  # argument is deprecated, if bool(html) is True you will get TypeError
        target: _Target | None = None,
        encoding: str | None = None,
        forbid_dtd: bool = False,
        forbid_entities: bool = True,
        forbid_external: bool = True,
    ) -> None: ...
    def defused_start_doctype_decl(self, name: str, sysid: str | None, pubid: str | None, has_internal_subset: bool) -> None: ...
    def defused_entity_decl(
        self,
        name: str,
        is_parameter_entity: bool,
        value: str | None,
        base: str | None,
        sysid: str,
        pubid: str | None,
        notation_name: str | None,
    ) -> None: ...
    def defused_unparsed_entity_decl(
        self, name: str, base: str | None, sysid: str, pubid: str | None, notation_name: str
    ) -> None: ...
    def defused_external_entity_ref_handler(
        self, context: str, base: str | None, sysid: str | None, pubid: str | None
    ) -> None: ...

XMLTreeBuilder = DefusedXMLParser
XMLParse = DefusedXMLParser
XMLParser = DefusedXMLParser

# wrapper to xml.etree.ElementTree.parse
def parse(
    source: _FileRead,
    parser: XMLParser | None = None,
    forbid_dtd: bool = False,
    forbid_entities: bool = True,
    forbid_external: bool = True,
) -> ElementTree: ...

# wrapper to xml.etree.ElementTree.iterparse
def iterparse(
    source: _FileRead,
    events: Sequence[str] | None = None,
    parser: XMLParser | None = None,
    forbid_dtd: bool = False,
    forbid_entities: bool = True,
    forbid_external: bool = True,
) -> _IterParseIterator: ...
def fromstring(
    text: str | ReadableBuffer, forbid_dtd: bool = False, forbid_entities: bool = True, forbid_external: bool = True
) -> Element: ...

XML = fromstring

__all__ = ["ParseError", "XML", "XMLParse", "XMLParser", "XMLTreeBuilder", "fromstring", "iterparse", "parse", "tostring"]
