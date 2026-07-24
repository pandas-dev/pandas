from _typeshed import SupportsRead
from typing import Final
from xml.dom.expatbuilder import ExpatBuilder as _ExpatBuilder, Namespaces as _Namespaces
from xml.dom.minidom import Document
from xml.dom.xmlbuilder import Options
from xml.parsers.expat import XMLParserType

__origin__: Final = "xml.dom.expatbuilder"

class DefusedExpatBuilder(_ExpatBuilder):
    forbid_dtd: bool
    forbid_entities: bool
    forbid_external: bool
    def __init__(
        self, options: Options | None = None, forbid_dtd: bool = False, forbid_entities: bool = True, forbid_external: bool = True
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
    def install(self, parser: XMLParserType) -> None: ...

class DefusedExpatBuilderNS(_Namespaces, DefusedExpatBuilder):
    def install(self, parser: XMLParserType) -> None: ...
    def reset(self) -> None: ...

def parse(
    file: str | SupportsRead[bytes | str],
    namespaces: bool = True,
    forbid_dtd: bool = False,
    forbid_entities: bool = True,
    forbid_external: bool = True,
) -> Document: ...
def parseString(
    string: str, namespaces: bool = True, forbid_dtd: bool = False, forbid_entities: bool = True, forbid_external: bool = True
) -> Document: ...
