from typing import Final
from xml.sax.expatreader import ExpatParser as _ExpatParser, _BoolType

__origin__: Final = "xml.sax.expatreader"

class DefusedExpatParser(_ExpatParser):
    forbid_dtd: bool
    forbid_entities: bool
    forbid_external: bool
    def __init__(
        self,
        namespaceHandling: _BoolType = 0,
        bufsize: int = 65516,
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
    def reset(self) -> None: ...

def create_parser(
    namespaceHandling: _BoolType = 0,
    bufsize: int = 65516,
    forbid_dtd: bool = False,
    forbid_entities: bool = True,
    forbid_external: bool = True,
) -> DefusedExpatParser: ...
