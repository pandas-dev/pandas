from _typeshed import ReadableBuffer, Unused
from typing import Final
from xml.sax import ErrorHandler as _ErrorHandler, _Source, xmlreader
from xml.sax.handler import _ContentHandlerProtocol

from .expatreader import DefusedExpatParser

__origin__: Final = "xml.sax"

def parse(
    source: xmlreader.InputSource | _Source,
    handler: _ContentHandlerProtocol,
    errorHandler: _ErrorHandler = ...,
    forbid_dtd: bool = False,
    forbid_entities: bool = True,
    forbid_external: bool = True,
) -> None: ...
def parseString(
    string: ReadableBuffer,
    handler: _ContentHandlerProtocol,
    errorHandler: _ErrorHandler = ...,
    forbid_dtd: bool = False,
    forbid_entities: bool = True,
    forbid_external: bool = True,
) -> None: ...
def make_parser(parser_list: Unused = []) -> DefusedExpatParser: ...
