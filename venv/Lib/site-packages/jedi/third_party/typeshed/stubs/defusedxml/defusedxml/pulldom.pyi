from typing import Final
from xml.dom.pulldom import DOMEventStream
from xml.sax import _SupportsReadClose
from xml.sax.xmlreader import XMLReader

__origin__: Final = "xml.dom.pulldom"

def parse(
    stream_or_string: str | _SupportsReadClose[str] | _SupportsReadClose[bytes],
    parser: XMLReader | None = None,
    bufsize: int | None = None,
    forbid_dtd: bool = False,
    forbid_entities: bool = True,
    forbid_external: bool = True,
) -> DOMEventStream: ...
def parseString(
    string: str,
    parser: XMLReader | None = None,
    forbid_dtd: bool = False,
    forbid_entities: bool = True,
    forbid_external: bool = True,
) -> DOMEventStream: ...
