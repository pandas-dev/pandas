from typing import Final

from .ElementTree import (
    XML as XML,
    ParseError as ParseError,
    XMLParse as XMLParse,
    XMLParser as XMLParser,
    XMLTreeBuilder as XMLTreeBuilder,
    fromstring as fromstring,
    iterparse as iterparse,
    parse as parse,
    tostring as tostring,
)

__origin__: Final = "xml.etree.cElementTree"
__all__ = ["ParseError", "XML", "XMLParse", "XMLParser", "XMLTreeBuilder", "fromstring", "iterparse", "parse", "tostring"]
