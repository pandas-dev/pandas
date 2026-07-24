from _typeshed import Incomplete
from re import Pattern
from typing import Final

# Can actually be imported from a mix of xml, lxml, et_xmlfile and defusedxml
# But et_xmlfile is untyped. openpyxl does not directly depend on lxml/lxml-stubs.
# And forcing a dependency on defusedxml felt overkill as it just wraps xml
# So for typing purposes, let's pretend xml is the only dependency.
# Prefer using protocols over these for parameters.
from xml.etree.ElementTree import Element as Element, QName as QName

from ._functions_overloads import (
    SubElement as SubElement,
    _HasTag,
    _HasText,
    fromstring as fromstring,
    iterparse as iterparse,
    tostring as tostring,
)

# from lxml.etree import xmlfile
# from et_xmlfile import xmlfile
xmlfile: Incomplete

NS_REGEX: Final[Pattern[str]]

def localname(node: _HasTag) -> str: ...
def whitespace(node: _HasText) -> None: ...
