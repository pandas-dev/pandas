import re
from xml.etree.ElementTree import Element

__all__ = ["to_html_string", "to_xhtml_string"]

RE_AMP: re.Pattern[str]

def to_html_string(element: Element) -> str: ...
def to_xhtml_string(element: Element) -> str: ...
