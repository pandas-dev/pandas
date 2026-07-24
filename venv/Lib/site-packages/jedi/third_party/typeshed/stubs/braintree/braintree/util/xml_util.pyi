from .generator import _XML

class XmlUtil:
    @staticmethod
    def xml_from_dict(dict: _XML) -> str: ...
    @staticmethod
    def dict_from_xml(xml: str | bytes) -> _XML: ...
