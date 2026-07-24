from typing import Any

from yaml.error import YAMLError
from yaml.nodes import Node

class SerializerError(YAMLError): ...

class Serializer:
    ANCHOR_TEMPLATE: Any
    use_encoding: Any
    use_explicit_start: Any
    use_explicit_end: Any
    use_version: Any
    use_tags: Any
    serialized_nodes: Any
    anchors: Any
    last_anchor_id: Any
    closed: Any
    def __init__(self, encoding=None, explicit_start=None, explicit_end=None, version=None, tags=None) -> None: ...
    def open(self) -> None: ...
    def close(self) -> None: ...
    def serialize(self, node: Node) -> None: ...
    def anchor_node(self, node): ...
    def generate_anchor(self, node): ...
    def serialize_node(self, node, parent, index): ...

__all__ = ["Serializer", "SerializerError"]
