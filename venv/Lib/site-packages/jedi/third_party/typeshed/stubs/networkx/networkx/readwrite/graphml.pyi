from _typeshed import Incomplete
from collections.abc import Generator
from typing import Final, Literal

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = [
    "write_graphml",
    "read_graphml",
    "generate_graphml",
    "write_graphml_xml",
    "write_graphml_lxml",
    "parse_graphml",
    "GraphMLWriter",
    "GraphMLReader",
]

def write_graphml_xml(
    G: Graph[_Node],
    path,
    encoding: str = "utf-8",
    prettyprint: bool = True,
    infer_numeric_types: bool = False,
    named_key_ids: bool = False,
    edge_id_from_attribute=None,
) -> None: ...
def write_graphml_lxml(
    G: Graph[_Node],
    path,
    encoding: str = "utf-8",
    prettyprint: bool = True,
    infer_numeric_types: bool = False,
    named_key_ids: bool = False,
    edge_id_from_attribute=None,
): ...
def generate_graphml(
    G: Graph[_Node], encoding: str = "utf-8", prettyprint: bool = True, named_key_ids: bool = False, edge_id_from_attribute=None
) -> Generator[Incomplete, Incomplete]: ...
@_dispatchable
def read_graphml(path, node_type=..., edge_key_type=..., force_multigraph: bool = False): ...
@_dispatchable
def parse_graphml(graphml_string, node_type=..., edge_key_type=..., force_multigraph: bool = False): ...

class GraphML:
    NS_GRAPHML: Final[str]
    NS_XSI: Final[str]
    NS_Y: Final[str]
    SCHEMALOCATION: Final[str]
    xml_type: Incomplete
    python_type: Incomplete
    def construct_types(self) -> None: ...
    convert_bool: Final[dict[Literal["true", "false", "0", 0, "1", 1], bool]]
    def get_xml_type(self, key): ...

class GraphMLWriter(GraphML):
    myElement: Incomplete
    infer_numeric_types: Incomplete
    prettyprint: Incomplete
    named_key_ids: Incomplete
    edge_id_from_attribute: Incomplete
    encoding: Incomplete
    xml: Incomplete
    keys: Incomplete
    attributes: Incomplete
    attribute_types: Incomplete
    def __init__(
        self,
        graph=None,
        encoding: str = "utf-8",
        prettyprint: bool = True,
        infer_numeric_types: bool = False,
        named_key_ids: bool = False,
        edge_id_from_attribute=None,
    ) -> None: ...
    def attr_type(self, name, scope, value): ...
    def get_key(self, name, attr_type, scope, default): ...
    def add_data(self, name, element_type, value, scope: str = "all", default=None): ...
    def add_attributes(self, scope, xml_obj, data, default) -> None: ...
    def add_nodes(self, G: Graph[_Node], graph_element) -> None: ...
    def add_edges(self, G: Graph[_Node], graph_element) -> None: ...
    def add_graph_element(self, G: Graph[_Node]) -> None: ...
    def add_graphs(self, graph_list) -> None: ...
    def dump(self, stream) -> None: ...
    def indent(self, elem, level: int = 0) -> None: ...

class IncrementalElement:
    xml: Incomplete
    prettyprint: Incomplete
    def __init__(self, xml, prettyprint) -> None: ...
    def append(self, element) -> None: ...

class GraphMLWriterLxml(GraphMLWriter):
    myElement: Incomplete
    named_key_ids: Incomplete
    edge_id_from_attribute: Incomplete
    infer_numeric_types: Incomplete
    xml: Incomplete
    keys: Incomplete
    attribute_types: Incomplete
    def __init__(
        self,
        path,
        graph=None,
        encoding: str = "utf-8",
        prettyprint: bool = True,
        infer_numeric_types: bool = False,
        named_key_ids: bool = False,
        edge_id_from_attribute=None,
    ) -> None: ...
    def add_graph_element(self, G: Graph[_Node]) -> None: ...
    def add_attributes(self, scope, xml_obj, data, default) -> None: ...
    def dump(self, stream=None) -> None: ...

write_graphml = write_graphml_lxml

class GraphMLReader(GraphML):
    node_type: Incomplete
    edge_key_type: Incomplete
    multigraph: Incomplete
    edge_ids: Incomplete
    def __init__(self, node_type=..., edge_key_type=..., force_multigraph: bool = False) -> None: ...
    xml: Incomplete
    def __call__(self, path=None, string=None) -> Generator[Incomplete]: ...
    def make_graph(self, graph_xml, graphml_keys, defaults, G=None): ...
    def add_node(self, G: Graph[_Node], node_xml, graphml_keys, defaults) -> None: ...
    def add_edge(self, G: Graph[_Node], edge_element, graphml_keys) -> None: ...
    def decode_data_elements(self, graphml_keys, obj_xml): ...
    def find_graphml_keys(self, graph_element): ...
