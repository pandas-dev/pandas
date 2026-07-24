from _typeshed import Incomplete
from collections.abc import Generator
from typing import ClassVar

__all__ = ["generate_network_text", "write_network_text"]

class BaseGlyphs:
    @classmethod
    def as_dict(cls) -> dict[str, str]: ...

class AsciiBaseGlyphs(BaseGlyphs):
    empty: ClassVar[str]
    newtree_last: ClassVar[str]
    newtree_mid: ClassVar[str]
    endof_forest: ClassVar[str]
    within_forest: ClassVar[str]
    within_tree: ClassVar[str]

class AsciiDirectedGlyphs(AsciiBaseGlyphs):
    last: ClassVar[str]
    mid: ClassVar[str]
    backedge: ClassVar[str]
    vertical_edge: ClassVar[str]

class AsciiUndirectedGlyphs(AsciiBaseGlyphs):
    last: ClassVar[str]
    mid: ClassVar[str]
    backedge: ClassVar[str]
    vertical_edge: ClassVar[str]

class UtfBaseGlyphs(BaseGlyphs):
    empty: ClassVar[str]
    newtree_last: ClassVar[str]
    newtree_mid: ClassVar[str]
    endof_forest: ClassVar[str]
    within_forest: ClassVar[str]
    within_tree: ClassVar[str]

class UtfDirectedGlyphs(UtfBaseGlyphs):
    last: ClassVar[str]
    mid: ClassVar[str]
    backedge: ClassVar[str]
    vertical_edge: ClassVar[str]

class UtfUndirectedGlyphs(UtfBaseGlyphs):
    last: ClassVar[str]
    mid: ClassVar[str]
    backedge: ClassVar[str]
    vertical_edge: ClassVar[str]

def generate_network_text(
    graph, with_labels: bool = True, sources=None, max_depth=None, ascii_only: bool = False, vertical_chains: bool = False
) -> Generator[Incomplete, None, Incomplete]: ...
def write_network_text(
    graph,
    path=None,
    with_labels: bool = True,
    sources=None,
    max_depth=None,
    ascii_only: bool = False,
    end: str = "\n",
    vertical_chains=False,
) -> None: ...
