from _typeshed import Incomplete

from antlr4.tree.ParseTreePattern import ParseTreePattern as ParseTreePattern
from antlr4.tree.Tree import ParseTree as ParseTree

class ParseTreeMatch:
    __slots__ = ("tree", "pattern", "labels", "mismatchedNode")
    tree: Incomplete
    pattern: Incomplete
    labels: Incomplete
    mismatchedNode: Incomplete
    def __init__(
        self, tree: ParseTree, pattern: ParseTreePattern, labels: dict[str, list[ParseTree]], mismatchedNode: ParseTree
    ) -> None: ...
    def get(self, label: str): ...
    def getAll(self, label: str): ...
    def succeeded(self): ...
