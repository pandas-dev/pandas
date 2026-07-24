from __future__ import annotations

from typing import TYPE_CHECKING, Final

if TYPE_CHECKING:
    from mypy.fixup import NodeFixer
    from mypy.nodes import MypyFile

# This is global mutable state. Don't add anything here unless there's a very
# good reason. This exists as a separate file to avoid method-level import in
# hot code in SymbolTableNode.node().


class ModulesState:
    def __init__(self) -> None:
        self.node_fixer: NodeFixer | None = None
        self.modules: dict[str, MypyFile] = {}


modules_state: Final = ModulesState()
