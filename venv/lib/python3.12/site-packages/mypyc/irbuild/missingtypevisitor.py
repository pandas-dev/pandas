from __future__ import annotations

from mypy.nodes import Expression, Node
from mypy.traverser import ExtendedTraverserVisitor
from mypy.types import AnyType, Type, TypeOfAny


class MissingTypesVisitor(ExtendedTraverserVisitor):
    """AST visitor that can be used to add any missing types as a generic AnyType."""

    def __init__(self, types: dict[Expression, Type]) -> None:
        super().__init__()
        self.types: dict[Expression, Type] = types

    def visit(self, o: Node) -> bool:
        if isinstance(o, Expression) and o not in self.types:
            self.types[o] = AnyType(TypeOfAny.special_form)

        # If returns True, will continue to nested nodes.
        return True
