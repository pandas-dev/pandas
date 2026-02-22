"""
Expression evaluation for Dataset.eval().

This module provides AST-based expression evaluation to support N-dimensional
arrays (N > 2), which pd.eval() doesn't support. See GitHub issue #11062.

We retain logical operator transformation ('and'/'or'/'not' to '&'/'|'/'~',
and chained comparisons) for consistency with query(), which still uses
pd.eval(). We don't migrate query() to this implementation because:
- query() typically works fine (expressions usually compare 1D coordinates)
- pd.eval() with numexpr is faster and well-tested for query's use case
"""

from __future__ import annotations

import ast
import builtins
from typing import Any

# Base namespace for eval expressions.
# We add common builtins back since we use an empty __builtins__ dict.
EVAL_BUILTINS: dict[str, Any] = {
    # Numeric/aggregation functions
    "abs": abs,
    "min": min,
    "max": max,
    "round": round,
    "len": len,
    "sum": sum,
    "pow": pow,
    "any": any,
    "all": all,
    # Type constructors
    "int": int,
    "float": float,
    "bool": bool,
    "str": str,
    "list": list,
    "tuple": tuple,
    "dict": dict,
    "set": set,
    "slice": slice,
    # Iteration helpers
    "range": range,
    "zip": zip,
    "enumerate": enumerate,
    "map": builtins.map,
    "filter": filter,
}


class LogicalOperatorTransformer(ast.NodeTransformer):
    """Transform operators for consistency with query().

    query() uses pd.eval() which transforms these operators automatically.
    We replicate that behavior here so syntax that works in query() also
    works in eval().

    Transformations:
    1. 'and'/'or'/'not' -> '&'/'|'/'~'
    2. 'a < b < c' -> '(a < b) & (b < c)'

    These constructs fail on arrays in standard Python because they call
    __bool__(), which is ambiguous for multi-element arrays.
    """

    def visit_BoolOp(self, node: ast.BoolOp) -> ast.AST:
        # Transform: a and b -> a & b, a or b -> a | b
        self.generic_visit(node)
        op: ast.BitAnd | ast.BitOr
        if isinstance(node.op, ast.And):
            op = ast.BitAnd()
        elif isinstance(node.op, ast.Or):
            op = ast.BitOr()
        else:
            return node

        # BoolOp can have multiple values: a and b and c
        # Transform to chained BinOp: (a & b) & c
        result = node.values[0]
        for value in node.values[1:]:
            result = ast.BinOp(left=result, op=op, right=value)
        return ast.fix_missing_locations(result)

    def visit_UnaryOp(self, node: ast.UnaryOp) -> ast.AST:
        # Transform: not a -> ~a
        self.generic_visit(node)
        if isinstance(node.op, ast.Not):
            return ast.fix_missing_locations(
                ast.UnaryOp(op=ast.Invert(), operand=node.operand)
            )
        return node

    def visit_Compare(self, node: ast.Compare) -> ast.AST:
        # Transform chained comparisons: 1 < x < 5 -> (1 < x) & (x < 5)
        # Python's chained comparisons use short-circuit evaluation at runtime,
        # which calls __bool__ on intermediate results. This fails for arrays.
        # We transform to bitwise AND which works element-wise.
        self.generic_visit(node)

        if len(node.ops) == 1:
            # Simple comparison, no transformation needed
            return node

        # Build individual comparisons and chain with BitAnd
        # For: a < b < c < d
        # We need: (a < b) & (b < c) & (c < d)
        comparisons = []
        left = node.left
        for op, comparator in zip(node.ops, node.comparators, strict=True):
            comp = ast.Compare(left=left, ops=[op], comparators=[comparator])
            comparisons.append(comp)
            left = comparator

        # Chain with BitAnd: (a < b) & (b < c) & ...
        result: ast.Compare | ast.BinOp = comparisons[0]
        for comp in comparisons[1:]:
            result = ast.BinOp(left=result, op=ast.BitAnd(), right=comp)
        return ast.fix_missing_locations(result)


def validate_expression(tree: ast.AST) -> None:
    """Validate that an AST doesn't contain patterns we don't support.

    These restrictions emulate pd.eval() behavior for consistency.
    """
    for node in ast.walk(tree):
        # Block lambda expressions (pd.eval: "Only named functions are supported")
        if isinstance(node, ast.Lambda):
            raise ValueError(
                "Lambda expressions are not allowed in eval(). "
                "Use direct operations on data variables instead."
            )
        # Block private/dunder attributes (consistent with pd.eval restrictions)
        if isinstance(node, ast.Attribute) and node.attr.startswith("_"):
            raise ValueError(
                f"Access to private attributes is not allowed: '{node.attr}'"
            )
