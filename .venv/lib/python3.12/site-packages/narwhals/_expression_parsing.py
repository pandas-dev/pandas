# Utilities for expression parsing
# Useful for backends which don't have any concept of expressions, such
# and pandas or PyArrow.
# ! Any change to this module will trigger the pyspark and pyspark-connect tests in CI
from __future__ import annotations

from enum import Enum, auto
from typing import TYPE_CHECKING, Any, Callable, Literal, cast

from narwhals._utils import zip_strict
from narwhals.dependencies import is_numpy_array_1d
from narwhals.exceptions import (
    InvalidIntoExprError,
    InvalidOperationError,
    MultiOutputExpressionError,
)

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence

    from typing_extensions import Never, TypeIs

    from narwhals._compliant import CompliantExpr, CompliantFrameT
    from narwhals._compliant.typing import (
        AliasNames,
        CompliantExprAny,
        CompliantFrameAny,
        CompliantNamespaceAny,
        EvalNames,
    )
    from narwhals.expr import Expr
    from narwhals.series import Series
    from narwhals.typing import IntoExpr, NonNestedLiteral, _1DArray


def is_expr(obj: Any) -> TypeIs[Expr]:
    """Check whether `obj` is a Narwhals Expr."""
    from narwhals.expr import Expr

    return isinstance(obj, Expr)


def is_series(obj: Any) -> TypeIs[Series[Any]]:
    """Check whether `obj` is a Narwhals Expr."""
    from narwhals.series import Series

    return isinstance(obj, Series)


def combine_evaluate_output_names(
    *exprs: CompliantExpr[CompliantFrameT, Any],
) -> EvalNames[CompliantFrameT]:
    # Follow left-hand-rule for naming. E.g. `nw.sum_horizontal(expr1, expr2)` takes the
    # first name of `expr1`.
    def evaluate_output_names(df: CompliantFrameT) -> Sequence[str]:
        return exprs[0]._evaluate_output_names(df)[:1]

    return evaluate_output_names


def combine_alias_output_names(*exprs: CompliantExprAny) -> AliasNames | None:
    # Follow left-hand-rule for naming. E.g. `nw.sum_horizontal(expr1.alias(alias), expr2)` takes the
    # aliasing function of `expr1` and apply it to the first output name of `expr1`.
    if exprs[0]._alias_output_names is None:
        return None

    def alias_output_names(names: Sequence[str]) -> Sequence[str]:
        return exprs[0]._alias_output_names(names)[:1]  # type: ignore[misc]

    return alias_output_names


def evaluate_output_names_and_aliases(
    expr: CompliantExprAny, df: CompliantFrameAny, exclude: Sequence[str]
) -> tuple[Sequence[str], Sequence[str]]:
    output_names = expr._evaluate_output_names(df)
    aliases = (
        output_names
        if expr._alias_output_names is None
        else expr._alias_output_names(output_names)
    )
    if exclude and expr._metadata.expansion_kind.is_multi_unnamed():
        output_names, aliases = zip_strict(
            *[
                (x, alias)
                for x, alias in zip_strict(output_names, aliases)
                if x not in exclude
            ]
        )
    return output_names, aliases


class ExprKind(Enum):
    """Describe which kind of expression we are dealing with."""

    LITERAL = auto()
    """e.g. `nw.lit(1)`"""

    AGGREGATION = auto()
    """Reduces to a single value, not affected by row order, e.g. `nw.col('a').mean()`"""

    ORDERABLE_AGGREGATION = auto()
    """Reduces to a single value, affected by row order, e.g. `nw.col('a').arg_max()`"""

    ELEMENTWISE = auto()
    """Preserves length, can operate without context for surrounding rows, e.g. `nw.col('a').abs()`."""

    ORDERABLE_WINDOW = auto()
    """Depends on the rows around it and on their order, e.g. `diff`."""

    WINDOW = auto()
    """Depends on the rows around it and possibly their order, e.g. `rank`."""

    FILTRATION = auto()
    """Changes length, not affected by row order, e.g. `drop_nulls`."""

    ORDERABLE_FILTRATION = auto()
    """Changes length, affected by row order, e.g. `tail`."""

    OVER = auto()
    """Results from calling `.over` on expression."""

    COL = auto()
    """Results from calling `nw.col`."""

    NTH = auto()
    """Results from calling `nw.nth`."""

    EXCLUDE = auto()
    """Results from calling `nw.exclude`."""

    ALL = auto()
    """Results from calling `nw.all`."""

    SELECTOR = auto()
    """Results from creating an expression with a selector."""

    WHEN_THEN = auto()
    """Results from `when/then expression`, possibly followed by `otherwise`."""

    SERIES = auto()
    """Results from converting a Series to Expr."""

    @property
    def is_orderable(self) -> bool:
        # Any operation which may be affected by `order_by`, such as `cum_sum`,
        # `diff`, `rank`, `arg_max`, ...
        return self in {
            ExprKind.ORDERABLE_WINDOW,
            ExprKind.WINDOW,
            ExprKind.ORDERABLE_AGGREGATION,
            ExprKind.ORDERABLE_FILTRATION,
        }

    @property
    def is_elementwise(self) -> bool:
        # Any operation which can operate on each row independently
        # of the rows around it, e.g. `abs(), __add__, sum_horizontal, ...`
        return self in {
            ExprKind.ALL,
            ExprKind.COL,
            ExprKind.ELEMENTWISE,
            ExprKind.EXCLUDE,
            ExprKind.LITERAL,
            ExprKind.NTH,
            ExprKind.SELECTOR,
            ExprKind.SERIES,
            ExprKind.WHEN_THEN,
        }

    @property
    def is_scalar_like(self) -> bool:
        return self in {
            ExprKind.AGGREGATION,
            ExprKind.LITERAL,
            ExprKind.ORDERABLE_AGGREGATION,
        }


def is_scalar_like(obj: CompliantExprAny) -> bool:
    return obj._metadata.is_scalar_like


class ExpansionKind(Enum):
    """Describe what kind of expansion the expression performs."""

    SINGLE = auto()
    """e.g. `nw.col('a'), nw.sum_horizontal(nw.all())`"""

    MULTI_NAMED = auto()
    """e.g. `nw.col('a', 'b')`"""

    MULTI_UNNAMED = auto()
    """e.g. `nw.all()`, nw.nth(0, 1)"""

    def is_multi_unnamed(self) -> bool:
        return self is ExpansionKind.MULTI_UNNAMED

    def is_multi_output(self) -> bool:
        return self in {ExpansionKind.MULTI_NAMED, ExpansionKind.MULTI_UNNAMED}

    def __and__(self, other: ExpansionKind) -> Literal[ExpansionKind.MULTI_UNNAMED]:
        if self is ExpansionKind.MULTI_UNNAMED and other is ExpansionKind.MULTI_UNNAMED:
            # e.g. nw.selectors.all() - nw.selectors.numeric().
            return ExpansionKind.MULTI_UNNAMED
        # Don't attempt anything more complex, keep it simple and raise in the face of ambiguity.
        msg = f"Unsupported ExpansionKind combination, got {self} and {other}, please report a bug."  # pragma: no cover
        raise AssertionError(msg)  # pragma: no cover


class ExprNode:
    """An operation to create or modify an expression.

    Parameters:
        kind: ExprKind of operation.
        name: Name of function, as defined in the compliant protocols.
        exprs: Expressifiable arguments to function.
        str_as_lit: Whether to interpret strings as literals when they
            are present in `exprs`.
        allow_multi_output: Whether to allow any of `exprs` to be multi-output.
        kwargs: Other (non-expressifiable) arguments to function.
    """

    def __init__(
        self,
        kind: ExprKind,
        name: str,
        /,
        *exprs: IntoExpr | NonNestedLiteral,
        str_as_lit: bool = False,
        allow_multi_output: bool = False,
        **kwargs: Any,
    ) -> None:
        self.kind: ExprKind = kind
        self.name: str = name
        self.exprs: Sequence[IntoExpr | NonNestedLiteral] = exprs
        self.kwargs: dict[str, Any] = kwargs
        self.str_as_lit: bool = str_as_lit
        self.allow_multi_output: bool = allow_multi_output

        # Cached methods.
        self._is_orderable_cached: bool | None = None
        self._is_elementwise_cached: bool | None = None

    def __repr__(self) -> str:
        if self.name == "col":
            names = ", ".join(str(x) for x in self.kwargs["names"])
            return f"col({names})"
        arg_str = []
        expr_repr = ", ".join(str(x) for x in self.exprs)
        kwargs_repr = ", ".join(f"{key}={value}" for key, value in self.kwargs.items())
        if self.exprs:
            arg_str.append(expr_repr)
        if self.kwargs:
            arg_str.append(kwargs_repr)
        return f"{self.name}({', '.join(arg_str)})"

    def as_dict(self) -> dict[str, Any]:  # pragma: no cover
        # Just for debugging.
        return {
            "kind": self.kind,
            "name": self.name,
            "exprs": self.exprs,
            "kwargs": self.kwargs,
            "str_as_lit": self.str_as_lit,
            "allow_multi_output": self.allow_multi_output,
        }

    def _with_kwargs(self, **kwargs: Any) -> ExprNode:
        return self.__class__(
            self.kind, self.name, *self.exprs, str_as_lit=self.str_as_lit, **kwargs
        )

    def _push_down_over_node_in_place(
        self, over_node: ExprNode, over_node_without_order_by: ExprNode
    ) -> None:
        exprs: list[IntoExpr | NonNestedLiteral] = []
        # Note: please keep this as a for-loop (rather than a list-comprehension)
        # so that pytest-cov highlights any uncovered branches.
        over_node_order_by = over_node.kwargs["order_by"]
        over_node_partition_by = over_node.kwargs["partition_by"]
        for expr in self.exprs:
            if not is_expr(expr):
                exprs.append(expr)
            elif over_node_order_by and any(
                expr_node.is_orderable() for expr_node in expr._nodes
            ):
                exprs.append(expr._with_over_node(over_node))
            elif over_node_partition_by and not all(
                expr_node.is_elementwise() for expr_node in expr._nodes
            ):
                exprs.append(expr._with_over_node(over_node_without_order_by))
            else:
                # If there's no `partition_by`, then `over_node_without_order_by` is a no-op.
                exprs.append(expr)
        self.exprs = exprs

    def is_orderable(self) -> bool:
        if self._is_orderable_cached is None:
            # Note: don't combine these if/then statements so that pytest-cov shows if
            # anything is uncovered.
            if self.kind.is_orderable:  # noqa: SIM114
                self._is_orderable_cached = True
            elif any(
                any(node.is_orderable() for node in expr._nodes)
                for expr in self.exprs
                if is_expr(expr)
            ):
                self._is_orderable_cached = True
            else:
                self._is_orderable_cached = False
        return self._is_orderable_cached

    def is_elementwise(self) -> bool:
        if self._is_elementwise_cached is None:
            # Note: don't combine these if/then statements so that pytest-cov shows if
            # anything is uncovered.
            if not self.kind.is_elementwise:  # noqa: SIM114
                self._is_elementwise_cached = False
            elif any(
                any(not node.is_elementwise() for node in expr._nodes)
                for expr in self.exprs
                if is_expr(expr)
            ):
                self._is_elementwise_cached = False
            else:
                self._is_elementwise_cached = True
        return self._is_elementwise_cached


class ExprMetadata:
    """Expression metadata.

    Parameters:
        expansion_kind: What kind of expansion the expression performs.
        has_windows: Whether it already contains window functions.
        is_elementwise: Whether it can operate row-by-row without context
            of the other rows around it.
        is_literal: Whether it is just a literal wrapped in an expression.
        is_scalar_like: Whether it is a literal or an aggregation.
        n_orderable_ops: The number of order-dependent operations. In the
            lazy case, this number must be `0` by the time the expression
            is evaluated.
        preserves_length: Whether the expression preserves the input length.
        current_node: The current ExprNode in the linked list.
        prev: Reference to the previous ExprMetadata in the linked list (None for root).
    """

    __slots__ = (
        "current_node",
        "expansion_kind",
        "has_windows",
        "is_elementwise",
        "is_literal",
        "is_scalar_like",
        "n_orderable_ops",
        "preserves_length",
        "prev",
    )

    def __init__(
        self,
        expansion_kind: ExpansionKind,
        *,
        has_windows: bool = False,
        n_orderable_ops: int = 0,
        preserves_length: bool = True,
        is_elementwise: bool = True,
        is_scalar_like: bool = False,
        is_literal: bool = False,
        current_node: ExprNode,
        prev: ExprMetadata | None = None,
    ) -> None:
        if is_literal:
            assert is_scalar_like  # noqa: S101  # debug assertion
        self.expansion_kind: ExpansionKind = expansion_kind
        self.has_windows: bool = has_windows
        self.n_orderable_ops: int = n_orderable_ops
        self.is_elementwise: bool = is_elementwise
        self.preserves_length: bool = preserves_length
        self.is_scalar_like: bool = is_scalar_like
        self.is_literal: bool = is_literal
        self.current_node: ExprNode = current_node
        self.prev: ExprMetadata | None = prev

    def __init_subclass__(cls, /, *args: Any, **kwds: Any) -> Never:  # pragma: no cover
        msg = f"Cannot subclass {cls.__name__!r}"
        raise TypeError(msg)

    def __repr__(self) -> str:  # pragma: no cover
        nodes = tuple(reversed(tuple(self.iter_nodes_reversed())))
        return (
            f"ExprMetadata(\n"
            f"  expansion_kind: {self.expansion_kind},\n"
            f"  has_windows: {self.has_windows},\n"
            f"  n_orderable_ops: {self.n_orderable_ops},\n"
            f"  is_elementwise: {self.is_elementwise},\n"
            f"  preserves_length: {self.preserves_length},\n"
            f"  is_scalar_like: {self.is_scalar_like},\n"
            f"  is_literal: {self.is_literal},\n"
            f"  nodes: {nodes},\n"
            ")"
        )

    def iter_nodes_reversed(self) -> Iterator[ExprNode]:
        """Iterate through all nodes from current to root."""
        current: ExprMetadata | None = self
        while current is not None:
            yield current.current_node
            current = current.prev

    @classmethod
    def from_node(
        cls, node: ExprNode, *compliant_exprs: CompliantExprAny
    ) -> ExprMetadata:
        return KIND_TO_METADATA_CONSTRUCTOR[node.kind](node, *compliant_exprs)

    def with_node(
        self,
        node: ExprNode,
        compliant_expr: CompliantExprAny,
        *compliant_expr_args: CompliantExprAny,
    ) -> ExprMetadata:
        return KIND_TO_METADATA_UPDATER[node.kind](
            self, node, compliant_expr, *compliant_expr_args
        )

    @classmethod
    def from_aggregation(cls, node: ExprNode) -> ExprMetadata:
        return cls(
            ExpansionKind.SINGLE,
            is_elementwise=False,
            preserves_length=False,
            is_scalar_like=True,
            current_node=node,
            prev=None,
        )

    @classmethod
    def from_literal(cls, node: ExprNode) -> ExprMetadata:
        return cls(
            ExpansionKind.SINGLE,
            is_elementwise=True,
            preserves_length=False,
            is_literal=True,
            is_scalar_like=True,
            current_node=node,
            prev=None,
        )

    @classmethod
    def from_series(cls, node: ExprNode) -> ExprMetadata:
        return cls(ExpansionKind.SINGLE, current_node=node, prev=None)

    @classmethod
    def from_col(cls, node: ExprNode) -> ExprMetadata:
        # e.g. `nw.col('a')`, `nw.nth(0)`
        return (
            cls(ExpansionKind.SINGLE, current_node=node, prev=None)
            if len(node.kwargs["names"]) == 1
            else cls.from_selector_multi_named(node)
        )

    @classmethod
    def from_nth(cls, node: ExprNode) -> ExprMetadata:
        return (
            cls(ExpansionKind.SINGLE, current_node=node, prev=None)
            if len(node.kwargs["indices"]) == 1
            else cls.from_selector_multi_named(node)
        )

    @classmethod
    def from_selector_multi_named(cls, node: ExprNode) -> ExprMetadata:
        # e.g. `nw.col('a', 'b')`
        return cls(ExpansionKind.MULTI_NAMED, current_node=node, prev=None)

    @classmethod
    def from_selector_multi_unnamed(cls, node: ExprNode) -> ExprMetadata:
        # e.g. `nw.all()`
        return cls(ExpansionKind.MULTI_UNNAMED, current_node=node, prev=None)

    @classmethod
    def from_elementwise(
        cls, node: ExprNode, *compliant_exprs: CompliantExprAny
    ) -> ExprMetadata:
        return combine_metadata(
            *compliant_exprs, to_single_output=True, current_node=node, prev=None
        )

    @property
    def is_filtration(self) -> bool:
        return not self.preserves_length and not self.is_scalar_like

    def with_aggregation(self, node: ExprNode, _ce: CompliantExprAny) -> ExprMetadata:
        if self.is_scalar_like:
            msg = "Can't apply aggregations to scalar-like expressions."
            raise InvalidOperationError(msg)
        return ExprMetadata(
            self.expansion_kind,
            has_windows=self.has_windows,
            n_orderable_ops=self.n_orderable_ops,
            preserves_length=False,
            is_elementwise=False,
            is_scalar_like=True,
            is_literal=False,
            current_node=node,
            prev=self,
        )

    def with_orderable_aggregation(
        self, node: ExprNode, _ce: CompliantExprAny
    ) -> ExprMetadata:
        # Deprecated, used only in stable.v1.
        if self.is_scalar_like:  # pragma: no cover
            msg = "Can't apply aggregations to scalar-like expressions."
            raise InvalidOperationError(msg)
        return ExprMetadata(
            self.expansion_kind,
            has_windows=self.has_windows,
            n_orderable_ops=self.n_orderable_ops + 1,
            preserves_length=False,
            is_elementwise=False,
            is_scalar_like=True,
            is_literal=False,
            current_node=node,
            prev=self,
        )

    def with_elementwise(
        self,
        node: ExprNode,
        compliant_expr: CompliantExprAny,
        *compliant_expr_args: CompliantExprAny,
    ) -> ExprMetadata:
        return combine_metadata(
            compliant_expr,
            *compliant_expr_args,
            to_single_output=False,
            current_node=node,
            prev=compliant_expr._metadata,
        )

    def with_window(self, node: ExprNode, _ce: CompliantExprAny) -> ExprMetadata:
        # Window function which may (but doesn't have to) be used with `over(order_by=...)`.
        if self.is_scalar_like:
            msg = "Can't apply window (e.g. `rank`) to scalar-like expression."
            raise InvalidOperationError(msg)
        return ExprMetadata(
            self.expansion_kind,
            has_windows=self.has_windows,
            # The function isn't order-dependent (but, users can still use `order_by` if they wish!),
            # so we don't increment `n_orderable_ops`.
            n_orderable_ops=self.n_orderable_ops,
            preserves_length=self.preserves_length,
            is_elementwise=False,
            is_scalar_like=False,
            is_literal=False,
            current_node=node,
            prev=self,
        )

    def with_orderable_window(
        self, node: ExprNode, _ce: CompliantExprAny
    ) -> ExprMetadata:
        # Window function which must be used with `over(order_by=...)`.
        if self.is_scalar_like:
            msg = "Can't apply orderable window (e.g. `diff`, `shift`) to scalar-like expression."
            raise InvalidOperationError(msg)
        return ExprMetadata(
            self.expansion_kind,
            has_windows=self.has_windows,
            n_orderable_ops=self.n_orderable_ops + 1,
            preserves_length=self.preserves_length,
            is_elementwise=False,
            is_scalar_like=False,
            is_literal=False,
            current_node=node,
            prev=self,
        )

    def with_ordered_over(self, node: ExprNode, _ce: CompliantExprAny) -> ExprMetadata:
        if self.has_windows:
            msg = "Cannot nest `over` statements."
            raise InvalidOperationError(msg)
        if self.is_elementwise or self.is_filtration:
            msg = (
                "Cannot use `over` on expressions which are elementwise\n"
                "(e.g. `abs`) or which change length (e.g. `drop_nulls`)."
            )
            raise InvalidOperationError(msg)
        n_orderable_ops = self.n_orderable_ops
        if (
            not n_orderable_ops
            and next(self.op_nodes_reversed()).kind is not ExprKind.WINDOW
        ):
            msg = (
                "Cannot use `order_by` in `over` on expression which isn't orderable.\n"
                "If your expression is orderable, then make sure that `over(order_by=...)`\n"
                "comes immediately after the order-dependent expression.\n\n"
                "Hint: instead of\n"
                "  - `(nw.col('price').diff() + 1).over(order_by='date')`\n"
                "write:\n"
                "  + `nw.col('price').diff().over(order_by='date') + 1`\n"
            )
            raise InvalidOperationError(msg)
        if next(self.op_nodes_reversed()).kind.is_orderable and n_orderable_ops > 0:
            n_orderable_ops -= 1
        return ExprMetadata(
            self.expansion_kind,
            has_windows=True,
            n_orderable_ops=n_orderable_ops,
            preserves_length=True,
            is_elementwise=False,
            is_scalar_like=False,
            is_literal=False,
            current_node=node,
            prev=self,
        )

    def with_partitioned_over(
        self, node: ExprNode, _ce: CompliantExprAny
    ) -> ExprMetadata:
        if self.has_windows:
            msg = "Cannot nest `over` statements."
            raise InvalidOperationError(msg)
        if self.is_elementwise or self.is_filtration:
            msg = (
                "Cannot use `over` on expressions which are elementwise\n"
                "(e.g. `abs`) or which change length (e.g. `drop_nulls`)."
            )
            raise InvalidOperationError(msg)
        return ExprMetadata(
            self.expansion_kind,
            has_windows=True,
            n_orderable_ops=self.n_orderable_ops,
            preserves_length=True,
            is_elementwise=False,
            is_scalar_like=False,
            is_literal=False,
            current_node=node,
            prev=self,
        )

    def with_over(self, node: ExprNode, _ce: CompliantExprAny) -> ExprMetadata:
        if node.kwargs["order_by"]:
            return self.with_ordered_over(node, _ce)
        if not node.kwargs["partition_by"]:  # pragma: no cover
            msg = "At least one of `partition_by` or `order_by` must be specified."
            raise InvalidOperationError(msg)
        return self.with_partitioned_over(node, _ce)

    def with_filtration(
        self, node: ExprNode, *compliant_exprs: CompliantExprAny
    ) -> ExprMetadata:
        if self.is_scalar_like:
            msg = "Can't apply filtration (e.g. `drop_nulls`) to scalar-like expression."
            raise InvalidOperationError(msg)
        result_has_windows = any(x._metadata.has_windows for x in compliant_exprs)
        result_n_orderable_ops = sum(x._metadata.n_orderable_ops for x in compliant_exprs)
        return ExprMetadata(
            self.expansion_kind,
            has_windows=result_has_windows,
            n_orderable_ops=result_n_orderable_ops,
            preserves_length=False,
            is_elementwise=False,
            is_scalar_like=False,
            is_literal=False,
            current_node=node,
            prev=self,
        )

    def with_orderable_filtration(
        self, node: ExprNode, _ce: CompliantExprAny
    ) -> ExprMetadata:
        if self.is_scalar_like:
            msg = "Can't apply filtration (e.g. `drop_nulls`) to scalar-like expression."
            raise InvalidOperationError(msg)
        return ExprMetadata(
            self.expansion_kind,
            has_windows=self.has_windows,
            n_orderable_ops=self.n_orderable_ops + 1,
            preserves_length=False,
            is_elementwise=False,
            is_scalar_like=False,
            is_literal=False,
            current_node=node,
            prev=self,
        )

    def op_nodes_reversed(self) -> Iterator[ExprNode]:
        for node in self.iter_nodes_reversed():
            if node.name.startswith(("name.", "alias")):
                # Skip nodes which only do aliasing.
                continue
            yield node


KIND_TO_METADATA_CONSTRUCTOR: dict[ExprKind, Callable[[ExprNode], ExprMetadata]] = {
    ExprKind.AGGREGATION: ExprMetadata.from_aggregation,
    ExprKind.ALL: ExprMetadata.from_selector_multi_unnamed,
    ExprKind.ELEMENTWISE: ExprMetadata.from_elementwise,
    ExprKind.EXCLUDE: ExprMetadata.from_selector_multi_unnamed,
    ExprKind.SERIES: ExprMetadata.from_series,
    ExprKind.COL: ExprMetadata.from_col,
    ExprKind.LITERAL: ExprMetadata.from_literal,
    ExprKind.NTH: ExprMetadata.from_nth,
    ExprKind.SELECTOR: ExprMetadata.from_selector_multi_unnamed,
}

KIND_TO_METADATA_UPDATER: dict[ExprKind, Callable[..., ExprMetadata]] = {
    ExprKind.AGGREGATION: ExprMetadata.with_aggregation,
    ExprKind.ELEMENTWISE: ExprMetadata.with_elementwise,
    ExprKind.FILTRATION: ExprMetadata.with_filtration,
    ExprKind.ORDERABLE_AGGREGATION: ExprMetadata.with_orderable_aggregation,
    ExprKind.ORDERABLE_FILTRATION: ExprMetadata.with_orderable_filtration,
    ExprKind.OVER: ExprMetadata.with_over,
    ExprKind.ORDERABLE_WINDOW: ExprMetadata.with_orderable_window,
    ExprKind.WINDOW: ExprMetadata.with_window,
}


def combine_metadata(
    *compliant_exprs: CompliantExprAny,
    to_single_output: bool,
    current_node: ExprNode,
    prev: ExprMetadata | None,
) -> ExprMetadata:
    """Combine metadata from `args`.

    Arguments:
        compliant_exprs: Expression arguments.
        to_single_output: Whether the result is always single-output, regardless
            of the inputs (e.g. `nw.sum_horizontal`).
        current_node: The current node being added.
        prev: ExprMetadata of previous node.
    """
    n_filtrations = 0
    result_expansion_kind = ExpansionKind.SINGLE
    result_has_windows = False
    result_n_orderable_ops = 0
    # result preserves length if at least one input does
    result_preserves_length = False
    # result is elementwise if all inputs are elementwise
    result_is_elementwise = True
    # result is scalar-like if all inputs are scalar-like
    result_is_scalar_like = True
    # result is literal if all inputs are literal
    result_is_literal = True

    for i, ce in enumerate(compliant_exprs):
        metadata = ce._metadata
        assert metadata is not None  # noqa: S101
        if metadata.expansion_kind.is_multi_output():
            expansion_kind = metadata.expansion_kind
            if not to_single_output:
                result_expansion_kind = (
                    result_expansion_kind & expansion_kind if i > 0 else expansion_kind
                )

        result_has_windows |= metadata.has_windows
        result_n_orderable_ops += metadata.n_orderable_ops
        result_preserves_length |= metadata.preserves_length
        result_is_elementwise &= metadata.is_elementwise
        result_is_scalar_like &= metadata.is_scalar_like
        result_is_literal &= metadata.is_literal
        n_filtrations += int(metadata.is_filtration)
    if n_filtrations > 1:
        msg = "Length-changing expressions can only be used in isolation, or followed by an aggregation"
        raise InvalidOperationError(msg)
    if result_preserves_length and n_filtrations:
        msg = "Cannot combine length-changing expressions with length-preserving ones or aggregations"
        raise InvalidOperationError(msg)
    return ExprMetadata(
        result_expansion_kind,
        has_windows=result_has_windows,
        n_orderable_ops=result_n_orderable_ops,
        preserves_length=result_preserves_length,
        is_elementwise=result_is_elementwise,
        is_scalar_like=result_is_scalar_like,
        is_literal=result_is_literal,
        current_node=current_node,
        prev=prev,
    )


def check_expressions_preserve_length(
    *args: CompliantExprAny, function_name: str
) -> None:
    # Raise if any argument in `args` isn't length-preserving.
    # For Series input, we don't raise (yet), we let such checks happen later,
    # as this function works lazily and so can't evaluate lengths.

    if not all(x._metadata.preserves_length for x in args):
        msg = f"Expressions which aggregate or change length cannot be passed to '{function_name}'."
        raise InvalidOperationError(msg)


def _parse_into_expr(
    arg: IntoExpr | NonNestedLiteral | _1DArray,
    *,
    str_as_lit: bool = False,
    backend: Any = None,
    allow_literal: bool = True,
) -> Expr:
    from narwhals.functions import col, lit, new_series

    if isinstance(arg, str) and not str_as_lit:
        return col(arg)
    if is_numpy_array_1d(arg):
        return new_series("", arg, backend=backend)._to_expr()
    if is_series(arg):
        return arg._to_expr()
    if is_expr(arg):
        return arg
    if not allow_literal:
        raise InvalidIntoExprError.from_invalid_type(type(arg))
    return lit(arg)


def evaluate_into_exprs(
    *exprs: IntoExpr | NonNestedLiteral | _1DArray,
    ns: CompliantNamespaceAny,
    str_as_lit: bool,
    allow_multi_output: bool,
) -> Iterator[CompliantExprAny]:
    for expr in exprs:
        ret = _parse_into_expr(
            expr, str_as_lit=str_as_lit, backend=ns._implementation
        )._to_compliant_expr(ns)
        if not allow_multi_output and ret._metadata.expansion_kind.is_multi_output():
            msg = "Multi-output expressions are not allowed in this context."
            raise MultiOutputExpressionError(msg)
        yield ret


def maybe_broadcast_ces(*compliant_exprs: CompliantExprAny) -> list[CompliantExprAny]:
    broadcast = any(not is_scalar_like(ce) for ce in compliant_exprs)
    results: list[CompliantExprAny] = []
    for compliant_expr in compliant_exprs:
        if broadcast and is_scalar_like(compliant_expr):
            _compliant_expr: CompliantExprAny = compliant_expr.broadcast()
            # Make sure to preserve metadata.
            _compliant_expr._opt_metadata = compliant_expr._metadata
            results.append(_compliant_expr)
        else:
            results.append(compliant_expr)
    return results


def evaluate_root_node(node: ExprNode, ns: CompliantNamespaceAny) -> CompliantExprAny:
    if node.name in {"col", "exclude"}:
        # There's too much potential for Sequence[str] vs str bugs, so we pass down
        # `names` positionally rather than as a sequence of strings.
        ce = getattr(ns, node.name)(*node.kwargs["names"])
        ces = []
    else:
        if "." in node.name:
            module, method = node.name.split(".")
            func = getattr(getattr(ns, module), method)
        else:
            func = getattr(ns, node.name)
        ces = maybe_broadcast_ces(
            *evaluate_into_exprs(
                *node.exprs,
                ns=ns,
                str_as_lit=node.str_as_lit,
                allow_multi_output=node.allow_multi_output,
            )
        )
        ce = cast("CompliantExprAny", func(*ces, **node.kwargs))
    md = ExprMetadata.from_node(node, *ces)
    ce._opt_metadata = md
    return ce


def evaluate_node(
    compliant_expr: CompliantExprAny, node: ExprNode, ns: CompliantNamespaceAny
) -> CompliantExprAny:
    md: ExprMetadata = compliant_expr._metadata
    compliant_expr, *compliant_expr_args = maybe_broadcast_ces(
        compliant_expr,
        *evaluate_into_exprs(
            *node.exprs,
            ns=ns,
            str_as_lit=node.str_as_lit,
            allow_multi_output=node.allow_multi_output,
        ),
    )
    md = md.with_node(node, compliant_expr, *compliant_expr_args)
    if "." in node.name:
        accessor, method = node.name.split(".")
        func = getattr(getattr(compliant_expr, accessor), method)
    else:
        func = getattr(compliant_expr, node.name)
    ret = cast("CompliantExprAny", func(*compliant_expr_args, **node.kwargs))
    ret._opt_metadata = md
    return ret


def evaluate_nodes(
    nodes: Sequence[ExprNode], ns: CompliantNamespaceAny
) -> CompliantExprAny:
    ce = evaluate_root_node(nodes[0], ns)
    for node in nodes[1:]:
        ce = evaluate_node(ce, node, ns)
    return ce
