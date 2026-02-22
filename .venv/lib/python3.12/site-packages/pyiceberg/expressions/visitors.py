# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.
import math
from abc import ABC, abstractmethod
from collections.abc import Callable
from functools import singledispatch
from typing import (
    Any,
    Generic,
    SupportsFloat,
    TypeVar,
)

from pyiceberg.conversions import from_bytes
from pyiceberg.expressions import (
    AlwaysFalse,
    AlwaysTrue,
    And,
    BooleanExpression,
    BoundEqualTo,
    BoundGreaterThan,
    BoundGreaterThanOrEqual,
    BoundIn,
    BoundIsNaN,
    BoundIsNull,
    BoundLessThan,
    BoundLessThanOrEqual,
    BoundLiteralPredicate,
    BoundNotEqualTo,
    BoundNotIn,
    BoundNotNaN,
    BoundNotNull,
    BoundNotStartsWith,
    BoundPredicate,
    BoundSetPredicate,
    BoundStartsWith,
    BoundTerm,
    BoundUnaryPredicate,
    Not,
    Or,
    UnboundPredicate,
)
from pyiceberg.manifest import DataFile, ManifestFile, PartitionFieldSummary
from pyiceberg.partitioning import UNPARTITIONED_PARTITION_SPEC, PartitionSpec
from pyiceberg.schema import Schema
from pyiceberg.typedef import EMPTY_DICT, L, LiteralValue, Record, StructProtocol
from pyiceberg.types import (
    DoubleType,
    FloatType,
    IcebergType,
    NestedField,
    PrimitiveType,
    StructType,
    TimestampType,
    TimestamptzType,
)
from pyiceberg.utils.datetime import micros_to_timestamp, micros_to_timestamptz

T = TypeVar("T")


class BooleanExpressionVisitor(Generic[T], ABC):
    @abstractmethod
    def visit_true(self) -> T:
        """Visit method for an AlwaysTrue boolean expression.

        Note: This visit method has no arguments since AlwaysTrue instances have no context.
        """

    @abstractmethod
    def visit_false(self) -> T:
        """Visit method for an AlwaysFalse boolean expression.

        Note: This visit method has no arguments since AlwaysFalse instances have no context.
        """

    @abstractmethod
    def visit_not(self, child_result: T) -> T:
        """Visit method for a Not boolean expression.

        Args:
            child_result (T): The result of visiting the child of the Not boolean expression.
        """

    @abstractmethod
    def visit_and(self, left_result: T, right_result: T) -> T:
        """Visit method for an And boolean expression.

        Args:
            left_result (T): The result of visiting the left side of the expression.
            right_result (T): The result of visiting the right side of the expression.
        """

    @abstractmethod
    def visit_or(self, left_result: T, right_result: T) -> T:
        """Visit method for an Or boolean expression.

        Args:
            left_result (T): The result of visiting the left side of the expression.
            right_result (T): The result of visiting the right side of the expression.
        """

    @abstractmethod
    def visit_unbound_predicate(self, predicate: UnboundPredicate) -> T:
        """Visit method for an unbound predicate in an expression tree.

        Args:
            predicate (UnboundPredicate): An instance of an UnboundPredicate.
        """

    @abstractmethod
    def visit_bound_predicate(self, predicate: BoundPredicate) -> T:
        """Visit method for a bound predicate in an expression tree.

        Args:
            predicate (BoundPredicate): An instance of a BoundPredicate.
        """


@singledispatch
def visit(obj: BooleanExpression, visitor: BooleanExpressionVisitor[T]) -> T:
    """Apply a boolean expression visitor to any point within an expression.

    The function traverses the expression in post-order fashion.

    Args:
        obj (BooleanExpression): An instance of a BooleanExpression.
        visitor (BooleanExpressionVisitor[T]): An instance of an implementation of the generic
            BooleanExpressionVisitor base class.

    Raises:
        NotImplementedError: If attempting to visit an unsupported expression.
    """
    raise NotImplementedError(f"Cannot visit unsupported expression: {obj}")


@visit.register(AlwaysTrue)
def _(_: AlwaysTrue, visitor: BooleanExpressionVisitor[T]) -> T:
    """Visit an AlwaysTrue boolean expression with a concrete BooleanExpressionVisitor."""
    return visitor.visit_true()


@visit.register(AlwaysFalse)
def _(_: AlwaysFalse, visitor: BooleanExpressionVisitor[T]) -> T:
    """Visit an AlwaysFalse boolean expression with a concrete BooleanExpressionVisitor."""
    return visitor.visit_false()


@visit.register(Not)
def _(obj: Not, visitor: BooleanExpressionVisitor[T]) -> T:
    """Visit a Not boolean expression with a concrete BooleanExpressionVisitor."""
    child_result: T = visit(obj.child, visitor=visitor)
    return visitor.visit_not(child_result=child_result)


@visit.register(And)
def _(obj: And, visitor: BooleanExpressionVisitor[T]) -> T:
    """Visit an And boolean expression with a concrete BooleanExpressionVisitor."""
    left_result: T = visit(obj.left, visitor=visitor)
    right_result: T = visit(obj.right, visitor=visitor)
    return visitor.visit_and(left_result=left_result, right_result=right_result)


@visit.register(UnboundPredicate)
def _(obj: UnboundPredicate, visitor: BooleanExpressionVisitor[T]) -> T:
    """Visit an unbound boolean expression with a concrete BooleanExpressionVisitor."""
    return visitor.visit_unbound_predicate(predicate=obj)


@visit.register(BoundPredicate)
def _(obj: BoundPredicate, visitor: BooleanExpressionVisitor[T]) -> T:
    """Visit a bound boolean expression with a concrete BooleanExpressionVisitor."""
    return visitor.visit_bound_predicate(predicate=obj)


@visit.register(Or)
def _(obj: Or, visitor: BooleanExpressionVisitor[T]) -> T:
    """Visit an Or boolean expression with a concrete BooleanExpressionVisitor."""
    left_result: T = visit(obj.left, visitor=visitor)
    right_result: T = visit(obj.right, visitor=visitor)
    return visitor.visit_or(left_result=left_result, right_result=right_result)


def bind(schema: Schema, expression: BooleanExpression, case_sensitive: bool) -> BooleanExpression:
    """Travers over an expression to bind the predicates to the schema.

    Args:
      schema (Schema): A schema to use when binding the expression.
      expression (BooleanExpression): An expression containing UnboundPredicates that can be bound.
      case_sensitive (bool): Whether to consider case when binding a reference to a field in a schema, defaults to True.

    Raises:
        TypeError: In the case a predicate is already bound.
    """
    return visit(expression, BindVisitor(schema, case_sensitive))


class BindVisitor(BooleanExpressionVisitor[BooleanExpression]):
    """Rewrites a boolean expression by replacing unbound references with references to fields in a struct schema.

    Args:
      schema (Schema): A schema to use when binding the expression.
      case_sensitive (bool): Whether to consider case when binding a reference to a field in a schema, defaults to True.

    Raises:
        TypeError: In the case a predicate is already bound.
    """

    schema: Schema
    case_sensitive: bool

    def __init__(self, schema: Schema, case_sensitive: bool) -> None:
        self.schema = schema
        self.case_sensitive = case_sensitive

    def visit_true(self) -> BooleanExpression:
        return AlwaysTrue()

    def visit_false(self) -> BooleanExpression:
        return AlwaysFalse()

    def visit_not(self, child_result: BooleanExpression) -> BooleanExpression:
        return Not(child=child_result)

    def visit_and(self, left_result: BooleanExpression, right_result: BooleanExpression) -> BooleanExpression:
        return And(left=left_result, right=right_result)

    def visit_or(self, left_result: BooleanExpression, right_result: BooleanExpression) -> BooleanExpression:
        return Or(left=left_result, right=right_result)

    def visit_unbound_predicate(self, predicate: UnboundPredicate) -> BooleanExpression:
        return predicate.bind(self.schema, case_sensitive=self.case_sensitive)

    def visit_bound_predicate(self, predicate: BoundPredicate) -> BooleanExpression:
        raise TypeError(f"Found already bound predicate: {predicate}")


class BoundBooleanExpressionVisitor(BooleanExpressionVisitor[T], ABC):
    @abstractmethod
    def visit_in(self, term: BoundTerm, literals: set[L]) -> T:
        """Visit a bound In predicate."""

    @abstractmethod
    def visit_not_in(self, term: BoundTerm, literals: set[L]) -> T:
        """Visit a bound NotIn predicate."""

    @abstractmethod
    def visit_is_nan(self, term: BoundTerm) -> T:
        """Visit a bound IsNan predicate."""

    @abstractmethod
    def visit_not_nan(self, term: BoundTerm) -> T:
        """Visit a bound NotNan predicate."""

    @abstractmethod
    def visit_is_null(self, term: BoundTerm) -> T:
        """Visit a bound IsNull predicate."""

    @abstractmethod
    def visit_not_null(self, term: BoundTerm) -> T:
        """Visit a bound NotNull predicate."""

    @abstractmethod
    def visit_equal(self, term: BoundTerm, literal: LiteralValue) -> T:
        """Visit a bound Equal predicate."""

    @abstractmethod
    def visit_not_equal(self, term: BoundTerm, literal: LiteralValue) -> T:
        """Visit a bound NotEqual predicate."""

    @abstractmethod
    def visit_greater_than_or_equal(self, term: BoundTerm, literal: LiteralValue) -> T:
        """Visit a bound GreaterThanOrEqual predicate."""

    @abstractmethod
    def visit_greater_than(self, term: BoundTerm, literal: LiteralValue) -> T:
        """Visit a bound GreaterThan predicate."""

    @abstractmethod
    def visit_less_than(self, term: BoundTerm, literal: LiteralValue) -> T:
        """Visit a bound LessThan predicate."""

    @abstractmethod
    def visit_less_than_or_equal(self, term: BoundTerm, literal: LiteralValue) -> T:
        """Visit a bound LessThanOrEqual predicate."""

    @abstractmethod
    def visit_true(self) -> T:
        """Visit a bound True predicate."""

    @abstractmethod
    def visit_false(self) -> T:
        """Visit a bound False predicate."""

    @abstractmethod
    def visit_not(self, child_result: T) -> T:
        """Visit a bound Not predicate."""

    @abstractmethod
    def visit_and(self, left_result: T, right_result: T) -> T:
        """Visit a bound And predicate."""

    @abstractmethod
    def visit_or(self, left_result: T, right_result: T) -> T:
        """Visit a bound Or predicate."""

    @abstractmethod
    def visit_starts_with(self, term: BoundTerm, literal: LiteralValue) -> T:
        """Visit bound StartsWith predicate."""

    @abstractmethod
    def visit_not_starts_with(self, term: BoundTerm, literal: LiteralValue) -> T:
        """Visit bound NotStartsWith predicate."""

    def visit_unbound_predicate(self, predicate: UnboundPredicate) -> T:
        """Visit an unbound predicate.

        Args:
            predicate (UnboundPredicate): An unbound predicate.
        Raises:
            TypeError: This always raises since an unbound predicate is not expected in a bound boolean expression.
        """
        raise TypeError(f"Not a bound predicate: {predicate}")

    def visit_bound_predicate(self, predicate: BoundPredicate) -> T:
        """Visit a bound predicate.

        Args:
            predicate (BoundPredicate): A bound predicate.
        """
        return visit_bound_predicate(predicate, self)


@singledispatch
def visit_bound_predicate(expr: BoundPredicate, _: BooleanExpressionVisitor[T]) -> T:
    raise TypeError(f"Unknown predicate: {expr}")


@visit_bound_predicate.register(BoundIn)
def _(expr: BoundIn, visitor: BoundBooleanExpressionVisitor[T]) -> T:
    return visitor.visit_in(term=expr.term, literals=expr.value_set)


@visit_bound_predicate.register(BoundNotIn)
def _(expr: BoundNotIn, visitor: BoundBooleanExpressionVisitor[T]) -> T:
    return visitor.visit_not_in(term=expr.term, literals=expr.value_set)


@visit_bound_predicate.register(BoundIsNaN)
def _(expr: BoundIsNaN, visitor: BoundBooleanExpressionVisitor[T]) -> T:
    return visitor.visit_is_nan(term=expr.term)


@visit_bound_predicate.register(BoundNotNaN)
def _(expr: BoundNotNaN, visitor: BoundBooleanExpressionVisitor[T]) -> T:
    return visitor.visit_not_nan(term=expr.term)


@visit_bound_predicate.register(BoundIsNull)
def _(expr: BoundIsNull, visitor: BoundBooleanExpressionVisitor[T]) -> T:
    return visitor.visit_is_null(term=expr.term)


@visit_bound_predicate.register(BoundNotNull)
def _(expr: BoundNotNull, visitor: BoundBooleanExpressionVisitor[T]) -> T:
    return visitor.visit_not_null(term=expr.term)


@visit_bound_predicate.register(BoundEqualTo)
def _(expr: BoundEqualTo, visitor: BoundBooleanExpressionVisitor[T]) -> T:
    return visitor.visit_equal(term=expr.term, literal=expr.literal)


@visit_bound_predicate.register(BoundNotEqualTo)
def _(expr: BoundNotEqualTo, visitor: BoundBooleanExpressionVisitor[T]) -> T:
    return visitor.visit_not_equal(term=expr.term, literal=expr.literal)


@visit_bound_predicate.register(BoundGreaterThanOrEqual)
def _(expr: BoundGreaterThanOrEqual, visitor: BoundBooleanExpressionVisitor[T]) -> T:
    """Visit a bound GreaterThanOrEqual predicate."""
    return visitor.visit_greater_than_or_equal(term=expr.term, literal=expr.literal)


@visit_bound_predicate.register(BoundGreaterThan)
def _(expr: BoundGreaterThan, visitor: BoundBooleanExpressionVisitor[T]) -> T:
    return visitor.visit_greater_than(term=expr.term, literal=expr.literal)


@visit_bound_predicate.register(BoundLessThan)
def _(expr: BoundLessThan, visitor: BoundBooleanExpressionVisitor[T]) -> T:
    return visitor.visit_less_than(term=expr.term, literal=expr.literal)


@visit_bound_predicate.register(BoundLessThanOrEqual)
def _(expr: BoundLessThanOrEqual, visitor: BoundBooleanExpressionVisitor[T]) -> T:
    return visitor.visit_less_than_or_equal(term=expr.term, literal=expr.literal)


@visit_bound_predicate.register(BoundStartsWith)
def _(expr: BoundStartsWith, visitor: BoundBooleanExpressionVisitor[T]) -> T:
    return visitor.visit_starts_with(term=expr.term, literal=expr.literal)


@visit_bound_predicate.register(BoundNotStartsWith)
def _(expr: BoundNotStartsWith, visitor: BoundBooleanExpressionVisitor[T]) -> T:
    return visitor.visit_not_starts_with(term=expr.term, literal=expr.literal)


def rewrite_not(expr: BooleanExpression) -> BooleanExpression:
    return visit(expr, _RewriteNotVisitor())


class _RewriteNotVisitor(BooleanExpressionVisitor[BooleanExpression]):
    """Inverts the negations."""

    def visit_true(self) -> BooleanExpression:
        return AlwaysTrue()

    def visit_false(self) -> BooleanExpression:
        return AlwaysFalse()

    def visit_not(self, child_result: BooleanExpression) -> BooleanExpression:
        return ~child_result

    def visit_and(self, left_result: BooleanExpression, right_result: BooleanExpression) -> BooleanExpression:
        return And(left=left_result, right=right_result)

    def visit_or(self, left_result: BooleanExpression, right_result: BooleanExpression) -> BooleanExpression:
        return Or(left=left_result, right=right_result)

    def visit_unbound_predicate(self, predicate: UnboundPredicate) -> BooleanExpression:
        return predicate

    def visit_bound_predicate(self, predicate: BoundPredicate) -> BooleanExpression:
        return predicate


def expression_evaluator(schema: Schema, unbound: BooleanExpression, case_sensitive: bool) -> Callable[[StructProtocol], bool]:
    return _ExpressionEvaluator(schema, unbound, case_sensitive).eval


class _ExpressionEvaluator(BoundBooleanExpressionVisitor[bool]):
    bound: BooleanExpression
    struct: StructProtocol

    def __init__(self, schema: Schema, unbound: BooleanExpression, case_sensitive: bool):
        self.bound = bind(schema, unbound, case_sensitive)

    def eval(self, struct: StructProtocol) -> bool:
        self.struct = struct
        return visit(self.bound, self)

    def visit_in(self, term: BoundTerm, literals: set[L]) -> bool:
        return term.eval(self.struct) in literals

    def visit_not_in(self, term: BoundTerm, literals: set[L]) -> bool:
        return term.eval(self.struct) not in literals

    def visit_is_nan(self, term: BoundTerm) -> bool:
        val = term.eval(self.struct)
        return val != val

    def visit_not_nan(self, term: BoundTerm) -> bool:
        val = term.eval(self.struct)
        return val == val

    def visit_is_null(self, term: BoundTerm) -> bool:
        return term.eval(self.struct) is None

    def visit_not_null(self, term: BoundTerm) -> bool:
        return term.eval(self.struct) is not None

    def visit_equal(self, term: BoundTerm, literal: LiteralValue) -> bool:
        return term.eval(self.struct) == literal.value

    def visit_not_equal(self, term: BoundTerm, literal: LiteralValue) -> bool:
        return term.eval(self.struct) != literal.value

    def visit_greater_than_or_equal(self, term: BoundTerm, literal: LiteralValue) -> bool:
        value = term.eval(self.struct)
        return value is not None and value >= literal.value

    def visit_greater_than(self, term: BoundTerm, literal: LiteralValue) -> bool:
        value = term.eval(self.struct)
        return value is not None and value > literal.value

    def visit_less_than(self, term: BoundTerm, literal: LiteralValue) -> bool:
        value = term.eval(self.struct)
        return value is not None and value < literal.value

    def visit_less_than_or_equal(self, term: BoundTerm, literal: LiteralValue) -> bool:
        value = term.eval(self.struct)
        return value is not None and value <= literal.value

    def visit_starts_with(self, term: BoundTerm, literal: LiteralValue) -> bool:
        eval_res = term.eval(self.struct)
        return eval_res is not None and str(eval_res).startswith(str(literal.value))

    def visit_not_starts_with(self, term: BoundTerm, literal: LiteralValue) -> bool:
        return not self.visit_starts_with(term, literal)

    def visit_true(self) -> bool:
        return True

    def visit_false(self) -> bool:
        return False

    def visit_not(self, child_result: bool) -> bool:
        return not child_result

    def visit_and(self, left_result: bool, right_result: bool) -> bool:
        return left_result and right_result

    def visit_or(self, left_result: bool, right_result: bool) -> bool:
        return left_result or right_result


ROWS_MIGHT_MATCH = True
ROWS_MUST_MATCH = True
ROWS_CANNOT_MATCH = False
ROWS_MIGHT_NOT_MATCH = False
IN_PREDICATE_LIMIT = 200


def _from_byte_buffer(field_type: IcebergType, val: bytes) -> Any:
    if not isinstance(field_type, PrimitiveType):
        raise ValueError(f"Expected a PrimitiveType, got: {type(field_type)}")
    return from_bytes(field_type, val)


class _ManifestEvalVisitor(BoundBooleanExpressionVisitor[bool]):
    partition_fields: list[PartitionFieldSummary]
    partition_filter: BooleanExpression

    def __init__(self, partition_struct_schema: Schema, partition_filter: BooleanExpression, case_sensitive: bool) -> None:
        self.partition_filter = bind(partition_struct_schema, rewrite_not(partition_filter), case_sensitive)

    def eval(self, manifest: ManifestFile) -> bool:
        if partitions := manifest.partitions:
            self.partition_fields = partitions
            return visit(self.partition_filter, self)

        # No partition information
        return ROWS_MIGHT_MATCH

    def visit_in(self, term: BoundTerm, literals: set[L]) -> bool:
        pos = term.ref().accessor.position
        field = self.partition_fields[pos]

        if field.lower_bound is None:
            return ROWS_CANNOT_MATCH

        if len(literals) > IN_PREDICATE_LIMIT:
            return ROWS_MIGHT_MATCH

        lower = _from_byte_buffer(term.ref().field.field_type, field.lower_bound)

        if all(lower > val for val in literals):
            return ROWS_CANNOT_MATCH

        if field.upper_bound is not None:
            upper = _from_byte_buffer(term.ref().field.field_type, field.upper_bound)
            if all(upper < val for val in literals):
                return ROWS_CANNOT_MATCH

        return ROWS_MIGHT_MATCH

    def visit_not_in(self, term: BoundTerm, literals: set[L]) -> bool:
        # because the bounds are not necessarily a min or max value, this cannot be answered using
        # them. notIn(col, {X, ...}) with (X, Y) doesn't guarantee that X is a value in col.
        return ROWS_MIGHT_MATCH

    def visit_is_nan(self, term: BoundTerm) -> bool:
        pos = term.ref().accessor.position
        field = self.partition_fields[pos]

        if field.contains_nan is False:
            return ROWS_CANNOT_MATCH

        return ROWS_MIGHT_MATCH

    def visit_not_nan(self, term: BoundTerm) -> bool:
        pos = term.ref().accessor.position
        field = self.partition_fields[pos]

        if field.contains_nan is True and field.contains_null is False and field.lower_bound is None:
            return ROWS_CANNOT_MATCH

        return ROWS_MIGHT_MATCH

    def visit_is_null(self, term: BoundTerm) -> bool:
        pos = term.ref().accessor.position

        if self.partition_fields[pos].contains_null is False:
            return ROWS_CANNOT_MATCH

        return ROWS_MIGHT_MATCH

    def visit_not_null(self, term: BoundTerm) -> bool:
        pos = term.ref().accessor.position

        # contains_null encodes whether at least one partition value is null,
        # lowerBound is null if all partition values are null
        all_null = self.partition_fields[pos].contains_null is True and self.partition_fields[pos].lower_bound is None

        if all_null and isinstance(term.ref().field.field_type, (DoubleType, FloatType)):
            # floating point types may include NaN values, which we check separately.
            # In case bounds don't include NaN value, contains_nan needs to be checked against.
            all_null = self.partition_fields[pos].contains_nan is False

        if all_null:
            return ROWS_CANNOT_MATCH

        return ROWS_MIGHT_MATCH

    def visit_equal(self, term: BoundTerm, literal: LiteralValue) -> bool:
        pos = term.ref().accessor.position
        field = self.partition_fields[pos]

        if field.lower_bound is None or field.upper_bound is None:
            # values are all null and literal cannot contain null
            return ROWS_CANNOT_MATCH

        lower = _from_byte_buffer(term.ref().field.field_type, field.lower_bound)

        if lower > literal.value:
            return ROWS_CANNOT_MATCH

        upper = _from_byte_buffer(term.ref().field.field_type, field.upper_bound)

        if literal.value > upper:
            return ROWS_CANNOT_MATCH

        return ROWS_MIGHT_MATCH

    def visit_not_equal(self, term: BoundTerm, literal: LiteralValue) -> bool:
        # because the bounds are not necessarily a min or max value, this cannot be answered using
        # them. notEq(col, X) with (X, Y) doesn't guarantee that X is a value in col.
        return ROWS_MIGHT_MATCH

    def visit_greater_than_or_equal(self, term: BoundTerm, literal: LiteralValue) -> bool:
        pos = term.ref().accessor.position
        field = self.partition_fields[pos]

        if field.upper_bound is None:
            return ROWS_CANNOT_MATCH

        upper = _from_byte_buffer(term.ref().field.field_type, field.upper_bound)

        if literal.value > upper:
            return ROWS_CANNOT_MATCH

        return ROWS_MIGHT_MATCH

    def visit_greater_than(self, term: BoundTerm, literal: LiteralValue) -> bool:
        pos = term.ref().accessor.position
        field = self.partition_fields[pos]

        if field.upper_bound is None:
            return ROWS_CANNOT_MATCH

        upper = _from_byte_buffer(term.ref().field.field_type, field.upper_bound)

        if literal.value >= upper:
            return ROWS_CANNOT_MATCH

        return ROWS_MIGHT_MATCH

    def visit_less_than(self, term: BoundTerm, literal: LiteralValue) -> bool:
        pos = term.ref().accessor.position
        field = self.partition_fields[pos]

        if field.lower_bound is None:
            return ROWS_CANNOT_MATCH

        lower = _from_byte_buffer(term.ref().field.field_type, field.lower_bound)

        if literal.value <= lower:
            return ROWS_CANNOT_MATCH

        return ROWS_MIGHT_MATCH

    def visit_less_than_or_equal(self, term: BoundTerm, literal: LiteralValue) -> bool:
        pos = term.ref().accessor.position
        field = self.partition_fields[pos]

        if field.lower_bound is None:
            return ROWS_CANNOT_MATCH

        lower = _from_byte_buffer(term.ref().field.field_type, field.lower_bound)

        if literal.value < lower:
            return ROWS_CANNOT_MATCH

        return ROWS_MIGHT_MATCH

    def visit_starts_with(self, term: BoundTerm, literal: LiteralValue) -> bool:
        pos = term.ref().accessor.position
        field = self.partition_fields[pos]
        prefix = str(literal.value)
        len_prefix = len(prefix)

        if field.lower_bound is None:
            return ROWS_CANNOT_MATCH

        lower = _from_byte_buffer(term.ref().field.field_type, field.lower_bound)
        # truncate lower bound so that its length is not greater than the length of prefix
        if lower is not None and lower[:len_prefix] > prefix:
            return ROWS_CANNOT_MATCH

        if field.upper_bound is None:
            return ROWS_CANNOT_MATCH

        upper = _from_byte_buffer(term.ref().field.field_type, field.upper_bound)
        # truncate upper bound so that its length is not greater than the length of prefix
        if upper is not None and upper[:len_prefix] < prefix:
            return ROWS_CANNOT_MATCH

        return ROWS_MIGHT_MATCH

    def visit_not_starts_with(self, term: BoundTerm, literal: LiteralValue) -> bool:
        pos = term.ref().accessor.position
        field = self.partition_fields[pos]
        prefix = str(literal.value)
        len_prefix = len(prefix)

        if field.contains_null or field.lower_bound is None or field.upper_bound is None:
            return ROWS_MIGHT_MATCH

        # not_starts_with will match unless all values must start with the prefix. This happens when
        # the lower and upper bounds both start with the prefix.
        lower = _from_byte_buffer(term.ref().field.field_type, field.lower_bound)
        upper = _from_byte_buffer(term.ref().field.field_type, field.upper_bound)

        if lower is not None and upper is not None:
            # if lower is shorter than the prefix then lower doesn't start with the prefix
            if len(lower) < len_prefix:
                return ROWS_MIGHT_MATCH

            if lower[:len_prefix] == prefix:
                # if upper is shorter than the prefix then upper can't start with the prefix
                if len(upper) < len_prefix:
                    return ROWS_MIGHT_MATCH

                if upper[:len_prefix] == prefix:
                    return ROWS_CANNOT_MATCH

        return ROWS_MIGHT_MATCH

    def visit_true(self) -> bool:
        return ROWS_MIGHT_MATCH

    def visit_false(self) -> bool:
        return ROWS_CANNOT_MATCH

    def visit_not(self, child_result: bool) -> bool:
        return not child_result

    def visit_and(self, left_result: bool, right_result: bool) -> bool:
        return left_result and right_result

    def visit_or(self, left_result: bool, right_result: bool) -> bool:
        return left_result or right_result


def manifest_evaluator(
    partition_spec: PartitionSpec, schema: Schema, partition_filter: BooleanExpression, case_sensitive: bool = True
) -> Callable[[ManifestFile], bool]:
    partition_type = partition_spec.partition_type(schema)
    partition_schema = Schema(*partition_type.fields)
    evaluator = _ManifestEvalVisitor(partition_schema, partition_filter, case_sensitive)
    return evaluator.eval


class ProjectionEvaluator(BooleanExpressionVisitor[BooleanExpression], ABC):
    schema: Schema
    spec: PartitionSpec
    case_sensitive: bool

    def __init__(self, schema: Schema, spec: PartitionSpec, case_sensitive: bool):
        self.schema = schema
        self.spec = spec
        self.case_sensitive = case_sensitive

    def project(self, expr: BooleanExpression) -> BooleanExpression:
        #  projections assume that there are no NOT nodes in the expression tree. to ensure that this
        #  is the case, the expression is rewritten to push all NOT nodes down to the expression
        #  leaf nodes.
        #  this is necessary to ensure that the default expression returned when a predicate can't be
        #  projected is correct.
        return visit(bind(self.schema, rewrite_not(expr), self.case_sensitive), self)

    def visit_true(self) -> BooleanExpression:
        return AlwaysTrue()

    def visit_false(self) -> BooleanExpression:
        return AlwaysFalse()

    def visit_not(self, child_result: BooleanExpression) -> BooleanExpression:
        raise ValueError(f"Cannot project not expression, should be rewritten: {child_result}")

    def visit_and(self, left_result: BooleanExpression, right_result: BooleanExpression) -> BooleanExpression:
        return And(left_result, right_result)

    def visit_or(self, left_result: BooleanExpression, right_result: BooleanExpression) -> BooleanExpression:
        return Or(left_result, right_result)

    def visit_unbound_predicate(self, predicate: UnboundPredicate) -> BooleanExpression:
        raise ValueError(f"Cannot project unbound predicate: {predicate}")


class InclusiveProjection(ProjectionEvaluator):
    def visit_bound_predicate(self, predicate: BoundPredicate) -> BooleanExpression:
        parts = self.spec.fields_by_source_id(predicate.term.ref().field.field_id)

        result: BooleanExpression = AlwaysTrue()
        for part in parts:
            # consider (d = 2019-01-01) with bucket(7, d) and bucket(5, d)
            # projections: b1 = bucket(7, '2019-01-01') = 5, b2 = bucket(5, '2019-01-01') = 0
            # any value where b1 != 5 or any value where b2 != 0 cannot be the '2019-01-01'
            #
            # similarly, if partitioning by day(ts) and hour(ts), the more restrictive
            # projection should be used. ts = 2019-01-01T01:00:00 produces day=2019-01-01 and
            # hour=2019-01-01-01. the value will be in 2019-01-01-01 and not in 2019-01-01-02.
            incl_projection = part.transform.project(name=part.name, pred=predicate)
            if incl_projection is not None:
                result = And(result, incl_projection)

        return result


def inclusive_projection(
    schema: Schema, spec: PartitionSpec, case_sensitive: bool = True
) -> Callable[[BooleanExpression], BooleanExpression]:
    return InclusiveProjection(schema, spec, case_sensitive).project


class _ColumnNameTranslator(BooleanExpressionVisitor[BooleanExpression]):
    """Converts the column names with the ones in the actual file.

    Args:
      file_schema (Schema): The schema of the file.
      case_sensitive (bool): Whether to consider case when binding a reference to a field in a schema, defaults to True.
      projected_field_values (Dict[int, Any]): Values for projected fields not present in the data file.

    Raises:
        TypeError: In the case of an UnboundPredicate.
        ValueError: When a column name cannot be found.
    """

    file_schema: Schema
    case_sensitive: bool
    projected_field_values: dict[int, Any]

    def __init__(self, file_schema: Schema, case_sensitive: bool, projected_field_values: dict[int, Any] = EMPTY_DICT) -> None:
        self.file_schema = file_schema
        self.case_sensitive = case_sensitive
        self.projected_field_values = projected_field_values

    def visit_true(self) -> BooleanExpression:
        return AlwaysTrue()

    def visit_false(self) -> BooleanExpression:
        return AlwaysFalse()

    def visit_not(self, child_result: BooleanExpression) -> BooleanExpression:
        return Not(child=child_result)

    def visit_and(self, left_result: BooleanExpression, right_result: BooleanExpression) -> BooleanExpression:
        return And(left=left_result, right=right_result)

    def visit_or(self, left_result: BooleanExpression, right_result: BooleanExpression) -> BooleanExpression:
        return Or(left=left_result, right=right_result)

    def visit_unbound_predicate(self, predicate: UnboundPredicate) -> BooleanExpression:
        raise TypeError(f"Expected Bound Predicate, got: {predicate.term}")

    def visit_bound_predicate(self, predicate: BoundPredicate) -> BooleanExpression:
        field = predicate.term.ref().field
        field_id = field.field_id
        file_column_name = self.file_schema.find_column_name(field_id)

        if file_column_name is None:
            # In the case of schema evolution or column projection, the field might not be present in the file schema.
            # we can use the projected value or the field's default value as a constant and evaluate it against the predicate
            pred: BooleanExpression
            if isinstance(predicate, BoundUnaryPredicate):
                pred = predicate.as_unbound(field.name)
            elif isinstance(predicate, BoundLiteralPredicate):
                pred = predicate.as_unbound(field.name, predicate.literal)
            elif isinstance(predicate, BoundSetPredicate):
                pred = predicate.as_unbound(field.name, predicate.literals)
            else:
                raise ValueError(f"Unsupported predicate: {predicate}")

            # In the order described by the "Column Projection" section of the Iceberg spec:
            # https://iceberg.apache.org/spec/#column-projection
            # Evaluate column projection first if it exists, otherwise default to the initial-default-value
            field_value = (
                self.projected_field_values[field_id] if field.field_id in self.projected_field_values else field.initial_default
            )
            return (
                AlwaysTrue()
                if expression_evaluator(Schema(field), pred, case_sensitive=self.case_sensitive)(Record(field_value))
                else AlwaysFalse()
            )

        if isinstance(predicate, BoundUnaryPredicate):
            return predicate.as_unbound(file_column_name)
        elif isinstance(predicate, BoundLiteralPredicate):
            return predicate.as_unbound(file_column_name, predicate.literal)
        elif isinstance(predicate, BoundSetPredicate):
            return predicate.as_unbound(file_column_name, predicate.literals)
        else:
            raise ValueError(f"Unsupported predicate: {predicate}")


def translate_column_names(
    expr: BooleanExpression, file_schema: Schema, case_sensitive: bool = True, projected_field_values: dict[int, Any] = EMPTY_DICT
) -> BooleanExpression:
    return visit(expr, _ColumnNameTranslator(file_schema, case_sensitive, projected_field_values))


class _ExpressionFieldIDs(BooleanExpressionVisitor[set[int]]):
    """Extracts the field IDs used in the BooleanExpression."""

    def visit_true(self) -> set[int]:
        return set()

    def visit_false(self) -> set[int]:
        return set()

    def visit_not(self, child_result: set[int]) -> set[int]:
        return child_result

    def visit_and(self, left_result: set[int], right_result: set[int]) -> set[int]:
        return left_result.union(right_result)

    def visit_or(self, left_result: set[int], right_result: set[int]) -> set[int]:
        return left_result.union(right_result)

    def visit_unbound_predicate(self, predicate: UnboundPredicate) -> set[int]:
        raise ValueError("Only works on bound records")

    def visit_bound_predicate(self, predicate: BoundPredicate) -> set[int]:
        return {predicate.term.ref().field.field_id}


def extract_field_ids(expr: BooleanExpression) -> set[int]:
    return visit(expr, _ExpressionFieldIDs())


class _RewriteToDNF(BooleanExpressionVisitor[tuple[BooleanExpression, ...]]):
    def visit_true(self) -> tuple[BooleanExpression, ...]:
        return (AlwaysTrue(),)

    def visit_false(self) -> tuple[BooleanExpression, ...]:
        return (AlwaysFalse(),)

    def visit_not(self, child_result: tuple[BooleanExpression, ...]) -> tuple[BooleanExpression, ...]:
        raise ValueError(f"Not expressions are not allowed: {child_result}")

    def visit_and(
        self, left_result: tuple[BooleanExpression, ...], right_result: tuple[BooleanExpression, ...]
    ) -> tuple[BooleanExpression, ...]:
        # Distributive law:
        # ((P OR Q) AND (R OR S)) AND (((P AND R) OR (P AND S)) OR ((Q AND R) OR ((Q AND S)))
        # A AND (B OR C) = (A AND B) OR (A AND C)
        # (A OR B) AND C = (A AND C) OR (B AND C)
        return tuple(And(le, re) for le in left_result for re in right_result)

    def visit_or(
        self, left_result: tuple[BooleanExpression, ...], right_result: tuple[BooleanExpression, ...]
    ) -> tuple[BooleanExpression, ...]:
        return left_result + right_result

    def visit_unbound_predicate(self, predicate: UnboundPredicate) -> tuple[BooleanExpression, ...]:
        return (predicate,)

    def visit_bound_predicate(self, predicate: BoundPredicate) -> tuple[BooleanExpression, ...]:
        return (predicate,)


def rewrite_to_dnf(expr: BooleanExpression) -> tuple[BooleanExpression, ...]:
    # Rewrites an arbitrary boolean expression to disjunctive normal form (DNF):
    # (A AND NOT(B) AND C) OR (NOT(D) AND E AND F) OR (G)
    expr_without_not = rewrite_not(expr)
    return visit(expr_without_not, _RewriteToDNF())


class ExpressionToPlainFormat(BoundBooleanExpressionVisitor[list[tuple[str, str, Any]]]):
    cast_int_to_date: bool

    def __init__(self, cast_int_to_date: bool = False) -> None:
        self.cast_int_to_date = cast_int_to_date

    def _cast_if_necessary(self, iceberg_type: IcebergType, literal: L | set[L]) -> L | set[L]:
        if self.cast_int_to_date:
            iceberg_type_class = type(iceberg_type)
            conversions = {TimestampType: micros_to_timestamp, TimestamptzType: micros_to_timestamptz}
            if iceberg_type_class in conversions:
                conversion_function = conversions[iceberg_type_class]
                if isinstance(literal, set):
                    return {conversion_function(lit) for lit in literal}  # type: ignore
                else:
                    return conversion_function(literal)  # type: ignore
        return literal

    def visit_in(self, term: BoundTerm, literals: set[L]) -> list[tuple[str, str, Any]]:
        field = term.ref().field
        return [(term.ref().field.name, "in", self._cast_if_necessary(field.field_type, literals))]

    def visit_not_in(self, term: BoundTerm, literals: set[L]) -> list[tuple[str, str, Any]]:
        field = term.ref().field
        return [(field.name, "not in", self._cast_if_necessary(field.field_type, literals))]

    def visit_is_nan(self, term: BoundTerm) -> list[tuple[str, str, Any]]:
        return [(term.ref().field.name, "==", float("nan"))]

    def visit_not_nan(self, term: BoundTerm) -> list[tuple[str, str, Any]]:
        return [(term.ref().field.name, "!=", float("nan"))]

    def visit_is_null(self, term: BoundTerm) -> list[tuple[str, str, Any]]:
        return [(term.ref().field.name, "==", None)]

    def visit_not_null(self, term: BoundTerm) -> list[tuple[str, str, Any]]:
        return [(term.ref().field.name, "!=", None)]

    def visit_equal(self, term: BoundTerm, literal: LiteralValue) -> list[tuple[str, str, Any]]:
        return [(term.ref().field.name, "==", self._cast_if_necessary(term.ref().field.field_type, literal.value))]

    def visit_not_equal(self, term: BoundTerm, literal: LiteralValue) -> list[tuple[str, str, Any]]:
        return [(term.ref().field.name, "!=", self._cast_if_necessary(term.ref().field.field_type, literal.value))]

    def visit_greater_than_or_equal(self, term: BoundTerm, literal: LiteralValue) -> list[tuple[str, str, Any]]:
        return [(term.ref().field.name, ">=", self._cast_if_necessary(term.ref().field.field_type, literal.value))]

    def visit_greater_than(self, term: BoundTerm, literal: LiteralValue) -> list[tuple[str, str, Any]]:
        return [(term.ref().field.name, ">", self._cast_if_necessary(term.ref().field.field_type, literal.value))]

    def visit_less_than(self, term: BoundTerm, literal: LiteralValue) -> list[tuple[str, str, Any]]:
        return [(term.ref().field.name, "<", self._cast_if_necessary(term.ref().field.field_type, literal.value))]

    def visit_less_than_or_equal(self, term: BoundTerm, literal: LiteralValue) -> list[tuple[str, str, Any]]:
        return [(term.ref().field.name, "<=", self._cast_if_necessary(term.ref().field.field_type, literal.value))]

    def visit_starts_with(self, term: BoundTerm, literal: LiteralValue) -> list[tuple[str, str, Any]]:
        return []

    def visit_not_starts_with(self, term: BoundTerm, literal: LiteralValue) -> list[tuple[str, str, Any]]:
        return []

    def visit_true(self) -> list[tuple[str, str, Any]]:
        return []  # Not supported

    def visit_false(self) -> list[tuple[str, str, Any]]:
        raise ValueError("Not supported: AlwaysFalse")

    def visit_not(self, child_result: list[tuple[str, str, Any]]) -> list[tuple[str, str, Any]]:
        raise ValueError(f"Not allowed: {child_result}")

    def visit_and(
        self, left_result: list[tuple[str, str, Any]], right_result: list[tuple[str, str, Any]]
    ) -> list[tuple[str, str, Any]]:
        return left_result + right_result

    def visit_or(
        self, left_result: list[tuple[str, str, Any]], right_result: list[tuple[str, str, Any]]
    ) -> list[tuple[str, str, Any]]:
        raise ValueError(f"Not allowed: {left_result} || {right_result}")


def expression_to_plain_format(
    expressions: tuple[BooleanExpression, ...], cast_int_to_datetime: bool = False
) -> list[list[tuple[str, str, Any]]]:
    """Format a Disjunctive Normal Form expression.

    These are the formats that the expression can be fed into:

    - https://arrow.apache.org/docs/python/generated/pyarrow.parquet.read_table.html
    - https://docs.dask.org/en/stable/generated/dask.dataframe.read_parquet.html

    Contrary to normal DNF that may contain Not expressions, but here they should have
    been rewritten. This can be done using ``rewrite_not(...)``.

    Keep in mind that this is only used for page skipping, and still needs to filter
    on a row level.

    Args:
        expressions: Expression in Disjunctive Normal Form.

    Returns:
        Formatter filter compatible with Dask and PyArrow.
    """
    # In the form of expr1 ∨ expr2 ∨ ... ∨ exprN
    visitor = ExpressionToPlainFormat(cast_int_to_datetime)
    return [visit(expression, visitor) for expression in expressions]


class _MetricsEvaluator(BoundBooleanExpressionVisitor[bool], ABC):
    value_counts: dict[int, int]
    null_counts: dict[int, int]
    nan_counts: dict[int, int]
    lower_bounds: dict[int, bytes]
    upper_bounds: dict[int, bytes]

    def visit_true(self) -> bool:
        # all rows match
        return ROWS_MIGHT_MATCH

    def visit_false(self) -> bool:
        # all rows fail
        return ROWS_CANNOT_MATCH

    def visit_not(self, child_result: bool) -> bool:
        raise ValueError(f"NOT should be rewritten: {child_result}")

    def visit_and(self, left_result: bool, right_result: bool) -> bool:
        return left_result and right_result

    def visit_or(self, left_result: bool, right_result: bool) -> bool:
        return left_result or right_result

    def _contains_nulls_only(self, field_id: int) -> bool:
        if (value_count := self.value_counts.get(field_id)) and (null_count := self.null_counts.get(field_id)):
            return value_count == null_count
        return False

    def _contains_nans_only(self, field_id: int) -> bool:
        if (nan_count := self.nan_counts.get(field_id)) and (value_count := self.value_counts.get(field_id)):
            return nan_count == value_count
        return False

    def _is_nan(self, val: Any) -> bool:
        try:
            return math.isnan(val)
        except TypeError:
            # In the case of None or other non-numeric types
            return False


class _InclusiveMetricsEvaluator(_MetricsEvaluator):
    struct: StructType
    expr: BooleanExpression

    def __init__(
        self, schema: Schema, expr: BooleanExpression, case_sensitive: bool = True, include_empty_files: bool = False
    ) -> None:
        self.struct = schema.as_struct()
        self.include_empty_files = include_empty_files
        self.expr = bind(schema, rewrite_not(expr), case_sensitive)

    def eval(self, file: DataFile) -> bool:
        """Test whether the file may contain records that match the expression."""
        if not self.include_empty_files and file.record_count == 0:
            return ROWS_CANNOT_MATCH

        if file.record_count < 0:
            # Older version don't correctly implement record count from avro file and thus
            # set record count -1 when importing avro tables to iceberg tables. This should
            # be updated once we implemented and set correct record count.
            return ROWS_MIGHT_MATCH

        self.value_counts = file.value_counts or EMPTY_DICT
        self.null_counts = file.null_value_counts or EMPTY_DICT
        self.nan_counts = file.nan_value_counts or EMPTY_DICT
        self.lower_bounds = file.lower_bounds or EMPTY_DICT
        self.upper_bounds = file.upper_bounds or EMPTY_DICT

        return visit(self.expr, self)

    def _may_contain_null(self, field_id: int) -> bool:
        return self.null_counts is None or (field_id in self.null_counts and self.null_counts.get(field_id) is not None)

    def _contains_nans_only(self, field_id: int) -> bool:
        if (nan_count := self.nan_counts.get(field_id)) and (value_count := self.value_counts.get(field_id)):
            return nan_count == value_count
        return False

    def visit_is_null(self, term: BoundTerm) -> bool:
        field_id = term.ref().field.field_id

        if self.null_counts.get(field_id) == 0:
            return ROWS_CANNOT_MATCH

        return ROWS_MIGHT_MATCH

    def visit_not_null(self, term: BoundTerm) -> bool:
        # no need to check whether the field is required because binding evaluates that case
        # if the column has no non-null values, the expression cannot match
        field_id = term.ref().field.field_id

        if self._contains_nulls_only(field_id):
            return ROWS_CANNOT_MATCH

        return ROWS_MIGHT_MATCH

    def visit_is_nan(self, term: BoundTerm) -> bool:
        field_id = term.ref().field.field_id

        if self.nan_counts.get(field_id) == 0:
            return ROWS_CANNOT_MATCH

        # when there's no nanCounts information, but we already know the column only contains null,
        # it's guaranteed that there's no NaN value
        if self._contains_nulls_only(field_id):
            return ROWS_CANNOT_MATCH

        return ROWS_MIGHT_MATCH

    def visit_not_nan(self, term: BoundTerm) -> bool:
        field_id = term.ref().field.field_id

        if self._contains_nans_only(field_id):
            return ROWS_CANNOT_MATCH

        return ROWS_MIGHT_MATCH

    def visit_less_than(self, term: BoundTerm, literal: LiteralValue) -> bool:
        field = term.ref().field
        field_id = field.field_id

        if self._contains_nulls_only(field_id) or self._contains_nans_only(field_id):
            return ROWS_CANNOT_MATCH

        if not isinstance(field.field_type, PrimitiveType):
            raise ValueError(f"Expected PrimitiveType: {field.field_type}")

        if lower_bound_bytes := self.lower_bounds.get(field_id):
            lower_bound = from_bytes(field.field_type, lower_bound_bytes)

            if self._is_nan(lower_bound):
                # NaN indicates unreliable bounds. See the InclusiveMetricsEvaluator docs for more.
                return ROWS_MIGHT_MATCH

            if lower_bound >= literal.value:
                return ROWS_CANNOT_MATCH

        return ROWS_MIGHT_MATCH

    def visit_less_than_or_equal(self, term: BoundTerm, literal: LiteralValue) -> bool:
        field = term.ref().field
        field_id = field.field_id

        if self._contains_nulls_only(field_id) or self._contains_nans_only(field_id):
            return ROWS_CANNOT_MATCH

        if not isinstance(field.field_type, PrimitiveType):
            raise ValueError(f"Expected PrimitiveType: {field.field_type}")

        if lower_bound_bytes := self.lower_bounds.get(field_id):
            lower_bound = from_bytes(field.field_type, lower_bound_bytes)
            if self._is_nan(lower_bound):
                # NaN indicates unreliable bounds. See the InclusiveMetricsEvaluator docs for more.
                return ROWS_MIGHT_MATCH

            if lower_bound > literal.value:
                return ROWS_CANNOT_MATCH

        return ROWS_MIGHT_MATCH

    def visit_greater_than(self, term: BoundTerm, literal: LiteralValue) -> bool:
        field = term.ref().field
        field_id = field.field_id

        if self._contains_nulls_only(field_id) or self._contains_nans_only(field_id):
            return ROWS_CANNOT_MATCH

        if not isinstance(field.field_type, PrimitiveType):
            raise ValueError(f"Expected PrimitiveType: {field.field_type}")

        if upper_bound_bytes := self.upper_bounds.get(field_id):
            upper_bound = from_bytes(field.field_type, upper_bound_bytes)
            if upper_bound <= literal.value:
                if self._is_nan(upper_bound):
                    # NaN indicates unreliable bounds. See the InclusiveMetricsEvaluator docs for more.
                    return ROWS_MIGHT_MATCH

                return ROWS_CANNOT_MATCH

        return ROWS_MIGHT_MATCH

    def visit_greater_than_or_equal(self, term: BoundTerm, literal: LiteralValue) -> bool:
        field = term.ref().field
        field_id = field.field_id

        if self._contains_nulls_only(field_id) or self._contains_nans_only(field_id):
            return ROWS_CANNOT_MATCH

        if not isinstance(field.field_type, PrimitiveType):
            raise ValueError(f"Expected PrimitiveType: {field.field_type}")

        if upper_bound_bytes := self.upper_bounds.get(field_id):
            upper_bound = from_bytes(field.field_type, upper_bound_bytes)
            if upper_bound < literal.value:
                if self._is_nan(upper_bound):
                    # NaN indicates unreliable bounds. See the InclusiveMetricsEvaluator docs for more.
                    return ROWS_MIGHT_MATCH

                return ROWS_CANNOT_MATCH

        return ROWS_MIGHT_MATCH

    def visit_equal(self, term: BoundTerm, literal: LiteralValue) -> bool:
        field = term.ref().field
        field_id = field.field_id

        if self._contains_nulls_only(field_id) or self._contains_nans_only(field_id):
            return ROWS_CANNOT_MATCH

        if not isinstance(field.field_type, PrimitiveType):
            raise ValueError(f"Expected PrimitiveType: {field.field_type}")

        if lower_bound_bytes := self.lower_bounds.get(field_id):
            lower_bound = from_bytes(field.field_type, lower_bound_bytes)
            if self._is_nan(lower_bound):
                # NaN indicates unreliable bounds. See the InclusiveMetricsEvaluator docs for more.
                return ROWS_MIGHT_MATCH

            if lower_bound > literal.value:
                return ROWS_CANNOT_MATCH

        if upper_bound_bytes := self.upper_bounds.get(field_id):
            upper_bound = from_bytes(field.field_type, upper_bound_bytes)
            if self._is_nan(upper_bound):
                # NaN indicates unreliable bounds. See the InclusiveMetricsEvaluator docs for more.
                return ROWS_MIGHT_MATCH

            if upper_bound < literal.value:
                return ROWS_CANNOT_MATCH

        return ROWS_MIGHT_MATCH

    def visit_not_equal(self, term: BoundTerm, literal: LiteralValue) -> bool:
        return ROWS_MIGHT_MATCH

    def visit_in(self, term: BoundTerm, literals: set[L]) -> bool:
        field = term.ref().field
        field_id = field.field_id

        if self._contains_nulls_only(field_id) or self._contains_nans_only(field_id):
            return ROWS_CANNOT_MATCH

        if len(literals) > IN_PREDICATE_LIMIT:
            # skip evaluating the predicate if the number of values is too big
            return ROWS_MIGHT_MATCH

        if not isinstance(field.field_type, PrimitiveType):
            raise ValueError(f"Expected PrimitiveType: {field.field_type}")

        if lower_bound_bytes := self.lower_bounds.get(field_id):
            lower_bound = from_bytes(field.field_type, lower_bound_bytes)
            if self._is_nan(lower_bound):
                # NaN indicates unreliable bounds. See the InclusiveMetricsEvaluator docs for more.
                return ROWS_MIGHT_MATCH

            literals = {lit for lit in literals if lower_bound <= lit}  # type: ignore[operator]
            if len(literals) == 0:
                return ROWS_CANNOT_MATCH

        if upper_bound_bytes := self.upper_bounds.get(field_id):
            upper_bound = from_bytes(field.field_type, upper_bound_bytes)
            # this is different from Java, here NaN is always larger
            if self._is_nan(upper_bound):
                return ROWS_MIGHT_MATCH

            literals = {lit for lit in literals if upper_bound >= lit}  # type: ignore[operator]
            if len(literals) == 0:
                return ROWS_CANNOT_MATCH

        return ROWS_MIGHT_MATCH

    def visit_not_in(self, term: BoundTerm, literals: set[L]) -> bool:
        # because the bounds are not necessarily a min or max value, this cannot be answered using
        # them. notIn(col, {X, ...}) with (X, Y) doesn't guarantee that X is a value in col.
        return ROWS_MIGHT_MATCH

    def visit_starts_with(self, term: BoundTerm, literal: LiteralValue) -> bool:
        field = term.ref().field
        field_id: int = field.field_id

        if self._contains_nulls_only(field_id):
            return ROWS_CANNOT_MATCH

        if not isinstance(field.field_type, PrimitiveType):
            raise ValueError(f"Expected PrimitiveType: {field.field_type}")

        prefix = str(literal.value)
        len_prefix = len(prefix)

        if lower_bound_bytes := self.lower_bounds.get(field_id):
            lower_bound = str(from_bytes(field.field_type, lower_bound_bytes))

            # truncate lower bound so that its length is not greater than the length of prefix
            if lower_bound and lower_bound[:len_prefix] > prefix:
                return ROWS_CANNOT_MATCH

        if upper_bound_bytes := self.upper_bounds.get(field_id):
            upper_bound = str(from_bytes(field.field_type, upper_bound_bytes))

            # truncate upper bound so that its length is not greater than the length of prefix
            if upper_bound is not None and upper_bound[:len_prefix] < prefix:
                return ROWS_CANNOT_MATCH

        return ROWS_MIGHT_MATCH

    def visit_not_starts_with(self, term: BoundTerm, literal: LiteralValue) -> bool:
        field = term.ref().field
        field_id: int = field.field_id

        if self._may_contain_null(field_id):
            return ROWS_MIGHT_MATCH

        if not isinstance(field.field_type, PrimitiveType):
            raise ValueError(f"Expected PrimitiveType: {field.field_type}")

        prefix = str(literal.value)
        len_prefix = len(prefix)

        # not_starts_with will match unless all values must start with the prefix. This happens when
        # the lower and upper bounds both start with the prefix.
        if (lower_bound_bytes := self.lower_bounds.get(field_id)) and (upper_bound_bytes := self.upper_bounds.get(field_id)):
            lower_bound = str(from_bytes(field.field_type, lower_bound_bytes))
            upper_bound = str(from_bytes(field.field_type, upper_bound_bytes))

            # if lower is shorter than the prefix then lower doesn't start with the prefix
            if len(lower_bound) < len_prefix:
                return ROWS_MIGHT_MATCH

            if lower_bound[:len_prefix] == prefix:
                # if upper is shorter than the prefix then upper can't start with the prefix
                if len(upper_bound) < len_prefix:
                    return ROWS_MIGHT_MATCH

                if upper_bound[:len_prefix] == prefix:
                    return ROWS_CANNOT_MATCH

        return ROWS_MIGHT_MATCH


def strict_projection(
    schema: Schema, spec: PartitionSpec, case_sensitive: bool = True
) -> Callable[[BooleanExpression], BooleanExpression]:
    return StrictProjection(schema, spec, case_sensitive).project


class StrictProjection(ProjectionEvaluator):
    def visit_bound_predicate(self, predicate: BoundPredicate) -> BooleanExpression:
        parts = self.spec.fields_by_source_id(predicate.term.ref().field.field_id)

        result: BooleanExpression = AlwaysFalse()
        for part in parts:
            # consider (ts > 2019-01-01T01:00:00) with day(ts) and hour(ts)
            # projections: d >= 2019-01-02 and h >= 2019-01-01-02 (note the inclusive bounds).
            # any timestamp where either projection predicate is true must match the original
            # predicate. For example, ts = 2019-01-01T03:00:00 matches the hour projection but not
            # the day, but does match the original predicate.
            strict_projection = part.transform.strict_project(name=part.name, pred=predicate)
            if strict_projection is not None:
                result = Or(result, strict_projection)

        return result


class _StrictMetricsEvaluator(_MetricsEvaluator):
    struct: StructType
    expr: BooleanExpression

    def __init__(
        self, schema: Schema, expr: BooleanExpression, case_sensitive: bool = True, include_empty_files: bool = False
    ) -> None:
        self.struct = schema.as_struct()
        self.include_empty_files = include_empty_files
        self.expr = bind(schema, rewrite_not(expr), case_sensitive)

    def eval(self, file: DataFile) -> bool:
        """Test whether all records within the file match the expression.

        Args:
            file: A data file

        Returns: false if the file may contain any row that doesn't match
                    the expression, true otherwise.
        """
        if file.record_count <= 0:
            # Older version don't correctly implement record count from avro file and thus
            # set record count -1 when importing avro tables to iceberg tables. This should
            # be updated once we implemented and set correct record count.
            return ROWS_MUST_MATCH

        self.value_counts = file.value_counts or EMPTY_DICT
        self.null_counts = file.null_value_counts or EMPTY_DICT
        self.nan_counts = file.nan_value_counts or EMPTY_DICT
        self.lower_bounds = file.lower_bounds or EMPTY_DICT
        self.upper_bounds = file.upper_bounds or EMPTY_DICT

        return visit(self.expr, self)

    def visit_is_null(self, term: BoundTerm) -> bool:
        # no need to check whether the field is required because binding evaluates that case
        # if the column has any non-null values, the expression does not match
        field_id = term.ref().field.field_id

        if self._contains_nulls_only(field_id):
            return ROWS_MUST_MATCH
        else:
            return ROWS_MIGHT_NOT_MATCH

    def visit_not_null(self, term: BoundTerm) -> bool:
        # no need to check whether the field is required because binding evaluates that case
        # if the column has any non-null values, the expression does not match
        field_id = term.ref().field.field_id

        if (null_count := self.null_counts.get(field_id)) is not None and null_count == 0:
            return ROWS_MUST_MATCH
        else:
            return ROWS_MIGHT_NOT_MATCH

    def visit_is_nan(self, term: BoundTerm) -> bool:
        field_id = term.ref().field.field_id

        if self._contains_nans_only(field_id):
            return ROWS_MUST_MATCH
        else:
            return ROWS_MIGHT_NOT_MATCH

    def visit_not_nan(self, term: BoundTerm) -> bool:
        field_id = term.ref().field.field_id

        if (nan_count := self.nan_counts.get(field_id)) is not None and nan_count == 0:
            return ROWS_MUST_MATCH

        if self._contains_nulls_only(field_id):
            return ROWS_MUST_MATCH

        return ROWS_MIGHT_NOT_MATCH

    def visit_less_than(self, term: BoundTerm, literal: LiteralValue) -> bool:
        # Rows must match when: <----------Min----Max---X------->

        field_id = term.ref().field.field_id

        if self._can_contain_nulls(field_id) or self._can_contain_nans(field_id):
            return ROWS_MIGHT_NOT_MATCH

        if upper_bytes := self.upper_bounds.get(field_id):
            field = self._get_field(field_id)
            upper = _from_byte_buffer(field.field_type, upper_bytes)

            if upper < literal.value:
                return ROWS_MUST_MATCH

        return ROWS_MIGHT_NOT_MATCH

    def visit_less_than_or_equal(self, term: BoundTerm, literal: LiteralValue) -> bool:
        # Rows must match when: <----------Min----Max---X------->

        field_id = term.ref().field.field_id

        if self._can_contain_nulls(field_id) or self._can_contain_nans(field_id):
            return ROWS_MIGHT_NOT_MATCH

        if upper_bytes := self.upper_bounds.get(field_id):
            field = self._get_field(field_id)
            upper = _from_byte_buffer(field.field_type, upper_bytes)

            if upper <= literal.value:
                return ROWS_MUST_MATCH

        return ROWS_MIGHT_NOT_MATCH

    def visit_greater_than(self, term: BoundTerm, literal: LiteralValue) -> bool:
        # Rows must match when: <-------X---Min----Max---------->

        field_id = term.ref().field.field_id

        if self._can_contain_nulls(field_id) or self._can_contain_nans(field_id):
            return ROWS_MIGHT_NOT_MATCH

        if lower_bytes := self.lower_bounds.get(field_id):
            field = self._get_field(field_id)
            lower = _from_byte_buffer(field.field_type, lower_bytes)

            if self._is_nan(lower):
                # NaN indicates unreliable bounds.
                # See the _StrictMetricsEvaluator docs for more.
                return ROWS_MIGHT_NOT_MATCH

            if lower > literal.value:
                return ROWS_MUST_MATCH

        return ROWS_MIGHT_NOT_MATCH

    def visit_greater_than_or_equal(self, term: BoundTerm, literal: LiteralValue) -> bool:
        # Rows must match when: <-------X---Min----Max---------->
        field_id = term.ref().field.field_id

        if self._can_contain_nulls(field_id) or self._can_contain_nans(field_id):
            return ROWS_MIGHT_NOT_MATCH

        if lower_bytes := self.lower_bounds.get(field_id):
            field = self._get_field(field_id)
            lower = _from_byte_buffer(field.field_type, lower_bytes)

            if self._is_nan(lower):
                # NaN indicates unreliable bounds.
                # See the _StrictMetricsEvaluator docs for more.
                return ROWS_MIGHT_NOT_MATCH

            if lower >= literal.value:
                return ROWS_MUST_MATCH

        return ROWS_MIGHT_NOT_MATCH

    def visit_equal(self, term: BoundTerm, literal: LiteralValue) -> bool:
        # Rows must match when Min == X == Max
        field_id = term.ref().field.field_id

        if self._can_contain_nulls(field_id) or self._can_contain_nans(field_id):
            return ROWS_MIGHT_NOT_MATCH

        if (lower_bytes := self.lower_bounds.get(field_id)) and (upper_bytes := self.upper_bounds.get(field_id)):
            field = self._get_field(field_id)
            lower = _from_byte_buffer(field.field_type, lower_bytes)
            upper = _from_byte_buffer(field.field_type, upper_bytes)

            if lower != literal.value or upper != literal.value:
                return ROWS_MIGHT_NOT_MATCH
            else:
                return ROWS_MUST_MATCH

        return ROWS_MIGHT_NOT_MATCH

    def visit_not_equal(self, term: BoundTerm, literal: LiteralValue) -> bool:
        # Rows must match when X < Min or Max < X because it is not in the range
        field_id = term.ref().field.field_id

        if self._can_contain_nulls(field_id) or self._can_contain_nans(field_id):
            return ROWS_MUST_MATCH

        field = self._get_field(field_id)

        if lower_bytes := self.lower_bounds.get(field_id):
            lower = _from_byte_buffer(field.field_type, lower_bytes)

            if self._is_nan(lower):
                # NaN indicates unreliable bounds.
                # See the _StrictMetricsEvaluator docs for more.
                return ROWS_MIGHT_NOT_MATCH

            if lower > literal.value:
                return ROWS_MUST_MATCH

        if upper_bytes := self.upper_bounds.get(field_id):
            upper = _from_byte_buffer(field.field_type, upper_bytes)

            if upper < literal.value:
                return ROWS_MUST_MATCH

        return ROWS_MIGHT_NOT_MATCH

    def visit_in(self, term: BoundTerm, literals: set[L]) -> bool:
        field_id = term.ref().field.field_id

        if self._can_contain_nulls(field_id) or self._can_contain_nans(field_id):
            return ROWS_MIGHT_NOT_MATCH

        field = self._get_field(field_id)

        if (lower_bytes := self.lower_bounds.get(field_id)) and (upper_bytes := self.upper_bounds.get(field_id)):
            # similar to the implementation in eq, first check if the lower bound is in the set
            lower = _from_byte_buffer(field.field_type, lower_bytes)
            if lower not in literals:
                return ROWS_MIGHT_NOT_MATCH

            # check if the upper bound is in the set
            upper = _from_byte_buffer(field.field_type, upper_bytes)
            if upper not in literals:
                return ROWS_MIGHT_NOT_MATCH

            # finally check if the lower bound and the upper bound are equal
            if lower != upper:
                return ROWS_MIGHT_NOT_MATCH

            # All values must be in the set if the lower bound and the upper bound are
            # in the set and are equal.
            return ROWS_MUST_MATCH

        return ROWS_MIGHT_NOT_MATCH

    def visit_not_in(self, term: BoundTerm, literals: set[L]) -> bool:
        field_id = term.ref().field.field_id

        if self._can_contain_nulls(field_id) or self._can_contain_nans(field_id):
            return ROWS_MUST_MATCH

        field = self._get_field(field_id)

        if lower_bytes := self.lower_bounds.get(field_id):
            lower = _from_byte_buffer(field.field_type, lower_bytes)

            if self._is_nan(lower):
                # NaN indicates unreliable bounds.
                # See the StrictMetricsEvaluator docs for more.
                return ROWS_MIGHT_NOT_MATCH

            literals = {val for val in literals if lower <= val}
            if len(literals) == 0:
                return ROWS_MUST_MATCH

        if upper_bytes := self.upper_bounds.get(field_id):
            upper = _from_byte_buffer(field.field_type, upper_bytes)

            literals = {val for val in literals if upper >= val}

            if len(literals) == 0:
                return ROWS_MUST_MATCH

        return ROWS_MIGHT_NOT_MATCH

    def visit_starts_with(self, term: BoundTerm, literal: LiteralValue) -> bool:
        return ROWS_MIGHT_NOT_MATCH

    def visit_not_starts_with(self, term: BoundTerm, literal: LiteralValue) -> bool:
        return ROWS_MIGHT_NOT_MATCH

    def _get_field(self, field_id: int) -> NestedField:
        field = self.struct.field(field_id=field_id)
        if field is None:
            raise ValueError(f"Cannot find field, might be nested or missing: {field_id}")

        return field

    def _can_contain_nulls(self, field_id: int) -> bool:
        return (null_count := self.null_counts.get(field_id)) is not None and null_count > 0

    def _can_contain_nans(self, field_id: int) -> bool:
        return (nan_count := self.nan_counts.get(field_id)) is not None and nan_count > 0


class ResidualVisitor(BoundBooleanExpressionVisitor[BooleanExpression], ABC):
    """Finds the residuals for an Expression the partitions in the given PartitionSpec.

    A residual expression is made by partially evaluating an expression using partition values.
    For example, if a table is partitioned by day(utc_timestamp) and is read with a filter expression
    utc_timestamp > a and utc_timestamp < b, then there are 4 possible residuals expressions
    for the partition data, d:


    1. If d > day(a) and d &lt; day(b), the residual is always true
    2. If d == day(a) and d != day(b), the residual is utc_timestamp > a
    3. if d == day(b) and d != day(a), the residual is utc_timestamp < b
    4. If d == day(a) == day(b), the residual is utc_timestamp > a and utc_timestamp < b
    Partition data is passed using StructLike. Residuals are returned by residualFor(StructLike).
    """

    schema: Schema
    spec: PartitionSpec
    case_sensitive: bool
    expr: BooleanExpression

    def __init__(self, schema: Schema, spec: PartitionSpec, case_sensitive: bool, expr: BooleanExpression) -> None:
        self.schema = schema
        self.spec = spec
        self.case_sensitive = case_sensitive
        self.expr = expr

    def eval(self, partition_data: Record) -> BooleanExpression:
        self.struct = partition_data
        return visit(self.expr, visitor=self)

    def visit_true(self) -> BooleanExpression:
        return AlwaysTrue()

    def visit_false(self) -> BooleanExpression:
        return AlwaysFalse()

    def visit_not(self, child_result: BooleanExpression) -> BooleanExpression:
        return Not(child_result)

    def visit_and(self, left_result: BooleanExpression, right_result: BooleanExpression) -> BooleanExpression:
        return And(left_result, right_result)

    def visit_or(self, left_result: BooleanExpression, right_result: BooleanExpression) -> BooleanExpression:
        return Or(left_result, right_result)

    def visit_is_null(self, term: BoundTerm) -> BooleanExpression:
        if term.eval(self.struct) is None:
            return AlwaysTrue()
        else:
            return AlwaysFalse()

    def visit_not_null(self, term: BoundTerm) -> BooleanExpression:
        if term.eval(self.struct) is not None:
            return AlwaysTrue()
        else:
            return AlwaysFalse()

    def visit_is_nan(self, term: BoundTerm) -> BooleanExpression:
        val = term.eval(self.struct)
        if isinstance(val, SupportsFloat) and math.isnan(val):
            return self.visit_true()
        else:
            return self.visit_false()

    def visit_not_nan(self, term: BoundTerm) -> BooleanExpression:
        val = term.eval(self.struct)
        if isinstance(val, SupportsFloat) and not math.isnan(val):
            return self.visit_true()
        else:
            return self.visit_false()

    def visit_less_than(self, term: BoundTerm, literal: LiteralValue) -> BooleanExpression:
        if term.eval(self.struct) < literal.value:
            return self.visit_true()
        else:
            return self.visit_false()

    def visit_less_than_or_equal(self, term: BoundTerm, literal: LiteralValue) -> BooleanExpression:
        if term.eval(self.struct) <= literal.value:
            return self.visit_true()
        else:
            return self.visit_false()

    def visit_greater_than(self, term: BoundTerm, literal: LiteralValue) -> BooleanExpression:
        if term.eval(self.struct) > literal.value:
            return self.visit_true()
        else:
            return self.visit_false()

    def visit_greater_than_or_equal(self, term: BoundTerm, literal: LiteralValue) -> BooleanExpression:
        if term.eval(self.struct) >= literal.value:
            return self.visit_true()
        else:
            return self.visit_false()

    def visit_equal(self, term: BoundTerm, literal: LiteralValue) -> BooleanExpression:
        if term.eval(self.struct) == literal.value:
            return self.visit_true()
        else:
            return self.visit_false()

    def visit_not_equal(self, term: BoundTerm, literal: LiteralValue) -> BooleanExpression:
        if term.eval(self.struct) != literal.value:
            return self.visit_true()
        else:
            return self.visit_false()

    def visit_in(self, term: BoundTerm, literals: set[L]) -> BooleanExpression:
        if term.eval(self.struct) in literals:
            return self.visit_true()
        else:
            return self.visit_false()

    def visit_not_in(self, term: BoundTerm, literals: set[L]) -> BooleanExpression:
        if term.eval(self.struct) not in literals:
            return self.visit_true()
        else:
            return self.visit_false()

    def visit_starts_with(self, term: BoundTerm, literal: LiteralValue) -> BooleanExpression:
        eval_res = term.eval(self.struct)
        if eval_res is not None and str(eval_res).startswith(str(literal.value)):
            return AlwaysTrue()
        else:
            return AlwaysFalse()

    def visit_not_starts_with(self, term: BoundTerm, literal: LiteralValue) -> BooleanExpression:
        if not self.visit_starts_with(term, literal):
            return AlwaysTrue()
        else:
            return AlwaysFalse()

    def visit_bound_predicate(self, predicate: BoundPredicate) -> BooleanExpression:
        """
        If there is no strict projection or if it evaluates to false, then return the predicate.

        Get the strict projection and inclusive projection of this predicate in partition data,
        then use them to determine whether to return the original predicate. The strict projection
        returns true iff the original predicate would have returned true, so the predicate can be
        eliminated if the strict projection evaluates to true. Similarly the inclusive projection
        returns false iff the original predicate would have returned false, so the predicate can
        also be eliminated if the inclusive projection evaluates to false.

        """
        parts = self.spec.fields_by_source_id(predicate.term.ref().field.field_id)
        if parts == []:
            return predicate

        def struct_to_schema(struct: StructType) -> Schema:
            return Schema(*struct.fields)

        for part in parts:
            strict_projection = part.transform.strict_project(part.name, predicate)
            strict_result = None

            if strict_projection is not None:
                bound = strict_projection.bind(
                    struct_to_schema(self.spec.partition_type(self.schema)), case_sensitive=self.case_sensitive
                )
                if isinstance(bound, BoundPredicate):
                    strict_result = super().visit_bound_predicate(bound)
                else:
                    # if the result is not a predicate, then it must be a constant like alwaysTrue or alwaysFalse
                    strict_result = bound

            if isinstance(strict_result, AlwaysTrue):
                return AlwaysTrue()

            inclusive_projection = part.transform.project(part.name, predicate)
            inclusive_result = None
            if inclusive_projection is not None:
                bound_inclusive = inclusive_projection.bind(
                    struct_to_schema(self.spec.partition_type(self.schema)), case_sensitive=self.case_sensitive
                )
                if isinstance(bound_inclusive, BoundPredicate):
                    # using predicate method specific to inclusive
                    inclusive_result = super().visit_bound_predicate(bound_inclusive)
                else:
                    # if the result is not a predicate, then it must be a constant like alwaysTrue or
                    # alwaysFalse
                    inclusive_result = bound_inclusive
            if isinstance(inclusive_result, AlwaysFalse):
                return AlwaysFalse()

        return predicate

    def visit_unbound_predicate(self, predicate: UnboundPredicate) -> BooleanExpression:
        bound = predicate.bind(self.schema, case_sensitive=self.case_sensitive)

        if isinstance(bound, BoundPredicate):
            bound_residual = self.visit_bound_predicate(predicate=bound)
            if not isinstance(bound_residual, (AlwaysFalse, AlwaysTrue)):
                # replace inclusive original unbound predicate
                return predicate

            # use the non-predicate residual (e.g. alwaysTrue)
            return bound_residual

        # if binding didn't result in a Predicate, return the expression
        return bound


class ResidualEvaluator(ResidualVisitor):
    def residual_for(self, partition_data: Record) -> BooleanExpression:
        return self.eval(partition_data)


class UnpartitionedResidualEvaluator(ResidualEvaluator):
    # Finds the residuals for an Expression the partitions in the given PartitionSpec
    def __init__(self, schema: Schema, expr: BooleanExpression):
        super().__init__(schema=schema, spec=UNPARTITIONED_PARTITION_SPEC, expr=expr, case_sensitive=False)
        self.expr = expr

    def residual_for(self, partition_data: Record) -> BooleanExpression:
        return self.expr


def residual_evaluator_of(
    spec: PartitionSpec, expr: BooleanExpression, case_sensitive: bool, schema: Schema
) -> ResidualEvaluator:
    return (
        UnpartitionedResidualEvaluator(schema=schema, expr=expr)
        if spec.is_unpartitioned()
        else ResidualEvaluator(spec=spec, expr=expr, schema=schema, case_sensitive=case_sensitive)
    )
